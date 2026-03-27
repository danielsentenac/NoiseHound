#!/usr/bin/env bash
#SBATCH -J nh-tower-wide
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-tower-wide-%j.out

# Extract a 240-second window (±120 s) around the catalog Omicron peak.
# Stages three 100-second raw files one at a time to limit peak disk usage.
# Channels: Hrec_hoft + all tower ASC ERR and CORR + DARM.
#
# Output: outputs/tower_glitch_probe_wide_1415578745.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"

CENTER_GPS=1415578745
WIN=120            # seconds on each side → 240 s total

export CENTER_GPS WIN
export STAGE_DIR="/tmp/nh_tower_wide_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Tower wide probe started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    CENTER_GPS = $CENTER_GPS   window = ±${WIN}s"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeriesDict

warnings.filterwarnings("ignore")

WORKDIR    = Path.home() / "NOISEHOUND"
STAGE_DIR  = Path(os.environ["STAGE_DIR"])
HPSS_BASE  = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
CENTER_GPS = float(os.environ["CENTER_GPS"])
WIN        = float(os.environ["WIN"])
T0         = CENTER_GPS - WIN
T1         = CENTER_GPS + WIN

CHANNELS = [
    "V1:Hrec_hoft_16384Hz",
    "V1:LSC_DARM_CORR",
    "V1:ASC_BS_TX_ERR",
    "V1:ASC_BS_TY_ERR",
    "V1:ASC_PR_TX_ERR",
    "V1:ASC_PR_TY_ERR",
    "V1:ASC_SR_TX_ERR",
    "V1:ASC_SR_TY_ERR",
    "V1:ASC_NI_TX_CORR",
    "V1:ASC_NI_TY_CORR",
    "V1:ASC_WI_TX_CORR",
    "V1:ASC_WI_TY_CORR",
    "V1:ASC_NE_TX_CORR",
    "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR",
    "V1:ASC_WE_TY_CORR",
]

# List all 100-second file blocks needed
file_starts = []
s = int(T0 // 100) * 100
while s <= int(T1 // 100) * 100:
    file_starts.append(s)
    s += 100
print(f"Files to stage: {file_starts}", flush=True)

TARGET_SR = 4096
rows = {}

for file_start in file_starts:
    subdir = str(file_start)[:4]
    fname  = f"V-raw-{file_start}-100.gwf"
    remote = f"{HPSS_BASE}/{subdir}/{fname}"
    local  = STAGE_DIR / fname

    print(f"\nStaging {fname} ...", flush=True)
    r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  FAILED: {r.stderr[:200]}", flush=True)
        continue
    print(f"  Staged ({local.stat().st_size/1e9:.2f} GB)", flush=True)

    seg_t0 = max(T0, float(file_start))
    seg_t1 = min(T1, float(file_start + 100))
    if seg_t0 >= seg_t1:
        local.unlink(missing_ok=True)
        continue

    print(f"  Reading t={seg_t0:.0f}–{seg_t1:.0f} ...", flush=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tsd = TimeSeriesDict.read(str(local), CHANNELS,
                                      start=seg_t0, end=seg_t1, nproc=1)
        for ch, ts in tsd.items():
            sr = ts.sample_rate.value
            if sr > TARGET_SR:
                ts = ts.resample(TARGET_SR)
            t = ts.times.value
            v = ts.value.astype(np.float32)
            chunk = pd.DataFrame({"gps": t, ch: v}).set_index("gps")
            rows[ch] = pd.concat([rows[ch], chunk]) if ch in rows else chunk
        print(f"  Read {len(tsd)} channels", flush=True)
    except Exception as e:
        print(f"  Error: {e}", flush=True)

    # Free disk immediately
    local.unlink(missing_ok=True)

if not rows:
    print("No data extracted!", flush=True)
    sys.exit(1)

merged = pd.concat(rows.values(), axis=1, join="outer").sort_index()
merged.index.name = "gps"
out = WORKDIR / "outputs" / "tower_glitch_probe_wide_1415578745.csv"
merged.reset_index().to_csv(out, index=False)
print(f"\nSaved {len(merged)} samples × {len(merged.columns)} channels → {out}", flush=True)
print(f"Time range: {merged.index.min():.3f} – {merged.index.max():.3f}", flush=True)

PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
