#!/usr/bin/env bash
#SBATCH -J nh-tower-probe
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-tower-probe-%j.out

# Extract a 30-second window around the 25-min glitch centered 12 s BEFORE
# the catalog Omicron peak (to capture the glitch onset).
# Channels: Hrec_hoft (reference) + all tower ASC ERR and CORR + DARM.
# Same raw file as the original probe (V-raw-1415578700-100.gwf).
#
# Output: outputs/tower_glitch_probe_1415578733.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"

# Center the window 12 s before the catalog GPS to catch glitch onset
CENTER_GPS=1415578733
WIN=15            # seconds on each side → 30 s total

export CENTER_GPS WIN
export STAGE_DIR="/tmp/nh_tower_probe_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Tower glitch probe started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    CENTER_GPS = $CENTER_GPS   window = ±${WIN}s"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeriesDict
from gwpy.io.gwf import iter_channel_names

warnings.filterwarnings("ignore")

WORKDIR    = Path.home() / "NOISEHOUND"
STAGE_DIR  = Path(os.environ["STAGE_DIR"])
HPSS_BASE  = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
CENTER_GPS = float(os.environ["CENTER_GPS"])
WIN        = float(os.environ["WIN"])
T0         = CENTER_GPS - WIN
T1         = CENTER_GPS + WIN

# ── Channels ─────────────────────────────────────────────────────────────────
CHANNELS = [
    # h(t) reference — the channel Omicron actually uses
    "V1:Hrec_hoft_16384Hz",
    # DARM
    "V1:LSC_DARM_CORR",
    "V1:LSC_DARM_ERR",
    # ASC error signals (angular sensors)
    "V1:ASC_BS_TX_ERR",
    "V1:ASC_BS_TY_ERR",
    "V1:ASC_PR_TX_ERR",
    "V1:ASC_PR_TY_ERR",
    "V1:ASC_SR_TX_ERR",
    "V1:ASC_SR_TY_ERR",
    # ASC correction signals (actuator outputs)
    "V1:ASC_NI_TX_CORR",
    "V1:ASC_NI_TY_CORR",
    "V1:ASC_WI_TX_CORR",
    "V1:ASC_WI_TY_CORR",
    "V1:ASC_NE_TX_CORR",
    "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR",
    "V1:ASC_WE_TY_CORR",
    # Sc_WI channels (non-causal reference from Step 2)
    "V1:Sc_WI_FF50HZ_P_ERR",
    "V1:Sc_WI_MIR_TX",
    "V1:Sc_WI_MIR_TY",
]

# ── Stage raw file ────────────────────────────────────────────────────────────
file_start = int(T0 // 100) * 100
subdir     = str(file_start)[:4]
fname      = f"V-raw-{file_start}-100.gwf"
remote     = f"{HPSS_BASE}/{subdir}/{fname}"
local      = STAGE_DIR / fname

print(f"Staging {fname} ...", flush=True)
r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
if r.returncode != 0:
    print(f"FAILED: {r.stderr[:200]}", flush=True)
    sys.exit(1)
print(f"  Staged ({local.stat().st_size/1e9:.2f} GB)", flush=True)

# Check if T1 spans a second file
if int(T1 // 100) * 100 != file_start:
    file_start2 = int(T1 // 100) * 100
    fname2      = f"V-raw-{file_start2}-100.gwf"
    remote2     = f"{HPSS_BASE}/{str(file_start2)[:4]}/{fname2}"
    local2      = STAGE_DIR / fname2
    print(f"Staging {fname2} ...", flush=True)
    r2 = subprocess.run(["rfcp", remote2, str(local2)], capture_output=True, text=True)
    if r2.returncode == 0:
        print(f"  Staged ({local2.stat().st_size/1e9:.2f} GB)", flush=True)
    else:
        print(f"  Second file failed (OK if within single 100s block)", flush=True)

all_files = sorted(STAGE_DIR.glob("V-raw-*.gwf"))

# ── Extract channels ──────────────────────────────────────────────────────────
# Downsample Hrec to 4096 Hz to keep file size manageable
TARGET_SR = 4096

rows = {}
for gwf in all_files:
    file_s = int(gwf.stem.split("-")[2])
    seg_t0 = max(T0, file_s)
    seg_t1 = min(T1, file_s + 100)
    if seg_t0 >= seg_t1:
        continue
    print(f"\nReading {gwf.name} (t={seg_t0:.0f}–{seg_t1:.0f}) ...", flush=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tsd = TimeSeriesDict.read(str(gwf), CHANNELS,
                                      start=seg_t0, end=seg_t1, nproc=1)
        for ch, ts in tsd.items():
            sr = ts.sample_rate.value
            if sr > TARGET_SR:
                ts = ts.resample(TARGET_SR)
            t  = ts.times.value
            v  = ts.value.astype(np.float32)
            chunk = pd.DataFrame({"gps": t, ch: v}).set_index("gps")
            if ch not in rows:
                rows[ch] = chunk
            else:
                rows[ch] = pd.concat([rows[ch], chunk])
        print(f"  Read {len(tsd)} channels", flush=True)
    except Exception as e:
        print(f"  Error: {e}", flush=True)

if not rows:
    print("No data extracted!", flush=True)
    sys.exit(1)

merged = pd.concat(rows.values(), axis=1, join="outer").sort_index()
merged.index.name = "gps"
out = WORKDIR / "outputs" / f"tower_glitch_probe_{int(CENTER_GPS)}.csv"
merged.reset_index().to_csv(out, index=False)
print(f"\nSaved {len(merged)} samples × {len(merged.columns)} channels → {out}", flush=True)
print(f"Time range: {merged.index.min():.3f} – {merged.index.max():.3f}", flush=True)

PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
