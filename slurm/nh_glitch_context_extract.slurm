#!/usr/bin/env bash
#SBATCH -J nh-glitch-ctx
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-glitch-ctx-%j.out

# Extract 30-55 Hz band-RMS of DARM_CORR + key ASC CORR/ERR channels
# one value per 100-second raw file, over ±30 min around the glitch.
# Files are staged one at a time to keep disk usage below ~6 GB.
#
# Usage:
#   sbatch slurm/nh_glitch_context_extract.slurm
#
# Output:
#   outputs/glitch_context_{GLITCH_GPS}.csv
#   One row per 100-second block: gps_mid + band-RMS per channel

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GLITCH_GPS="${GLITCH_GPS:-1415578745}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
WIN_SEC=1800    # ±30 min
export GLITCH_GPS WIN_SEC
export STAGE_DIR="/tmp/nh_glitch_ctx_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Glitch context extract started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    GLITCH_GPS = $GLITCH_GPS   window = ±${WIN_SEC}s"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from gwpy.timeseries import TimeSeriesDict

warnings.filterwarnings("ignore")

WORKDIR    = Path.home() / "NOISEHOUND"
STAGE_DIR  = Path(os.environ["STAGE_DIR"])
HPSS_BASE  = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
GLITCH_GPS = float(os.environ["GLITCH_GPS"])
WIN_SEC    = float(os.environ.get("WIN_SEC", 1800))

T0_WINDOW = GLITCH_GPS - WIN_SEC
T1_WINDOW = GLITCH_GPS + WIN_SEC

# ── Channels to extract ───────────────────────────────────────────────────────
CHANNELS = [
    "V1:LSC_DARM_CORR",
    "V1:ASC_NI_TX_CORR",
    "V1:ASC_NI_TY_CORR",
    "V1:ASC_WI_TX_CORR",
    "V1:ASC_WI_TY_CORR",
    "V1:ASC_NE_TX_CORR",
    "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR",
    "V1:ASC_WE_TY_CORR",
    "V1:ASC_SR_TX_ERR",
    "V1:ASC_SR_TY_ERR",
    "V1:ASC_BS_TX_ERR",
    "V1:ASC_PR_TX_ERR",
]

# ── Bandpass ──────────────────────────────────────────────────────────────────
def band_rms(v, fs, flo=30, fhi=55, order=4):
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    bp  = sosfiltfilt(sos, v)
    return float(np.sqrt(np.mean(bp**2)))

# ── List of 100-second file starts to process ─────────────────────────────────
file_starts = []
s = int(T0_WINDOW // 100) * 100
while s < T1_WINDOW:
    file_starts.append(s)
    s += 100
print(f"Files to process: {len(file_starts)}", flush=True)

# ── Process one file at a time ────────────────────────────────────────────────
records = []

for start in file_starts:
    end   = start + 100
    subdir = str(start)[:4]
    fname  = f"V-raw-{start}-100.gwf"
    remote = f"{HPSS_BASE}/{subdir}/{fname}"
    local  = STAGE_DIR / fname

    # Stage
    print(f"\n  Staging {fname} ...", flush=True)
    r = subprocess.run(["rfcp", remote, str(local)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    FAILED: {r.stderr[:100]}", flush=True)
        continue
    print(f"    Staged ({local.stat().st_size/1e9:.2f} GB)", flush=True)

    # Read channels (clamp segment to file bounds and window bounds)
    seg_t0 = max(T0_WINDOW, start)
    seg_t1 = min(T1_WINDOW, end)
    row = {"gps_start": seg_t0, "gps_end": seg_t1,
           "gps_mid": (seg_t0 + seg_t1) / 2}

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tsd = TimeSeriesDict.read(str(local), CHANNELS,
                                      start=seg_t0, end=seg_t1, nproc=1)
        for ch, ts in tsd.items():
            v  = ts.value.astype(float)
            fs = float(ts.sample_rate.value)
            if len(v) < 100:
                continue
            rms = band_rms(v, fs)
            row[ch] = rms
        print(f"    Computed RMS for {len(tsd)} channels", flush=True)
    except Exception as e:
        print(f"    Read error: {e}", flush=True)

    records.append(row)

    # Delete staged file immediately to free disk
    local.unlink(missing_ok=True)

# ── Save ─────────────────────────────────────────────────────────────────────
if records:
    out_df = pd.DataFrame(records).sort_values("gps_mid")
    out = WORKDIR / "outputs" / f"glitch_context_{int(GLITCH_GPS)}.csv"
    out_df.to_csv(out, index=False)
    print(f"\nSaved {len(out_df)} rows → {out}", flush=True)
    print(out_df[["gps_mid"] + [c for c in out_df.columns if "DARM" in c or "NI" in c]].to_string(), flush=True)
else:
    print("No records — check HPSS access.", flush=True)
    sys.exit(1)

PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
