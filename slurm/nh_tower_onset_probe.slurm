#!/usr/bin/env bash
#SBATCH -J nh-tower-onset
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-tower-onset-%j.out

# Extract all tower ERR/CORR channels + Hrec_hoft around a pre-Christmas 2025 glitch.
#
# Step 1 result: GPS 1446106739.098 (2025-11-02 08:18:59 UTC, SNR=499, f=40.4 Hz, O4b)
# is selected as the reference pre-Christmas 2025 event.
#
# Usage (default: 2025-11-02 event):
#   sbatch slurm/nh_tower_onset_probe.slurm
#
# Override glitch GPS:
#   sbatch --export=GLITCH_GPS=1446106739 slurm/nh_tower_onset_probe.slurm
#
# Window: ±60 s around catalog GPS (two 100-s raw frames staged).
#
# Channels extracted:
#   ASC: NI/WI/NE/WE (CORR), BS/PR (ERR+CORR), SR (B1p_DCP+ERR+OUT), DIFFp (ERR+CORR)
#   IB:  Sc_IB_MAR_{TX,TY}_CORR (10 kHz → downsampled)
#   DET: Sa_OB_F0_{TY,X}_CORR_500Hz
#   Ref: V1:Hrec_hoft_16384Hz, V1:LSC_DARM_CORR
#
# Output:
#   outputs/tower_onset_probe_{GLITCH_GPS}.csv   — all channels, 120 s window

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
export GLITCH_GPS="${GLITCH_GPS:-1446106739}"
export HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4c/raw"
export STAGE_DIR="/tmp/nh_tower_onset_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Tower onset probe started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    GLITCH_GPS = $GLITCH_GPS"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess, os, re
from pathlib import Path
import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeries, TimeSeriesDict
from gwpy.io.gwf import iter_channel_names

warnings.filterwarnings("ignore")

WORKDIR    = Path.home() / "NOISEHOUND"
STAGE_DIR  = Path(os.environ["STAGE_DIR"])
HPSS_BASE  = os.environ.get("HPSS_BASE",
             "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4c/raw")
GLITCH_GPS = float(os.environ["GLITCH_GPS"])

WIN_BEFORE = 60.0   # seconds before catalog GPS
WIN_AFTER  = 60.0   # seconds after  catalog GPS
T0 = GLITCH_GPS - WIN_BEFORE
T1 = GLITCH_GPS + WIN_AFTER

def gps_to_dir(gps):
    """Return HPSS subdirectory: first 4 digits of GPS integer."""
    return str(int(gps))[:4]

def find_raw_files(t_start, t_end, stage_dir):
    """Stage all 100-s raw GWF frames that overlap [t_start, t_end]."""
    block_start = int(t_start // 100) * 100
    block_end   = int(t_end   // 100) * 100
    staged = []
    for start in range(block_start, block_end + 1, 100):
        subdir = gps_to_dir(start)
        fname  = f"V-raw-{start}-100.gwf"
        remote = f"{HPSS_BASE}/{subdir}/{fname}"
        local  = stage_dir / fname
        if local.exists():
            staged.append(local)
            continue
        print(f"  Staging {remote} ...", flush=True)
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
        if r.returncode == 0:
            staged.append(local)
            print(f"  OK ({local.stat().st_size/1e9:.2f} GB)", flush=True)
        else:
            print(f"  FAILED: {r.stderr[:200]}", flush=True)
    return staged

# ── Stage frames covering T0..T1 ──────────────────────────────────────────────
all_files = find_raw_files(T0, T1, STAGE_DIR)
if not all_files:
    print("No raw files staged — aborting.", flush=True)
    sys.exit(1)

# ── Build channel selection ────────────────────────────────────────────────────
print(f"\nListing channels from {all_files[0].name} ...", flush=True)
all_ch = list(iter_channel_names(str(all_files[0])))
print(f"Total channels in raw frame: {len(all_ch)}", flush=True)

# Exact channel names we want (fall back gracefully if absent)
WANT = [
    # NI/WI: ASC_NI/WI_TX/TY_CORR = identically zero in O4c (ASC DOF matrix disabled for them).
    # Actual angular drive goes through the Sc_ suspension namespace at 10 kHz.
    "V1:Sc_NI_MAR_TX_CORR", "V1:Sc_NI_MAR_TY_CORR",
    "V1:Sc_WI_MAR_TX_CORR", "V1:Sc_WI_MAR_TY_CORR",
    # NE, WE: ASC_*_CORR still active at 2000 Hz (confirmed)
    "V1:ASC_NE_TX_CORR", "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR", "V1:ASC_WE_TY_CORR",
    # BS: ERR + CORR
    "V1:ASC_BS_TX_ERR",  "V1:ASC_BS_TY_ERR",
    "V1:ASC_BS_TX_CORR", "V1:ASC_BS_TY_CORR",
    # PR: ERR + CORR
    "V1:ASC_PR_TX_ERR",  "V1:ASC_PR_TY_ERR",
    "V1:ASC_PR_TX_CORR", "V1:ASC_PR_TY_CORR",
    # SR: B1p_DCP = raw optical sensor (2000 Hz) — direct SR angular readout before loop filter.
    # Step-2 ranked ASC_SR_TY_DCP_mag_B1_mag is the 10 Hz envelope of this signal.
    # Also keep ERR (post input-matrix) and OUT (actuator correction).
    "V1:ASC_SR_TX_B1p_DCP", "V1:ASC_SR_TY_B1p_DCP",
    "V1:ASC_SR_TX_ERR",     "V1:ASC_SR_TY_ERR",
    "V1:ASC_SR_TX_OUT",     "V1:ASC_SR_TY_OUT",
    # IB: marionette CORR at 10 kHz (downsampled to TARGET_SR=4096 Hz)
    "V1:Sc_IB_MAR_TX_CORR", "V1:Sc_IB_MAR_TY_CORR",
    # DET = OB (Output Bench): no Sc_DET_MAR_* (confirmed). Sa_OB_F0 at 500 Hz.
    "V1:Sa_OB_F0_TY_CORR_500Hz",
    "V1:Sa_OB_F0_X_CORR_500Hz",
    # DIFFp = differential angular DOF (NE−WE combination); Step-2 rank 19; all 2000 Hz.
    "V1:ASC_DIFFp_TX_ERR",  "V1:ASC_DIFFp_TY_ERR",
    "V1:ASC_DIFFp_TX_CORR", "V1:ASC_DIFFp_TY_CORR",
    # Reference: h(t) at 16384 Hz and DARM correction
    "V1:Hrec_hoft_16384Hz",
    "V1:LSC_DARM_CORR",
    # Step-2 top ranked channels (raw equivalents of 1 Hz trend products):
    #   rank 1: LSC_DARM_IMC_LINE_mag_100Hz_mean  → IMC 100 Hz calibration line amplitude in DARM
    #   rank 2: LSC_DARM_PSTAB0_COUPLING_100Hz_mean → PSTAB0 coupling at 100 Hz in DARM
    #   rank 6: LSC_DARM_SIB1_LINE_mag_100Hz_mean → SIB1 line amplitude in DARM
    #   rank 7: Sc_WI_FF50HZ_P_ERR_mean → WI 50 Hz feedforward phase error (first non-DARM)
    "V1:LSC_DARM_IMC_LINE_mag_100Hz",
    "V1:LSC_DARM_PSTAB0_COUPLING_100Hz",
    "V1:LSC_DARM_SIB1_LINE_mag_100Hz",
    "V1:Sc_WI_FF50HZ_P_ERR",
]

ch_set    = set(all_ch)
selected  = [c for c in WANT if c in ch_set]
missing   = [c for c in WANT if c not in ch_set]
print(f"\nSelected {len(selected)} / {len(WANT)} channels:", flush=True)
for c in selected:
    print(f"  {c}", flush=True)
if missing:
    print(f"\nMissing (will be absent from CSV):", flush=True)
    for c in missing:
        print(f"  {c}", flush=True)

if not selected:
    print("No matching channels — aborting.", flush=True)
    sys.exit(1)

# ── Extract time series from all staged frames ─────────────────────────────────
print(f"\nExtracting {WIN_BEFORE+WIN_AFTER:.0f} s window: GPS {T0:.1f} – {T1:.1f}", flush=True)

rows = {}
for gwf in sorted(all_files):
    parts      = gwf.stem.split("-")   # V-raw-{start}-100
    file_start = int(parts[2])
    file_end   = file_start + 100
    seg_t0 = max(T0, file_start)
    seg_t1 = min(T1, file_end)
    if seg_t0 >= seg_t1:
        continue
    print(f"\n  Reading {gwf.name} (segment {seg_t0:.1f}–{seg_t1:.1f}) ...", flush=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tsd = TimeSeriesDict.read(str(gwf), selected,
                                      start=seg_t0, end=seg_t1, nproc=1)
        for ch, ts in tsd.items():
            rows[ch] = rows[ch].append(ts) if ch in rows else ts
        print(f"    Read {len(tsd)} channels", flush=True)
    except Exception as e:
        print(f"    Error: {e}", flush=True)

if not rows:
    print("No data extracted!", flush=True)
    sys.exit(1)

# ── Resample and save ─────────────────────────────────────────────────────────
TARGET_SR = 4096   # Hz – sufficient for 150 Hz bandpass; downsamples 16384 Hz h(t)

dfs = []
for ch, ts in rows.items():
    sr = ts.sample_rate.value
    if sr > TARGET_SR:
        ts = ts.resample(TARGET_SR)
    dfs.append(
        pd.DataFrame({"gps": ts.times.value, ch: ts.value.astype(np.float32)})
        .set_index("gps")
    )

merged = pd.concat(dfs, axis=1, join="outer").sort_index()
merged.index.name = "gps"
(WORKDIR / "outputs").mkdir(exist_ok=True)
out = WORKDIR / "outputs" / f"tower_onset_probe_{int(GLITCH_GPS)}.csv"
merged.reset_index().to_csv(out, index=False)
print(f"\nSaved {len(merged)} samples × {len(merged.columns)} channels → {out}",
      flush=True)
print(f"Time range: {merged.index.min():.3f} – {merged.index.max():.3f}", flush=True)

PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
