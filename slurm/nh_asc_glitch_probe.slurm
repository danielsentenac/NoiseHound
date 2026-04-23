#!/usr/bin/env bash
#SBATCH -J nh-asc-probe
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-asc-probe-%j.out

# Extract ASC and LSC error signals around one 25-min glitch from raw data.
#
# Usage (default: Nov 14 2024 glitch at GPS 1415579346):
#   sbatch slurm/nh_asc_glitch_probe.slurm
#
# Override glitch GPS time:
#   sbatch --export=GLITCH_GPS=1415579346 slurm/nh_asc_glitch_probe.slurm
#
# Output:
#   outputs/asc_glitch_probe_{GLITCH_GPS}.csv   — all ASC/LSC channels, 10 s window
#   outputs/asc_glitch_channels_{GLITCH_GPS}.txt — full channel list from raw file

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GLITCH_GPS="${GLITCH_GPS:-1415579346}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
export STAGE_DIR="/tmp/nh_asc_probe_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== ASC glitch probe started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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
HPSS_BASE  = os.environ.get("HPSS_BASE", "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw")
GLITCH_GPS = float(os.environ["GLITCH_GPS"])

# ── Window: 5 s before glitch, 5 s after ─────────────────────────────────────
WIN_BEFORE = 5.0   # seconds
WIN_AFTER  = 5.0
T0 = GLITCH_GPS - WIN_BEFORE
T1 = GLITCH_GPS + WIN_AFTER

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

def gps_to_dir(gps):
    """Return the HPSS subdirectory name: first 4 digits of GPS."""
    return str(int(gps))[:4]

def find_raw_file(gps_t, stage_dir):
    """
    Find and stage the raw GWF file(s) covering gps_t.
    Raw files are 100 s: V-raw-{start}-100.gwf where start = floor(gps/100)*100
    """
    subdir = gps_to_dir(gps_t)
    # Try the 100-s block containing gps_t and also the one before (overlap safety)
    candidates = [
        int(gps_t // 100) * 100,
        int(gps_t // 100) * 100 - 100,
    ]
    staged = []
    for start in candidates:
        fname = f"V-raw-{start}-100.gwf"
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
            print(f"  → OK ({local.stat().st_size/1e9:.2f} GB)", flush=True)
        else:
            print(f"  → FAILED: {r.stderr[:120]}", flush=True)
    return staged

# ── Stage files covering T0..T1 ───────────────────────────────────────────────
files_t0 = find_raw_file(T0, STAGE_DIR)
files_t1 = find_raw_file(T1, STAGE_DIR)
all_files = list(dict.fromkeys(files_t0 + files_t1))   # deduplicate, preserve order

if not all_files:
    print("No raw files staged — aborting.", flush=True)
    sys.exit(1)

# ── Get channel list ──────────────────────────────────────────────────────────
print(f"\nListing channels from {all_files[0].name} ...", flush=True)
all_channels = list(iter_channel_names(str(all_files[0])))
print(f"Total channels in raw file: {len(all_channels)}", flush=True)

# Save full list
ch_out = WORKDIR / "outputs" / f"asc_glitch_channels_{int(GLITCH_GPS)}.txt"
with open(ch_out, "w") as f:
    f.write("\n".join(all_channels))
print(f"Channel list saved → {ch_out}", flush=True)

# ── Select ASC and LSC error/sensor channels ─────────────────────────────────
# Mirrors: NE, WE, NI, WI, PR, BS, SR, IB, DET
# Want: ASC *_TX* *_TY* error/sensor; LSC DARM/MICH/PRCL/SRCL/CARM error signals
MIRRORS = ["NE", "WE", "NI", "WI", "PR", "BS", "SR", "IB", "DET"]

patterns = (
    [re.compile(rf"V1:ASC_{m}_(TX|TY)", re.IGNORECASE) for m in MIRRORS]
    + [re.compile(r"V1:LSC_(DARM|MICH|PRCL|SRCL|CARM|DARM_ERR|B1_DC)", re.IGNORECASE)]
    + [re.compile(r"V1:h_", re.IGNORECASE)]
)

selected = [c for c in all_channels
            if any(p.search(c) for p in patterns)]

print(f"\nSelected channels ({len(selected)}):", flush=True)
for c in sorted(selected):
    print(f"  {c}", flush=True)

if not selected:
    print("\nNo matching channels found! Saving full ASC/LSC list only.", flush=True)
    # Fallback: just show all with ASC or LSC in name
    asc_lsc = [c for c in all_channels if "ASC" in c or "LSC" in c]
    print(f"  Total ASC/LSC channels: {len(asc_lsc)}", flush=True)
    for c in sorted(asc_lsc)[:50]:
        print(f"    {c}", flush=True)
    sys.exit(0)

# ── Extract time series from all staged files ─────────────────────────────────
print(f"\nExtracting {WIN_BEFORE+WIN_AFTER:.0f} s window around GPS {GLITCH_GPS:.1f} ...", flush=True)
print(f"  T0 = {T0:.1f}  T1 = {T1:.1f}", flush=True)

rows = {}   # channel -> TimeSeries

for gwf in all_files:
    print(f"\n  Reading {gwf.name} ...", flush=True)
    # Clamp to file bounds
    file_start = int(gwf.stem.split("-")[2])
    file_end   = file_start + 100
    seg_t0 = max(T0, file_start)
    seg_t1 = min(T1, file_end)
    if seg_t0 >= seg_t1:
        continue
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tsd = TimeSeriesDict.read(str(gwf), selected,
                                      start=seg_t0, end=seg_t1, nproc=1)
        for ch, ts in tsd.items():
            if ch not in rows:
                rows[ch] = ts
            else:
                rows[ch] = rows[ch].append(ts)
        print(f"    Read {len(tsd)} channels", flush=True)
    except Exception as e:
        print(f"    Error reading {gwf.name}: {e}", flush=True)

if not rows:
    print("No data extracted!", flush=True)
    sys.exit(1)

# ── Resample to common grid and save ─────────────────────────────────────────
# Find max sample rate among selected channels (for info)
rates = {ch: ts.sample_rate.value for ch, ts in rows.items()}
print(f"\nSample rates found: {set(rates.values())}", flush=True)

# Save as wide CSV: time column + one column per channel
# Downsample high-rate channels to 4096 Hz for manageable file size
TARGET_SR = 4096   # Hz

dfs = []
for ch, ts in rows.items():
    sr = ts.sample_rate.value
    if sr > TARGET_SR:
        ts = ts.resample(TARGET_SR)
    t = ts.times.value   # GPS seconds
    v = ts.value
    dfs.append(pd.DataFrame({
        "gps": t,
        ch:    v.astype(np.float32),
    }).set_index("gps"))

if dfs:
    merged = pd.concat(dfs, axis=1, join="outer").sort_index()
    merged.index.name = "gps"
    out = WORKDIR / "outputs" / f"asc_glitch_probe_{int(GLITCH_GPS)}.csv"
    merged.reset_index().to_csv(out, index=False)
    print(f"\nSaved {len(merged)} samples × {len(merged.columns)} channels → {out}", flush=True)
    print(f"Time range: {merged.index.min():.3f} – {merged.index.max():.3f}", flush=True)

PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
