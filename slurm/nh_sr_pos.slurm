#!/usr/bin/env bash
#SBATCH -J nh-sr-pos
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-sr-pos-%j.out

# Extract SR marionette setpoint channels from trend GWF
# for the disappearance epoch: Oct 2025 – Apr 2026.
# Channels: SAT_SR_MAR_TY_SET / SAT_SR_MAR_TX_SET
# (marionette angular setpoints — directly encode the LN3 vs LN3_ALIGNED step)
# Bins to 1-h averages, saves outputs/sr_pos_epoch.csv
# Expected: TX steps from ~+6 to ~-757 at Christmas shutdown.

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_sr_pos_$$"
mkdir -p "$STAGE_DIR"

echo "=== SR position extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

WORKDIR   = Path.home() / "NOISEHOUND"
STAGE_DIR = Path("/tmp/nh_sr_pos_$$".replace("$$", __import__("os").environ.get("SLURM_JOB_ID", "local")))
STAGE_DIR.mkdir(exist_ok=True)

HPSS_BASE     = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
GPS_EPOCH     = pd.Timestamp("1980-01-06", tz="UTC")
GPS_START     = 1443657600   # 2025-10-01
GPS_END       = 1459468800   # 2026-04-01
BIN_S         = 3600         # 1-h bins

CHANNELS = ["V1:SAT_SR_MAR_TY_SET", "V1:SAT_SR_MAR_TX_SET"]

from gwpy.timeseries import TimeSeries

def stage(gps_day: int) -> Path | None:
    year   = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  rfcp failed {gps_day}: {r.stderr[:80]}", flush=True)
            return None
    return local

day_start = (GPS_START // 86400) * 86400
day_end   = (GPS_END   // 86400) * 86400

records = []
for gps_day in range(day_start, day_end + 86400, 86400):
    gwf = stage(gps_day)
    if gwf is None:
        continue
    t0 = max(GPS_START, gps_day)
    t1 = min(GPS_END,   gps_day + 86400)
    row = {"t0": t0}
    ok  = False
    for ch in CHANNELS:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts = TimeSeries.read(str(gwf), ch, start=t0, end=t1)
            t_arr = np.asarray(ts.times.value, float)
            v_arr = np.asarray(ts.value, float)
            bins  = (t_arr // BIN_S).astype(int) * BIN_S
            df_ch = (pd.DataFrame({"gps_bin": bins, "val": v_arr})
                       .groupby("gps_bin")["val"].mean())
            key = "sr_mar_ty_set" if "TY" in ch else "sr_mar_tx_set"
            for gps_bin, val in df_ch.items():
                records.append({"gps_bin": int(gps_bin), key: val})
            ok = True
        except Exception as e:
            print(f"  {ch} @ {gps_day}: {e}", flush=True)
    if ok:
        print(f"  day {gps_day} OK", flush=True)
    local = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if local.exists():
        local.unlink()

if not records:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = pd.DataFrame(records)
df = df.groupby("gps_bin").mean().reset_index().sort_values("gps_bin")
out = WORKDIR / "outputs" / "sr_pos_epoch.csv"
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins → {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
