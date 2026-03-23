#!/usr/bin/env bash
#SBATCH -J nh-sr-o4
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-sr-o4-%j.out

# Extract SR alignment and angular oscillation channels for full O4 (2023-2026).
# Parameterized: pass GPS_START and GPS_END as env vars or use defaults.
#
# Channels:
#   V1:SAT_SR_MAR_TX_SET          — SR marionette TX setpoint (alignment state)
#   V1:ASC_SR_TY_B1p_mag_10Hz_FS_mean — SR TY angular amplitude on B1p
#
# Output: outputs/sr_o4_${GPS_START}_${GPS_END}.csv
# Columns: gps_bin, sr_mar_tx_set, asc_sr_ty_b1p_mag
#
# Submit 4 parallel jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_sr_o4_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_sr_o4_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_sr_o4_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_sr_o4_extract.slurm

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_sr_o4_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== SR O4 extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    GPS range: $GPS_START – $GPS_END"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << PYEOF
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeries

warnings.filterwarnings("ignore")

WORKDIR   = Path.home() / "NOISEHOUND"
STAGE_DIR = Path(os.environ["STAGE_DIR"])
STAGE_DIR.mkdir(exist_ok=True)

HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
GPS_START = int(os.environ["GPS_START"])
GPS_END   = int(os.environ["GPS_END"])
BIN_S     = 3600   # 1-h bins

CHANNELS = {
    "V1:SAT_SR_MAR_TX_SET":              "sr_mar_tx_set",
    "V1:ASC_SR_TY_B1p_mag_10Hz_FS_mean": "asc_sr_ty_b1p_mag",
}

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
    ok = False
    row_data = {}
    for ch, key in CHANNELS.items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts = TimeSeries.read(str(gwf), ch, start=t0, end=t1)
            t_arr = np.asarray(ts.times.value, float)
            v_arr = np.asarray(ts.value, float)
            bins  = (t_arr // BIN_S).astype(int) * BIN_S
            df_ch = (pd.DataFrame({"gps_bin": bins, "val": v_arr})
                       .groupby("gps_bin")["val"].mean())
            for gps_bin, val in df_ch.items():
                records.append({"gps_bin": int(gps_bin), key: val})
            ok = True
        except Exception as e:
            print(f"  {ch} @ {gps_day}: {e}", flush=True)
    if ok:
        print(f"  day {gps_day} OK", flush=True)
    gwf.unlink(missing_ok=True)

if not records:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = pd.DataFrame(records)
df = df.groupby("gps_bin").mean().reset_index().sort_values("gps_bin")
out = WORKDIR / "outputs" / f"sr_o4_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
