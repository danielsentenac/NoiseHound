#!/usr/bin/env bash
#SBATCH -J nh-lock-o4
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 04:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-lock-o4-%j.out

# Extract V1:META_ITF_LOCK_index for one O4 quarter.
# Pass GPS_START and GPS_END via --export.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_lock_o4_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_lock_o4_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_lock_o4_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_lock_o4_extract.slurm
#
# Output: outputs/lock_o4_${GPS_START}_${GPS_END}.csv
# Columns: gps_bin, lock_mean

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_lock_o4_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Lock O4 extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    GPS range: $GPS_START – $GPS_END"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
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

CHANNEL = "V1:META_ITF_LOCK_index"


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
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = TimeSeries.read(str(gwf), CHANNEL, start=t0, end=t1)
        t_arr = np.asarray(ts.times.value, float)
        v_arr = np.asarray(ts.value, float)
        bins  = (t_arr // BIN_S).astype(int) * BIN_S
        df_ch = (pd.DataFrame({"gps_bin": bins, "val": v_arr})
                   .groupby("gps_bin")["val"].mean())
        for gps_bin, val in df_ch.items():
            records.append({"gps_bin": int(gps_bin), "lock_mean": val})
        print(f"  day {gps_day} OK", flush=True)
    except Exception as e:
        print(f"  {CHANNEL} @ {gps_day}: {e}", flush=True)
    gwf.unlink(missing_ok=True)

if not records:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = pd.DataFrame(records)
df = df.groupby("gps_bin").mean().reset_index().sort_values("gps_bin")
out = WORKDIR / "outputs" / f"lock_o4_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
