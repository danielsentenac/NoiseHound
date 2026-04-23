#!/usr/bin/env bash
#SBATCH -J nh-hrec-range-bns
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-hrec-range-bns-%j.out

# Extract hourly Hrec_Range_BNS from Virgo trend GWF files over a chosen GPS range.
#
# Pass GPS_START / GPS_END with --export, for example:
#   sbatch --export=GPS_START=1451779200,GPS_END=1455667200 slurm/nh_hrec_range_bns_extract.slurm
#
# Output:
#   outputs/hrec_range_bns_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1451779200}"
GPS_END="${GPS_END:-1455667200}"
export STAGE_DIR="/tmp/nh_hrec_range_bns_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Hrec_Range_BNS extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    GPS range: $GPS_START – $GPS_END"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import os
import subprocess
import sys
import warnings
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
BIN_S     = 3600
CHANNEL   = "V1:Hrec_Range_BNS"


def stage(gps_day: int) -> Path | None:
    year   = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  rfcp failed {gps_day}: {r.stderr[:120]}", flush=True)
            return None
    return local


def read_hourly(gwf: Path, channel: str, start: int, end: int) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ts = TimeSeries.read(str(gwf), channel, start=start, end=end)
    t_arr = np.asarray(ts.times.value, float)
    v_arr = np.asarray(ts.value, float)
    bins  = (t_arr // BIN_S).astype(int) * BIN_S
    return (
        pd.DataFrame({"gps_bin": bins, "hrec_range_bns": v_arr})
        .groupby("gps_bin", as_index=False)["hrec_range_bns"]
        .mean()
    )


day_start = (GPS_START // 86400) * 86400
day_end   = (GPS_END   // 86400) * 86400

records = []
for gps_day in range(day_start, day_end + 86400, 86400):
    gwf = stage(gps_day)
    if gwf is None:
        continue
    t0 = max(GPS_START, gps_day)
    t1 = min(GPS_END, gps_day + 86400)
    if t1 <= t0:
        gwf.unlink(missing_ok=True)
        continue
    try:
        df_ch = read_hourly(gwf, CHANNEL, t0, t1)
        records.append(df_ch)
        print(f"  day {gps_day} OK ({len(df_ch)} bins)", flush=True)
    except Exception as e:
        print(f"  {CHANNEL} @ {gps_day}: failed ({e})", flush=True)
    gwf.unlink(missing_ok=True)

if not records:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = (
    pd.concat(records, ignore_index=True)
    .drop_duplicates(subset=["gps_bin"])
    .sort_values("gps_bin")
)

out = WORKDIR / "outputs" / f"hrec_range_bns_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)

print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
