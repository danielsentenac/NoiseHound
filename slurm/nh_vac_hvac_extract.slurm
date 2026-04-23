#!/usr/bin/env bash
#SBATCH -J nh-vac-hvac
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-vac-hvac-%j.out

# Extract VAC_VALVEBIGNI_CRYO_TEMP and HVAC_INJ_TE_SATUR_CORR_COLD for full O4.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_vac_hvac_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_vac_hvac_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_vac_hvac_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_vac_hvac_extract.slurm
#
# Output: outputs/vac_hvac_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_vac_hvac_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== VAC/HVAC extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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

CHANNELS = [
    ("V1:VAC_VALVEBIGNI_CRYO_TEMP",       "vac_ni_cryo_temp"),
    ("V1:HVAC_INJ_TE_SATUR_CORR_COLD",    "hvac_inj_te_cold"),
]

def stage(gps_day):
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

# Build 1-h bin grid
bins = list(range(GPS_START, GPS_END, BIN_S))
rows = []

day_start = (GPS_START // 86400) * 86400
day_end   = ((GPS_END + 86399) // 86400) * 86400

for gps_day in range(day_start, day_end, 86400):
    gwf = stage(gps_day)
    if gwf is None:
        continue
    day_bins = [b for b in bins if gps_day <= b < gps_day + 86400]
    if not day_bins:
        gwf.unlink(missing_ok=True)
        continue
    for col, key in CHANNELS:
        try:
            ts = TimeSeries.read(str(gwf), col,
                                 start=gps_day, end=gps_day + 86400,
                                 nproc=1)
            arr = np.asarray(ts.value, float)
            t0  = float(ts.t0.value)
            dt  = float(ts.dt.value)
            for b in day_bins:
                i0 = max(0, int((b - t0) / dt))
                i1 = min(len(arr), int((b + BIN_S - t0) / dt))
                v  = np.nanmean(arr[i0:i1]) if i1 > i0 else np.nan
                rows.append({"gps_bin": b, "channel": key, "value": v})
        except Exception as e:
            for b in day_bins:
                rows.append({"gps_bin": b, "channel": key, "value": np.nan})
    gwf.unlink(missing_ok=True)
    print(f"  day {gps_day} done", flush=True)

if not rows:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = pd.DataFrame(rows)
pivot = df.pivot_table(index="gps_bin", columns="channel", values="value", aggfunc="mean")
pivot.reset_index(inplace=True)
pivot.columns.name = None

out = WORKDIR / "outputs" / f"vac_hvac_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
pivot.to_csv(out, index=False)
print(f"Saved {len(pivot)} rows → {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
