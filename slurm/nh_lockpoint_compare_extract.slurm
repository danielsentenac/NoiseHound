#!/usr/bin/env bash
#SBATCH -J nh-lockpoint
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-lockpoint-%j.out

# Extract a narrow set of lock-point observables used for the SR Christmas 2025
# comparison (Oct 2025 / Dec 2025 / Jan 2026).
#
# Intended observables:
#   - B1p DC level
#   - B4 112 MHz sideband level
#   - MICH setpoint
#   - SRCL setpoint
#
# Pass GPS_START / GPS_END with --export, for example:
#   sbatch --export=GPS_START=1444435200,GPS_END=1444867200 slurm/nh_lockpoint_compare_extract.slurm
#   sbatch --export=GPS_START=1449705600,GPS_END=1450051200 slurm/nh_lockpoint_compare_extract.slurm
#   sbatch --export=GPS_START=1452038400,GPS_END=1452297600 slurm/nh_lockpoint_compare_extract.slurm
#   sbatch --export=GPS_START=1452556800,GPS_END=1454112000 slurm/nh_lockpoint_compare_extract.slurm
#
# Output: outputs/lockpoint_compare_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1444435200}"
GPS_END="${GPS_END:-1444867200}"
export STAGE_DIR="/tmp/nh_lockpoint_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Lock-point comparison extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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

CHANNEL_GROUPS = [
    ("b1p_dc_trig_p", [
        "V1:LSC_B1p_DC_TRIG_p_mean",
        "V1:LSC_B1p_DC_TRIG_p",
    ]),
    ("b1p_dc", [
        "V1:LSC_B1p_DC_mean",
        "V1:LSC_B1p_DC",
    ]),
    ("det_b1p_dc", [
        "V1:DET_B1p_DC_mean",
        "V1:DET_B1p_DC",
    ]),
    ("b4_112mhz_mag", [
        "V1:LSC_B4_112MHz_MAG_mean",
        "V1:LSC_B4_112MHz_MAG",
    ]),
    ("b4_112mhz_mag_norm", [
        "V1:LSC_B4_112MHz_MAG_NORM_mean",
        "V1:LSC_B4_112MHz_MAG_NORM",
    ]),
    ("sprb_b4_112mhz_mag", [
        "V1:SPRB_B4_112MHz_mag_mean",
        "V1:SPRB_B4_112MHz_mag",
        "V1:SPRB_B4_PD1_112MHz_mag_mean",
    ]),
    ("sprb_b4_qd2_ul_112mhz", [
        "V1:SPRB_B4_QD2_UL_112MHz_mag_FS_mean",
    ]),
    ("mich_set", [
        "V1:LSC_MICH_SET_mean",
        "V1:LSC_MICH_SET",
    ]),
    ("mich_set_tot", [
        "V1:LSC_MICH_SET_TOT_mean",
        "V1:LSC_MICH_SET_TOT",
    ]),
    ("srcl_set", [
        "V1:LSC_SRCL_SET_mean",
        "V1:LSC_SRCL_SET",
    ]),
    ("srcl_set_tot", [
        "V1:LSC_SRCL_SET_TOT_mean",
        "V1:LSC_SRCL_SET_TOT",
    ]),
]


def stage(gps_day: int) -> Path | None:
    year   = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
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
    return (pd.DataFrame({"gps_bin": bins, "val": v_arr})
              .groupby("gps_bin", as_index=False)["val"].mean())


day_start = (GPS_START // 86400) * 86400
day_end   = (GPS_END   // 86400) * 86400

records = []
selected = {}

for gps_day in range(day_start, day_end + 86400, 86400):
    gwf = stage(gps_day)
    if gwf is None:
        continue
    t0 = max(GPS_START, gps_day)
    t1 = min(GPS_END, gps_day + 86400)
    any_ok = False
    for key, candidates in CHANNEL_GROUPS:
        last_err = None
        for ch in candidates:
            try:
                df_ch = read_hourly(gwf, ch, t0, t1)
                for gps_bin, val in df_ch.itertuples(index=False):
                    records.append({"gps_bin": int(gps_bin), key: val})
                selected.setdefault(key, ch)
                any_ok = True
                break
            except Exception as e:
                last_err = e
        else:
            print(f"  {key} @ {gps_day}: all candidates failed ({last_err})", flush=True)
    if any_ok:
        print(f"  day {gps_day} OK", flush=True)
    gwf.unlink(missing_ok=True)

if not records:
    print("No data extracted!", flush=True)
    sys.exit(1)

df = (pd.DataFrame(records)
        .groupby("gps_bin", as_index=False)
        .mean()
        .sort_values("gps_bin"))

out = WORKDIR / "outputs" / f"lockpoint_compare_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)

print(f"Saved {len(df)} bins -> {out}", flush=True)
print("Selected channels:", flush=True)
for key, ch in selected.items():
    print(f"  {key}: {ch}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
