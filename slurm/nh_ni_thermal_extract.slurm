#!/usr/bin/env bash
#SBATCH -J nh-ni-thermal
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-ni-thermal-%j.out

# Extract NI tower thermal + angular motion channels for one O4 quarter.
# Relevant for investigating the NI thermal oscillation → SR glitch coupling.
#
# Pass GPS_START and GPS_END via --export.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_ni_thermal_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_ni_thermal_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_ni_thermal_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_ni_thermal_extract.slurm
#
# Output: outputs/ni_thermal_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_ni_thermal_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== NI thermal extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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

# (channel_name, output_column_key)
CHANNELS = [
    # Bottom tower temperature — primary NI thermal oscillation proxy
    ("V1:INF_NI_BOTTOM_TE1",                    "ni_bottom_te1"),
    # Mirror coil & RH temperatures
    ("V1:INF_NI_MIR_COIL_UL_TE",               "ni_mir_coil_ul_te"),
    ("V1:INF_TCS_NI_RH_TE",                     "ni_rh_te"),
    # Environmental temperatures at 3 tower floors
    ("V1:ENV_NI_F0_TE1",                        "ni_f0_te1"),
    ("V1:ENV_NI_F0_TE2",                        "ni_f0_te2"),
    ("V1:ENV_NI_F4_TE1",                        "ni_f4_te1"),
    ("V1:ENV_NI_F4_TE2",                        "ni_f4_te2"),
    ("V1:ENV_NI_F7_TE1",                        "ni_f7_te1"),
    ("V1:ENV_NI_F7_TE2",                        "ni_f7_te2"),
    # NI laser body temperatures near the etalon thermal actuators
    ("V1:TCS_NI_TE_CO2Laser",                   "tcs_ni_co2laser_te"),
    ("V1:TCS_NI_TE_AUXLaser",                   "tcs_ni_auxlaser_te"),
    # CO2 laser power (thermal actuator on NI mirror)
    ("V1:TCS_NI_CO2_PWRLAS_mean",               "ni_co2_pwrlas"),
    ("V1:TCS_NI_CO2_PWRIN_mean",                "ni_co2_pwrin"),
    ("V1:TCS_NI_CO2_PWROUT_mean",               "ni_co2_pwrout"),
    # Ring heater (etalon control)
    ("V1:LSC_Etalon_NI_RH_OUT_mean",            "ni_rh_out"),
    ("V1:LSC_Etalon_NI_RH_SET_mean",            "ni_rh_set"),
    # NI angular motion amplitude (TX axis, 10 Hz FS)
    ("V1:ASC_NI_TX_LF_B7_mag_10Hz_FS_mean",    "ni_tx_mag"),
]


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
    for ch, key in CHANNELS:
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
out = WORKDIR / "outputs" / f"ni_thermal_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
