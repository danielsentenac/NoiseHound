#!/usr/bin/env bash
#SBATCH -J nh-sr-pwr
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-sr-pwr-%j.out

# Extract all SR power/electrical channels for the disappearance epoch Oct 2025 – Apr 2026.
#
# Channel families:
#   Sc_SR_MAR_PSDM_PWR_mean   — marionette mirror PSD total optical power
#   Sc_SR_MAR_PSDT_PWR_mean   — marionette top PSD total optical power
#   Sc_SR_MIR_PSDF_PWR_mean   — mirror front PSD total optical power
#   Sc_SR_MIR_PSDI_PWR_mean   — mirror inside PSD total optical power
#   Sc_SR_MIR_VOUT_{UL,UR,DL,DR}_rms  — mirror coil voltage outputs (4 corners)
#   TCS_SR_RH3_CURR / _VOLT   — ring heater 3 current & voltage
#   TCS_SR_RH4_CURR / _VOLT   — ring heater 4 current & voltage
#   TCS_SR_CHROCC_I_OUT_rms   — CO2 chrono-occupation output current
#   SAT_SR_mtr{A,B,C}_I / _P  — motor coil current & power (3 motors)
#   ASC_SR_TX_B1p_mag_10Hz_FS_mean  — SR TX angular oscillation on B1p
#   ASC_SR_TX_B4_112MHz_mag_10Hz_FS_mean — SR TX on B4 (112 MHz)
#
# Output: outputs/sr_pwr_epoch.csv
# Columns: gps_bin + one column per channel (snake_case key)

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_sr_pwr_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"
export STAGE_DIR

echo "=== SR power channel extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

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
GPS_START = 1443657600   # 2025-10-01
GPS_END   = 1459468800   # 2026-04-01
BIN_S     = 3600         # 1-h bins

# (channel_name, output_column_key)
CHANNELS = [
    # PSD optical power on local control sensors
    ("V1:Sc_SR_MAR_PSDM_PWR_mean",              "psdm_pwr"),
    ("V1:Sc_SR_MAR_PSDT_PWR_mean",              "psdt_pwr"),
    ("V1:Sc_SR_MIR_PSDF_PWR_mean",              "psdf_pwr"),
    ("V1:Sc_SR_MIR_PSDI_PWR_mean",              "psdi_pwr"),
    # Mirror coil voltage outputs (4 corners — angular state encoded here)
    ("V1:Sc_SR_MIR_VOUT_UL_rms",                "mir_vout_ul"),
    ("V1:Sc_SR_MIR_VOUT_UR_rms",                "mir_vout_ur"),
    ("V1:Sc_SR_MIR_VOUT_DL_rms",                "mir_vout_dl"),
    ("V1:Sc_SR_MIR_VOUT_DR_rms",                "mir_vout_dr"),
    # TCS ring heaters (thermal control of SR mirror)
    ("V1:TCS_SR_RH3_CURR",                       "tcs_rh3_curr"),
    ("V1:TCS_SR_RH3_VOLT",                       "tcs_rh3_volt"),
    ("V1:TCS_SR_RH4_CURR",                       "tcs_rh4_curr"),
    ("V1:TCS_SR_RH4_VOLT",                       "tcs_rh4_volt"),
    ("V1:TCS_SR_CHROCC_I_OUT_rms",               "tcs_chrocc_i_out"),
    # SAT box motor currents & power (3 motors driving the marionette)
    ("V1:SAT_SR_mtrA_I",                         "mtr_a_curr"),
    ("V1:SAT_SR_mtrA_P",                         "mtr_a_pwr"),
    ("V1:SAT_SR_mtrB_I",                         "mtr_b_curr"),
    ("V1:SAT_SR_mtrB_P",                         "mtr_b_pwr"),
    ("V1:SAT_SR_mtrC_I",                         "mtr_c_curr"),
    ("V1:SAT_SR_mtrC_P",                         "mtr_c_pwr"),
    # ASC angular oscillation amplitude (TX axis, already have TY from sr_pos)
    ("V1:ASC_SR_TX_B1p_mag_10Hz_FS_mean",        "asc_sr_tx_b1p_mag"),
    ("V1:ASC_SR_TX_B4_112MHz_mag_10Hz_FS_mean",  "asc_sr_tx_b4_mag"),
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
out = WORKDIR / "outputs" / "sr_pwr_epoch.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins → {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
