#!/usr/bin/env bash
#SBATCH -J nh-logbook-ch
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-logbook-ch-%j.out

# Extract logbook-referenced channels and LSC DARM calibration lines for full O4.
# Channels sourced from:
#   - Virgo logbooks #60143, #62992, #67414, #68210 (NI/WI thermal correlators)
#   - EXCAVATor per-event ranking (SNEB geophones, CEB electrical)
#   - Manual selection: LSC DARM IMC/PSTAB lines (sensitivity sanity check)
#
# Channels already covered by previous jobs (ni_thermal, tower_thermal) are excluded.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_logbook_channels_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_logbook_channels_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_logbook_channels_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_logbook_channels_extract.slurm
#
# Output: outputs/logbook_channels_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_logbook_ch_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Logbook channels extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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
    # ── HWS bench thermistors ─────────────────────────────────────────────
    ("V1:TCS_HWS_NI_TE1_mean",               "tcs_hws_ni_te1"),
    ("V1:TCS_HWS_NI_TE2_mean",               "tcs_hws_ni_te2"),
    ("V1:TCS_HWS_WI_TE1_mean",               "tcs_hws_wi_te1"),
    ("V1:TCS_HWS_WI_TE2_mean",               "tcs_hws_wi_te2"),
    ("V1:TCS_HWS_NE_TE1_mean",               "tcs_hws_ne_te1"),   # control (remote)

    # ── CO2 bench / laser temperatures ───────────────────────────────────
    ("V1:ENV_TCS_CO2_NI_TE",                 "env_co2_ni_te"),
    ("V1:TCS_NI_TE_CO2Laser",                "tcs_ni_co2laser_te"),
    ("V1:ENV_TCS_CO2_WI_TE",                 "env_co2_wi_te"),
    ("V1:TCS_WI_TE_CO2Laser",                "tcs_wi_co2laser_te"),

    # ── CEB ambient temperature ───────────────────────────────────────────
    # Logbook #60143, #66142: partial glitch-rate correlator
    ("V1:ENV_CEB_N_TE",                      "env_ceb_n_te"),

    # ── Ring heater ERR / IN (set+out already extracted) ─────────────────
    ("V1:LSC_Etalon_NI_RH_IN_mean",          "rh_ni_in"),
    ("V1:LSC_Etalon_NI_RH_ERR_mean",         "rh_ni_err"),
    ("V1:LSC_Etalon_WI_RH_ERR_mean",         "rh_wi_err"),

    # ── Mains electrical — UPS voltage & current ─────────────────────────
    # Logbook suspicion: mains coupling to glitch rate
    ("V1:ENV_NEB_UPS_VOLT_R_mean",           "env_neb_ups_volt_r"),
    ("V1:ENV_CEB_UPS_VOLT_R_mean",           "env_ceb_ups_volt_r"),
    ("V1:ENV_WEB_UPS_VOLT_R_mean",           "env_web_ups_volt_r"),
    ("V1:ENV_MCB_IPS_CURR_T_mean",           "env_mcb_ips_curr_t"),
    # EXCAVATor per-event rank 15: CEB electrical coupling
    ("V1:ENV_CEB_UPS_CURR_R_mean",           "env_ceb_ups_curr_r"),

    # ── SNEB / SWEB suspended bench geophones ────────────────────────────
    # EXCAVATor per-event rank 5/20: seismic at SNEB correlates with glitches
    ("V1:SBE_SNEB_GEO_GRWE_raw_mean",        "sbe_sneb_geo_we"),
    ("V1:SBE_SNEB_GEO_GRNS_raw_mean",        "sbe_sneb_geo_ns"),
    ("V1:SBE_SWEB_GEO_GRWE_raw_mean",        "sbe_sweb_geo_we"),   # control
    ("V1:SBE_SWEB_GEO_GRNS_raw_mean",        "sbe_sweb_geo_ns"),   # control

    # ── LSC DARM calibration / noise lines ───────────────────────────────
    # IMC line: inversely proportional to sensitivity (sanity check)
    ("V1:LSC_DARM_IMC_LINE_mag_100Hz_mean",  "lsc_darm_imc_mag"),
    ("V1:LSC_DARM_IMC_LINE_Q_100Hz_mean",    "lsc_darm_imc_q"),
    # PSTAB0 coupling: power stabilisation noise coupling to DARM
    ("V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean", "lsc_darm_pstab0_coupling"),
    ("V1:LSC_DARM_PSTAB0_I_FS_mean",         "lsc_darm_pstab0_i_fs"),
    ("V1:LSC_DARM_PSTAB0_I_100Hz_mean",      "lsc_darm_pstab0_i_100hz"),
    # SIB1 line: squeezing injection bench calibration line
    ("V1:LSC_DARM_SIB1_LINE_mag_100Hz_mean", "lsc_darm_sib1_mag"),
    # WI feed-forward 50 Hz error and phase
    ("V1:Sc_WI_FF50HZ_P_ERR_mean",           "sc_wi_ff50hz_p_err"),
    ("V1:Sc_WI_FF50HZ_PHASE_mean",           "sc_wi_ff50hz_phase"),

    # ── SQB1 camera — squeezing beam waist ───────────────────────────────
    ("V1:SQB1_Cam1_FitWaistY_mean",          "sqb1_cam1_waist_y"),
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
out = WORKDIR / "outputs" / f"logbook_channels_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
