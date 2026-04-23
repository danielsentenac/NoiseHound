#!/usr/bin/env bash
#SBATCH -J nh-tower-thermal
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-tower-thermal-%j.out

# Extract thermal/TCS channels for all towers (WI, WE, NE, BS, PR, SR, IB, DET)
# for one O4 quarter. Mirrors the structure of nh_ni_thermal_extract.slurm.
#
# Pass GPS_START and GPS_END via --export.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_tower_thermal_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_tower_thermal_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_tower_thermal_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_tower_thermal_extract.slurm
#
# Output: outputs/tower_thermal_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_tower_thermal_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Tower thermal extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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
    # ── WI ──────────────────────────────────────────────────────────────────
    ("V1:INF_WI_BOTTOM_TE1",                    "wi_bottom_te1"),
    ("V1:INF_WI_MIR_COIL_DR_TE",               "wi_mir_coil_dr_te"),
    ("V1:INF_TCS_WI_RH_TE",                     "wi_rh_te"),
    ("V1:ENV_WI_F0_TE1",                         "wi_f0_te1"),
    ("V1:ENV_WI_F0_TE2",                         "wi_f0_te2"),
    ("V1:ENV_WI_F4_TE1",                         "wi_f4_te1"),
    ("V1:ENV_WI_F4_TE2",                         "wi_f4_te2"),
    ("V1:ENV_WI_F7_TE1",                         "wi_f7_te1"),
    ("V1:ENV_WI_F7_TE2",                         "wi_f7_te2"),
    ("V1:TCS_WI_CO2_PWRLAS_mean",               "wi_co2_pwrlas"),
    ("V1:TCS_WI_CO2_PWRIN_mean",                "wi_co2_pwrin"),
    ("V1:LSC_Etalon_WI_RH_OUT_mean",            "wi_rh_out"),
    ("V1:LSC_Etalon_WI_RH_SET_mean",            "wi_rh_set"),

    # ── WE ──────────────────────────────────────────────────────────────────
    ("V1:INF_WE_MIR_COIL_UL_TE",               "we_mir_coil_ul_te"),
    ("V1:INF_TCS_WE_RH_TE",                     "we_rh_te"),
    ("V1:ENV_WE_F0_TE1",                         "we_f0_te1"),
    ("V1:ENV_WE_F0_TE2",                         "we_f0_te2"),
    ("V1:ENV_WE_F4_TE1",                         "we_f4_te1"),
    ("V1:ENV_WE_F4_TE2",                         "we_f4_te2"),
    ("V1:ENV_WE_F7_TE1",                         "we_f7_te1"),
    ("V1:ENV_WE_F7_TE2",                         "we_f7_te2"),
    ("V1:TCS_WE_RH_SET",                         "we_rh_set"),

    # ── NE ──────────────────────────────────────────────────────────────────
    ("V1:INF_NE_MIR_COIL_UL_TE",               "ne_mir_coil_ul_te"),
    ("V1:INF_TCS_NE_RH_TE",                     "ne_rh_te"),
    ("V1:ENV_NE_F0_TE1",                         "ne_f0_te1"),
    ("V1:ENV_NE_F0_TE2",                         "ne_f0_te2"),
    ("V1:ENV_NE_F4_TE1",                         "ne_f4_te1"),
    ("V1:ENV_NE_F4_TE2",                         "ne_f4_te2"),
    ("V1:ENV_NE_F7_TE1",                         "ne_f7_te1"),
    ("V1:ENV_NE_F7_TE2",                         "ne_f7_te2"),
    ("V1:TCS_NE_RH_SET",                         "ne_rh_set"),

    # ── BS ──────────────────────────────────────────────────────────────────
    ("V1:ENV_BS_F0_TE1",                         "bs_f0_te1"),
    ("V1:ENV_BS_F0_TE2",                         "bs_f0_te2"),
    ("V1:ENV_BS_F4_TE1",                         "bs_f4_te1"),
    ("V1:ENV_BS_F4_TE2",                         "bs_f4_te2"),
    ("V1:ENV_BS_F7_TE2",                         "bs_f7_te2"),

    # ── PR ──────────────────────────────────────────────────────────────────
    ("V1:INF_PR_MIR_COIL_UL_TE",               "pr_mir_coil_ul_te"),
    ("V1:INF_TCS_PR_RH_TE",                     "pr_rh_te"),
    ("V1:ENV_PR_F0_TE1",                         "pr_f0_te1"),
    ("V1:ENV_PR_F0_TE2",                         "pr_f0_te2"),
    ("V1:ENV_PR_F4_TE1",                         "pr_f4_te1"),
    ("V1:ENV_PR_F4_TE2",                         "pr_f4_te2"),
    ("V1:ENV_PR_F7_TE1",                         "pr_f7_te1"),
    ("V1:ENV_PR_F7_TE2",                         "pr_f7_te2"),

    # ── SR ──────────────────────────────────────────────────────────────────
    ("V1:INF_SR_MIR_COIL_UL_TE",               "sr_mir_coil_ul_te"),
    ("V1:INF_TCS_SR_RH_TE",                     "sr_rh_te"),
    ("V1:ENV_SR_F0_TE1",                         "sr_f0_te1"),
    ("V1:ENV_SR_F0_TE2",                         "sr_f0_te2"),
    ("V1:ENV_SR_F4_TE1",                         "sr_f4_te1"),
    ("V1:ENV_SR_F4_TE2",                         "sr_f4_te2"),
    ("V1:ENV_SR_F7_TE1",                         "sr_f7_te1"),
    ("V1:ENV_SR_F7_TE2",                         "sr_f7_te2"),
    ("V1:TCS_SR_RH_SET",                         "sr_rh_set"),

    # ── IB ──────────────────────────────────────────────────────────────────
    ("V1:ENV_IB_BENCH_TE1",                      "ib_bench_te1"),
    ("V1:ENV_IB_BENCH_TE2",                      "ib_bench_te2"),
    ("V1:ENV_IB_BOX_TE1",                        "ib_box_te1"),
    ("V1:ENV_IB_BOX_TE2",                        "ib_box_te2"),
    ("V1:ENV_IB_F0_TE1",                         "ib_f0_te1"),
    ("V1:ENV_IB_F0_TE2",                         "ib_f0_te2"),
    ("V1:ENV_IB_F4_TE1",                         "ib_f4_te1"),
    ("V1:ENV_IB_F4_TE2",                         "ib_f4_te2"),
    ("V1:ENV_IB_F7_TE1",                         "ib_f7_te1"),
    ("V1:ENV_IB_F7_TE2",                         "ib_f7_te2"),

    # ── DET ─────────────────────────────────────────────────────────────────
    ("V1:HVAC_DET_TE_IN",                        "det_te_in"),
    ("V1:HVAC_DET_TE_OUT",                       "det_te_out"),
    ("V1:HVAC_DET_TE_SET",                       "det_te_set"),
    ("V1:HVAC_DET_DLAB_TE",                      "det_dlab_te"),
    ("V1:HVAC_DET_SAS_TE_IN",                    "det_sas_te_in"),
    ("V1:HVAC_DET_SAS_TE_OUT",                   "det_sas_te_out"),
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
out = WORKDIR / "outputs" / f"tower_thermal_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
