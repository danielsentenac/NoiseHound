#!/usr/bin/env bash
#SBATCH -J nh-snr-corr-ch2
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-snr-corr-ch2-%j.out

# Extract second batch of top SNR-correlated physical channels for full O4.
# Selected from snr_channel_correlations.csv — new candidates:
#   - CAL permline amplitudes (direct sensitivity proxies, rho~-0.84 to -0.90)
#   - SQB green PSD channels (squeezing-related, rho~0.93-0.96)
#   - SPRB humidity + pressure (rho~-0.89 to -0.92)
#   - BsX QN (rho=+0.85, lag=-56d)
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_snr_corr_channels2_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_snr_corr_channels2_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_snr_corr_channels2_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_snr_corr_channels2_extract.slurm
#
# Output: outputs/snr_corr_channels2_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_snr_corr_ch2_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== SNR corr channels2 extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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
    # ── CAL permanent calibration line amplitudes ─────────────────────────
    # Inversely proportional to sensitivity → should anti-correlate with SNR
    # (rho ~ -0.84 to -0.90)
    ("V1:CAL_NE_MIR_permline1_ampl",         "cal_ne_mir_permline1"),
    ("V1:CAL_WE_MIR_permline1_ampl",         "cal_we_mir_permline1"),
    ("V1:CAL_BS_MIR_permline1_ampl",         "cal_bs_mir_permline1"),
    ("V1:CAL_SR_MIR_permline1_ampl",         "cal_sr_mir_permline1"),
    ("V1:CAL_WE_PCAL_permline1_ampl",        "cal_we_pcal_permline1"),
    ("V1:CAL_NE_PCAL_permline1_ampl",        "cal_ne_pcal_permline1"),
    ("V1:CAL_WE_PCAL_permline0_ampl",        "cal_we_pcal_permline0"),

    # ── SQB green reference photodetectors (rho ~ +0.93-0.96) ────────────
    # Squeezing-related channels; large lag may indicate infrastructure changes
    ("V1:SQB2_GR_PSD2_X1_DC_mean",          "sqb2_gr_psd2_x1_dc"),
    ("V1:SQB1_GR_PSD1_SUM_mean",            "sqb1_gr_psd1_sum"),
    ("V1:SQB2_GR_PSD2_SUM_mean",            "sqb2_gr_psd2_sum"),
    ("V1:SQB1_GR_PSD1_X1_DC_mean",          "sqb1_gr_psd1_x1_dc"),

    # ── SPRB environment (rho ~ -0.89 to -0.92) ──────────────────────────
    ("V1:SPRB_DBOX_LeftDown_ext_humidity",   "sprb_leftdown_humidity"),
    ("V1:SPRB_DBOX_LeftUp_ext_humidity",     "sprb_leftup_humidity"),
    ("V1:SPRB_B4_QD2_Pressure",              "sprb_b4_pressure"),
    ("V1:SPRB_DBOX_LeftDown_mezz1_temp0",    "sprb_leftdown_temp"),

    # ── BsX homodyne other (rho=+0.85, lag=-56d) ─────────────────────────
    ("V1:BsX_QN_v_mean",                     "bsx_qn_v"),

    # ── SPRB LC LVDT (rho ~ -0.93, lag=-14d) — SR bench suspension ───────
    ("V1:SPRB_LC_LVDT_BL_V_out_q_mean",     "sprb_lvdt_bl_v_q"),
    ("V1:SPRB_LC_LVDT_FL_V_out_q_mean",     "sprb_lvdt_fl_v_q"),
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
out = WORKDIR / "outputs" / f"snr_corr_channels2_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
