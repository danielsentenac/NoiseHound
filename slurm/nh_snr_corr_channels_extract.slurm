#!/usr/bin/env bash
#SBATCH -J nh-snr-corr-ch
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-snr-corr-ch-%j.out

# Extract the top SNR-correlated physical channels for full O4 timeseries.
# Selected from snr_channel_correlations.csv — lag=0 physical channels
# plus the key VAC air/vacuum pressure channels.
#
# Submit 4 jobs covering full O4:
#   sbatch --export=GPS_START=1364774400,GPS_END=1388534400 slurm/nh_snr_corr_channels_extract.slurm
#   sbatch --export=GPS_START=1388534400,GPS_END=1406937600 slurm/nh_snr_corr_channels_extract.slurm
#   sbatch --export=GPS_START=1406937600,GPS_END=1430524800 slurm/nh_snr_corr_channels_extract.slurm
#   sbatch --export=GPS_START=1430524800,GPS_END=1459468800 slurm/nh_snr_corr_channels_extract.slurm
#
# Output: outputs/snr_corr_channels_${GPS_START}_${GPS_END}.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:-1364774400}"
GPS_END="${GPS_END:-1388534400}"
export STAGE_DIR="/tmp/nh_snr_corr_ch_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== SNR corr channels extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
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
    # ── BS homodyne PZT (rho=-0.87, lag=0) ──────────────────────────────
    ("V1:BsX_PZT_DHC_mean",                        "bsx_pzt_dhc"),
    ("V1:BsX_PZT_DH_mean",                         "bsx_pzt_dh"),
    ("V1:BsX_X_CORR_mean",                         "bsx_x_corr"),

    # ── SPRB B4 photodetector signals (rho~0.77-0.79, lag=0) ────────────
    ("V1:SPRB_LC_B4_QD2_TX_err_mean",              "sprb_b4_tx_err"),
    ("V1:SPRB_B4_QD2_UL_112MHz_mag_FS_mean",       "sprb_b4_ul_112mhz"),
    ("V1:SPRB_B4_QD2_UL_12MHz_mag_FS_mean",        "sprb_b4_ul_12mhz"),
    ("V1:SPRB_B4_QD2_Sum_mean",                    "sprb_b4_sum"),

    # ── SPRB longitudinal/angular control (rho=-0.75, lag=0) ────────────
    ("V1:SPRB_LC_Y_mean",                          "sprb_lc_y"),
    ("V1:SPRB_LC_Y_err_mean",                      "sprb_lc_y_err"),

    # ── CEB electronics platform temperature (rho=+0.77, lag=0) ─────────
    ("V1:CEB_DBOX_PLATFORM_mezz0_temp0",           "ceb_platform_temp"),

    # ── SR F0 accelerometer (rho=+0.64, lag=0) ──────────────────────────
    ("V1:Sa_SR_F0_ACC_H2_FB_mean",                 "sa_sr_f0_acc_h2"),

    # ── VAC air pressure — tube atmosphere (rho~0.80, lag=0) ────────────
    ("V1:VAC_TUBELAL3000N_AIRPRESSURE1",           "vac_tube3000n_airpres"),
    ("V1:VAC_TUBELAL3000W_AIRPRESSURE1",           "vac_tube3000w_airpres"),
    ("V1:VAC_TUBELAL1200N_AIRPRESSURE1",           "vac_tube1200n_airpres"),
    ("V1:VAC_TUBELAL600W_AIRPRESSURE1",            "vac_tube600w_airpres"),

    # ── VAC vacuum pressure — tube & link (rho~0.77-0.80, lag~0) ────────
    ("V1:VAC_LINK_PR_P21HR",                       "vac_link_pr_p21"),
    ("V1:VAC_TUBELAL3000W_P21HR",                  "vac_tube3000w_p21"),
    ("V1:VAC_TUBELAL600W_P21HR",                   "vac_tube600w_p21"),
    ("V1:VAC_CRYOLINK_DET_P21HR",                  "vac_cryolink_det_p21"),
    ("V1:VAC_CRYOLINK_IB_P21HR",                   "vac_cryolink_ib_p21"),

    # ── VAC tiltmeter PR (rho=-0.76, lag=+14d) ──────────────────────────
    ("V1:VAC_TILTMETER_PR1",                       "vac_tiltmeter_pr1"),
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
out = WORKDIR / "outputs" / f"snr_corr_channels_{GPS_START}_{GPS_END}.csv"
out.parent.mkdir(exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {len(df)} bins -> {out}", flush=True)
print(f"Columns: {list(df.columns)}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
