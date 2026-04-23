#!/usr/bin/env python3
"""
Plot top SNR-correlated physical channels (full O4 timeseries).

Channels grouped into logical panels, each with SNR trend overlaid on
a twin y-axis.  Correlation coefficients from snr_channel_correlations.csv
annotated on each panel.
"""

import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator

# ── paths ─────────────────────────────────────────────────────────────────────
WORKDIR    = __import__("pathlib").Path(__file__).parent.parent
DATA_FILE   = WORKDIR / "outputs" / "snr_corr_channels_full_o4.csv"
SNR_FILE    = WORKDIR / "outputs" / "snr_trend_daily.csv"
GLITCH_FILE = WORKDIR / "data"    / "full_25min_glitches_ER16-O4b.csv"
LOCK_FILE   = WORKDIR / "outputs" / "itf_lock_full_o4.csv"
CORR_FILE   = WORKDIR / "outputs" / "snr_channel_correlations.csv"
OUT_FILE    = WORKDIR / "usecases" / "25-minute-glitch" / "snr_corr_channels_timeseries.png"

GPS_EPOCH  = pd.Timestamp("1980-01-06", tz="UTC")

# ── load data ─────────────────────────────────────────────────────────────────
df     = pd.read_csv(DATA_FILE)
lock   = pd.read_csv(LOCK_FILE)
glitch = pd.read_csv(GLITCH_FILE)

df["time"]     = GPS_EPOCH + pd.to_timedelta(df["gps_bin"],   unit="s")
lock["time"]   = GPS_EPOCH + pd.to_timedelta(lock["gps_bin"], unit="s")
glitch["time"] = GPS_EPOCH + pd.to_timedelta(glitch["time"],  unit="s")

# load correlation table for annotation
corr = {}
try:
    cdf = pd.read_csv(CORR_FILE).set_index("channel")
    corr = cdf["rho"].to_dict()
except Exception:
    pass

# ── lock state spans ──────────────────────────────────────────────────────────
STATE_META = {
    134: ("LN3",         0), 135: ("LN3",         0),
    144: ("LN3_ALIGNED", 1), 145: ("LN3_ALIGNED", 1),
    149: ("LN3_SQZ",     2), 150: ("LN3_SQZ",     2),
    1:   ("SCIENCE",     3), 2:   ("SCIENCE",     3),
}
STATE_IDS = sorted({k: v for k, v in STATE_META.items() if v[1] < 4}, key=lambda k: STATE_META[k][1])
_cmap = plt.cm.turbo(np.linspace(0.08, 0.95, 8))

def add_lock_spans(ax):
    prev = None
    for _, row in lock.iterrows():
        sid = int(row["state_id"]) if "state_id" in lock.columns else None
        if sid is None:
            continue
        t = row["time"]
        if prev is not None:
            ps, pt = prev
            if ps in STATE_META:
                idx = STATE_META[ps][1]
                col = _cmap[min(idx, 7)]
                ax.axvspan(pt, t, alpha=0.08, color=col, linewidth=0)
        prev = (sid, t)

# ── panel definitions ─────────────────────────────────────────────────────────
# Each entry: (panel_title, [(col_key, label, channel_name_for_rho), ...])
PANELS = [
    ("BsX homodyne PZT  [ρ ≈ −0.87 with SNR trend]", [
        ("bsx_pzt_dhc", "BsX_PZT_DHC",  "V1:BsX_PZT_DHC_mean"),
        ("bsx_pzt_dh",  "BsX_PZT_DH",   "V1:BsX_PZT_DH_mean"),
        ("bsx_x_corr",  "BsX_X_CORR",   "V1:BsX_X_CORR_mean"),
    ]),
    ("SPRB B4 QD2 photodetector  [ρ ≈ +0.77–0.79]", [
        ("sprb_b4_sum",       "SPRB B4 Sum",        "V1:SPRB_B4_QD2_Sum_mean"),
        ("sprb_b4_tx_err",    "SPRB B4 TX err",     "V1:SPRB_LC_B4_QD2_TX_err_mean"),
        ("sprb_b4_ul_112mhz", "SPRB B4 UL 112 MHz", "V1:SPRB_B4_QD2_UL_112MHz_mag_FS_mean"),
        ("sprb_b4_ul_12mhz",  "SPRB B4 UL 12 MHz",  "V1:SPRB_B4_QD2_UL_12MHz_mag_FS_mean"),
    ]),
    ("SPRB LC_Y / CEB temp / SR F0 acc", [
        ("sprb_lc_y",       "SPRB LC_Y",       "V1:SPRB_LC_Y_mean"),
        ("sprb_lc_y_err",   "SPRB LC_Y err",   "V1:SPRB_LC_Y_err_mean"),
        ("ceb_platform_temp","CEB platform T°", "V1:CEB_DBOX_PLATFORM_mezz0_temp0"),
        ("sa_sr_f0_acc_h2", "SR F0 ACC H2",    "V1:Sa_SR_F0_ACC_H2_FB_mean"),
    ]),
    ("VAC air pressure — tube atmosphere  [ρ ≈ +0.80, seasonal]", [
        ("vac_tube3000n_airpres", "Tube 3000N air",  "V1:VAC_TUBELAL3000N_AIRPRESSURE1"),
        ("vac_tube3000w_airpres", "Tube 3000W air",  "V1:VAC_TUBELAL3000W_AIRPRESSURE1"),
        ("vac_tube1200n_airpres", "Tube 1200N air",  "V1:VAC_TUBELAL1200N_AIRPRESSURE1"),
        ("vac_tube600w_airpres",  "Tube 600W air",   "V1:VAC_TUBELAL600W_AIRPRESSURE1"),
    ]),
    ("VAC vacuum pressure — P21HR  [ρ ≈ +0.77–0.80]", [
        ("vac_link_pr_p21",      "Link PR P21HR",      "V1:VAC_LINK_PR_P21HR"),
        ("vac_tube3000w_p21",    "Tube 3000W P21HR",   "V1:VAC_TUBELAL3000W_P21HR"),
        ("vac_tube600w_p21",     "Tube 600W P21HR",    "V1:VAC_TUBELAL600W_P21HR"),
        ("vac_cryolink_det_p21", "CryoLink DET P21HR", "V1:VAC_CRYOLINK_DET_P21HR"),
        ("vac_cryolink_ib_p21",  "CryoLink IB P21HR",  "V1:VAC_CRYOLINK_IB_P21HR"),
    ]),
]

N_PANELS = len(PANELS)

# ── figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(N_PANELS, 1, figsize=(18, 4 * N_PANELS),
                         sharex=True)
fig.subplots_adjust(hspace=0.08, top=0.96, bottom=0.06, left=0.08, right=0.90)

colors6 = plt.cm.tab10(np.linspace(0, 0.9, 10))

for ax, (title, chan_list) in zip(axes, PANELS):
    add_lock_spans(ax)

    # plot channels (normalise each to [0,1] range for comparability)
    for i, (col, lbl, vch) in enumerate(chan_list):
        if col not in df.columns:
            continue
        y = df[col].values.copy().astype(float)
        y[~np.isfinite(y)] = np.nan
        # clip extreme outliers (>10 MAD from median) before normalising
        med = np.nanmedian(y)
        mad = np.nanmedian(np.abs(y - med))
        if mad > 0:
            y = np.where(np.abs(y - med) > 10 * mad, np.nan, y)
            med = np.nanmedian(y)
            mad = np.nanmedian(np.abs(y - med))
            yn = (y - med) / (6 * mad)   # ≈ ±1 range for most data
        else:
            yn = y - med
        # 24h rolling median smooth for hourly data (reduces lock/unlock spikes)
        yn_s = pd.Series(yn).rolling(24, center=True, min_periods=3).median().values
        rho = corr.get(vch, np.nan)
        rho_str = f"ρ={rho:+.3f}" if np.isfinite(rho) else ""
        label = f"{vch}  {rho_str}"
        ax.plot(df["time"], yn_s, lw=0.8, alpha=0.80, color=colors6[i], label=label)

    ax.set_ylabel("norm. value\n(robust σ units)", fontsize=8)
    ax.set_title(title, fontsize=9, loc="left", pad=3)
    ax.tick_params(labelsize=8)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(loc="upper left", fontsize=7, ncol=2, framealpha=0.6)

    # Glitch amplitude scatter on twin axis
    ax2 = ax.twinx()
    ax2.scatter(glitch["time"], glitch["amplitude"], s=1, color="black", alpha=0.25,
                linewidths=0, zorder=5, label="glitch amplitude")
    ax2.set_yscale("log")
    ax2.set_ylabel("glitch amplitude (strain)", fontsize=8, color="black")
    ax2.tick_params(labelsize=8, colors="black")
    ax2.spines["right"].set_color("black")

axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
axes[-1].tick_params(axis="x", labelsize=8, rotation=30)

fig.suptitle("Top SNR-correlated physical channels — full O4 vs glitch amplitude (strain)",
             fontsize=11, fontweight="bold", y=0.99)

OUT_FILE.parent.mkdir(exist_ok=True)
fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
print(f"Saved → {OUT_FILE}")
