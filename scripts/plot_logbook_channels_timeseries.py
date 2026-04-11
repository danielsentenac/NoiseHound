#!/usr/bin/env python3
"""
Plot logbook-referenced channels vs SNR glitch scatter — full O4.

Panels grouped by physical theme; each panel has the individual glitch SNR
scatter on a twin right axis.  Channels not present in the merged CSV are
silently skipped.
"""

import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
WORKDIR     = Path(__file__).parent.parent
DATA_FILE   = WORKDIR / "outputs" / "logbook_channels_full_o4.csv"
GLITCH_FILE = WORKDIR / "data"    / "full_25min_glitches_ER16-O4b.csv"
LOCK_FILE   = WORKDIR / "outputs" / "itf_lock_full_o4.csv"
OUT_FILE    = WORKDIR / "usecases" / "25-minute-glitch" / "logbook_channels_timeseries.png"

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

# ── load ──────────────────────────────────────────────────────────────────────
df     = pd.read_csv(DATA_FILE)
lock   = pd.read_csv(LOCK_FILE)
glitch = pd.read_csv(GLITCH_FILE)

df["time"]     = GPS_EPOCH + pd.to_timedelta(df["gps_bin"],   unit="s")
lock["time"]   = GPS_EPOCH + pd.to_timedelta(lock["gps_bin"], unit="s")
glitch["time"] = GPS_EPOCH + pd.to_timedelta(glitch["time"],  unit="s")

# ── lock state spans ──────────────────────────────────────────────────────────
STATE_META = {
    134: 0, 135: 0,   # LN3
    144: 1, 145: 1,   # LN3_ALIGNED
    149: 2, 150: 2,   # LN3_SQZ
    1:   3, 2:   3,   # SCIENCE
}
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
                ax.axvspan(pt, t, alpha=0.08,
                           color=_cmap[min(STATE_META[ps], 7)], linewidth=0)
        prev = (sid, t)

# ── panel definitions ─────────────────────────────────────────────────────────
# (title, [(col_key, full_channel_name), ...])
PANELS = [
    ("LSC DARM calibration lines — sensitivity sanity check", [
        ("lsc_darm_imc_mag",        "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean"),
        ("lsc_darm_imc_q",          "V1:LSC_DARM_IMC_LINE_Q_100Hz_mean"),
        ("lsc_darm_sib1_mag",       "V1:LSC_DARM_SIB1_LINE_mag_100Hz_mean"),
        ("lsc_darm_pstab0_coupling","V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean"),
        ("lsc_darm_pstab0_i_fs",    "V1:LSC_DARM_PSTAB0_I_FS_mean"),
        ("lsc_darm_pstab0_i_100hz", "V1:LSC_DARM_PSTAB0_I_100Hz_mean"),
        ("sc_wi_ff50hz_p_err",      "V1:Sc_WI_FF50HZ_P_ERR_mean"),
        ("sc_wi_ff50hz_phase",      "V1:Sc_WI_FF50HZ_PHASE_mean"),
        ("sqb1_cam1_waist_y",       "V1:SQB1_Cam1_FitWaistY_mean"),
    ]),
    ("HWS bench thermistors — NI / WI / NE (control)", [
        ("tcs_hws_ni_te1", "V1:TCS_HWS_NI_TE1_mean"),
        ("tcs_hws_ni_te2", "V1:TCS_HWS_NI_TE2_mean"),
        ("tcs_hws_wi_te1", "V1:TCS_HWS_WI_TE1_mean"),
        ("tcs_hws_wi_te2", "V1:TCS_HWS_WI_TE2_mean"),
        ("tcs_hws_ne_te1", "V1:TCS_HWS_NE_TE1_mean"),
    ]),
    ("CO2 bench / laser temperatures — NI & WI", [
        ("env_co2_ni_te",      "V1:ENV_TCS_CO2_NI_TE"),
        ("tcs_ni_co2laser_te", "V1:TCS_NI_TE_CO2Laser"),
        ("env_co2_wi_te",      "V1:ENV_TCS_CO2_WI_TE"),
        ("tcs_wi_co2laser_te", "V1:TCS_WI_TE_CO2Laser"),
        ("env_ceb_n_te",       "V1:ENV_CEB_N_TE"),
    ]),
    ("NI / WI ring heater control — IN / ERR signals", [
        ("rh_ni_in",  "V1:LSC_Etalon_NI_RH_IN_mean"),
        ("rh_ni_err", "V1:LSC_Etalon_NI_RH_ERR_mean"),
        ("rh_wi_err", "V1:LSC_Etalon_WI_RH_ERR_mean"),
    ]),
    ("Mains electrical — UPS voltage & current", [
        ("env_neb_ups_volt_r", "V1:ENV_NEB_UPS_VOLT_R_mean"),
        ("env_ceb_ups_volt_r", "V1:ENV_CEB_UPS_VOLT_R_mean"),
        ("env_web_ups_volt_r", "V1:ENV_WEB_UPS_VOLT_R_mean"),
        ("env_mcb_ips_curr_t", "V1:ENV_MCB_IPS_CURR_T_mean"),
        ("env_ceb_ups_curr_r", "V1:ENV_CEB_UPS_CURR_R_mean"),
    ]),
    ("SNEB / SWEB suspended bench geophones", [
        ("sbe_sneb_geo_we", "V1:SBE_SNEB_GEO_GRWE_raw_mean"),
        ("sbe_sneb_geo_ns", "V1:SBE_SNEB_GEO_GRNS_raw_mean"),
        ("sbe_sweb_geo_we", "V1:SBE_SWEB_GEO_GRWE_raw_mean"),
        ("sbe_sweb_geo_ns", "V1:SBE_SWEB_GEO_GRNS_raw_mean"),
    ]),
]

# ── figure ────────────────────────────────────────────────────────────────────
N = len(PANELS)
fig, axes = plt.subplots(N, 1, figsize=(18, 4 * N), sharex=True)
fig.subplots_adjust(hspace=0.08, top=0.96, bottom=0.05, left=0.08, right=0.90)
colors = plt.cm.tab10(np.linspace(0, 0.9, 10))

for ax, (title, chan_list) in zip(axes, PANELS):
    add_lock_spans(ax)
    plotted = 0
    for i, (col, fullname) in enumerate(chan_list):
        if col not in df.columns:
            print(f"  skipping missing column: {col}")
            continue
        y = df[col].values.copy().astype(float)
        y[~np.isfinite(y)] = np.nan
        # robust normalise + clip outliers
        med = np.nanmedian(y)
        mad = np.nanmedian(np.abs(y - med))
        if mad > 0:
            y = np.where(np.abs(y - med) > 10 * mad, np.nan, y)
            med = np.nanmedian(y)
            mad = np.nanmedian(np.abs(y - med))
            yn = (y - med) / (6 * mad)
        else:
            yn = y - med
        yn_s = pd.Series(yn).rolling(24, center=True, min_periods=3).median().values
        ax.plot(df["time"], yn_s, lw=0.9, alpha=0.85,
                color=colors[i % 10], label=fullname)
        plotted += 1

    ax.set_ylabel("norm. value\n(robust σ units)", fontsize=8)
    ax.set_title(title, fontsize=9, loc="left", pad=3)
    ax.tick_params(labelsize=8)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    if plotted:
        ax.legend(loc="upper left", fontsize=7, ncol=2, framealpha=0.6)
    else:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                ha="center", va="center", color="gray")

    # Glitch amplitude scatter on twin axis
    ax2 = ax.twinx()
    ax2.scatter(glitch["time"], glitch["amplitude"], s=1, color="black",
                alpha=0.25, linewidths=0, zorder=5)
    ax2.set_yscale("log")
    ax2.set_ylabel("glitch amplitude (strain)", fontsize=8)
    ax2.tick_params(labelsize=8)

axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
axes[-1].tick_params(axis="x", labelsize=8, rotation=30)

fig.suptitle("Logbook-referenced channels — full O4 vs glitch amplitude (strain)",
             fontsize=11, fontweight="bold", y=0.99)

OUT_FILE.parent.mkdir(exist_ok=True)
fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
print(f"Saved → {OUT_FILE}")
