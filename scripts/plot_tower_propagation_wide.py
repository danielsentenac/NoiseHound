#!/usr/bin/env python3
"""
Tower propagation — wide ±120 s window.

Same stacked layout as plot_tower_propagation.py but using the 240-second
probe (tower_glitch_probe_wide_1415578745.csv).  Four figures:
  1. 240 s overview  — shows full background before and calm-down after
  2.  60 s zoom      — ±30 s around catalog GPS
  3.  15 s zoom      — ±7.5 s around catalog GPS
  4.   6 s zoom      — ±3 s around catalog GPS

Baseline for RMS normalisation: t = -110 to -90 s (well before the glitch).
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR          = Path(__file__).parent.parent
GLITCH_GPS_EXACT = 1415578745.894531
CSV              = WORKDIR / "outputs" / "tower_glitch_probe_wide_1415578745.csv"
OUT_DIR          = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH        = pd.Timestamp("1980-01-06", tz="UTC")

FLO, FHI = 20, 150

# Baseline: quiet stretch well before the glitch
BASELINE_T0 = -110.0
BASELINE_T1 =  -90.0

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV).sort_values("gps").set_index("gps")
t_gps = df.index.values
t_rel = t_gps - GLITCH_GPS_EXACT
print(f"Probe window: t_rel = {t_rel.min():.2f} to {t_rel.max():.2f} s")

# ── Helpers ───────────────────────────────────────────────────────────────────
def get_fs(col):
    idx = df[col].dropna().index
    return round(1.0 / np.median(np.diff(idx.values[:200])))

def process(col):
    if col not in df.columns:
        return None
    v = df[col].fillna(0).values.astype(float)
    fs = get_fs(col)
    sos = butter(4, [FLO, FHI], btype="band", fs=fs, output="sos")
    bp  = sosfiltfilt(sos, v)
    mask_bl = (t_rel >= BASELINE_T0) & (t_rel <= BASELINE_T1)
    rms_bl  = np.std(bp[mask_bl]) if mask_bl.any() else np.std(bp)
    norm    = bp / rms_bl if rms_bl > 0 else bp
    return t_rel, norm

# ── Channel definitions ───────────────────────────────────────────────────────
CHANNELS = [
    ("h(t)",        "V1:Hrec_hoft_16384Hz", "black",          "ref"),
    ("DARM corr",   "V1:LSC_DARM_CORR",     "dimgray",        "ref"),
    ("BS TX err",   "V1:ASC_BS_TX_ERR",     "tab:blue",       "err"),
    ("BS TY err",   "V1:ASC_BS_TY_ERR",     "cornflowerblue", "err"),
    ("PR TX err",   "V1:ASC_PR_TX_ERR",     "tab:green",      "err"),
    ("PR TY err",   "V1:ASC_PR_TY_ERR",     "limegreen",      "err"),
    ("SR TX err",   "V1:ASC_SR_TX_ERR",     "tab:red",        "err"),
    ("SR TY err",   "V1:ASC_SR_TY_ERR",     "salmon",         "err"),
    ("NI TX corr",  "V1:ASC_NI_TX_CORR",    "tab:blue",       "corr"),
    ("NI TY corr",  "V1:ASC_NI_TY_CORR",    "cornflowerblue", "corr"),
    ("WI TX corr",  "V1:ASC_WI_TX_CORR",    "tab:orange",     "corr"),
    ("WI TY corr",  "V1:ASC_WI_TY_CORR",    "goldenrod",      "corr"),
    ("NE TX corr",  "V1:ASC_NE_TX_CORR",    "tab:green",      "corr"),
    ("NE TY corr",  "V1:ASC_NE_TY_CORR",    "limegreen",      "corr"),
    ("WE TX corr",  "V1:ASC_WE_TX_CORR",    "tab:red",        "corr"),
    ("WE TY corr",  "V1:ASC_WE_TY_CORR",    "salmon",         "corr"),
]

# ── Process all channels ──────────────────────────────────────────────────────
traces = {}
for label, col, color, grp in CHANNELS:
    r = process(col)
    if r is not None:
        traces[label] = (color, grp, r[0], r[1])
        print(f"  {label:18s}  fs={get_fs(col)} Hz")

# ── Stacked plot ──────────────────────────────────────────────────────────────
glitch_utc = GPS_EPOCH + pd.to_timedelta(GLITCH_GPS_EXACT, unit="s")

def make_stacked_plot(win_lo, win_hi, suffix, xtick_step):
    chans = [(l, c, g) for l, col, c, g in CHANNELS if l in traces]
    n = len(chans)
    fig, axes = plt.subplots(n, 1, figsize=(18, n * 1.1),
                             sharex=True, sharey=False)
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.0, top=0.95, bottom=0.06,
                        left=0.13, right=0.91)

    group_colors = {"ref": "#e8e8e8", "err": "#fff3e0", "corr": "#e8f5e9"}

    for i, (label, color, grp) in enumerate(chans):
        ax = axes[i]
        _, _, t_arr, norm = traces[label]
        mask = (t_arr >= win_lo) & (t_arr <= win_hi)
        data_win = norm[mask]
        ax.plot(t_arr[mask], data_win, color=color, lw=0.8, alpha=0.9)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.7)
        # Shade the baseline region only if it overlaps the current window
        bl_lo = max(BASELINE_T0, win_lo)
        bl_hi = min(BASELINE_T1, win_hi)
        if bl_lo < bl_hi:
            ax.axvspan(bl_lo, bl_hi, color="lightblue", alpha=0.15, zorder=0)
        ax.set_ylabel(label, fontsize=7.5, rotation=0, ha="right",
                      va="center", labelpad=5)
        ax.tick_params(axis="y", labelsize=6)
        ax.grid(axis="x", ls=":", alpha=0.3)
        ax.set_facecolor(group_colors.get(grp, "white"))
        ax.axhline(0, color="gray", lw=0.4, ls="-", alpha=0.4)

        # Peak amplitude label (max |norm| in the window, in units of baseline RMS)
        peak = np.nanmax(np.abs(data_win)) if len(data_win) > 0 else 0
        ax.text(1.002, 0.5, f"peak={peak:.0f}σ", transform=ax.transAxes,
                fontsize=6.5, color=color, va="center", ha="left", clip_on=False)

        if i == 0 or chans[i-1][2] != grp:
            group_name = {"ref": "Reference", "err": "ERR (sensors)",
                          "corr": "CORR (actuators)"}[grp]
            ax.text(-0.135, 0.5, group_name, transform=ax.transAxes,
                    fontsize=7, color="gray", rotation=90,
                    ha="center", va="center", style="italic")

    xt = np.arange(np.ceil(win_lo / xtick_step) * xtick_step,
                   win_hi + xtick_step * 0.01, xtick_step)
    for ax in axes:
        ax.set_xticks(xt)
    axes[-1].set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
    axes[-1].set_xlabel(
        f"Time relative to catalog GPS {GLITCH_GPS_EXACT:.4f}  (s)", fontsize=9)

    fig.suptitle(
        f"Tower propagation — ERR (sensors) then CORR (actuators)  |  "
        f"GPS {GLITCH_GPS_EXACT:.4f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})\n"
        f"Bandpass {FLO}–{FHI} Hz  |  Normalised to RMS in [{BASELINE_T0:.0f}, {BASELINE_T1:.0f}] s  |  "
        f"Blue band = baseline window  |  Crimson = catalog GPS",
        fontsize=9, fontweight="bold")

    out = OUT_DIR / f"tower_propagation_wide_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

OUT_DIR.mkdir(parents=True, exist_ok=True)
make_stacked_plot(-120.0, 120.0, "240s", xtick_step=20)
make_stacked_plot( -30.0,  30.0, "60s",  xtick_step=5)
make_stacked_plot(  -7.5,   7.5, "15s",  xtick_step=1)
make_stacked_plot(  -3.0,   3.0, "6s",   xtick_step=0.5)
