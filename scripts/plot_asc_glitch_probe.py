#!/usr/bin/env python3
"""
Plot ASC/LSC error and correction signals around one 25-min glitch
to identify which mirror is affected first.

Two figures:
  1. Overview: 2-second window, all channels
  2. Zoom:     200 ms window around the glitch peak, all channels
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR    = Path(__file__).parent.parent
GLITCH_GPS = 1415578745   # Nov 14 2024 00:09:05 UTC
CSV        = WORKDIR / "outputs" / f"asc_glitch_probe_{GLITCH_GPS}.csv"
OUT_DIR    = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH  = pd.Timestamp("1980-01-06", tz="UTC")

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV)
df = df.sort_values("gps").set_index("gps")

glitch_t = 1415578745.894531   # catalog GPS of glitch peak

# ── Channel groups ────────────────────────────────────────────────────────────
# (label, channel, color)
LSC_CHANS = [
    ("DARM corr", "V1:LSC_DARM_CORR", "black"),
    ("DARM",      "V1:LSC_DARM",      "dimgray"),
    ("MICH",      "V1:LSC_MICH_ERR",  "tab:blue"),
    ("PRCL",      "V1:LSC_PRCL_ERR",  "tab:green"),
    ("SRCL",      "V1:LSC_SRCL_ERR",  "tab:red"),
    ("CARM",      "V1:LSC_CARM_ERR",  "tab:purple"),
]

ASC_ERR_CHANS = [
    ("BS TX err", "V1:ASC_BS_TX_ERR", "tab:blue"),
    ("BS TY err", "V1:ASC_BS_TY_ERR", "tab:cyan"),
    ("PR TX err", "V1:ASC_PR_TX_ERR", "tab:orange"),
    ("PR TY err", "V1:ASC_PR_TY_ERR", "gold"),
    ("SR TX err", "V1:ASC_SR_TX_ERR", "tab:red"),
    ("SR TY err", "V1:ASC_SR_TY_ERR", "salmon"),
]

ASC_CORR_CHANS = [
    ("NI TX corr", "V1:ASC_NI_TX_CORR", "tab:blue"),
    ("NI TY corr", "V1:ASC_NI_TY_CORR", "steelblue"),
    ("WI TX corr", "V1:ASC_WI_TX_CORR", "tab:orange"),
    ("WI TY corr", "V1:ASC_WI_TY_CORR", "goldenrod"),
    ("NE TX corr", "V1:ASC_NE_TX_CORR", "tab:green"),
    ("NE TY corr", "V1:ASC_NE_TY_CORR", "limegreen"),
    ("WE TX corr", "V1:ASC_WE_TX_CORR", "tab:red"),
    ("WE TY corr", "V1:ASC_WE_TY_CORR", "salmon"),
]

def normalise(series, robust=True):
    """Centre and normalise a series to unit RMS (or robust sigma)."""
    v = series.dropna().values
    if len(v) == 0:
        return series
    med = np.median(v)
    if robust:
        mad = np.median(np.abs(v - med))
        sigma = mad * 1.4826 if mad > 0 else 1.0
    else:
        sigma = v.std() if v.std() > 0 else 1.0
    return (series - med) / sigma

def plot_panel(ax, groups, df_win, title, glitch_t, lw=0.8, alpha=0.9):
    """Plot normalised channels on one axes."""
    offset = 0
    yticks, ylabels = [], []
    for label, col, color in groups:
        if col not in df_win.columns:
            continue
        t = df_win.index.values - glitch_t   # seconds relative to glitch
        y = normalise(df_win[col])
        if y.isna().all():
            continue
        ax.plot(t, y + offset, lw=lw, color=color, alpha=alpha, label=label)
        yticks.append(offset)
        ylabels.append(label)
        offset += 6   # vertical spacing
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8, label="glitch")
    ax.set_title(title, fontsize=9, loc="left")
    ax.grid(axis="x", ls=":", alpha=0.4)
    ax.tick_params(axis="x", labelsize=8)

def make_figure(win_before, win_after, suffix):
    t0 = glitch_t - win_before
    t1 = glitch_t + win_after
    df_win = df.loc[t0:t1].copy()

    n_panels = 3
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 4 * n_panels), sharex=True)
    fig.subplots_adjust(hspace=0.10, top=0.93, bottom=0.07, left=0.16, right=0.97)

    plot_panel(axes[0], LSC_CHANS,      df_win,
               "LSC length error signals (DARM / MICH / PRCL / SRCL / CARM)",
               glitch_t)

    plot_panel(axes[1], ASC_ERR_CHANS,  df_win,
               "ASC angular error signals — BS, PR, SR",
               glitch_t)

    plot_panel(axes[2], ASC_CORR_CHANS, df_win,
               "ASC angular correction output — NI, WI, NE, WE  (no direct ERR available)",
               glitch_t)

    axes[-1].set_xlabel("Time relative to glitch (s)", fontsize=9)

    glitch_utc = GPS_EPOCH + pd.to_timedelta(glitch_t, unit="s")
    fig.suptitle(
        f"25-min glitch propagation — GPS {glitch_t}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})\n"
        f"Window: −{win_before}s to +{win_after}s   |   All traces: normalised (zero-mean, unit robust-σ) + offset",
        fontsize=10, fontweight="bold", y=0.98)

    out = OUT_DIR / f"asc_glitch_probe_{GLITCH_GPS}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Figure 1: 2-second overview
make_figure(win_before=1.0, win_after=1.0, suffix="2s")

# Figure 2: 200 ms zoom
make_figure(win_before=0.10, win_after=0.10, suffix="200ms")

# Figure 3: 50 ms ultra-zoom to see propagation delay
make_figure(win_before=0.025, win_after=0.025, suffix="50ms")
