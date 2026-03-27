#!/usr/bin/env python3
"""
25-minute glitch context plot.

Panel 1 — DARM_CORR 30-55 Hz band-RMS over ±30 min.
           Normalised to the pre-glitch baseline (first 5 minutes of window).
           Shows when the glitch starts, peaks, and ends.

Panel 2 — ASC CORR 30-55 Hz band-RMS (mirror corrections), same normalisation.
           If a mirror correction is elevated during the DARM glitch,
           it is part of the angular drive.

Panel 3 — ASC ERR 30-55 Hz band-RMS (SR, BS, PR error signals).

Each panel has a dashed line at 1.0 (baseline level) and a vertical
line at the catalog GPS (t = 0).
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR          = Path(__file__).parent.parent
GLITCH_GPS       = 1415578745
GLITCH_GPS_EXACT = 1415578745.894531
CSV              = WORKDIR / "outputs" / f"glitch_context_{GLITCH_GPS}.csv"
OUT_DIR          = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH        = pd.Timestamp("1980-01-06", tz="UTC")

# Baseline: first 5 minutes of the extracted window
BASELINE_MIN = -30    # minutes from glitch
BASELINE_MAX = -20    # minutes — "quiet" part before glitch

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV).sort_values("gps_mid")
t_rel = (df["gps_mid"].values - GLITCH_GPS_EXACT) / 60.0   # minutes

print(f"Loaded {len(df)} rows, t_rel range: {t_rel.min():.1f} to {t_rel.max():.1f} min")

# ── Channel groups ────────────────────────────────────────────────────────────
DARM_COL = "V1:LSC_DARM_CORR"

ASC_CORR = [
    ("NI TX", "V1:ASC_NI_TX_CORR", "tab:blue"),
    ("NI TY", "V1:ASC_NI_TY_CORR", "steelblue"),
    ("WI TX", "V1:ASC_WI_TX_CORR", "tab:orange"),
    ("WI TY", "V1:ASC_WI_TY_CORR", "goldenrod"),
    ("NE TX", "V1:ASC_NE_TX_CORR", "tab:green"),
    ("NE TY", "V1:ASC_NE_TY_CORR", "limegreen"),
    ("WE TX", "V1:ASC_WE_TX_CORR", "tab:red"),
    ("WE TY", "V1:ASC_WE_TY_CORR", "salmon"),
]

ASC_ERR = [
    ("SR TX", "V1:ASC_SR_TX_ERR", "tab:red"),
    ("SR TY", "V1:ASC_SR_TY_ERR", "salmon"),
    ("BS TX", "V1:ASC_BS_TX_ERR", "tab:blue"),
    ("PR TX", "V1:ASC_PR_TX_ERR", "tab:green"),
]

# ── Normalise to baseline ────────────────────────────────────────────────────
def norm(col):
    if col not in df.columns:
        return None
    v = df[col].values.astype(float)
    mask_bl = (t_rel >= BASELINE_MIN) & (t_rel <= BASELINE_MAX)
    bl = np.nanmedian(v[mask_bl]) if mask_bl.any() else np.nanmedian(v)
    return v / bl if bl > 0 else v

darm_norm = norm(DARM_COL)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
fig.subplots_adjust(hspace=0.08, top=0.91, bottom=0.08, left=0.10, right=0.97)

def plot_panel(ax, label, col_or_norm, color="black", lw=2.0):
    if isinstance(col_or_norm, str):
        v = norm(col_or_norm)
    else:
        v = col_or_norm
    if v is None:
        return
    ax.plot(t_rel, v, color=color, lw=lw, alpha=0.85, label=label)

# Panel 1: DARM
ax1 = axes[0]
if darm_norm is not None:
    ax1.plot(t_rel, darm_norm, color="black", lw=2, alpha=0.9,
             label="DARM_CORR band-RMS")
ax1.axhline(1.0, color="black", lw=0.8, ls=":", alpha=0.5)
ax1.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8, label="catalog GPS")
ax1.set_ylabel("Band-RMS\n(× baseline)", fontsize=9)
ax1.set_title(
    "DARM_CORR — 30–55 Hz band-RMS (1 value per 100 s block)\n"
    "1.0 = baseline level (pre-glitch median)", fontsize=9, loc="left")
ax1.legend(fontsize=8, loc="upper right")
ax1.grid(axis="x", ls=":", alpha=0.4)

# Panel 2: ASC CORR
ax2 = axes[1]
for label, col, color in ASC_CORR:
    v = norm(col)
    if v is None:
        continue
    ax2.plot(t_rel, v, color=color, lw=1.5, alpha=0.8, label=label)
ax2.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
ax2.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
ax2.set_ylabel("Band-RMS\n(× baseline)", fontsize=9)
ax2.set_title(
    "ASC mirror correction outputs — 30–55 Hz band-RMS\n"
    "Elevated during glitch → correction is driving or reacting at 40 Hz",
    fontsize=9, loc="left")
ax2.legend(fontsize=8, loc="upper right", ncol=2,
           handlelength=1.2, handletextpad=0.4, borderpad=0.4)
ax2.grid(axis="x", ls=":", alpha=0.4)

# Panel 3: ASC ERR
ax3 = axes[2]
for label, col, color in ASC_ERR:
    v = norm(col)
    if v is None:
        continue
    ax3.plot(t_rel, v, color=color, lw=1.5, alpha=0.8, label=label)
ax3.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
ax3.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
ax3.set_ylabel("Band-RMS\n(× baseline)", fontsize=9)
ax3.set_title(
    "ASC error signals (SR, BS, PR) — 30–55 Hz band-RMS\n"
    "Elevated → mirror sees a 40 Hz angular disturbance",
    fontsize=9, loc="left")
ax3.legend(fontsize=8, loc="upper right", ncol=2,
           handlelength=1.2, handletextpad=0.4, borderpad=0.4)
ax3.grid(axis="x", ls=":", alpha=0.4)

# Shared x-axis
xt = np.arange(-30, 31, 5)
for ax in axes:
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
axes[-1].set_xlabel(
    f"Time relative to catalog GPS {GLITCH_GPS_EXACT:.4f} (minutes)", fontsize=9)

glitch_utc = GPS_EPOCH + pd.to_timedelta(GLITCH_GPS_EXACT, unit="s")
fig.suptitle(
    f"25-min glitch — 30-minute context  |  "
    f"GPS {GLITCH_GPS_EXACT:.4f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})\n"
    f"30–55 Hz band-RMS, 1 value per 100 s block — shows glitch onset, duration, end",
    fontsize=10, fontweight="bold")

out = OUT_DIR / f"glitch_context_{GLITCH_GPS}.png"
OUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved → {out}")
