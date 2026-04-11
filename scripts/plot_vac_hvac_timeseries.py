#!/usr/bin/env python3
"""
Plot VAC_VALVEBIGNI_CRYO_TEMP and HVAC_INJ_TE_SATUR_CORR_COLD
vs glitch SNR scatter — full O4.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path

WORKDIR     = Path(__file__).parent.parent
GLITCH_FILE = WORKDIR / "data"    / "full_25min_glitches_ER16-O4b.csv"
OUT_FILE    = WORKDIR / "usecases" / "25-minute-glitch" / "vac_hvac_timeseries.png"

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

# ── Merge quarterly CSVs ──────────────────────────────────────────────────────
parts = sorted((WORKDIR / "outputs").glob("vac_hvac_*_*.csv"))
if not parts:
    raise FileNotFoundError("No vac_hvac_*.csv found in outputs/")

df = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
df = df.sort_values("gps_bin").drop_duplicates("gps_bin")
df["time"] = GPS_EPOCH + pd.to_timedelta(df["gps_bin"], unit="s")

glitch = pd.read_csv(GLITCH_FILE)
glitch["time"] = GPS_EPOCH + pd.to_timedelta(glitch["time"], unit="s")

print(f"Loaded {len(df)} hourly bins from {len(parts)} files")
print(f"Columns: {[c for c in df.columns if c != 'gps_bin']}")

# ── SNR dip boundaries ────────────────────────────────────────────────────────
DIP1_DN_S = GPS_EPOCH + pd.to_timedelta(1415318400, unit="s")
DIP1_DN_E = GPS_EPOCH + pd.to_timedelta(1416528000, unit="s")
DIP2_DN_S = GPS_EPOCH + pd.to_timedelta(1429228800, unit="s")
DIP2_DN_E = GPS_EPOCH + pd.to_timedelta(1431043200, unit="s")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)
fig.subplots_adjust(hspace=0.08, top=0.94, bottom=0.07, left=0.07, right=0.90)

CHANNELS = [
    ("vac_ni_cryo_temp",  "V1:VAC_VALVEBIGNI_CRYO_TEMP",    "NI cryo valve temp",    "#1f77b4"),
    ("hvac_inj_te_cold",  "V1:HVAC_INJ_TE_SATUR_CORR_COLD", "Inj hall HVAC cold TE", "#d62728"),
]

for ax, (col, fullname, label, color) in zip(axes, CHANNELS):
    # dip shading
    ax.axvspan(DIP1_DN_S, DIP1_DN_E, color="#d62728", alpha=0.12, linewidth=0, label="SNR DOWN Dip1")
    ax.axvspan(DIP2_DN_S, DIP2_DN_E, color="#9467bd", alpha=0.12, linewidth=0, label="SNR DOWN Dip2")

    if col in df.columns:
        y = df[col].values.copy().astype(float)
        y[~np.isfinite(y)] = np.nan
        # 24h rolling median to smooth
        ys = pd.Series(y).rolling(24, center=True, min_periods=3).median().values
        ax.plot(df["time"], ys, lw=1.2, color=color, alpha=0.9, label=fullname, zorder=3)
        ax.set_ylabel(f"{label}\n(hourly mean)", fontsize=9)
    else:
        ax.text(0.5, 0.5, f"column '{col}' not found", transform=ax.transAxes,
                ha="center", va="center", color="gray")

    # Glitch amplitude scatter on twin axis
    ax2 = ax.twinx()
    ax2.scatter(glitch["time"], glitch["amplitude"], s=1, color="black",
                alpha=0.20, linewidths=0, zorder=2)
    ax2.set_yscale("log")
    ax2.set_ylabel("glitch amplitude (strain)", fontsize=8, color="gray")
    ax2.tick_params(labelsize=8, labelcolor="gray")

    ax.legend(loc="upper left", fontsize=8, framealpha=0.7)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", ls="--", alpha=0.3)
    ax.set_title(fullname, fontsize=9, loc="left")

axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
axes[-1].tick_params(axis="x", labelsize=8, rotation=30)

fig.suptitle("VAC NI cryo valve temp & Inj hall HVAC cold — full O4 vs glitch amplitude (strain)\n"
             "Red shading = SNR Dip1 (Nov–Dec 2024)   Purple = SNR Dip2 (Apr–May 2025)",
             fontsize=10, fontweight="bold", y=0.98)

OUT_FILE.parent.mkdir(exist_ok=True)
fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
print(f"Saved → {OUT_FILE}")
