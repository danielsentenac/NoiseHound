"""
analyze_glitch_vs_sr_tx.py — Step 8: glitch properties vs SR alignment.

For each glitch in the O4 catalog, looks up the contemporaneous
SAT_SR_MAR_TX_SET value and ASC_SR_TY_B1p_mag_10Hz_FS_mean.
Plots 3 glitch properties as a function of SR TX marionette setpoint:
  1. Center frequency  (fstart+fend)/2
  2. Bandwidth         fend-fstart
  3. Amplitude (SNR from catalog, and SR B1p mag from trend)

Display: binned medians + IQR bands over raw scatter (2D histogram).

Usage:
    python scripts/analyze_glitch_vs_sr_tx.py \\
        --glitches  data/full_25min_glitches_ER16-O4b.csv \\
        --sr-csv    outputs/sr_o4_merged.csv \\
        --output    usecases/25-minute-glitch/glitch_vs_sr_tx.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")


def lookup(gps: np.ndarray, sr: pd.DataFrame, col: str) -> np.ndarray:
    idx = np.searchsorted(sr["gps_bin"].values, gps, side="right") - 1
    idx = np.clip(idx, 0, len(sr) - 1)
    return sr[col].values[idx]


def binned_stats(x: np.ndarray, y: np.ndarray, n_bins: int = 40):
    """Return bin centres, medians, 25th and 75th percentiles."""
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    edges = np.percentile(x, np.linspace(0, 100, n_bins + 1))
    edges = np.unique(edges)
    centres, med, p25, p75, counts = [], [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (x >= lo) & (x < hi)
        if sel.sum() < 5:
            continue
        centres.append((lo + hi) / 2)
        med.append(np.median(y[sel]))
        p25.append(np.percentile(y[sel], 25))
        p75.append(np.percentile(y[sel], 75))
        counts.append(sel.sum())
    return (np.array(centres), np.array(med),
            np.array(p25), np.array(p75), np.array(counts))


def panel(ax, x, y, xlabel, ylabel, title, color, x_clip=None):
    if x_clip is not None:
        mask = (x >= x_clip[0]) & (x <= x_clip[1]) & np.isfinite(y)
        x, y = x[mask], y[mask]

    # 2D histogram background
    h, xe, ye = np.histogram2d(x, y, bins=60)
    ax.pcolormesh(xe, ye, h.T, cmap="Blues", alpha=0.6, rasterized=True)

    # Binned median + IQR
    cx, med, p25, p75, _ = binned_stats(x, y, n_bins=40)
    ax.plot(cx, med, color=color, lw=2, label="median")
    ax.fill_between(cx, p25, p75, color=color, alpha=0.25, label="IQR")

    ax.set_xlabel("V1:SAT_SR_MAR_TX_SET [arb]", fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glitches", required=True)
    ap.add_argument("--sr-csv",   required=True)
    ap.add_argument("--output",   default="usecases/25-minute-glitch/glitch_vs_sr_tx.png")
    args = ap.parse_args()

    # ── Load ────────────────────────────────────────────────────────────────
    gl = pd.read_csv(args.glitches).sort_values("time").reset_index(drop=True)
    sr = pd.read_csv(args.sr_csv).sort_values("gps_bin").reset_index(drop=True)
    sr["sr_mar_tx_set"] = sr["sr_mar_tx_set"].ffill().bfill()
    print(f"Glitches: {len(gl)}   SR bins: {len(sr)}")

    # ── Compute properties ──────────────────────────────────────────────────
    gps = gl["time"].values
    tx       = lookup(gps, sr, "sr_mar_tx_set")
    b1p_mag  = lookup(gps, sr, "asc_sr_ty_b1p_mag") if "asc_sr_ty_b1p_mag" in sr.columns else np.full(len(gps), np.nan)
    center_f = (gl["fstart"].values + gl["fend"].values) / 2
    bandwidth= gl["fend"].values - gl["fstart"].values
    snr      = gl["snr"].values

    # Restrict to main cluster (exclude extreme outliers)
    x_clip = (np.percentile(tx[np.isfinite(tx)], 0.5),
              np.percentile(tx[np.isfinite(tx)], 99.5))
    print(f"TX range (0.5–99.5 pct): {x_clip[0]:.1f} – {x_clip[1]:.1f}")

    # ── Plot ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle(
        "25-min glitch properties vs SR marionette TX setpoint (V1:SAT_SR_MAR_TX_SET) — full O4",
        fontsize=11, y=1.01
    )
    gs = GridSpec(1, 3, figure=fig, wspace=0.38)

    panel(fig.add_subplot(gs[0, 0]), tx, center_f,
          "V1:SAT_SR_MAR_TX_SET [arb]", "Center frequency [Hz]",
          "Glitch center frequency", "tab:blue", x_clip)

    panel(fig.add_subplot(gs[0, 1]), tx, bandwidth,
          "V1:SAT_SR_MAR_TX_SET [arb]", "Bandwidth (fend−fstart) [Hz]",
          "Glitch bandwidth", "tab:orange", x_clip)

    valid_b1p = np.isfinite(b1p_mag) & (b1p_mag > 0)
    if valid_b1p.sum() > 100:
        panel(fig.add_subplot(gs[0, 2]), tx, np.log10(b1p_mag + 1e-12),
              "V1:SAT_SR_MAR_TX_SET [arb]", "log₁₀(ASC_SR_TY_B1p_mag) [arb]",
              "SR B1p angular amplitude", "tab:purple", x_clip)
    else:
        panel(fig.add_subplot(gs[0, 2]), tx, snr,
              "V1:SAT_SR_MAR_TX_SET [arb]", "SNR",
              "Glitch SNR", "tab:purple", x_clip)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
