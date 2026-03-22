"""
analyze_glitch_vs_sr_tx.py — Step 8: glitch properties vs SR alignment.

Correlates 25-min glitch properties (frequency range, cadence, amplitude)
against the SR marionette TX setpoint (SAT_SR_MAR_TX_SET) over full O4.

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

# Glitches closer than this are considered consecutive for period computation
MAX_PERIOD_S = 60 * 60   # 1h — ignore gaps longer than this


def gps_to_dt(gps: np.ndarray) -> pd.DatetimeIndex:
    return pd.to_datetime(gps, unit="s", origin=GPS_EPOCH, utc=True)


def load_sr(sr_csv: str) -> pd.DataFrame:
    dfs = []
    for p in Path(sr_csv).parent.glob(Path(sr_csv).name):
        dfs.append(pd.read_csv(p))
    if not dfs:
        dfs.append(pd.read_csv(sr_csv))
    df = pd.concat(dfs).sort_values("gps_bin").reset_index(drop=True)
    # Forward-fill SAT_SR_MAR_TX_SET (piecewise constant between shutdowns)
    df["sr_mar_tx_set"] = df["sr_mar_tx_set"].ffill().bfill()
    return df


def lookup_sr_tx(glitch_gps: np.ndarray, sr: pd.DataFrame) -> np.ndarray:
    """For each glitch GPS time, return the nearest sr_mar_tx_set value."""
    idx = np.searchsorted(sr["gps_bin"].values, glitch_gps, side="right") - 1
    idx = np.clip(idx, 0, len(sr) - 1)
    return sr["sr_mar_tx_set"].values[idx]


def lookup_b1p_mag(glitch_gps: np.ndarray, sr: pd.DataFrame) -> np.ndarray:
    """For each glitch GPS time, return the nearest asc_sr_ty_b1p_mag value."""
    if "asc_sr_ty_b1p_mag" not in sr.columns:
        return np.full(len(glitch_gps), np.nan)
    idx = np.searchsorted(sr["gps_bin"].values, glitch_gps, side="right") - 1
    idx = np.clip(idx, 0, len(sr) - 1)
    return sr["asc_sr_ty_b1p_mag"].values[idx]


def compute_inter_glitch_period(gps: np.ndarray) -> np.ndarray:
    """
    Compute inter-glitch interval for each glitch.
    Returns NaN for glitches that are isolated (gap > MAX_PERIOD_S).
    """
    gps = np.sort(gps)
    dt = np.diff(gps)
    # Assign each interval to the later glitch
    periods = np.full(len(gps), np.nan)
    for i, d in enumerate(dt):
        if d < MAX_PERIOD_S:
            periods[i + 1] = d
    return periods


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glitches", required=True)
    ap.add_argument("--sr-csv",   required=True,
                    help="Merged SR O4 CSV (outputs/sr_o4_merged.csv) "
                         "or glob pattern outputs/sr_o4_*.csv")
    ap.add_argument("--output",   default="usecases/25-minute-glitch/glitch_vs_sr_tx.png")
    args = ap.parse_args()

    # ── Load data ──────────────────────────────────────────────────────────
    glitches = pd.read_csv(args.glitches)
    glitches = glitches.sort_values("time").reset_index(drop=True)

    # Handle glob pattern for sr-csv
    sr_files = sorted(Path(".").glob(args.sr_csv)) if "*" in args.sr_csv else [Path(args.sr_csv)]
    if not sr_files:
        sr_files = [Path(args.sr_csv)]
    sr = pd.concat([pd.read_csv(f) for f in sr_files]).sort_values("gps_bin").reset_index(drop=True)
    sr["sr_mar_tx_set"] = sr["sr_mar_tx_set"].ffill().bfill()
    print(f"Glitches: {len(glitches)}")
    print(f"SR bins:  {len(sr)}  ({sr.gps_bin.min()} – {sr.gps_bin.max()})")

    # ── Compute glitch properties ──────────────────────────────────────────
    gps = glitches["time"].values

    sr_tx       = lookup_sr_tx(gps, sr)
    b1p_mag     = lookup_b1p_mag(gps, sr)
    center_freq = (glitches["fstart"].values + glitches["fend"].values) / 2
    bandwidth   = glitches["fend"].values - glitches["fstart"].values
    period_s    = compute_inter_glitch_period(gps)
    period_min  = period_s / 60.0

    # Filter: only use well-defined periods and non-NaN TX
    valid_period = np.isfinite(period_min) & np.isfinite(sr_tx)
    valid_all    = np.isfinite(sr_tx)

    # ── Plot ────────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        "25-min glitch properties vs SR marionette TX setpoint (V1:SAT_SR_MAR_TX_SET) — full O4",
        fontsize=11, y=1.002
    )
    gs = GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    def scatter(ax, x, y, mask, xlabel, ylabel, title, color="tab:blue", alpha=0.15):
        ax.scatter(x[mask], y[mask], s=4, alpha=alpha, color=color, rasterized=True)
        # Running median
        order = np.argsort(x[mask])
        xo, yo = x[mask][order], y[mask][order]
        if len(xo) > 50:
            from scipy.ndimage import uniform_filter1d
            win = max(1, len(xo) // 80)
            ax.plot(xo, uniform_filter1d(yo.astype(float), win), color="tab:red", lw=1.5, label="running median")
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.tick_params(labelsize=8)

    # Panel 1: center frequency vs SR TX
    scatter(fig.add_subplot(gs[0, 0]),
            sr_tx, center_freq, valid_all,
            "SAT_SR_MAR_TX_SET [arb]", "Center frequency [Hz]",
            "Glitch center frequency vs SR TX setpoint")

    # Panel 2: bandwidth vs SR TX
    scatter(fig.add_subplot(gs[0, 1]),
            sr_tx, bandwidth, valid_all,
            "SAT_SR_MAR_TX_SET [arb]", "Bandwidth (fend−fstart) [Hz]",
            "Glitch bandwidth vs SR TX setpoint", color="tab:orange")

    # Panel 3: inter-glitch period vs SR TX
    ax3 = fig.add_subplot(gs[1, 0])
    scatter(ax3, sr_tx, period_min, valid_period,
            "SAT_SR_MAR_TX_SET [arb]", "Inter-glitch period [min]",
            "Glitch cadence vs SR TX setpoint", color="tab:green")
    ax3.axhline(25, color="gray", lw=0.8, ls="--", label="25 min")
    ax3.set_ylim(0, 60)
    ax3.legend(fontsize=7)

    # Panel 4: B1p amplitude vs SR TX (if available)
    ax4 = fig.add_subplot(gs[1, 1])
    valid_b1p = valid_all & np.isfinite(b1p_mag) & (b1p_mag > 0)
    if valid_b1p.sum() > 10:
        scatter(ax4, sr_tx, b1p_mag, valid_b1p,
                "SAT_SR_MAR_TX_SET [arb]", "ASC_SR_TY_B1p_mag [arb]",
                "SR B1p angular amplitude vs SR TX setpoint", color="tab:purple")
    else:
        ax4.text(0.5, 0.5, "B1p mag data not available",
                 ha="center", va="center", transform=ax4.transAxes, fontsize=9)
        ax4.set_title("SR B1p angular amplitude vs SR TX setpoint", fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 1])
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
