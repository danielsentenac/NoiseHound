"""
analyze_sr_steps.py — Option 2: glitch SNR response to discrete SR alignment steps.

Detects step events in V1:SAT_SR_MAR_TX_SET, then for each step measures
mean glitch SNR in windows before and after. Plots ΔSNR vs ΔTX — free of
the slow O4-long time confound.

Usage:
    python scripts/analyze_sr_steps.py \\
        --glitches  data/full_25min_glitches_ER16-O4b.csv \\
        --sr-csv    outputs/sr_o4_merged.csv \\
        --output    usecases/25-minute-glitch/sr_steps_vs_snr.png
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
WINDOW_H  = 72    # hours before/after step to average SNR
MIN_GLITCHES = 5  # minimum glitches required in each window


def detect_steps(tx: np.ndarray, gps: np.ndarray,
                 threshold: float = 15.0, min_gap_h: int = 48) -> list[dict]:
    """
    Detect persistent level changes in TX using a sliding-window median comparison.
    A step at index i means: median(tx[i-6:i]) vs median(tx[i:i+6]) differs by > threshold.
    min_gap_h prevents double-counting steps within the same relock event.
    """
    steps = []
    half = 6
    last_step_idx = -min_gap_h
    for i in range(half, len(tx) - half):
        before = np.nanmedian(tx[max(0, i-half):i])
        after  = np.nanmedian(tx[i:i+half])
        delta  = after - before
        if abs(delta) > threshold and (i - last_step_idx) > min_gap_h:
            steps.append({
                "idx": i,
                "gps": int(gps[i]),
                "tx_before": before,
                "tx_after":  after,
                "delta_tx":  delta,
            })
            last_step_idx = i
    return steps


def mean_snr_window(gps_step: int, gl: pd.DataFrame,
                    window_s: int, side: str) -> float | None:
    if side == "before":
        mask = (gl["time"] >= gps_step - window_s) & (gl["time"] < gps_step)
    else:
        mask = (gl["time"] > gps_step) & (gl["time"] <= gps_step + window_s)
    sub = gl[mask]
    if len(sub) < MIN_GLITCHES:
        return None
    return float(sub["snr"].median())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glitches", required=True)
    ap.add_argument("--sr-csv",   required=True)
    ap.add_argument("--output",   default="usecases/25-minute-glitch/sr_steps_vs_snr.png")
    ap.add_argument("--threshold", type=float, default=15.0,
                    help="Min TX change (arb) to count as a step")
    args = ap.parse_args()

    gl = pd.read_csv(args.glitches).sort_values("time").reset_index(drop=True)
    sr = pd.read_csv(args.sr_csv).sort_values("gps_bin").reset_index(drop=True)
    sr["sr_mar_tx_set"] = sr["sr_mar_tx_set"].ffill().bfill()

    tx  = sr["sr_mar_tx_set"].values
    gps = sr["gps_bin"].values
    window_s = WINDOW_H * 3600

    steps = detect_steps(tx, gps, threshold=args.threshold)
    print(f"Detected {len(steps)} steps (threshold={args.threshold} arb)")

    # For each step compute ΔSNR
    records = []
    for s in steps:
        snr_before = mean_snr_window(s["gps"], gl, window_s, "before")
        snr_after  = mean_snr_window(s["gps"], gl, window_s, "after")
        if snr_before is None or snr_after is None:
            continue
        s["snr_before"] = snr_before
        s["snr_after"]  = snr_after
        s["delta_snr"]  = snr_after - snr_before
        s["dt"] = GPS_EPOCH + pd.to_timedelta(s["gps"], unit="s")
        records.append(s)

    df = pd.DataFrame(records)
    print(f"Steps with enough glitches in both windows: {len(df)}")
    if df.empty:
        print("No usable steps found.")
        return

    # ── Plot ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"SNR response to discrete SR alignment steps (±{WINDOW_H}h window, "
        f"threshold={args.threshold} arb)",
        fontsize=11, y=1.02
    )

    # Panel 1: ΔSNR vs ΔTX, coloured by time
    ax = axes[0]
    year_frac = 1980 + df["gps"] / (365.25 * 86400)
    sc = ax.scatter(df["delta_tx"], df["delta_snr"],
                    c=year_frac, cmap="plasma", s=60, zorder=3)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axvline(0, color="gray", lw=0.8, ls="--")
    # Linear fit
    m, b = np.polyfit(df["delta_tx"], df["delta_snr"], 1)
    xfit = np.linspace(df["delta_tx"].min(), df["delta_tx"].max(), 100)
    ax.plot(xfit, m*xfit + b, color="red", lw=1.5,
            label=f"fit: slope={m:.2f} SNR/arb")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Year", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_xlabel("ΔTX (after − before step) [arb]", fontsize=9)
    ax.set_ylabel("ΔSNR (median after − before) [arb]", fontsize=9)
    ax.set_title("ΔSNR vs ΔTX at each alignment step", fontsize=9)
    ax.legend(fontsize=8)

    # Panel 2: SNR before vs SNR after each step
    ax2 = axes[1]
    ax2.scatter(df["snr_before"], df["snr_after"],
                c=year_frac, cmap="plasma", s=60, zorder=3)
    lim = [min(df["snr_before"].min(), df["snr_after"].min()) - 10,
           max(df["snr_before"].max(), df["snr_after"].max()) + 10]
    ax2.plot(lim, lim, color="gray", lw=0.8, ls="--", label="no change")
    ax2.set_xlim(lim); ax2.set_ylim(lim)
    ax2.set_xlabel("Median SNR before step", fontsize=9)
    ax2.set_ylabel("Median SNR after step", fontsize=9)
    ax2.set_title("SNR before vs after alignment step", fontsize=9)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
