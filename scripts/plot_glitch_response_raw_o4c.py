#!/usr/bin/env python3
"""
Plot raw ASC mirror CORR signals around a glitch with no signal transformation.

- No envelope
- No filtering
- No normalization
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"
TOWERS = [
    ("NI", ["V1:ASC_NI_TX_CORR", "V1:ASC_NI_TY_CORR"], ["tab:blue", "steelblue"]),
    ("WI", ["V1:ASC_WI_TX_CORR", "V1:ASC_WI_TY_CORR"], ["tab:orange", "goldenrod"]),
    ("NE", ["V1:ASC_NE_TX_CORR", "V1:ASC_NE_TY_CORR"], ["tab:green", "limegreen"]),
    ("WE", ["V1:ASC_WE_TX_CORR", "V1:ASC_WE_TY_CORR"], ["tab:red", "salmon"]),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Raw glitch response (no transforms) for ASC CORR channels.")
    p.add_argument("--event-gps", type=float, required=True, help="Glitch peak GPS.")
    p.add_argument("--win", type=float, default=1.0, help="Half-window in seconds (default 1.0 => 2s total).")
    p.add_argument(
        "--csv",
        default=None,
        help="Optional CSV path. Default: outputs/asc_err_corr_recovery_<event>.csv",
    )
    p.add_argument("--out", default=None, help="Optional output PNG path.")
    return p.parse_args()


def choose_csv(event_id: int, explicit_csv: str | None) -> Path:
    if explicit_csv:
        p = Path(explicit_csv)
        if not p.exists():
            raise FileNotFoundError(f"CSV not found: {p}")
        return p
    c1 = OUTPUTS / f"asc_err_corr_recovery_{event_id}.csv"
    c2 = OUTPUTS / f"asc_err_corr_shortlist_{event_id}.csv"
    if c1.exists():
        return c1
    if c2.exists():
        return c2
    raise FileNotFoundError(f"No cached CSV found for event {event_id}: expected {c1.name} or {c2.name}")


def load_series(df: pd.DataFrame, channel: str, event_gps: float) -> tuple[np.ndarray, np.ndarray]:
    d = df[df["channel"] == channel].sort_values("gps")
    t = d["gps"].to_numpy(dtype=float) - event_gps
    v = d["value"].to_numpy(dtype=float)
    return t, v


def apply_robust_ylim(ax: plt.Axes, ys: list[np.ndarray]) -> None:
    vals = [y[np.isfinite(y)] for y in ys if y.size > 0]
    if not vals:
        return
    arr = np.concatenate(vals)
    if arr.size == 0:
        return
    lo = float(np.percentile(arr, 1))
    hi = float(np.percentile(arr, 99))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return
    if hi <= lo:
        span = max(1.0, abs(hi) * 0.1 + 1e-9)
        ax.set_ylim(lo - span, hi + span)
        return
    pad = 0.12 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)


def main() -> int:
    args = parse_args()
    event = float(args.event_gps)
    event_id = int(event)
    csv_path = choose_csv(event_id, args.csv)

    df = pd.read_csv(csv_path)
    have = set(df["channel"])
    win = float(args.win)

    nrows = 1 + len(TOWERS)
    fig, axes = plt.subplots(nrows, 1, figsize=(16, 2.8 * nrows), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.92, bottom=0.08, left=0.10, right=0.98)

    # Panel 1: Hrec raw
    ax0 = axes[0]
    if REF_CHANNEL in have:
        t, v = load_series(df, REF_CHANNEL, event)
        m = (t >= -win) & (t <= win)
        ax0.plot(t[m], v[m], color="black", lw=0.8, label=REF_CHANNEL)
        apply_robust_ylim(ax0, [v[m]])
    ax0.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax0.set_ylabel("Raw")
    ax0.set_title("Hrec raw", loc="left", fontsize=10)
    ax0.grid(axis="x", ls=":", alpha=0.35)
    if REF_CHANNEL in have:
        ax0.legend(fontsize=8, loc="upper right")

    # Tower panels
    for i, (tower, chs, cols) in enumerate(TOWERS, start=1):
        ax = axes[i]
        panel_ys: list[np.ndarray] = []
        for ch, col in zip(chs, cols):
            if ch not in have:
                continue
            t, v = load_series(df, ch, event)
            m = (t >= -win) & (t <= win)
            y = v[m]
            panel_ys.append(y)
            ax.plot(t[m], y, color=col, lw=1.2, alpha=0.9, label=ch.replace("V1:ASC_", ""))
        apply_robust_ylim(ax, panel_ys)
        ax.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
        ax.set_ylabel("Raw")
        ax.set_title(f"{tower} raw CORR", loc="left", fontsize=10)
        ax.grid(axis="x", ls=":", alpha=0.35)
        if panel_ys:
            ax.legend(fontsize=8, loc="upper right")

    for ax in axes:
        ax.set_xlim(-win, win)
    axes[-1].set_xlabel(f"Time relative to glitch GPS {event:.6f} (s)")

    suffix = f"{int(2 * win)}s" if win >= 1 else f"{int(2 * win * 1000)}ms"
    out_png = Path(args.out) if args.out else OUT_DIR / f"glitch_response_raw_{event_id}_{suffix}.png"
    fig.suptitle(
        f"Raw glitch response (no envelope/filter/normalization) | event {event:.6f} | window ±{win:.3f}s",
        fontsize=11,
        fontweight="bold",
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Using CSV: {csv_path}")
    print(f"Saved -> {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
