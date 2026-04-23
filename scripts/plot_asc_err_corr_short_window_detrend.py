#!/usr/bin/env python3
"""
Short-window ASC ERR/CORR plot around glitch onset with mirror detrending.

- Uses cached raw extraction from outputs/asc_err_corr_recovery_<event>.csv/json
- No frequency cutoff on mirror channels
- Removes pre-glitch linear drift per channel to reveal onset transient
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"
SHORTLIST = [
    "V1:ASC_NI_TX_CORR",
    "V1:ASC_NI_TY_CORR",
    "V1:ASC_WI_TX_CORR",
    "V1:ASC_WI_TY_CORR",
    "V1:ASC_NE_TX_CORR",
    "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR",
    "V1:ASC_WE_TY_CORR",
    "V1:ASC_BS_TX_ERR",
    "V1:ASC_BS_TY_ERR",
    "V1:ASC_BS_TX_CORR",
    "V1:ASC_BS_TY_CORR",
    "V1:ASC_PR_TX_ERR",
    "V1:ASC_PR_TY_ERR",
    "V1:ASC_PR_TX_CORR",
    "V1:ASC_PR_TY_CORR",
    "V1:ASC_SR_TX_ERR",
    "V1:ASC_SR_TY_ERR",
    "V1:ASC_SR_DOF_TX_CORR",
    "V1:ASC_SR_DOF_TY_CORR",
]

ERR_CH = [c for c in SHORTLIST if c.endswith("_ERR")]
CORR_CH = [c for c in SHORTLIST if c.endswith("_CORR")]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Short-window ASC ERR/CORR plot with pre-glitch detrending.")
    p.add_argument("--event-gps", type=float, required=True, help="Glitch peak GPS.")
    p.add_argument("--win-lo", type=float, default=-0.1, help="Window start relative to glitch (s).")
    p.add_argument("--win-hi", type=float, default=0.5, help="Window end relative to glitch (s).")
    p.add_argument(
        "--fit-hi",
        type=float,
        default=-0.01,
        help="Upper bound of pre-glitch fit window for drift removal (s).",
    )
    p.add_argument(
        "--norm-lo",
        type=float,
        default=-5.0,
        help="Lower bound of normalization baseline (s, default -5.0).",
    )
    p.add_argument(
        "--norm-hi",
        type=float,
        default=-1.0,
        help="Upper bound of normalization baseline (s, default -1.0).",
    )
    p.add_argument(
        "--hrec-band",
        action="store_true",
        help="Apply 30-55 Hz bandpass on HREC (recommended for onset visibility).",
    )
    return p.parse_args()


def color_for(channel: str) -> str:
    if "_NI_" in channel:
        return "#1f77b4"
    if "_WI_" in channel:
        return "#ff7f0e"
    if "_NE_" in channel:
        return "#2ca02c"
    if "_WE_" in channel:
        return "#d62728"
    if "_BS_" in channel:
        return "#9467bd"
    if "_PR_" in channel:
        return "#8c564b"
    if "_SR_" in channel:
        return "#e377c2"
    return "black"


def short_label(channel: str) -> str:
    return channel.replace("V1:ASC_", "")


def series(df: pd.DataFrame, channel: str, event_peak: float) -> tuple[np.ndarray, np.ndarray]:
    d = df[df["channel"] == channel].sort_values("gps")
    t = d["gps"].to_numpy(dtype=float) - event_peak
    v = d["value"].to_numpy(dtype=float)
    return t, v


def hrec_bandpass(t: np.ndarray, v: np.ndarray) -> np.ndarray:
    if len(v) < 5:
        return v
    dt = np.median(np.diff(t))
    if not np.isfinite(dt) or dt <= 0:
        return v
    fs = 1.0 / dt
    if fs <= 2 * 55:
        return v
    sos = butter(4, [30, 55], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def detrend_prewindow(
    t: np.ndarray,
    v: np.ndarray,
    win_lo: float,
    fit_hi: float,
    norm_lo: float,
    norm_hi: float,
) -> np.ndarray:
    fit = (t >= win_lo) & (t <= fit_hi)
    if np.sum(fit) >= 2:
        p1, p0 = np.polyfit(t[fit], v[fit], deg=1)
        trend = p1 * t + p0
        y = v - trend
    else:
        y = v - np.median(v[t < 0]) if np.any(t < 0) else v - np.median(v)
    norm = (t >= norm_lo) & (t <= norm_hi)
    if not np.any(norm):
        norm = fit
    if not np.any(norm):
        norm = t < 0
    if not np.any(norm):
        norm = np.ones_like(t, dtype=bool)
    sig = np.std(y[norm])
    if not np.isfinite(sig) or sig == 0:
        sig = 1.0
    return y / sig


def z_norm_pre(t: np.ndarray, v: np.ndarray, base_lo: float, base_hi: float) -> np.ndarray:
    base = (t >= base_lo) & (t <= base_hi)
    if not np.any(base):
        base = t < 0
    if not np.any(base):
        base = np.ones_like(t, dtype=bool)
    mu = np.mean(v[base])
    sig = np.std(v[base])
    if not np.isfinite(sig) or sig == 0:
        sig = 1.0
    return (v - mu) / sig


def robust_ylim(values: list[np.ndarray], floor: float = 2.5) -> tuple[float, float]:
    if not values:
        return (-floor, floor)
    arr = np.concatenate([v[np.isfinite(v)] for v in values if v.size > 0])
    if arr.size == 0:
        return (-floor, floor)
    a = float(np.percentile(np.abs(arr), 99.0))
    lim = max(floor, a)
    return (-lim, lim)


def plot_window(
    df: pd.DataFrame,
    meta: dict,
    event_peak: float,
    win_lo: float,
    win_hi: float,
    fit_hi: float,
    norm_lo: float,
    norm_hi: float,
    hrec_band: bool,
    out_png: Path,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.92, bottom=0.08, left=0.09, right=0.98)

    ax0, ax1, ax2 = axes

    if REF_CHANNEL in meta.get("found", []):
        t, v = series(df, REF_CHANNEL, event_peak)
        if hrec_band:
            v = hrec_bandpass(t, v)
        y = z_norm_pre(t, v, norm_lo, norm_hi)
        m = (t >= win_lo) & (t <= win_hi)
        ax0.plot(t[m], y[m], color="black", lw=0.9)
    ax0.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    h_title = "V1:Hrec_hoft_16384Hz (30-55 Hz bandpass)" if hrec_band else "V1:Hrec_hoft_16384Hz (raw)"
    ax0.set_title(f"{h_title}, z-normalized on baseline", loc="left", fontsize=10)
    ax0.set_ylabel("Hrec z")
    ax0.grid(axis="x", ls=":", alpha=0.35)

    err_found = [c for c in ERR_CH if c in meta.get("found", [])]
    err_vals: list[np.ndarray] = []
    for ch in err_found:
        t, v = series(df, ch, event_peak)
        y = detrend_prewindow(t, v, win_lo, fit_hi, norm_lo, norm_hi)
        m = (t >= win_lo) & (t <= win_hi)
        err_vals.append(y[m])
        ax1.plot(t[m], y[m], lw=1.1, alpha=0.9, color=color_for(ch), label=short_label(ch))
    ax1.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax1.set_title("ASC ERR (no frequency filter, pre-glitch linear detrend)", loc="left", fontsize=10)
    ax1.set_ylabel("detrended z")
    ax1.grid(axis="x", ls=":", alpha=0.35)
    y0, y1 = robust_ylim(err_vals, floor=2.5)
    ax1.set_ylim(y0, y1)
    if err_found:
        ax1.legend(fontsize=8, ncol=3, loc="upper right")

    corr_found = [c for c in CORR_CH if c in meta.get("found", [])]
    corr_vals: list[np.ndarray] = []
    for ch in corr_found:
        t, v = series(df, ch, event_peak)
        y = detrend_prewindow(t, v, win_lo, fit_hi, norm_lo, norm_hi)
        m = (t >= win_lo) & (t <= win_hi)
        corr_vals.append(y[m])
        ax2.plot(t[m], y[m], lw=1.0, alpha=0.85, color=color_for(ch), label=short_label(ch))
    ax2.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax2.set_title("ASC CORR (no frequency filter, pre-glitch linear detrend)", loc="left", fontsize=10)
    ax2.set_ylabel("detrended z")
    ax2.set_xlabel("Time relative to glitch peak GPS (s)")
    ax2.grid(axis="x", ls=":", alpha=0.35)
    y0, y1 = robust_ylim(corr_vals, floor=2.5)
    ax2.set_ylim(y0, y1)
    if corr_found:
        ax2.legend(fontsize=7, ncol=4, loc="upper right")

    for ax in axes:
        ax.set_xlim(win_lo, win_hi)

    fig.suptitle(
        f"ASC ERR/CORR onset around glitch peak {event_peak:.6f} | window {win_lo:+.3f}..{win_hi:+.3f} s",
        fontsize=12,
        fontweight="bold",
    )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    event_peak = float(args.event_gps)
    event_id = int(event_peak)

    csv_path = OUTPUTS / f"asc_err_corr_recovery_{event_id}.csv"
    meta_path = OUTPUTS / f"asc_err_corr_recovery_{event_id}.json"
    if not csv_path.exists() or not meta_path.exists():
        raise FileNotFoundError(
            f"Missing cached inputs for event {event_id}: {csv_path.name} and/or {meta_path.name}"
        )

    df = pd.read_csv(csv_path)
    meta = json.loads(meta_path.read_text())
    tag = f"m{int(abs(args.win_lo) * 1000)}ms_p{int(args.win_hi * 1000)}ms"
    if args.hrec_band:
        tag += "_hrecbp"
    out_png = OUT_DIR / f"asc_err_corr_shortlist_{event_id}_{tag}_detrend.png"
    plot_window(
        df=df,
        meta=meta,
        event_peak=event_peak,
        win_lo=args.win_lo,
        win_hi=args.win_hi,
        fit_hi=args.fit_hi,
        norm_lo=args.norm_lo,
        norm_hi=args.norm_hi,
        hrec_band=args.hrec_band,
        out_png=out_png,
    )
    print(f"Saved plot: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
