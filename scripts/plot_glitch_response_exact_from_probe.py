#!/usr/bin/env python3
"""
Exact reproduction of old glitch_response style from a wide probe CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import butter, hilbert, sosfiltfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

DARM_COL = "V1:LSC_DARM_CORR"
ASC_TOWERS = [
    ("NI", [("V1:ASC_NI_TX_CORR", "tab:blue"), ("V1:ASC_NI_TY_CORR", "steelblue")]),
    ("WI", [("V1:ASC_WI_TX_CORR", "tab:orange"), ("V1:ASC_WI_TY_CORR", "goldenrod")]),
    ("NE", [("V1:ASC_NE_TX_CORR", "tab:green"), ("V1:ASC_NE_TY_CORR", "limegreen")]),
    ("WE", [("V1:ASC_WE_TX_CORR", "tab:red"), ("V1:ASC_WE_TY_CORR", "salmon")]),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reproduce old glitch_response plot from probe CSV.")
    p.add_argument("--csv", required=True, help="Wide probe CSV with gps column and channel columns.")
    p.add_argument("--event-gps", type=float, required=True, help="Catalog glitch GPS (used for t=0).")
    p.add_argument("--win", type=float, default=1.0, help="Half-window in seconds (default 1.0 => 2s).")
    p.add_argument("--flo", type=float, default=30.0, help="Bandpass low cut (Hz).")
    p.add_argument("--fhi", type=float, default=55.0, help="Bandpass high cut (Hz).")
    p.add_argument("--env-smooth-ms", type=float, default=10.0, help="Envelope smoothing in ms.")
    p.add_argument("--baseline-start", type=float, default=-5.0, help="Baseline start (s) rel. to glitch.")
    p.add_argument("--baseline-end", type=float, default=-3.0, help="Baseline end (s) rel. to glitch.")
    p.add_argument(
        "--center-on",
        choices=["catalog", "darm-peak"],
        default="catalog",
        help="Time origin for the vertical dashed line: catalog GPS or local DARM envelope peak.",
    )
    p.add_argument(
        "--center-search-ms",
        type=float,
        default=300.0,
        help="When --center-on darm-peak, search local DARM envelope peak in ±center-search-ms around catalog.",
    )
    p.add_argument(
        "--skip-flat",
        action="store_true",
        help="Skip channels that are flat/constant in raw data (e.g. NI/WI zero lines).",
    )
    p.add_argument(
        "--flat-std-thr",
        type=float,
        default=1e-20,
        help="Std threshold below which a channel is considered flat when --skip-flat is used.",
    )
    p.add_argument("--out", default=None, help="Optional output PNG path.")
    return p.parse_args()


def bandpass(v: np.ndarray, fs: float, flo: float = 30.0, fhi: float = 55.0, order: int = 4) -> np.ndarray:
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def envelope(v: np.ndarray, fs: float, smooth_ms: float) -> np.ndarray:
    env = np.abs(hilbert(v))
    n = max(1, int(round(smooth_ms * 1e-3 * fs)))
    if n > 1:
        env = np.convolve(env, np.ones(n) / n, mode="same")
    return env


def get_fs_from_times(t: np.ndarray) -> float:
    if t.size < 3:
        return 1.0
    diffs = np.diff(t[: min(200, t.size)])
    dt = np.median(diffs)
    if not np.isfinite(dt) or dt <= 0:
        return 1.0
    return float(round(1.0 / dt))


def channel_series(df: pd.DataFrame, col: str) -> tuple[np.ndarray, np.ndarray]:
    d = df[["gps", col]].dropna()
    return d["gps"].to_numpy(dtype=float), d[col].to_numpy(dtype=float)


def main() -> int:
    args = parse_args()
    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path).sort_values("gps")
    print(f"Loaded {df.shape[0]} samples × {df.shape[1]} channels from {csv_path}")

    if DARM_COL not in df.columns:
        raise RuntimeError(f"Missing required column: {DARM_COL}")

    glitch_t_catalog = float(args.event_gps)
    win = float(args.win)

    # DARM
    t_darm, darm_raw = channel_series(df, DARM_COL)
    t_rel_darm_catalog = t_darm - glitch_t_catalog
    darm_fs = get_fs_from_times(t_darm)
    darm_bp = bandpass(darm_raw, darm_fs, args.flo, args.fhi)
    darm_env = envelope(darm_bp, darm_fs, args.env_smooth_ms)

    mask_bl = (t_rel_darm_catalog >= args.baseline_start) & (t_rel_darm_catalog <= args.baseline_end)
    if not np.any(mask_bl):
        mask_bl = t_rel_darm_catalog < 0
    if not np.any(mask_bl):
        mask_bl = np.ones_like(t_rel_darm_catalog, dtype=bool)

    darm_baseline = darm_env[mask_bl].mean()
    if not np.isfinite(darm_baseline) or darm_baseline <= 0:
        darm_baseline = 1.0
    darm_env_norm = darm_env / darm_baseline

    darm_rms_bl = np.std(darm_bp[mask_bl])
    if not np.isfinite(darm_rms_bl) or darm_rms_bl == 0:
        darm_rms_bl = 1.0
    darm_bp_norm = darm_bp / darm_rms_bl

    # Time origin: catalog or local DARM peak
    if args.center_on == "darm-peak":
        search = args.center_search_ms * 1e-3
        m_center = (t_rel_darm_catalog >= -search) & (t_rel_darm_catalog <= search)
        if np.any(m_center):
            t_center = t_darm[m_center][np.argmax(darm_env_norm[m_center])]
        else:
            t_center = glitch_t_catalog
    else:
        t_center = glitch_t_catalog
    t_rel_darm = t_darm - t_center

    # ASC tower traces
    asc_towers_traces: list[tuple[str, list[tuple[str, str, np.ndarray, np.ndarray]]]] = []
    skipped_flat: list[str] = []
    for tower_label, chans in ASC_TOWERS:
        tower_traces = []
        for col, color in chans:
            if col not in df.columns:
                continue
            t, v = channel_series(df, col)
            if t.size < 3:
                continue
            if args.skip_flat and np.nanstd(v) <= args.flat_std_thr:
                skipped_flat.append(col)
                continue
            fs = get_fs_from_times(t)
            bp = bandpass(v, fs, args.flo, args.fhi)
            env = envelope(bp, fs, args.env_smooth_ms)
            t_rel_catalog = t - glitch_t_catalog
            bl = (t_rel_catalog >= args.baseline_start) & (t_rel_catalog <= args.baseline_end)
            if not np.any(bl):
                bl = t_rel_catalog < 0
            if not np.any(bl):
                bl = np.ones_like(t_rel_catalog, dtype=bool)
            bl_mean = env[bl].mean() if np.any(bl) else 1.0
            if not np.isfinite(bl_mean) or bl_mean <= 0:
                bl_mean = 1.0
            env_norm = env / bl_mean
            tower_traces.append((col, color, t - t_center, env_norm))
        if tower_traces:
            asc_towers_traces.append((tower_label, tower_traces))

    # Plot (exact old layout style)
    use_ms = win <= 0.5
    scale = 1e3 if use_ms else 1.0
    unit = "ms" if use_ms else "s"
    suffix = f"{int(win * 2 * 1000)}ms" if use_ms else f"{int(win * 2)}s"

    mask_win_darm = (t_rel_darm >= -win) & (t_rel_darm <= win)
    t_plot_darm = t_rel_darm[mask_win_darm] * scale

    n_rows = 2 + len(asc_towers_traces)
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 3 * n_rows), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.13, right=0.97)

    ax1 = axes[0]
    ax1.plot(
        t_plot_darm,
        darm_bp_norm[mask_win_darm],
        color="black",
        lw=0.8,
        alpha=0.9,
        label=f"{DARM_COL} (30–55 Hz)",
    )
    ax1.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax1.set_ylabel("Amplitude\n(pre-glitch RMS)", fontsize=9)
    ax1.set_title(f"{DARM_COL} — bandpassed waveform", fontsize=9, loc="left")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(axis="x", ls=":", alpha=0.4)

    ax2 = axes[1]
    ax2.plot(
        t_plot_darm,
        darm_env_norm[mask_win_darm],
        color="black",
        lw=1.5,
        alpha=0.9,
        label=f"{DARM_COL} envelope",
    )
    ax2.axhline(1.0, color="black", lw=0.8, ls=":", alpha=0.5, label="baseline = 1")
    ax2.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax2.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
    ax2.set_title(f"{DARM_COL} — Hilbert envelope  (1.0 = pre-glitch level)", fontsize=9, loc="left")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(axis="x", ls=":", alpha=0.4)

    for i, (tower_label, tower_traces) in enumerate(asc_towers_traces):
        ax = axes[2 + i]
        for col, color, t_arr, env_norm in tower_traces:
            asc_mask = (t_arr >= -win) & (t_arr <= win)
            ax.plot(t_arr[asc_mask] * scale, env_norm[asc_mask], color=color, lw=1.2, alpha=0.9, label=col)
        ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
        ax.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
        ax.set_title(f"{tower_label} — ASC CORR Hilbert envelope", fontsize=9, loc="left")
        ax.legend(fontsize=8, loc="upper right", handlelength=1.5)
        ax.grid(axis="x", ls=":", alpha=0.4)

    xtick_step = 50.0 if use_ms else (0.2 if win <= 1.0 else 1.0)
    xt = np.arange(-win * scale, win * scale + xtick_step * 0.01, xtick_step)
    for ax in axes:
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
    center_label = "catalog GPS" if args.center_on == "catalog" else "DARM local peak"
    axes[-1].set_xlabel(f"Time relative to {center_label} ({unit})", fontsize=9)

    glitch_utc = GPS_EPOCH + pd.to_timedelta(glitch_t_catalog, unit="s")
    shift_ms = (t_center - glitch_t_catalog) * 1e3
    fig.suptitle(
        "25-min glitch — did DARM glitch drive ASC correction responses?\n"
        f"GPS {glitch_t_catalog:.4f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})  "
        f"Window: ±{win * 1e3:.0f} ms   Bandpass {args.flo:.0f}-{args.fhi:.0f} Hz   "
        f"Envelope smoothed {args.env_smooth_ms:.0f} ms   center={args.center_on} ({shift_ms:+.1f} ms)",
        fontsize=10,
        fontweight="bold",
    )
    if skipped_flat:
        fig.text(
            0.13,
            0.012,
            "Skipped flat channels: " + ", ".join(skipped_flat),
            fontsize=8,
            color="dimgray",
        )

    out = Path(args.out) if args.out else OUT_DIR / f"glitch_response_{int(glitch_t_catalog)}_{suffix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
