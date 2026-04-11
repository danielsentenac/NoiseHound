#!/usr/bin/env python3
"""
Superpose reference glitch and mirror glitches on a 1 ms grid.

Reference channel defaults to V1:LSC_DARM_CORR (Hrec_hoft is not present in the probe CSV).
Mirror channels default to NE/WE TX/TY CORR.
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
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Superpose glitch response on 1 ms grid.")
    p.add_argument("--event-gps", type=float, required=True, help="Catalog glitch GPS around which data are extracted.")
    p.add_argument(
        "--csv",
        default=None,
        help="Wide probe CSV path. Default: outputs/asc_glitch_probe_<event>.csv",
    )
    p.add_argument("--ref-channel", default="V1:LSC_DARM_CORR", help="Reference glitch channel.")
    p.add_argument(
        "--mirror-channels",
        nargs="+",
        default=[
            "V1:ASC_NE_TX_CORR",
            "V1:ASC_NE_TY_CORR",
            "V1:ASC_WE_TX_CORR",
            "V1:ASC_WE_TY_CORR",
        ],
        help="Mirror channels to superpose.",
    )
    p.add_argument("--flo", type=float, default=30.0, help="Bandpass low cut (Hz).")
    p.add_argument("--fhi", type=float, default=55.0, help="Bandpass high cut (Hz).")
    p.add_argument("--env-smooth-ms", type=float, default=3.0, help="Envelope smoothing in ms (small to preserve timing).")
    p.add_argument("--grid-ms", type=float, default=1.0, help="Interpolation grid in ms (default 1.0).")
    p.add_argument("--plot-ms-left", type=float, default=120.0, help="Left plot window in ms after alignment.")
    p.add_argument("--plot-ms-right", type=float, default=220.0, help="Right plot window in ms after alignment.")
    p.add_argument("--align-search-ms", type=float, default=300.0, help="Search half-window in ms for reference local peak.")
    p.add_argument(
        "--align-mode",
        choices=["catalog", "local-peak"],
        default="catalog",
        help="Alignment origin: catalog GPS (default) or local ref-envelope peak.",
    )
    p.add_argument("--lag-search-left-ms", type=float, default=60.0, help="Left lag-search window in ms.")
    p.add_argument("--lag-search-right-ms", type=float, default=180.0, help="Right lag-search window in ms.")
    p.add_argument(
        "--wave-norm-ms",
        type=float,
        default=20.0,
        help="Normalize each waveform by its abs-peak within ±wave-norm-ms around alignment origin. Set <=0 to disable.",
    )
    p.add_argument(
        "--baseline-start",
        type=float,
        default=-5.0,
        help="Baseline start in seconds relative to catalog GPS (for normalization).",
    )
    p.add_argument(
        "--baseline-end",
        type=float,
        default=-3.0,
        help="Baseline end in seconds relative to catalog GPS (for normalization).",
    )
    p.add_argument("--out", default=None, help="Output PNG path.")
    return p.parse_args()


def default_csv(event_id: int) -> Path:
    return OUTPUTS / f"asc_glitch_probe_{event_id}.csv"


def get_fs(df: pd.DataFrame, col: str) -> float:
    idx = df[col].dropna().index.to_numpy(dtype=float)
    if idx.size < 3:
        return 1.0
    dt = np.median(np.diff(idx[: min(200, idx.size)]))
    if not np.isfinite(dt) or dt <= 0:
        return 1.0
    return float(1.0 / dt)


def bandpass(v: np.ndarray, fs: float, flo: float, fhi: float, order: int = 4) -> np.ndarray:
    if v.size < 8 or fs <= 2.2 * fhi:
        return v.copy()
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def envelope(v: np.ndarray, fs: float, smooth_ms: float) -> np.ndarray:
    env = np.abs(hilbert(v))
    n = max(1, int(round((smooth_ms * 1e-3) * fs)))
    if n > 1:
        env = np.convolve(env, np.ones(n) / n, mode="same")
    return env


def baseline_mask(t_rel: np.ndarray, start: float, end: float) -> np.ndarray:
    m = (t_rel >= start) & (t_rel <= end)
    if not np.any(m):
        m = t_rel < 0
    if not np.any(m):
        m = np.ones_like(t_rel, dtype=bool)
    return m


def short(ch: str) -> str:
    return ch.replace("V1:ASC_", "")


def robust_abs_limit(arrays: list[np.ndarray], floor: float = 0.8) -> float:
    vals = [a[np.isfinite(a)] for a in arrays if a.size > 0]
    if not vals:
        return floor
    z = np.concatenate(vals)
    if z.size == 0:
        return floor
    lim = float(np.percentile(np.abs(z), 99.0))
    return max(floor, lim)


def main() -> int:
    args = parse_args()
    event = float(args.event_gps)
    event_id = int(event)
    csv_path = Path(args.csv) if args.csv else default_csv(event_id)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path).sort_values("gps").set_index("gps")
    if args.ref_channel not in df.columns:
        raise RuntimeError(f"Reference channel missing: {args.ref_channel}")

    channels = [args.ref_channel] + [c for c in args.mirror_channels if c in df.columns]
    if len(channels) <= 1:
        raise RuntimeError("No mirror channels found in CSV.")

    t_ref = df.index.to_numpy(dtype=float)
    t_rel = t_ref - event

    # Build 1 ms grid with enough margin to remain valid after optional local-peak alignment.
    left = -(args.plot_ms_left + args.align_search_ms + 20.0) * 1e-3
    right = (args.plot_ms_right + args.align_search_ms + 20.0) * 1e-3
    step = args.grid_ms * 1e-3
    t_grid = np.arange(left, right + 0.5 * step, step)

    waves: dict[str, np.ndarray] = {}
    envs: dict[str, np.ndarray] = {}
    used: list[str] = []
    skipped: list[str] = []
    for ch in channels:
        if ch not in df.columns:
            skipped.append(ch)
            continue
        s = df[ch].astype(float)
        valid = s.dropna()
        if valid.empty or np.nanstd(valid.to_numpy()) == 0:
            skipped.append(ch)
            continue
        t = valid.index.to_numpy(dtype=float)
        v = valid.to_numpy(dtype=float)
        tr = t - event
        fs = get_fs(df, ch)
        vb = bandpass(v, fs, args.flo, args.fhi)
        env = envelope(vb, fs, args.env_smooth_ms)
        bl = baseline_mask(tr, args.baseline_start, args.baseline_end)
        rms = np.std(vb[bl]) if np.any(bl) else np.std(vb)
        if not np.isfinite(rms) or rms == 0:
            rms = 1.0
        vb_norm = vb / rms
        env_bl = np.mean(env[bl]) if np.any(bl) else np.mean(env)
        if not np.isfinite(env_bl) or env_bl <= 0:
            env_bl = 1.0
        env_norm = env / env_bl
        waves[ch] = np.interp(t_grid, tr, vb_norm)
        envs[ch] = np.interp(t_grid, tr, env_norm)
        used.append(ch)

    if args.ref_channel not in used:
        raise RuntimeError(f"Reference channel is unavailable/constant after filtering: {args.ref_channel}")

    # Alignment origin
    ref_env = envs[args.ref_channel]
    if args.align_mode == "local-peak":
        align_m = (t_grid >= -args.align_search_ms * 1e-3) & (t_grid <= args.align_search_ms * 1e-3)
        i_ref_peak = np.argmax(ref_env[align_m])
        t_ref_peak = t_grid[align_m][i_ref_peak]
    else:
        t_ref_peak = 0.0
    t_al = t_grid - t_ref_peak

    # Estimate per-channel lag from envelope peak within lag-search window after alignment.
    lag_m = (t_al >= -args.lag_search_left_ms * 1e-3) & (t_al <= args.lag_search_right_ms * 1e-3)
    lags_ms: dict[str, float] = {}
    for ch in used:
        i = np.argmax(envs[ch][lag_m])
        tpk = t_al[lag_m][i]
        lags_ms[ch] = float(tpk * 1e3)

    # Plot only requested display window after alignment.
    plot_m = (t_al >= -args.plot_ms_left * 1e-3) & (t_al <= args.plot_ms_right * 1e-3)
    x_ms = t_al[plot_m] * 1e3

    # Optional local peak normalization for waveform overlays to avoid flattening.
    waves_plot: dict[str, np.ndarray] = {}
    if args.wave_norm_ms > 0:
        m_norm = (t_al >= -args.wave_norm_ms * 1e-3) & (t_al <= args.wave_norm_ms * 1e-3)
    else:
        m_norm = plot_m
    for ch in used:
        y = waves[ch].copy()
        amp = float(np.max(np.abs(y[m_norm]))) if np.any(m_norm) else float(np.max(np.abs(y[plot_m])))
        if np.isfinite(amp) and amp > 0:
            y = y / amp
        waves_plot[ch] = y

    fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
    fig.subplots_adjust(hspace=0.10, top=0.90, bottom=0.09, left=0.10, right=0.98)
    axw, axe = axes

    # Waveform overlay
    axw.plot(x_ms, waves_plot[args.ref_channel][plot_m], color="black", lw=1.6, label=f"{short(args.ref_channel)} (ref)")
    color_cycle = ["tab:green", "limegreen", "tab:red", "salmon", "tab:blue", "steelblue", "tab:orange", "goldenrod"]
    for i, ch in enumerate([c for c in used if c != args.ref_channel]):
        col = color_cycle[i % len(color_cycle)]
        axw.plot(x_ms, waves_plot[ch][plot_m], color=col, lw=1.1, alpha=0.95, label=f"{short(ch)}")
    axw.axvline(0.0, color="crimson", lw=1.5, ls="--", alpha=0.85, label="ref peak")
    wave_lim = robust_abs_limit([waves_plot[ch][plot_m] for ch in used], floor=0.8)
    axw.set_ylim(-1.15 * wave_lim, 1.15 * wave_lim)
    axw.set_ylabel("Bandpassed waveform\n(local-peak normalized)")
    axw.set_title("Waveform superposition (30–55 Hz), aligned on reference", loc="left", fontsize=10)
    axw.grid(axis="x", ls=":", alpha=0.35)
    axw.legend(fontsize=8, ncol=3, loc="upper right")

    # Envelope overlay with lag annotations
    axe.plot(x_ms, envs[args.ref_channel][plot_m], color="black", lw=2.0, label=f"{short(args.ref_channel)} lag = 0.0 ms")
    for i, ch in enumerate([c for c in used if c != args.ref_channel]):
        col = color_cycle[i % len(color_cycle)]
        lag = lags_ms[ch]
        axe.plot(x_ms, envs[ch][plot_m], color=col, lw=1.4, alpha=0.95, label=f"{short(ch)} lag = {lag:+.1f} ms")
        axe.axvline(lag, color=col, lw=1.0, ls=":", alpha=0.45)
    axe.axvline(0.0, color="crimson", lw=1.5, ls="--", alpha=0.85)
    axe.axhline(1.0, color="gray", lw=0.9, ls=":", alpha=0.6)
    axe.set_ylabel("Envelope\n(× pre-glitch)")
    axe.set_title(f"Hilbert envelope superposition on {args.grid_ms:.0f} ms grid", loc="left", fontsize=10)
    env_vals = np.concatenate([envs[ch][plot_m] for ch in used if np.any(plot_m)])
    if env_vals.size > 0:
        lo = float(np.percentile(env_vals, 1.0))
        hi = float(np.percentile(env_vals, 99.0))
        lo = min(lo, 1.0)
        hi = max(hi, 1.0)
        if hi <= lo:
            hi = lo + 0.2
        pad = 0.10 * (hi - lo)
        axe.set_ylim(lo - pad, hi + pad)
    axe.grid(axis="x", ls=":", alpha=0.35)
    axe.legend(fontsize=8, ncol=2, loc="upper right")

    for ax in axes:
        ax.set_xlim(-args.plot_ms_left, args.plot_ms_right)
    axes[-1].set_xlabel("Time relative to reference peak (ms)")

    utc = GPS_EPOCH + pd.to_timedelta(event, unit="s")
    fig.suptitle(
        f"Reference/mirror glitch superposition at 1 ms precision\n"
        f"Catalog GPS {event:.6f} ({utc.strftime('%Y-%m-%d %H:%M:%S UTC')}), "
        f"ref={short(args.ref_channel)}, align mode={args.align_mode}, align shift={t_ref_peak*1e3:+.2f} ms",
        fontsize=11,
        fontweight="bold",
    )
    if skipped:
        fig.text(
            0.10,
            0.015,
            "Skipped missing/constant channels: " + ", ".join(short(c) for c in skipped),
            fontsize=8,
            color="dimgray",
        )

    out = Path(args.out) if args.out else OUT_DIR / f"glitch_superpose_1ms_{event_id}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)

    print(f"Loaded: {csv_path}")
    print(f"Saved -> {out}")
    print("Lags (ms):")
    for ch in [args.ref_channel] + [c for c in used if c != args.ref_channel]:
        lag = 0.0 if ch == args.ref_channel else lags_ms[ch]
        print(f"  {short(ch):24s} {lag:+7.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
