#!/usr/bin/env python3
"""
Reproduce the old glitch_response style for current ASC ERR/CORR shortlist.

Panels:
1) Hrec bandpassed waveform (30-55 Hz), normalized by pre-glitch RMS
2) Hrec Hilbert envelope, normalized to pre-glitch baseline (=1)
3+) One panel per tower with ASC ERR/CORR Hilbert envelopes (normalized)

Default output is a 2 s window (±1 s around glitch peak).
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

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"

TOWERS: list[tuple[str, list[str]]] = [
    ("NI", ["V1:ASC_NI_TX_CORR", "V1:ASC_NI_TY_CORR"]),
    ("WI", ["V1:ASC_WI_TX_CORR", "V1:ASC_WI_TY_CORR"]),
    ("NE", ["V1:ASC_NE_TX_CORR", "V1:ASC_NE_TY_CORR"]),
    ("WE", ["V1:ASC_WE_TX_CORR", "V1:ASC_WE_TY_CORR"]),
    ("BS", ["V1:ASC_BS_TX_ERR", "V1:ASC_BS_TY_ERR", "V1:ASC_BS_TX_CORR", "V1:ASC_BS_TY_CORR"]),
    ("PR", ["V1:ASC_PR_TX_ERR", "V1:ASC_PR_TY_ERR", "V1:ASC_PR_TX_CORR", "V1:ASC_PR_TY_CORR"]),
    ("SR", ["V1:ASC_SR_TX_ERR", "V1:ASC_SR_TY_ERR", "V1:ASC_SR_DOF_TX_CORR", "V1:ASC_SR_DOF_TY_CORR"]),
]

TOWER_COLORS = {
    "NI": "tab:blue",
    "WI": "tab:orange",
    "NE": "tab:green",
    "WE": "tab:red",
    "BS": "tab:purple",
    "PR": "tab:brown",
    "SR": "tab:pink",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Old-style glitch-response plot for ASC ERR/CORR.")
    p.add_argument("--event-gps", type=float, required=True, help="Glitch peak GPS.")
    p.add_argument("--win", type=float, default=1.0, help="Half-window in seconds (default 1.0 => 2s total).")
    p.add_argument("--flo", type=float, default=30.0, help="Bandpass low cut (Hz).")
    p.add_argument("--fhi", type=float, default=55.0, help="Bandpass high cut (Hz).")
    p.add_argument("--env-smooth-ms", type=float, default=10.0, help="Envelope smoothing window (ms).")
    p.add_argument("--baseline-start", type=float, default=-1.0, help="Baseline start rel. to glitch (s).")
    p.add_argument("--baseline-end", type=float, default=-0.2, help="Baseline end rel. to glitch (s).")
    p.add_argument("--csv", default=None, help="Optional CSV path; default outputs/asc_err_corr_recovery_<event>.csv")
    p.add_argument("--out", default=None, help="Optional output PNG path.")
    return p.parse_args()


def short(ch: str) -> str:
    return ch.replace("V1:ASC_", "")


def infer_fs(t: np.ndarray) -> float:
    if t.size < 3:
        return 1.0
    dt = np.median(np.diff(t))
    if not np.isfinite(dt) or dt <= 0:
        return 1.0
    return float(1.0 / dt)


def bandpass(v: np.ndarray, fs: float, flo: float, fhi: float, order: int = 4) -> np.ndarray:
    if fs <= 2.2 * fhi or v.size < 8:
        return v.copy()
    sos = butter(order, [flo, fhi], btype="bandpass", fs=fs, output="sos")
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
        m = (t_rel >= -0.5) & (t_rel <= -0.1)
    if not np.any(m):
        m = t_rel < 0
    if not np.any(m):
        m = np.ones_like(t_rel, dtype=bool)
    return m


def robust_top(y: np.ndarray, floor: float = 2.0) -> float:
    if y.size == 0:
        return floor
    z = y[np.isfinite(y)]
    if z.size == 0:
        return floor
    return max(floor, float(np.percentile(z, 99.5) * 1.15))


def load_channel(df: pd.DataFrame, ch: str) -> tuple[np.ndarray, np.ndarray]:
    d = df[df["channel"] == ch].sort_values("gps")
    t = d["gps"].to_numpy(dtype=float)
    v = d["value"].to_numpy(dtype=float)
    return t, v


def plot_response(args: argparse.Namespace) -> Path:
    event = float(args.event_gps)
    event_id = int(event)
    csv_path = Path(args.csv) if args.csv else OUTPUTS / f"asc_err_corr_recovery_{event_id}.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)

    if REF_CHANNEL not in set(df["channel"]):
        raise RuntimeError(f"Missing {REF_CHANNEL} in {csv_path.name}")

    win = float(args.win)
    use_ms = win <= 0.5
    scale = 1e3 if use_ms else 1.0
    unit = "ms" if use_ms else "s"
    suffix = f"{int(win*2*1000)}ms" if use_ms else f"{int(win*2)}s"

    t_ref, v_ref = load_channel(df, REF_CHANNEL)
    t_rel_ref = t_ref - event
    fs_ref = infer_fs(t_ref)
    ref_bp = bandpass(v_ref, fs_ref, args.flo, args.fhi)
    ref_env = envelope(ref_bp, fs_ref, args.env_smooth_ms)
    bl_ref = baseline_mask(t_rel_ref, args.baseline_start, args.baseline_end)
    ref_rms = np.std(ref_bp[bl_ref]) if np.any(bl_ref) else np.std(ref_bp)
    if not np.isfinite(ref_rms) or ref_rms == 0:
        ref_rms = 1.0
    ref_bp_norm = ref_bp / ref_rms
    ref_env_bl = np.mean(ref_env[bl_ref]) if np.any(bl_ref) else np.mean(ref_env)
    if not np.isfinite(ref_env_bl) or ref_env_bl <= 0:
        ref_env_bl = 1.0
    ref_env_norm = ref_env / ref_env_bl

    tower_traces: list[tuple[str, list[tuple[str, str, str, np.ndarray, np.ndarray]]]] = []
    dead_channels: list[str] = []
    present = set(df["channel"])
    for tower, channels in TOWERS:
        traces = []
        for ch in channels:
            if ch not in present:
                continue
            t, v = load_channel(df, ch)
            if np.nanstd(v) < 1e-30:
                dead_channels.append(ch)
                continue
            t_rel = t - event
            fs = infer_fs(t)
            bp = bandpass(v, fs, args.flo, args.fhi)
            env = envelope(bp, fs, args.env_smooth_ms)
            bl = baseline_mask(t_rel, args.baseline_start, args.baseline_end)
            bl_mean = np.mean(env[bl]) if np.any(bl) else np.mean(env)
            if not np.isfinite(bl_mean) or bl_mean <= 0:
                bl_mean = 1.0
            env_norm = env / bl_mean
            ls = "--" if ch.endswith("_ERR") else "-"
            traces.append((ch, TOWER_COLORS[tower], ls, t_rel, env_norm))
        tower_traces.append((tower, traces))

    tower_traces = [(tw, tr) for tw, tr in tower_traces if tr]
    n_rows = 2 + len(tower_traces)
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 3 * n_rows), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.13, right=0.97)

    mask_ref = (t_rel_ref >= -win) & (t_rel_ref <= win)
    tplot_ref = t_rel_ref[mask_ref] * scale

    ax = axes[0]
    ax.plot(
        tplot_ref,
        ref_bp_norm[mask_ref],
        color="black",
        lw=0.8,
        alpha=0.9,
        label=f"{REF_CHANNEL} ({args.flo:.0f}-{args.fhi:.0f} Hz)",
    )
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax.set_ylabel("Amplitude\n(pre-glitch RMS)", fontsize=9)
    ax.set_title(f"{REF_CHANNEL} — bandpassed waveform", fontsize=9, loc="left")
    yabs = np.percentile(np.abs(ref_bp_norm[mask_ref]), 99.5) if np.any(mask_ref) else 2.0
    ax.set_ylim(-max(2.0, 1.2 * yabs), max(2.0, 1.2 * yabs))
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="x", ls=":", alpha=0.4)

    ax = axes[1]
    ax.plot(tplot_ref, ref_env_norm[mask_ref], color="black", lw=1.5, alpha=0.9, label=f"{REF_CHANNEL} envelope")
    ax.axhline(1.0, color="black", lw=0.8, ls=":", alpha=0.5, label="baseline = 1")
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
    ax.set_title(f"{REF_CHANNEL} — Hilbert envelope", fontsize=9, loc="left")
    ax.set_ylim(0.0, robust_top(ref_env_norm[mask_ref], floor=2.0))
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="x", ls=":", alpha=0.4)

    for i, (tower, traces) in enumerate(tower_traces):
        ax = axes[2 + i]
        panel_vals = []
        tower_mean = None
        for ch, color, ls, t_rel, env_norm in traces:
            m = (t_rel >= -win) & (t_rel <= win)
            if not np.any(m):
                continue
            x = t_rel[m] * scale
            y = env_norm[m]
            panel_vals.append(y)
            ax.plot(x, y, color=color, lw=1.2, ls=ls, alpha=0.9, label=short(ch))
            if tower_mean is None:
                tower_mean = y.copy()
            else:
                tower_mean += y
        if tower_mean is not None and panel_vals:
            tower_mean = tower_mean / len(panel_vals)
            ax.plot(x, tower_mean, color="black", lw=2.0, alpha=0.85, label=f"{tower} mean")
        ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
        ax.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
        ax.set_title(f"{tower} — ASC ERR/CORR Hilbert envelope", fontsize=9, loc="left")
        if panel_vals:
            ax.set_ylim(0.0, robust_top(np.concatenate(panel_vals), floor=2.0))
            ax.legend(fontsize=8, loc="upper right", ncol=2, handlelength=1.5)
        ax.grid(axis="x", ls=":", alpha=0.4)

    xtick_step = 0.2 if win >= 1.0 else (0.1 if win >= 0.5 else 50.0 / scale)
    xt = np.arange(-win * scale, win * scale + xtick_step * scale * 0.01, xtick_step * scale)
    for ax in axes:
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
    axes[-1].set_xlabel(f"Time relative to GPS {event:.6f} ({unit})", fontsize=9)

    utc = GPS_EPOCH + pd.to_timedelta(event, unit="s")
    fig.suptitle(
        f"ASC ERR/CORR response to glitch around Hrec\n"
        f"GPS {event:.6f} ({utc.strftime('%Y-%m-%d %H:%M:%S UTC')})   "
        f"Window: ±{win*1e3:.0f} ms   Bandpass {args.flo:.0f}-{args.fhi:.0f} Hz   "
        f"Envelope smoothing {args.env_smooth_ms:.0f} ms",
        fontsize=10,
        fontweight="bold",
    )
    if dead_channels:
        fig.text(
            0.13,
            0.012,
            "Skipped constant-zero channels: " + ", ".join(short(c) for c in dead_channels),
            fontsize=8,
            color="dimgray",
        )

    out = Path(args.out) if args.out else OUT_DIR / f"glitch_response_{event_id}_{suffix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    args = parse_args()
    out = plot_response(args)
    print(f"Saved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
