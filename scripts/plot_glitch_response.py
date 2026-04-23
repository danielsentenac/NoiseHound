#!/usr/bin/env python3
"""
Generate the 2-second glitch-response cross-check plot.

Default behaviour: plot the underlying bandpassed time series directly,
without any Hilbert/envelope step. Use ``--mode envelope`` only if the
older envelope-style view is explicitly needed.
"""
import argparse
import os
import warnings; warnings.filterwarnings("ignore")

from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, hilbert, sosfiltfilt


WORKDIR = Path(__file__).parent.parent
OUT_DIR_DEFAULT = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

DARM_COL = ("DARM corr", "V1:LSC_DARM_CORR", "black")

ASC_GROUPS = [
    ("NI correction outputs", [
        ("NI TX", "V1:ASC_NI_TX_CORR", "tab:blue"),
        ("NI TY", "V1:ASC_NI_TY_CORR", "steelblue"),
    ]),
    ("WI correction outputs", [
        ("WI TX", "V1:ASC_WI_TX_CORR", "tab:orange"),
        ("WI TY", "V1:ASC_WI_TY_CORR", "goldenrod"),
    ]),
    ("NE correction outputs", [
        ("NE TX", "V1:ASC_NE_TX_CORR", "tab:green"),
        ("NE TY", "V1:ASC_NE_TY_CORR", "limegreen"),
    ]),
    ("WE correction outputs", [
        ("WE TX", "V1:ASC_WE_TX_CORR", "tab:red"),
        ("WE TY", "V1:ASC_WE_TY_CORR", "salmon"),
    ]),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glitch-gps", type=int, default=1415579346,
                        help="Integer GPS used in the CSV/output filename.")
    parser.add_argument("--glitch-gps-exact", type=float, default=1415579346.0,
                        help="Reference GPS for t=0 on the x-axis.")
    parser.add_argument("--csv", type=Path, default=None,
                        help="Input CSV. Defaults to outputs/asc_glitch_probe_{glitch_gps}.csv")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT,
                        help="Directory for the rendered PNG.")
    parser.add_argument("--suffix", default="2s_exact",
                        help="Suffix used in glitch_response_{glitch_gps}_{suffix}.png")
    parser.add_argument("--window", type=float, default=1.0,
                        help="Half-width of the plotted time window in seconds.")
    parser.add_argument("--mode", choices=("waveform", "envelope"), default="waveform",
                        help="Plot bandpassed time series directly or the legacy Hilbert envelope.")
    parser.add_argument("--flo", type=float, default=30.0,
                        help="Bandpass low edge in Hz.")
    parser.add_argument("--fhi", type=float, default=55.0,
                        help="Bandpass high edge in Hz.")
    parser.add_argument("--order", type=int, default=4,
                        help="Butterworth bandpass order.")
    parser.add_argument("--baseline-start", type=float, default=-5.0,
                        help="Baseline window start relative to glitch GPS in seconds.")
    parser.add_argument("--baseline-end", type=float, default=-3.0,
                        help="Baseline window end relative to glitch GPS in seconds.")
    parser.add_argument("--env-smooth-ms", type=float, default=10.0,
                        help="Envelope smoothing window in ms when --mode envelope is used.")
    return parser.parse_args()


def default_csv_path(glitch_gps):
    return WORKDIR / "outputs" / f"asc_glitch_probe_{glitch_gps}.csv"


def load_probe_csv(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing probe CSV: {path}\n"
            f"Stage/extract it first, e.g. with slurm/nh_asc_glitch_probe.slurm "
            f"for the requested GLITCH_GPS."
        )
    df = pd.read_csv(path).sort_values("gps").set_index("gps")
    if df.empty:
        raise RuntimeError(f"Probe CSV is empty: {path}")
    return df


def get_fs(index):
    if len(index) < 2:
        return 1.0
    diffs = np.diff(index[: min(len(index), 200)])
    median_dt = np.median(diffs)
    if not np.isfinite(median_dt) or median_dt <= 0:
        return 1.0
    return round(1.0 / median_dt)


def bandpass(values, fs, flo, fhi, order):
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, values)


def smoothed_envelope(values, fs, smooth_ms):
    env = np.abs(hilbert(values))
    n_smooth = max(1, int(round(smooth_ms * 1e-3 * fs)))
    if n_smooth > 1:
        env = np.convolve(env, np.ones(n_smooth) / n_smooth, mode="same")
    return env


def baseline_mask(t_rel, start, end):
    return (t_rel >= start) & (t_rel <= end)


def safe_rms(values):
    rms = np.std(values)
    return rms if np.isfinite(rms) and rms > 0 else 1.0


def safe_mean(values):
    mean = np.mean(values)
    return mean if np.isfinite(mean) and mean > 0 else 1.0


def prepare_trace(df, channel, args):
    if channel not in df.columns:
        return None
    series = df[channel].dropna()
    if len(series) < max(8, args.order * 4):
        return None

    times = series.index.values.astype(float)
    values = series.values.astype(float)
    fs = get_fs(times)
    bp = bandpass(values, fs, args.flo, args.fhi, args.order)
    t_rel = times - args.glitch_gps_exact
    bl_mask = baseline_mask(t_rel, args.baseline_start, args.baseline_end)

    if args.mode == "envelope":
        y = smoothed_envelope(bp, fs, args.env_smooth_ms)
        denom = safe_mean(y[bl_mask]) if bl_mask.any() else safe_mean(y)
        y = y / denom
    else:
        denom = safe_rms(bp[bl_mask]) if bl_mask.any() else safe_rms(bp)
        y = bp / denom

    return {"t_rel": t_rel, "y": y, "fs": fs}


def window_config(window):
    use_ms = window <= 0.5
    if use_ms:
        scale = 1e3
        unit = "ms"
        if window <= 0.05:
            step = 5
        elif window <= 0.2:
            step = 10
        else:
            step = 50
        formatter = lambda value: f"{value:+.0f}"
    else:
        scale = 1.0
        unit = "s"
        if window <= 1.0:
            step = 0.2
        elif window <= 2.0:
            step = 0.5
        else:
            step = 1.0
        formatter = lambda value: f"{value:+.1f}"
    return use_ms, scale, unit, step, formatter


def plot_group(ax, title, traces, args, scale, step, formatter):
    any_trace = False
    for label, channel, color in traces:
        trace = prepare_trace(args.df, channel, args)
        if trace is None:
            continue
        mask = (trace["t_rel"] >= -args.window) & (trace["t_rel"] <= args.window)
        if mask.sum() < 2:
            continue
        x = trace["t_rel"][mask] * scale
        y = trace["y"][mask]
        ax.plot(x, y, color=color, lw=1.0, alpha=0.9, label=label)
        any_trace = True

    baseline = 1.0 if args.mode == "envelope" else 0.0
    ax.axhline(baseline, color="gray", lw=0.8, ls=":", alpha=0.5)
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.85, label="catalog GPS")
    ax.set_title(title, fontsize=9, loc="left")
    ax.grid(axis="x", ls=":", alpha=0.4)
    ax.tick_params(axis="y", labelsize=8)
    if any_trace:
        ax.legend(fontsize=8, loc="upper right", ncol=2,
                  handlelength=1.2, handletextpad=0.4, borderpad=0.4)
    else:
        ax.text(0.5, 0.5, "No data in CSV for this panel",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="dimgray")

    xticks = np.arange(-args.window * scale, args.window * scale + step * 0.01, step)
    ax.set_xticks(xticks)
    ax.set_xticklabels([formatter(value) for value in xticks], fontsize=8)


def main():
    args = parse_args()
    args.csv = args.csv or default_csv_path(args.glitch_gps)
    args.df = load_probe_csv(args.csv)

    use_ms, scale, unit, step, formatter = window_config(args.window)
    del use_ms

    # 6-panel layout: DARM waveform | DARM envelope | 4 x ASC group envelopes
    n_panels = 2 + len(ASC_GROUPS)
    fig, axes = plt.subplots(n_panels, 1, figsize=(16, 12), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.91, bottom=0.08, left=0.10, right=0.97)

    # Panel 1 — DARM bandpassed waveform
    waveform_args = type("A", (), vars(args) | {"mode": "waveform"})()
    plot_group(axes[0], "V1:LSC_DARM_CORR — bandpassed waveform",
               [DARM_COL], waveform_args, scale, step, formatter)
    axes[0].set_ylabel("Amplitude\n(x pre-glitch RMS)", fontsize=9)

    # Panel 2 — DARM Hilbert envelope
    envelope_args = type("A", (), vars(args) | {"mode": "envelope"})()
    plot_group(axes[1],
               "V1:LSC_DARM_CORR — Hilbert envelope  (1.0 = pre-glitch level)",
               [DARM_COL], envelope_args, scale, step, formatter)
    axes[1].set_ylabel("Envelope\n(x pre-glitch mean)", fontsize=9)

    # Panels 3-6 — ASC group envelopes
    for ax, (group_title, traces) in zip(axes[2:], ASC_GROUPS):
        plot_group(ax, f"{group_title} — ASC CORR Hilbert envelope",
                   traces, envelope_args, scale, step, formatter)
        ax.set_ylabel("Envelope\n(x pre-glitch mean)", fontsize=9)

    axes[-1].set_xlabel(
        f"Time relative to catalog GPS {args.glitch_gps_exact:.4f}  ({unit})",
        fontsize=9,
    )

    glitch_utc = GPS_EPOCH + pd.to_timedelta(args.glitch_gps_exact, unit="s")
    fig.suptitle(
        f"25-min glitch \u2014 did DARM glitch drive ASC correction responses?\n"
        f"GPS {args.glitch_gps_exact:.4f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})  "
        f"Window: \u00b1{args.window * scale:.0f} {unit}   "
        f"Bandpass {args.flo:.0f}\u2013{args.fhi:.0f} Hz   "
        f"Envelope smoothed {args.env_smooth_ms:.0f} ms",
        fontsize=10, fontweight="bold",
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / f"glitch_response_{args.glitch_gps}_{args.suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Loaded {args.df.shape[0]} samples x {args.df.shape[1]} columns from {args.csv}")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
