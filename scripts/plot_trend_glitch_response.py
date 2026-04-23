#!/usr/bin/env python3
"""
Plot a single 25-minute glitch from trend data.

Default layout:
- optional q-transform panel on top from a high-rate raw channel
- one subplot per trend channel
- raw values or MAD-normalized excess for the trend channels
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


WORKDIR = Path(__file__).parent.parent
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

TREND_CHANNELS = [
    ("V1:LSC_DARM_CORR_mean", "black"),
    ("V1:ASC_SR_TX_ERR_mean", "tab:red"),
    ("V1:ASC_SR_TY_ERR_mean", "salmon"),
    ("V1:ASC_PR_TX_ERR_mean", "tab:orange"),
    ("V1:ASC_PR_TY_ERR_mean", "gold"),
    ("V1:ASC_BS_TX_ERR_mean", "tab:blue"),
    ("V1:ASC_BS_TY_ERR_mean", "tab:cyan"),
    ("V1:ASC_NI_TX_LF_B7_I_10Hz_mean", "tab:blue"),
    ("V1:ASC_NI_TY_LF_B7_I_10Hz_mean", "steelblue"),
    ("V1:ASC_WI_TX_LF_B8_I_10Hz_mean", "tab:orange"),
    ("V1:ASC_WI_TY_LF_B8_I_10Hz_mean", "goldenrod"),
    ("V1:ASC_NE_TX_CORR_mean", "tab:green"),
    ("V1:ASC_NE_TY_CORR_mean", "limegreen"),
    ("V1:ASC_WE_TX_CORR_mean", "tab:red"),
    ("V1:ASC_WE_TY_CORR_mean", "lightsalmon"),
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--glitch-gps", type=float, required=True,
                        help="Exact glitch GPS used as t = 0.")
    parser.add_argument("--window", type=float, default=30.0,
                        help="Half-width of the plotted window in seconds.")
    parser.add_argument("--trend-gwf", type=Path, default=None,
                        help="Local daily trend GWF to read directly.")
    parser.add_argument("--csv", type=Path, default=None,
                        help="Pre-extracted trend CSV with gps + channel columns.")
    parser.add_argument("--save-csv", type=Path, default=None,
                        help="Optional path to save the extracted trend window as CSV.")
    parser.add_argument("--trend-mode", choices=("raw", "mad-excess"), default="raw",
                        help="How to render the trend channels.")
    parser.add_argument("--q-gwf", type=Path, default=None,
                        help="Optional high-rate GWF for the top q-transform panel.")
    parser.add_argument("--q-channel", default="V1:LSC_DARM_CORR",
                        help="High-rate channel for the q-transform.")
    parser.add_argument("--q-fmin", type=float, default=20.0,
                        help="Minimum frequency for the q-transform in Hz.")
    parser.add_argument("--q-fmax", type=float, default=80.0,
                        help="Maximum frequency for the q-transform in Hz.")
    parser.add_argument("--q-qmin", type=float, default=4.0,
                        help="Minimum Q for the q-transform.")
    parser.add_argument("--q-qmax", type=float, default=64.0,
                        help="Maximum Q for the q-transform.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output PNG path.")
    return parser.parse_args()


def gps_to_utc_string(gps: float) -> str:
    try:
        from astropy.time import Time
        return Time(gps, format="gps", scale="utc").strftime("%Y-%m-%d %H:%M:%S UTC")
    except Exception:
        return (GPS_EPOCH + pd.to_timedelta(gps, unit="s")).strftime("%Y-%m-%d %H:%M:%S UTC")


def read_trend_gwf_window(gwf_path: Path, glitch_gps: float, window: float) -> pd.DataFrame:
    from gwpy.timeseries import TimeSeries

    t0 = glitch_gps - window
    t1 = glitch_gps + window
    frames = []

    for channel, _color in TREND_CHANNELS:
        try:
            ts = TimeSeries.read(str(gwf_path), channel, start=t0, end=t1)
        except Exception:
            continue
        frames.append(
            pd.DataFrame(
                {"gps": np.asarray(ts.times.value, float), channel: np.asarray(ts.value, float)}
            ).set_index("gps")
        )

    if not frames:
        raise RuntimeError(f"No requested trend channels could be read from {gwf_path}")

    merged = pd.concat(frames, axis=1, join="outer").sort_index()
    merged.index.name = "gps"
    return merged.reset_index()


def compute_qtransform(q_gwf: Path, q_channel: str, glitch_gps: float, window: float,
                       fmin: float, fmax: float, qmin: float, qmax: float):
    from gwpy.timeseries import TimeSeries

    pad = 2.0
    t0 = glitch_gps - window - pad
    t1 = glitch_gps + window + pad
    ts = TimeSeries.read(str(q_gwf), q_channel, start=t0, end=t1)
    qgram = ts.q_transform(
        outseg=(glitch_gps - window, glitch_gps + window),
        frange=(fmin, fmax),
        qrange=(qmin, qmax),
        whiten=True,
        logf=False,
    )

    times = np.asarray(qgram.xindex.value, dtype=float) - glitch_gps
    freqs = np.asarray(qgram.yindex.value, dtype=float)
    values = np.asarray(qgram.value, dtype=float)
    if values.shape == (len(times), len(freqs)):
        values = values.T
    if values.shape != (len(freqs), len(times)):
        raise RuntimeError(
            f"Unexpected q-transform shape {values.shape} for {len(freqs)} freqs x {len(times)} times"
        )
    return times, freqs, values


def plot_q_panel(ax, args):
    times, freqs, values = compute_qtransform(
        args.q_gwf,
        args.q_channel,
        args.glitch_gps,
        args.window,
        args.q_fmin,
        args.q_fmax,
        args.q_qmin,
        args.q_qmax,
    )

    finite = values[np.isfinite(values)]
    if finite.size:
        vmin = np.nanpercentile(finite, 15)
        vmax = np.nanpercentile(finite, 99.5)
    else:
        vmin, vmax = None, None

    ax.imshow(
        values,
        aspect="auto",
        origin="lower",
        extent=[times[0], times[-1], freqs[0], freqs[-1]],
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    ax.axvline(0, color="crimson", lw=1.3, ls="--", alpha=0.85)
    ax.set_ylabel("Hz", fontsize=8)
    ax.set_title(f"{args.q_channel} q-transform", fontsize=8, loc="left")
    ax.grid(False)


def mad_normalized_excess(values: np.ndarray) -> np.ndarray:
    center = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - center))
    if not np.isfinite(mad) or mad <= 0:
        spread = np.nanstd(values)
        if not np.isfinite(spread) or spread <= 0:
            return np.zeros_like(values)
        return (values - center) / spread
    return (values - center) / mad


def transform_trend_values(values: np.ndarray, trend_mode: str) -> np.ndarray:
    if trend_mode == "mad-excess":
        return mad_normalized_excess(values)
    return values


def plot_trend_channel(
    ax,
    channel: str,
    color: str,
    df: pd.DataFrame,
    glitch_gps: float,
    trend_mode: str,
    y_limits=None,
):
    if channel not in df.columns:
        return False

    sub = df[["gps", channel]].dropna().copy()
    if sub.empty:
        return False

    t_rel = sub["gps"].to_numpy(dtype=float) - glitch_gps
    values = transform_trend_values(sub[channel].to_numpy(dtype=float), trend_mode)

    ax.plot(
        t_rel,
        values,
        color=color,
        lw=0.9,
        marker="o",
        ms=2.2,
        label=channel,
    )
    ax.axvline(0, color="crimson", lw=1.1, ls="--", alpha=0.85)
    if trend_mode == "mad-excess":
        ax.axhline(0, color="0.45", lw=0.8, ls=":", alpha=0.7)
        if y_limits is not None:
            ax.set_ylim(*y_limits)
    ax.grid(axis="x", ls=":", alpha=0.35)
    ax.tick_params(axis="both", labelsize=7)
    ax.legend(
        fontsize=6.8,
        loc="upper right",
        framealpha=0.8,
        handlelength=1.3,
        borderpad=0.25,
        handletextpad=0.4,
    )
    ax.margins(x=0.01)
    return True


def main():
    args = parse_args()
    if not args.csv and not args.trend_gwf:
        raise SystemExit("Provide either --csv or --trend-gwf")

    if args.csv:
        df = pd.read_csv(args.csv).sort_values("gps")
    else:
        df = read_trend_gwf_window(args.trend_gwf, args.glitch_gps, args.window)
        if args.save_csv:
            args.save_csv.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(args.save_csv, index=False)

    mask = (df["gps"] >= args.glitch_gps - args.window) & (df["gps"] <= args.glitch_gps + args.window)
    df = df.loc[mask].copy().sort_values("gps")
    if df.empty:
        raise RuntimeError("No data in the requested time window")

    available = [(channel, color) for channel, color in TREND_CHANNELS if channel in df.columns and df[channel].notna().any()]
    if not available:
        raise RuntimeError("None of the requested trend channels are populated in the selected window")

    common_trend_limits = None
    if args.trend_mode == "mad-excess":
        transformed = []
        for channel, _color in available:
            values = df[channel].dropna().to_numpy(dtype=float)
            if values.size:
                transformed.append(transform_trend_values(values, args.trend_mode))
        finite_chunks = [vals[np.isfinite(vals)] for vals in transformed if np.isfinite(vals).any()]
        if finite_chunks:
            finite = np.concatenate(finite_chunks)
            limit = max(3.0, float(np.nanpercentile(np.abs(finite), 99.0)))
            limit = np.ceil(limit * 2.0) / 2.0
        else:
            limit = 3.0
        common_trend_limits = (-limit, limit)

    n_rows = len(available) + (1 if args.q_gwf else 0)
    height_ratios = ([2.4] if args.q_gwf else []) + [1.0] * len(available)
    fig_height = max(10, 2.2 + 0.95 * n_rows)
    fig, axes = plt.subplots(
        n_rows,
        1,
        figsize=(16, fig_height),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios},
    )
    if n_rows == 1:
        axes = [axes]
    else:
        axes = list(axes)
    fig.subplots_adjust(hspace=0.06, top=0.95, bottom=0.06, left=0.10, right=0.98)

    ax_index = 0
    if args.q_gwf:
        plot_q_panel(axes[ax_index], args)
        axes[ax_index].tick_params(axis="x", labelbottom=False)
        ax_index += 1

    for idx, (channel, color) in enumerate(available):
        ax = axes[ax_index + idx]
        plotted = plot_trend_channel(
            ax,
            channel,
            color,
            df,
            args.glitch_gps,
            args.trend_mode,
            y_limits=common_trend_limits,
        )
        if not plotted:
            ax.text(0.5, 0.5, f"No data for {channel}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="dimgray")

    xticks = np.arange(-args.window, args.window + 5.0, 5.0)
    for ax in axes:
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{value:+.0f}" for value in xticks], fontsize=7)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False)
    axes[-1].set_xlabel(f"Time relative to glitch GPS {args.glitch_gps:.6f} (s)", fontsize=9)

    glitch_utc = gps_to_utc_string(args.glitch_gps)
    top_note = f"GPS {args.glitch_gps:.6f} ({glitch_utc})  |  Window +/-{args.window:.0f} s"
    if args.trend_mode == "mad-excess":
        top_note += "  |  trend: (x - median) / MAD"
    fig.suptitle(top_note, fontsize=10, fontweight="bold")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")
    print(f"Rows: {len(df)}  Columns: {len(df.columns)}")
    if args.save_csv:
        print(f"Saved CSV -> {args.save_csv}")


if __name__ == "__main__":
    main()
