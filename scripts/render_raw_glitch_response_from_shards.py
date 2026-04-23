#!/usr/bin/env python3
"""
Render a raw glitch interval from pre-extracted per-file shard bundles.

The lower timeseries panels are built from shard data, while the q-transform
panel can still be computed directly from the raw GWF list.
"""
import argparse
import gzip
import os
import pickle
import time
import warnings; warnings.filterwarnings("ignore")

from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_raw_glitch_response import (
    RAW_CHANNELS,
    estimate_sample_rate,
    format_elapsed,
    gps_to_utc_string,
    log_progress,
    resolve_glitch_times,
    set_robust_ylim,
    thin_for_plot,
    transform_panel,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--glitch-gps", type=float, default=None,
                     help="Exact glitch GPS used as t = 0 (pair with --window).")
    grp.add_argument("--gps-start", type=float, default=None,
                     help="Start of GPS interval (pair with --gps-end).")
    parser.add_argument("--gps-end", type=float, default=None,
                        help="End of GPS interval (required when --gps-start is used).")
    parser.add_argument("--window", type=float, default=60.0,
                        help="Half-width of the plotted window in seconds (used with --glitch-gps).")
    parser.add_argument("--shard", nargs="+", type=Path, required=True,
                        help="One or more shard files produced by extract_raw_interval_shard.py.")
    parser.add_argument("--gwf", nargs="+", type=Path, default=None,
                        help="Optional raw GWF files used for the q-transform panel.")
    parser.add_argument("--q-channel", default="V1:Hrec_hoft_16384Hz",
                        help="High-rate channel for the q-transform.")
    parser.add_argument("--q-fmin", type=float, default=20.0,
                        help="Minimum frequency for the q-transform in Hz.")
    parser.add_argument("--q-fmax", type=float, default=80.0,
                        help="Maximum frequency for the q-transform in Hz.")
    parser.add_argument("--q-qmin", type=float, default=4.0,
                        help="Minimum Q for the q-transform.")
    parser.add_argument("--q-qmax", type=float, default=64.0,
                        help="Maximum Q for the q-transform.")
    parser.add_argument("--q-max-rate", type=float, default=4096.0,
                        help="Maximum sample rate used for the q-transform.")
    parser.add_argument("--panel-mode", choices=("raw", "uniform-bandpass", "mixed-env", "all-envelope"),
                        default="raw",
                        help="How to render the lower timeseries panels.")
    parser.add_argument("--bandpass-low", type=float, default=None,
                        help="Optional low cutoff frequency in Hz for the timeseries panels.")
    parser.add_argument("--bandpass-high", type=float, default=None,
                        help="Optional high cutoff frequency in Hz for the timeseries panels.")
    parser.add_argument("--bandpass-order", type=int, default=4,
                        help="Butterworth band-pass order for the timeseries panels.")
    parser.add_argument("--smooth-s", type=float, default=0.0,
                        help="Optional smoothing scale in seconds for direct waveform panels.")
    parser.add_argument("--hrec-bandpass-low", type=float, default=35.0,
                        help="Low cutoff in Hz for the Hrec waveform in mixed mode.")
    parser.add_argument("--hrec-bandpass-high", type=float, default=50.0,
                        help="High cutoff in Hz for the Hrec waveform in mixed mode.")
    parser.add_argument("--aux-bandpass-low", type=float, default=30.0,
                        help="Low cutoff in Hz for non-Hrec channels in mixed mode.")
    parser.add_argument("--aux-bandpass-high", type=float, default=60.0,
                        help="High cutoff in Hz for non-Hrec channels in mixed mode.")
    parser.add_argument("--aux-envelope-smooth-s", type=float, default=0.25,
                        help="Smoothing scale in seconds for envelope-excess panels.")
    parser.add_argument("--baseline-end-s", type=float, default=-10.0,
                        help="End of the pre-glitch baseline window in seconds.")
    parser.add_argument("--glitch-times", nargs="+", default=None,
                        help="GPS times of glitches to mark with vertical dashed lines.")
    parser.add_argument("--max-plot-points", type=int, default=5000,
                        help="Maximum points per channel kept for display.")
    parser.add_argument("--filter-guard", type=float, default=30.0,
                        help="Extra seconds loaded on each side for filter guard band (default 30s).")
    parser.add_argument("--channel-list", type=Path, default=None,
                        help="Text file with one 'CHANNEL COLOR' entry per line overriding the default channel set.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output PNG path.")
    args = parser.parse_args()
    if args.gps_start is not None:
        if args.gps_end is None:
            parser.error("--gps-end is required when --gps-start is used")
        if args.gps_end <= args.gps_start:
            parser.error("--gps-end must be greater than --gps-start")
        args.window = (args.gps_end - args.gps_start) / 2.0
        args.glitch_gps = args.gps_start + args.window
        args._interval_mode = True
    else:
        args._interval_mode = False
        args.gps_start = args.glitch_gps - args.window
        args.gps_end = args.glitch_gps + args.window
    # Extend load range by filter_guard on each side to avoid edge transients.
    # The display window (args.window / gps_start / gps_end) is kept unchanged.
    args._load_start = args.gps_start - args.filter_guard
    args._load_end = args.gps_end + args.filter_guard
    return args


def load_shards(shard_paths: list[Path], gps_start: float, gps_end: float):
    shard_paths = sorted(shard_paths)
    if not shard_paths:
        raise RuntimeError("No shard files provided")

    channel_chunks = {}
    for idx, path in enumerate(shard_paths, start=1):
        if not path.exists():
            raise FileNotFoundError(f"Missing shard file: {path}")
        started = time.perf_counter()
        log_progress(f"Loading shard {idx}/{len(shard_paths)}: {path}")
        with gzip.open(path, "rb") as fp:
            payload = pickle.load(fp)
        for channel, series in payload.get("channels", {}).items():
            times = np.asarray(series["times"], dtype=np.float64)
            values = np.asarray(series["values"], dtype=np.float32)
            mask = (times >= gps_start) & (times <= gps_end)
            if not np.any(mask):
                continue
            channel_chunks.setdefault(channel, []).append((times[mask], values[mask]))
        log_progress(
            f"Loaded shard {idx}/{len(shard_paths)} in "
            f"{format_elapsed(time.perf_counter() - started)}"
        )

    series_map = {}
    for channel, chunks in channel_chunks.items():
        times = np.concatenate([chunk[0] for chunk in chunks])
        values = np.concatenate([chunk[1] for chunk in chunks])
        order = np.argsort(times, kind="mergesort")
        times = times[order]
        values = values[order]
        if len(times) > 1:
            keep = np.ones(len(times), dtype=bool)
            keep[1:] = times[1:] > times[:-1]
            times = times[keep]
            values = values[keep]
        series_map[channel] = (times, values)

    if not series_map:
        raise RuntimeError("No populated channels were found in the shard files")

    log_progress(f"Built merged shard map with {len(series_map)} channels")
    return series_map


def compute_qtransform_from_series(series_map, q_channel: str, glitch_gps: float, window: float,
                                   fmin: float, fmax: float, qmin: float, qmax: float,
                                   q_max_rate: float):
    from gwpy.timeseries import TimeSeries

    if q_channel not in series_map:
        raise RuntimeError(f"Q-transform channel not present in shards: {q_channel}")

    times, values = series_map[q_channel]
    sample_rate = estimate_sample_rate(times)
    if sample_rate is None:
        raise RuntimeError(f"Could not estimate sample rate for q-transform channel {q_channel}")

    log_progress(
        f"Preparing q-transform from shard data for {q_channel} at estimated "
        f"{sample_rate:.1f} Hz over {len(values)} samples"
    )
    ts = TimeSeries(values, sample_rate=sample_rate, t0=times[0])
    if sample_rate > q_max_rate:
        log_progress(f"Resampling q-transform channel from {sample_rate:.1f} Hz to {q_max_rate:.1f} Hz")
        ts = ts.resample(q_max_rate)

    q_started = time.perf_counter()
    log_progress(
        f"Computing q-transform from shards for {2.0 * window:.0f}s window, "
        f"{fmin:.0f}-{fmax:.0f} Hz, Q {qmin:g}-{qmax:g}"
    )
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
    log_progress(
        f"Computed shard q-transform in {format_elapsed(time.perf_counter() - q_started)} "
        f"with shape {values.shape[0]}x{values.shape[1]}"
    )
    return times, freqs, values


def plot_q_panel_from_source(ax, args, series_map, glitch_times=()):
    if args.gwf:
        from plot_raw_glitch_response import compute_qtransform

        times, freqs, values = compute_qtransform(
            args.gwf,
            args.q_channel,
            args.glitch_gps,
            args.window,
            args.q_fmin,
            args.q_fmax,
            args.q_qmin,
            args.q_qmax,
            args.q_max_rate,
        )
    else:
        times, freqs, values = compute_qtransform_from_series(
            series_map,
            args.q_channel,
            args.glitch_gps,
            args.window,
            args.q_fmin,
            args.q_fmax,
            args.q_qmin,
            args.q_qmax,
            args.q_max_rate,
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
    if not getattr(args, "_interval_mode", False):
        ax.axvline(0, color="crimson", lw=1.3, ls="--", alpha=0.85)
    for gt in glitch_times:
        ax.axvline(gt - args.glitch_gps, color="crimson", lw=1.1, ls="--", alpha=0.85)
    ax.set_ylabel("Hz", fontsize=8)
    ax.set_title(f"{args.q_channel} q-transform", fontsize=8, loc="left")
    ax.grid(False)


def plot_channel_from_series(ax, channel: str, color: str, series_map, glitch_gps: float,
                             max_plot_points: int, args, glitch_times=()):
    if channel not in series_map:
        return False

    times, values = series_map[channel]
    if len(times) == 0:
        return False

    values, note = transform_panel(channel, times, values, glitch_gps, args)
    # Trim to display window after filtering to discard guard-band transients.
    display_mask = (times >= args.gps_start) & (times <= args.gps_end)
    times = times[display_mask]
    values = values[display_mask]
    if len(times) == 0:
        return False
    y_for_limits = values.copy()
    t_rel = times - glitch_gps
    t_rel, values = thin_for_plot(t_rel, values, max_plot_points)

    ax.plot(t_rel, values, color=color, lw=0.7, label=channel)
    if not getattr(args, "_interval_mode", False):
        ax.axvline(0, color="crimson", lw=1.1, ls="--", alpha=0.85)
    for gt in glitch_times:
        ax.axvline(gt - glitch_gps, color="crimson", lw=1.1, ls="--", alpha=0.85)
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
    if args.panel_mode in ("mixed-env", "all-envelope") and (
        args.panel_mode == "all-envelope" or channel != "V1:Hrec_hoft_16384Hz"
    ):
        set_robust_ylim(ax, y_for_limits, lower_q=1.0, upper_q=99.0)
    elif args.panel_mode == "uniform-bandpass":
        set_robust_ylim(ax, y_for_limits, lower_q=0.5, upper_q=99.5)
    if note:
        ax.text(0.01, 0.88, note, transform=ax.transAxes, fontsize=6.5, color="0.35")
    return True


def main():
    args = parse_args()
    render_started = time.perf_counter()
    log_progress(
        f"Starting shard render for GPS {args.gps_start:.0f}-{args.gps_end:.0f} "
        f"in {args.panel_mode} mode using {len(args.shard)} shard files"
    )

    series_map = load_shards(args.shard, args._load_start, args._load_end)
    glitch_times = resolve_glitch_times(args.glitch_times, args.gps_start, args.gps_end)
    log_progress(f"Loaded {len(glitch_times)} glitch markers in window")

    if args.channel_list is not None:
        channel_list = []
        for line in args.channel_list.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            channel_list.append((parts[0], parts[1] if len(parts) > 1 else "black"))
    else:
        channel_list = RAW_CHANNELS
    available = [(ch, col) for ch, col in channel_list if ch in series_map and len(series_map[ch][0]) > 0]
    if not available:
        raise RuntimeError("None of the requested raw channels are populated in the shard inputs")

    show_q = bool(args.gwf) or args.q_channel in series_map
    n_rows = len(available) + (1 if show_q else 0)
    height_ratios = ([2.4] if show_q else []) + [1.0] * len(available)
    fig_height = max(10, 2.2 + 0.90 * n_rows)
    q_note = " plus q-transform panel" if show_q else ""
    log_progress(f"Preparing figure with {len(available)} timeseries panels{q_note}")
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

    ts_axes = axes
    if show_q:
        plot_q_panel_from_source(axes[0], args, series_map, glitch_times=glitch_times)
        axes[0].tick_params(axis="x", labelbottom=False)
        ts_axes = axes[1:]

    for idx, (channel, color) in enumerate(available):
        log_progress(f"Plotting panel {idx + 1}/{len(available)}: {channel}")
        ax = ts_axes[idx]
        plotted = plot_channel_from_series(
            ax,
            channel,
            color,
            series_map,
            args.glitch_gps,
            args.max_plot_points,
            args,
            glitch_times=glitch_times,
        )
        if not plotted:
            ax.text(0.5, 0.5, f"No data for {channel}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="dimgray")

    total_span = 2.0 * args.window
    for step in (1, 2, 5, 10, 20, 30, 60, 120, 300, 600, 900, 1800):
        if total_span / step <= 25:
            tick_step = step
            break
    else:
        tick_step = int(round(total_span / 20.0 / 100.0) * 100) or 600
    xticks = np.arange(-args.window, args.window + tick_step, tick_step)
    for ax in axes:
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{v:+.0f}" for v in xticks], fontsize=7)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", labelbottom=False)
    if args._interval_mode:
        axes[-1].set_xlabel(
            f"Time since GPS {args.gps_start:.0f} (s)  "
            f"[{args.gps_start:.0f} – {args.gps_end:.0f}]",
            fontsize=9,
        )
    else:
        axes[-1].set_xlabel(
            f"Time relative to GPS {args.glitch_gps:.6f} (s)", fontsize=9
        )

    mode_note = ""
    if args.panel_mode == "uniform-bandpass" and args.bandpass_low is not None and args.bandpass_high is not None:
        mode_note = f"  |  band-pass {args.bandpass_low:.0f}-{args.bandpass_high:.0f} Hz"
        if args.smooth_s > 0:
            mode_note += f"  |  smooth {int(round(args.smooth_s * 1000))} ms"
    elif args.panel_mode == "mixed-env":
        mode_note = (
            f"  |  Hrec {args.hrec_bandpass_low:.0f}-{args.hrec_bandpass_high:.0f} Hz"
            f"  |  others {args.aux_bandpass_low:.0f}-{args.aux_bandpass_high:.0f} Hz env excess"
        )
    elif args.panel_mode == "all-envelope":
        mode_note = (
            f"  |  all channels {args.aux_bandpass_low:.0f}-{args.aux_bandpass_high:.0f} Hz"
            f" env SNR  |  smooth {int(round(args.aux_envelope_smooth_s * 1000))} ms"
        )
    if args._interval_mode:
        start_utc = gps_to_utc_string(args.gps_start)
        end_utc = gps_to_utc_string(args.gps_end)
        top_note = (
            f"GPS {args.gps_start:.0f} – {args.gps_end:.0f}"
            f"  ({start_utc} … {end_utc})"
            f"  |  {total_span:.0f} s  |  raw data from shards{mode_note}"
        )
    else:
        glitch_utc = gps_to_utc_string(args.glitch_gps)
        top_note = (
            f"GPS {args.glitch_gps:.6f} ({glitch_utc})"
            f"  |  Window +/-{args.window:.0f} s  |  raw data from shards{mode_note}"
        )
    fig.suptitle(top_note, fontsize=9, fontweight="bold")

    log_progress(f"Saving figure to {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log_progress(f"Saved -> {args.out}")
    log_progress(
        f"Rendered {len(available)} channels from {len(args.shard)} shards in "
        f"{format_elapsed(time.perf_counter() - render_started)}"
    )


if __name__ == "__main__":
    main()
