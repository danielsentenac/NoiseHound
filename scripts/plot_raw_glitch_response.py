#!/usr/bin/env python3
"""
Plot a raw glitch window with one subplot per channel.

Default channels:
- V1:Hrec_hoft_16384Hz
- the same ASC/LSC/LF channels as the trend-response figure,
  but using their raw channel names (no trend suffix)
"""
import argparse
import os
import time
import warnings; warnings.filterwarnings("ignore")

from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKDIR = Path(__file__).parent.parent
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

RAW_CHANNELS = [
    ("V1:Hrec_hoft_16384Hz", "black"),
    ("V1:LSC_DARM_ERR", "dimgray"),
    ("V1:LSC_DARM_CORR", "tab:blue"),
    ("V1:LSC_MICH_CORR", "tab:purple"),
    ("V1:LSC_PRCL_CORR", "mediumorchid"),
    ("V1:LSC_SRCL_CORR", "plum"),
    ("V1:LSC_CARM_CORR", "indigo"),
    ("V1:LSC_BS_CORR", "tab:brown"),
    ("V1:LSC_PR_CORR", "sienna"),
    ("V1:LSC_SR_CORR", "peru"),
    ("V1:LSC_NI_CORR", "tab:olive"),
    ("V1:LSC_WI_CORR", "yellowgreen"),
    ("V1:LSC_NE_CORR", "tab:green"),
    ("V1:LSC_WE_CORR", "limegreen"),
]


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, seconds = divmod(seconds, 60.0)
    if minutes < 60:
        return f"{int(minutes)}m {seconds:.1f}s"
    hours, minutes = divmod(minutes, 60.0)
    return f"{int(hours)}h {int(minutes)}m {seconds:.0f}s"


def log_progress(message: str) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{stamp}] {message}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    # Time window: either glitch-gps + window  OR  gps-start + gps-end
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--glitch-gps", type=float, default=None,
                     help="Exact glitch GPS used as t = 0 (pair with --window).")
    grp.add_argument("--gps-start", type=float, default=None,
                     help="Start of GPS interval (pair with --gps-end).")
    parser.add_argument("--gps-end", type=float, default=None,
                        help="End of GPS interval (required when --gps-start is used).")
    parser.add_argument("--window", type=float, default=60.0,
                        help="Half-width of the plotted window in seconds (used with --glitch-gps).")
    parser.add_argument("--gwf", nargs="+", type=Path, default=None,
                        help="One or more local raw GWF files.")
    parser.add_argument("--csv", type=Path, default=None,
                        help="Optional pre-extracted raw CSV with gps + channel columns.")
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
                        help="GPS times of glitches to mark with vertical dashed lines. "
                             "Pass either explicit GPS values or a single CSV file path "
                             "(column 'gps_time' or 'peak_time' or first numeric column).")
    parser.add_argument("--save-csv", type=Path, default=None,
                        help="Optional path to save the extracted window as CSV.")
    parser.add_argument("--max-plot-points", type=int, default=5000,
                        help="Maximum points per channel kept for display.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output PNG path.")
    args = parser.parse_args()
    # Resolve glitch_gps / window from interval if --gps-start used
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
    return args


def resolve_glitch_times(glitch_times_arg, gps_start: float, gps_end: float):
    """Return list of GPS floats within [gps_start, gps_end] from --glitch-times."""
    if not glitch_times_arg:
        return []
    # Single argument that looks like a file path
    if len(glitch_times_arg) == 1:
        p = Path(glitch_times_arg[0])
        if p.exists():
            df = pd.read_csv(p)
            for col in ("gps_time", "peak_time", "gps"):
                if col in df.columns:
                    times = df[col].dropna().astype(float).tolist()
                    break
            else:
                # fall back to first numeric column
                numeric = df.select_dtypes(include=[float, int])
                if numeric.empty:
                    return []
                times = numeric.iloc[:, 0].dropna().astype(float).tolist()
            return [t for t in times if gps_start <= t <= gps_end]
    # Otherwise treat all values as GPS floats
    times = []
    for v in glitch_times_arg:
        try:
            times.append(float(v))
        except ValueError:
            pass
    return [t for t in times if gps_start <= t <= gps_end]


def gps_to_utc_string(gps: float) -> str:
    try:
        from astropy.time import Time
        return Time(gps, format="gps", scale="utc").strftime("%Y-%m-%d %H:%M:%S UTC")
    except Exception:
        return (GPS_EPOCH + pd.to_timedelta(gps, unit="s")).strftime("%Y-%m-%d %H:%M:%S UTC")


def thin_for_plot(times: np.ndarray, values: np.ndarray, max_points: int):
    if len(times) <= max_points:
        return times, values
    step = max(1, int(np.ceil(len(times) / max_points)))
    trim = (len(times) // step) * step
    if trim == 0:
        return times[:1], values[:1]
    abs_chunks = np.abs(values[:trim]).reshape(-1, step)
    peak_offsets = np.argmax(abs_chunks, axis=1)
    chunk_starts = np.arange(len(peak_offsets)) * step
    idx = chunk_starts + peak_offsets
    return times[idx], values[idx]


def estimate_sample_rate(times: np.ndarray):
    if times.size < 3:
        return None
    dt = np.diff(times)
    dt = dt[np.isfinite(dt) & (dt > 0)]
    if dt.size == 0:
        return None
    median_dt = float(np.median(dt))
    if median_dt <= 0:
        return None
    return 1.0 / median_dt


def maybe_bandpass(times: np.ndarray, values: np.ndarray, low, high, order: int):
    if low is None or high is None:
        return values, False, estimate_sample_rate(times)

    sample_rate = estimate_sample_rate(times)
    if sample_rate is None:
        return values, False, None

    nyquist = 0.5 * sample_rate
    if low <= 0 or high <= low or high >= nyquist:
        return values, False, sample_rate

    from scipy.signal import butter, sosfiltfilt

    sos = butter(order, [low, high], btype="bandpass", fs=sample_rate, output="sos")
    return sosfiltfilt(sos, values), True, sample_rate


def smooth_series(values: np.ndarray, sample_rate, smooth_s: float) -> np.ndarray:
    if sample_rate is None or smooth_s <= 0:
        return values
    width = max(1, int(round(sample_rate * smooth_s)))
    if width <= 1:
        return values
    kernel = np.ones(width, dtype=float) / width
    return np.convolve(values, kernel, mode="same")


def set_robust_ylim(ax, values: np.ndarray, lower_q: float = 1.0, upper_q: float = 99.0):
    finite = values[np.isfinite(values)]
    if finite.size < 10:
        return

    lo = float(np.nanpercentile(finite, lower_q))
    hi = float(np.nanpercentile(finite, upper_q))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return

    if hi <= lo:
        center = float(np.nanmedian(finite))
        spread = float(np.nanstd(finite))
        if not np.isfinite(spread) or spread <= 0:
            spread = max(abs(center) * 0.1, 1.0)
        ax.set_ylim(center - spread, center + spread)
        return

    pad = 0.08 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)


def transform_panel(channel: str, times: np.ndarray, values: np.ndarray,
                    glitch_gps: float, args):
    t_rel = times - glitch_gps

    if args.panel_mode == "uniform-bandpass":
        out, filtered, sample_rate = maybe_bandpass(
            times, values, args.bandpass_low, args.bandpass_high, args.bandpass_order
        )
        note = None
        if filtered and args.bandpass_low is not None and args.bandpass_high is not None:
            note = f"{args.bandpass_low:.0f}-{args.bandpass_high:.0f} Hz"
        if args.smooth_s > 0:
            out = smooth_series(out, sample_rate, args.smooth_s)
            smooth_ms = int(round(args.smooth_s * 1000))
            note = f"{note} + {smooth_ms} ms smooth" if note else f"{smooth_ms} ms smooth"
        return out, note

    if args.panel_mode == "mixed-env":
        if channel == "V1:Hrec_hoft_16384Hz":
            out, filtered, _fs = maybe_bandpass(
                times, values, args.hrec_bandpass_low, args.hrec_bandpass_high, args.bandpass_order
            )
            note = "waveform"
            if filtered:
                note = f"{args.hrec_bandpass_low:.0f}-{args.hrec_bandpass_high:.0f} Hz waveform"
            return out, note

        bandpassed, filtered, sample_rate = maybe_bandpass(
            times, values, args.aux_bandpass_low, args.aux_bandpass_high, args.bandpass_order
        )
        if not filtered:
            return values, "raw"

        from scipy.signal import hilbert

        envelope = np.abs(hilbert(bandpassed))
        envelope = smooth_series(envelope, sample_rate, args.aux_envelope_smooth_s)
        baseline_mask = t_rel <= args.baseline_end_s
        baseline = envelope[baseline_mask]
        baseline = baseline[np.isfinite(baseline)]
        if baseline.size == 0:
            baseline = envelope[np.isfinite(envelope)]
        baseline_median = float(np.nanmedian(baseline)) if baseline.size else 0.0
        baseline_std = float(np.nanstd(baseline)) if baseline.size else 1.0
        if not np.isfinite(baseline_std) or baseline_std <= 0:
            baseline_std = 1.0
        return (envelope - baseline_median) / baseline_std, (
            f"{args.aux_bandpass_low:.0f}-{args.aux_bandpass_high:.0f} Hz env SNR"
        )

    if args.panel_mode == "all-envelope":
        bandpassed, filtered, sample_rate = maybe_bandpass(
            times, values, args.aux_bandpass_low, args.aux_bandpass_high, args.bandpass_order
        )
        if not filtered:
            return values, "raw"
        from scipy.signal import hilbert
        envelope = np.abs(hilbert(bandpassed))
        envelope = smooth_series(envelope, sample_rate, args.aux_envelope_smooth_s)
        baseline_mask = t_rel <= args.baseline_end_s
        baseline = envelope[baseline_mask]
        baseline = baseline[np.isfinite(baseline)]
        if baseline.size == 0:
            baseline = envelope[np.isfinite(envelope)]
        baseline_median = float(np.nanmedian(baseline)) if baseline.size else 0.0
        baseline_std = float(np.nanstd(baseline)) if baseline.size else 1.0
        if not np.isfinite(baseline_std) or baseline_std <= 0:
            baseline_std = 1.0
        return (envelope - baseline_median) / baseline_std, (
            f"{args.aux_bandpass_low:.0f}-{args.aux_bandpass_high:.0f} Hz env SNR"
        )

    return values, None


def compute_qtransform(gwf_paths, q_channel: str, glitch_gps: float, window: float,
                       fmin: float, fmax: float, qmin: float, qmax: float, q_max_rate: float):
    from gwpy.timeseries import TimeSeries

    pad = 2.0
    t0 = glitch_gps - window - pad
    t1 = glitch_gps + window + pad
    read_started = time.perf_counter()
    log_progress(
        f"Reading q-transform channel {q_channel} from {len(gwf_paths)} GWF files "
        f"over GPS {t0:.0f}-{t1:.0f}"
    )
    ts = TimeSeries.read([str(p) for p in gwf_paths], q_channel, start=t0, end=t1)
    log_progress(f"Loaded q-transform channel in {format_elapsed(time.perf_counter() - read_started)}")
    try:
        sample_rate = float(ts.sample_rate.value)
    except Exception:
        sample_rate = None
    if sample_rate and sample_rate > q_max_rate:
        log_progress(f"Resampling q-transform channel from {sample_rate:.1f} Hz to {q_max_rate:.1f} Hz")
        ts = ts.resample(q_max_rate)
        sample_rate = q_max_rate
    q_started = time.perf_counter()
    log_progress(
        f"Computing q-transform for {2.0 * window:.0f}s window, "
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
        f"Computed q-transform in {format_elapsed(time.perf_counter() - q_started)} "
        f"with shape {values.shape[0]}x{values.shape[1]}"
    )
    return times, freqs, values


def plot_q_panel(ax, args, glitch_times=()):
    if not args.gwf:
        raise RuntimeError("Q-transform requires --gwf")

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
    for gt in glitch_times:
        ax.axvline(gt - args.glitch_gps, color="white", lw=0.9, ls="--", alpha=0.7)
    ax.set_ylabel("Hz", fontsize=8)
    ax.set_title(f"{args.q_channel} q-transform", fontsize=8, loc="left")
    ax.grid(False)


def extract_raw_window(gwf_paths, glitch_gps: float, window: float):
    from gwpy.io.gwf import iter_channel_names
    from gwpy.timeseries import TimeSeries

    t0 = glitch_gps - window
    t1 = glitch_gps + window
    gwf_paths = [str(path) for path in gwf_paths]
    channel_names = set(iter_channel_names(gwf_paths[0]))

    selected = []
    for channel, _color in RAW_CHANNELS:
        if channel not in channel_names:
            base = channel.split(":", 1)[1]
            close = sorted(ch for ch in channel_names if base.upper() in ch.upper())
            if close:
                print(f"Missing {channel}; close matches: {', '.join(close[:5])}")
            else:
                print(f"Missing {channel}")
            continue
        selected.append(channel)

    if not selected:
        raise RuntimeError("None of the requested raw channels could be read")

    read_started = time.perf_counter()
    log_progress(
        f"Reading {len(selected)} raw channels from {len(gwf_paths)} GWF files "
        f"over GPS {t0:.0f}-{t1:.0f}"
    )

    rows = {}
    for idx, channel in enumerate(selected, start=1):
        channel_started = time.perf_counter()
        log_progress(f"Reading raw channel {idx}/{len(selected)}: {channel}")
        try:
            ts = TimeSeries.read(gwf_paths, channel, start=t0, end=t1)
        except Exception as exc:
            log_progress(
                f"Read failed for {channel} after {format_elapsed(time.perf_counter() - channel_started)}: "
                f"{exc.__class__.__name__}: {exc}"
            )
            continue
        sample_count = len(ts)
        rows[channel] = pd.DataFrame(
            {"gps": np.asarray(ts.times.value, float), channel: np.asarray(ts.value, float)}
        ).set_index("gps")
        log_progress(
            f"Loaded raw channel {idx}/{len(selected)}: {channel} in "
            f"{format_elapsed(time.perf_counter() - channel_started)} "
            f"with {sample_count} samples"
        )

    if not rows:
        raise RuntimeError("None of the requested raw channels are populated")

    log_progress(f"Loaded raw channel data in {format_elapsed(time.perf_counter() - read_started)}")
    merged = pd.concat(rows.values(), axis=1, join="outer").sort_index()
    merged.index.name = "gps"
    log_progress(
        f"Merged {len(rows)} populated channels into a dataframe with "
        f"{len(merged)} rows and {len(merged.columns)} columns"
    )
    return merged.reset_index()

def plot_channel(ax, channel: str, color: str, df: pd.DataFrame, glitch_gps: float,
                 max_plot_points: int, args, glitch_times=()):
    if channel not in df.columns:
        return False

    sub = df[["gps", channel]].dropna().copy()
    if sub.empty:
        return False

    times = sub["gps"].to_numpy(dtype=float)
    values = sub[channel].to_numpy(dtype=float)
    values, note = transform_panel(channel, times, values, glitch_gps, args)
    y_for_limits = values.copy()
    t_rel = times - glitch_gps
    t_rel, values = thin_for_plot(t_rel, values, max_plot_points)

    ax.plot(t_rel, values, color=color, lw=0.7, label=channel)
    ax.axvline(0, color="crimson", lw=1.1, ls="--", alpha=0.85)
    for gt in glitch_times:
        ax.axvline(gt - glitch_gps, color="gray", lw=0.8, ls="--", alpha=0.6)
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
        ax.text(0.01, 0.88, note,
                transform=ax.transAxes, fontsize=6.5, color="0.35")
    return True


def main():
    args = parse_args()
    if not args.csv and not args.gwf:
        raise SystemExit("Provide --csv or --gwf")

    render_started = time.perf_counter()
    input_desc = f"{len(args.gwf)} GWF files" if args.gwf else f"CSV {args.csv}"
    log_progress(
        f"Starting render for GPS {args.gps_start:.0f}-{args.gps_end:.0f} "
        f"in {args.panel_mode} mode using {input_desc}"
    )

    if args.csv:
        csv_started = time.perf_counter()
        log_progress(f"Loading CSV {args.csv}")
        df = pd.read_csv(args.csv).sort_values("gps")
        log_progress(f"Loaded CSV with {len(df)} rows in {format_elapsed(time.perf_counter() - csv_started)}")
    else:
        df = extract_raw_window(args.gwf, args.glitch_gps, args.window)

    mask = (df["gps"] >= args.glitch_gps - args.window) & (df["gps"] <= args.glitch_gps + args.window)
    df = df.loc[mask].copy().sort_values("gps")
    if df.empty:
        raise RuntimeError("No data in the requested time window")
    log_progress(f"Windowed dataframe has {len(df)} rows and {len(df.columns) - 1} data channels")

    if args.save_csv:
        args.save_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.save_csv, index=False)
        log_progress(f"Saved extracted window CSV to {args.save_csv}")

    log_progress("Resolving glitch markers")
    glitch_times = resolve_glitch_times(args.glitch_times, args.gps_start, args.gps_end)
    log_progress(f"Loaded {len(glitch_times)} glitch markers in window")

    available = [(channel, color) for channel, color in RAW_CHANNELS if channel in df.columns and df[channel].notna().any()]
    if not available:
        raise RuntimeError("None of the requested raw channels are populated in the selected window")

    show_q = bool(args.gwf)
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
        plot_q_panel(axes[0], args, glitch_times=glitch_times)
        axes[0].tick_params(axis="x", labelbottom=False)
        ts_axes = axes[1:]

    for idx, (channel, color) in enumerate(available):
        log_progress(f"Plotting panel {idx + 1}/{len(available)}: {channel}")
        ax = ts_axes[idx]
        plotted = plot_channel(
            ax,
            channel,
            color,
            df,
            args.glitch_gps,
            args.max_plot_points,
            args,
            glitch_times=glitch_times,
        )
        if not plotted:
            ax.text(0.5, 0.5, f"No data for {channel}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="dimgray")

    # Adaptive xtick spacing so we never exceed ~25 ticks across the plot
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
            f"  |  {total_span:.0f} s  |  raw data{mode_note}"
        )
    else:
        glitch_utc = gps_to_utc_string(args.glitch_gps)
        top_note = (
            f"GPS {args.glitch_gps:.6f} ({glitch_utc})"
            f"  |  Window +/-{args.window:.0f} s  |  raw data{mode_note}"
        )
    fig.suptitle(top_note, fontsize=9, fontweight="bold")

    log_progress(f"Saving figure to {args.out}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log_progress(f"Saved -> {args.out}")
    log_progress(
        f"Rows: {len(df)}  Columns: {len(df.columns)}  "
        f"Total runtime: {format_elapsed(time.perf_counter() - render_started)}"
    )
    if args.save_csv:
        log_progress(f"Saved CSV -> {args.save_csv}")


if __name__ == "__main__":
    main()
