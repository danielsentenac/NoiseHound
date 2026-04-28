#!/usr/bin/env python3
"""
Stack inter-glitch intervals for watchlist channels aligned to the next glitch (τ=0).
For each channel plots individual intervals (faint) and their mean, so systematic
amplitude build-up/decay before each glitch becomes visible.
"""
import argparse
import csv
import os
import re
import warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import numpy as np
from pathlib import Path
from datetime import datetime, timezone
from scipy.ndimage import uniform_filter1d


BP_LOW, BP_HIGH = 37.0, 44.0
MIN_RATE_FOR_BANDPASS = 100
ENVELOPE_SMOOTH_S = 1.0   # Hilbert envelope smoothing (same as scan)
AMP_SMOOTH_S = 30.0        # Additional smoothing to reveal slow modulation


def log(msg):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{stamp}] {msg}", flush=True)


def find_gwf_files(gwf_dir, t0, t1):
    files = []
    for p in sorted(Path(gwf_dir).glob("V-raw-*.gwf")):
        m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", p.name)
        if not m:
            continue
        fs, fd = float(m.group(1)), float(m.group(2))
        if fs < t1 and fs + fd > t0:
            files.append(p)
    return files


def read_and_process_file(gwf_file, channels, ft0, ft1, rows_by_channel):
    """Read one GWF file, immediately compute 1 Hz amplitude proxy per channel.
    Returns dict: channel -> (tau_1hz_array, amp_1hz_array) relative to ft1.
    Never accumulates raw high-rate data across files — keeps memory at O(one file).
    """
    from gwpy.timeseries import TimeSeriesDict
    try:
        tsd = dict(TimeSeriesDict.read([str(gwf_file)], channels, start=ft0, end=ft1))
    except Exception as exc:
        log(f"  Warning: {gwf_file.name}: {exc}")
        return {}

    result = {}
    for ch, row in rows_by_channel.items():
        if ch not in tsd:
            continue
        ts     = tsd[ch]
        times  = np.asarray(ts.times.value, dtype=np.float64)
        values = np.asarray(ts.value,       dtype=np.float64)
        rate   = int(row["rate"])
        method = row["method"]

        amp = amplitude_proxy(values, rate, method)

        # Downsample to 1 Hz immediately — discard raw
        tau_1hz = np.arange(np.ceil(times[0]), times[-1], 1.0)
        amp_1hz = np.interp(tau_1hz, times, amp)
        result[ch] = (tau_1hz, amp_1hz)

    return result


def amplitude_proxy(values, rate, method):
    """Slowly-varying amplitude proxy at the channel's native rate."""
    values = np.asarray(values, dtype=np.float64)
    nyq = rate / 2.0

    if method == "bp_envelope" and rate >= MIN_RATE_FOR_BANDPASS and BP_HIGH < nyq:
        from scipy.signal import butter, sosfiltfilt, hilbert
        sos = butter(4, [BP_LOW / nyq, BP_HIGH / nyq], btype="bandpass", output="sos")
        env = np.abs(hilbert(sosfiltfilt(sos, values)))
        # 1s smooth (same as scan), then 30s smooth to expose slow modulation
        hw1 = max(1, int(rate * ENVELOPE_SMOOTH_S / 2))
        env = uniform_filter1d(env, size=2 * hw1 + 1, mode="reflect")
        hw2 = max(1, int(rate * AMP_SMOOTH_S / 2))
        return uniform_filter1d(env, size=2 * hw2 + 1, mode="reflect")
    else:
        # DC / slow channels: smoothed absolute deviation from median
        absdev = np.abs(values - np.median(values))
        hw = max(1, int(rate * AMP_SMOOTH_S / 2))
        return uniform_filter1d(absdev, size=2 * hw + 1, mode="reflect")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates",   type=Path, required=True,
                   help="Watchlist CSV (channel, rate, method, ...)")
    p.add_argument("--gwf-dir",      type=Path, required=True,
                   help="Directory containing V-raw-*.gwf files")
    p.add_argument("--glitch-times", type=str,  required=True,
                   help="Comma-separated GPS times of consecutive glitches")
    p.add_argument("--output",       type=Path, default=None)
    p.add_argument("--checkpoint-dir", type=Path, default=None,
                   help="Directory for per-interval .npy checkpoints (enables resume)")
    p.add_argument("--plot-only",    action="store_true",
                   help="Skip data reading, load all checkpoints and plot")
    p.add_argument("--scan-window",  type=float, default=15.0,
                   help="Highlight this many seconds before glitch [default 15]")
    p.add_argument("--smooth",       type=float, default=30.0,
                   help="Amplitude smoothing window in seconds [default 30]")
    return p.parse_args()


def main():
    args = parse_args()
    global AMP_SMOOTH_S
    AMP_SMOOTH_S = args.smooth

    # Load candidates — preserve CSV order
    rows = []
    with open(args.candidates) as f:
        for row in csv.DictReader(f):
            if not row.get("error"):
                rows.append(row)
    log(f"Loaded {len(rows)} candidates")

    glitch_times = [float(t) for t in args.glitch_times.split(",")]
    intervals = list(zip(glitch_times[:-1], glitch_times[1:]))
    log(f"{len(intervals)} inter-glitch intervals: " +
        ", ".join(f"{t1 - t0:.0f}s" for t0, t1 in intervals))

    channels = [r["channel"] for r in rows]
    rows_by_channel = {r["channel"]: r for r in rows}

    out = args.output or Path("interglitch_profile.png")
    ckpt_dir = args.checkpoint_dir or out.parent
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    # Per-channel list of (tau_array_1Hz, amp_array_1Hz) for each interval
    channel_traces = {r["channel"]: [] for r in rows}

    if not args.plot_only:
        for i, (t_prev, t_next) in enumerate(intervals):
            ckpt_file = ckpt_dir / f"interglitch_ckpt_{i}.npy"
            if ckpt_file.exists():
                log(f"Interval {i + 1}: checkpoint found, skipping data read")
                continue

            log(f"Interval {i + 1}: {t_prev:.0f} → {t_next:.0f}")
            gwf_files = find_gwf_files(args.gwf_dir, t_prev, t_next)
            if not gwf_files:
                log("  No GWF files — skipping")
                continue
            log(f"  {len(gwf_files)} GWF files (processing one at a time) ...")

            # Accumulate 1 Hz chunks per channel — never hold full interval in RAM
            interval_chunks = {ch: ([], []) for ch in channels}

            for gwf_f in gwf_files:
                m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", gwf_f.name)
                if m:
                    fs, fd = float(m.group(1)), float(m.group(2))
                    ft0, ft1 = max(t_prev, fs), min(t_next, fs + fd)
                else:
                    ft0, ft1 = t_prev, t_next

                file_result = read_and_process_file(
                    gwf_f, channels, ft0, ft1, rows_by_channel)

                for ch, (tau_1hz, amp_1hz) in file_result.items():
                    interval_chunks[ch][0].append(tau_1hz)
                    interval_chunks[ch][1].append(amp_1hz)

            log("  Done. Assembling and saving checkpoint ...")
            interval_data = {}
            for ch in channels:
                tau_parts, amp_parts = interval_chunks[ch]
                if not tau_parts:
                    continue
                tau_cat = np.concatenate(tau_parts)
                amp_cat = np.concatenate(amp_parts)
                order = np.argsort(tau_cat)
                tau_cat, amp_cat = tau_cat[order], amp_cat[order]
                tau_rel = tau_cat - t_next
                interval_data[ch] = (tau_rel, amp_cat)

            np.save(ckpt_file, interval_data)
            log(f"  Checkpoint saved: {ckpt_file}")

    # Load all checkpoints and merge into channel_traces
    loaded = 0
    for i in range(len(intervals)):
        ckpt_file = ckpt_dir / f"interglitch_ckpt_{i}.npy"
        if not ckpt_file.exists():
            log(f"Interval {i + 1}: no checkpoint found — skipped in plot")
            continue
        interval_data = np.load(ckpt_file, allow_pickle=True).item()
        for ch in channels:
            if ch in interval_data:
                channel_traces[ch].append(interval_data[ch])
        loaded += 1
    log(f"Loaded {loaded} checkpoint(s) for plotting")

    # --- Plot ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    n = len(rows)
    panel_h = 1.8
    fig, axes = plt.subplots(n, 1, figsize=(15, n * panel_h + 1.5), sharex=False)
    if n == 1:
        axes = [axes]

    fig.suptitle(
        f"Inter-glitch amplitude profile — {len(intervals)} intervals  "
        f"(shaded: last {args.scan_window:.0f} s before glitch)",
        fontsize=11, y=1.002,
    )

    tab10 = plt.cm.tab10.colors

    for idx, (ax, row) in enumerate(zip(axes, rows)):
        ch     = row["channel"]
        rate   = int(row["rate"])
        method = row["method"]
        color  = tab10[idx % len(tab10)]
        traces = channel_traces[ch]

        if not traces:
            ax.text(0.5, 0.5, f"{ch}  [no data]", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8)
            ax.set_xlim(-1800, 0)
            continue

        # Build common tau grid from the longest interval
        tau_min = min(tr[0][0] for tr in traces)
        tau_grid = np.arange(np.ceil(tau_min), 0.5, 1.0)

        stacked = []
        for tau_rel, amp in traces:
            amp_interp = np.interp(tau_grid, tau_rel, amp,
                                   left=np.nan, right=np.nan)
            stacked.append(amp_interp)
            ax.plot(tau_grid, amp_interp, color="#cccccc", lw=0.7,
                    alpha=0.8, zorder=1)

        mean_amp = np.nanmean(stacked, axis=0)
        ax.plot(tau_grid, mean_amp, color=color, lw=1.4, zorder=2)

        ax.axvspan(-args.scan_window, 0, color=color, alpha=0.08, zorder=0)
        ax.axvline(0, color="black", lw=0.8, alpha=0.5, zorder=3)

        # Y label: normalise axis to median of mean so relative changes are visible
        med = np.nanmedian(mean_amp)
        if med > 0:
            ax2 = ax.twinx()
            ax2.plot(tau_grid, mean_amp / med, alpha=0)  # invisible, just for scale
            ax2.set_ylabel("/ median", fontsize=5, color="grey")
            ax2.tick_params(labelsize=5, colors="grey")

        label = f"#{idx + 1}  {ch.replace('V1:', '')}  {rate} Hz  [{method}]"
        ax.text(0.01, 0.92, label, transform=ax.transAxes, fontsize=7.5,
                va="top", color=color, fontweight="bold")
        ax.set_xlim(tau_min, 0)
        ax.tick_params(labelsize=7)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3))
        ax.grid(True, lw=0.3, alpha=0.4)

    axes[-1].set_xlabel("Time relative to next glitch  τ  [s]", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 1.0])

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved: {out}  ({out.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
