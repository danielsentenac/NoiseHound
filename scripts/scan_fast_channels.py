#!/usr/bin/env python3
"""
Scan auxiliary channels for correlation with a glitch in Hrec_hoft.

For each channel in a batch:
  1. Read ±(window + guard) seconds around the glitch GPS time
  2. Bandpass 37–44 Hz (requires sample rate >= 100 Hz; low-rate channels fall
     back to absolute-deviation SNR)
  3. Hilbert envelope
  4. Smooth with 1 s sliding mean
  5. SNR = (value - median) / MAD over the full loaded window
  6. Report peak SNR within ±window seconds of glitch

Channel selection: all channels with actual sample rate >= min_rate Hz, read
from GWF frame metadata via frameCPP (falls back to name-suffix heuristic if
frameCPP is unavailable).
Use --channel-offset / --n-channels to process in batches.
If --plot-dir is given, a stacked-panel PNG is saved for each batch.
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


BP_LOW = 37.0
BP_HIGH = 44.0
SMOOTH_S = 1.0
MIN_RATE_FOR_BANDPASS = 100   # Hz


def log(msg):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{stamp}] {msg}", flush=True)


def parse_rate_from_name(name: str):
    m = re.search(r"_(\d+)Hz$", name)
    return int(m.group(1)) if m else None


def find_gwf_files(gwf_dir: Path, t0: float, t1: float) -> list[Path]:
    """Return GWF files from gwf_dir whose time range overlaps [t0, t1]."""
    files = []
    for p in sorted(gwf_dir.glob("V-raw-*.gwf")):
        m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", p.name)
        if not m:
            continue
        fs, fd = float(m.group(1)), float(m.group(2))
        if fs < t1 and fs + fd > t0:
            files.append(p)
    return files


def get_all_channels_with_rate(gwf_path: Path, min_rate: int) -> list[tuple[str, int]]:
    """Return (name, rate) for every time-series channel with rate >= min_rate.

    Uses lalframe TOC + FrFileQueryChanVectorLength to read actual sample rates
    from GWF frame metadata without loading any signal data.
    Falls back to name-suffix heuristic if lalframe is unavailable.
    """
    log(f"Reading channel metadata from {gwf_path.name} ...")

    try:
        import lalframe
        import os

        # Suppress C-level XLAL error messages for channels that aren't
        # time-series (frequency series, etc.) — they're harmless.
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        saved_stderr = os.dup(2)
        os.dup2(devnull_fd, 2)
        try:
            frfile  = lalframe.FrFileOpenURL(str(gwf_path))
            ufrfile = lalframe.FrameUFrFileOpen(str(gwf_path), "r")
            toc     = lalframe.FrameUFrTOCRead(ufrfile)
            frame_dt = lalframe.FrFileQueryDt(frfile, 0)  # frame duration in seconds

            n_adc  = lalframe.FrameUFrTOCQueryAdcN(toc)
            n_proc = lalframe.FrameUFrTOCQueryProcN(toc)
            n_sim  = lalframe.FrameUFrTOCQuerySimN(toc)
            total  = n_adc + n_proc + n_sim

            channels = []
            for getter, n, ctype in [
                (lalframe.FrameUFrTOCQueryAdcName,  n_adc,  lalframe.ADC_CHAN),
                (lalframe.FrameUFrTOCQueryProcName, n_proc, lalframe.PROC_CHAN),
                (lalframe.FrameUFrTOCQuerySimName,  n_sim,  lalframe.SIM_CHAN),
            ]:
                for i in range(n):
                    name = getter(toc, i)
                    try:
                        vlen = lalframe.FrFileQueryChanVectorLength(frfile, name, ctype)
                        if vlen > 0:
                            rate = int(round(vlen / frame_dt))
                            if rate >= min_rate:
                                channels.append((name, rate))
                    except Exception:
                        pass
        finally:
            os.dup2(saved_stderr, 2)
            os.close(devnull_fd)
            os.close(saved_stderr)

        channels.sort(key=lambda x: x[0])
        log(f"  {total} channels total, {len(channels)} with rate >= {min_rate} Hz (lalframe metadata)")
        return channels

    except ImportError:
        log("  lalframe not available — falling back to name-suffix heuristic")
        from gwpy.io.gwf import iter_channel_names
        all_names = list(iter_channel_names(str(gwf_path)))
        log(f"  {len(all_names)} channels total")
        fast = [
            (name, rate)
            for name in all_names
            if (rate := parse_rate_from_name(name)) is not None and rate >= min_rate
        ]
        fast.sort(key=lambda x: x[0])
        log(f"  {len(fast)} channels with rate >= {min_rate} Hz (name-suffix fallback)")
        return fast


def snr_normalize(arr: np.ndarray) -> np.ndarray:
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    if mad < 1e-30:
        return np.zeros_like(arr)
    return (arr - med) / mad


def boxcar_smooth(arr: np.ndarray, half_width: int) -> np.ndarray:
    if half_width <= 0:
        return arr
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(arr.astype(np.float64), size=2 * half_width + 1, mode="reflect")


def read_batch(
    gwf_files: list[Path],
    channels: list[str],
    t0: float,
    t1: float,
) -> dict:
    """Read all channels in one TimeSeriesDict pass. Returns dict channel→TimeSeries.

    When multiple GWF files are needed, gwpy fails if a channel has incompatible
    units across files (unit mismatch on append). We avoid this by reading each
    file separately and concatenating manually with unit override.
    """
    from gwpy.timeseries import TimeSeriesDict

    files = [str(f) for f in gwf_files]

    # Fast path: single file, no concatenation needed
    if len(files) == 1:
        return dict(TimeSeriesDict.read(files, channels, start=t0, end=t1))

    # Read each GWF file independently, clipping the time range to that file's
    # actual coverage so gwpy doesn't raise on partial-span requests
    parts = []
    for gwf_f, f in zip(gwf_files, files):
        m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", gwf_f.name)
        if m:
            fs, fd = float(m.group(1)), float(m.group(2))
            ft0, ft1 = max(t0, fs), min(t1, fs + fd)
        else:
            ft0, ft1 = t0, t1
        try:
            part = TimeSeriesDict.read([f], channels, start=ft0, end=ft1)
            parts.append(dict(part))
        except Exception as exc:
            log(f"  Warning: failed to read {gwf_f.name}: {exc}")
            parts.append({})

    # Manually concatenate segments, bypassing gwpy unit check on append
    result = {}
    for ch in channels:
        segments = [p[ch] for p in parts if ch in p]
        if not segments:
            result[ch] = {"_chan_error": "not found in any GWF file"}
            continue
        ts = segments[0]
        for seg in segments[1:]:
            try:
                ts = ts.append(seg, gap="pad", inplace=False)
            except Exception:
                # Unit mismatch: strip unit from the incoming segment and retry
                try:
                    seg.override_unit(ts.unit)
                    ts = ts.append(seg, gap="pad", inplace=False)
                except Exception:
                    pass  # keep whatever we have from the first file
        result[ch] = ts
    return result


def process_channel_from_ts(
    ts,
    channel: str,
    rate: int,
    gps_center: float,
    display_window: float,
) -> dict:
    """
    Analyse one TimeSeries object already loaded in memory.
    Returns a result dict.  On success includes:
      peak_snr, peak_time_rel, method
      _t  (relative time array, display window only)
      _snr (SNR array matching _t)
      _raw (raw values array matching _t)
    On failure includes:
      error
    """
    from scipy.signal import butter, sosfiltfilt, hilbert

    times = np.asarray(ts.times.value, dtype=np.float64)
    values = np.asarray(ts.value, dtype=np.float64)

    disp_mask = (times >= gps_center - display_window) & (times <= gps_center + display_window)
    if disp_mask.sum() < 3:
        return {"channel": channel, "rate": rate, "error": "too few samples in display window"}

    t_rel_disp = times[disp_mask] - gps_center
    raw_disp = values[disp_mask]

    # Search for the peak only in the pre-glitch half (t < 0 = causal priors).
    pre_mask = t_rel_disp < 0.0
    if pre_mask.sum() < 2:
        return {"channel": channel, "rate": rate, "error": "too few pre-glitch samples"}

    nyq = rate / 2.0
    if rate < MIN_RATE_FOR_BANDPASS or BP_HIGH >= nyq:
        # Cannot bandpass at 37-44 Hz — use absolute-deviation SNR
        snr_full = snr_normalize(np.abs(values - np.median(values)))
        snr_disp = snr_full[disp_mask]
        peak_idx = int(np.argmax(snr_disp[pre_mask]))
        peak_idx = int(np.where(pre_mask)[0][peak_idx])
        return {
            "channel": channel,
            "rate": rate,
            "peak_snr": float(snr_disp[peak_idx]),
            "peak_time_rel": float(t_rel_disp[peak_idx]),
            "method": "absdev",
            "_t": t_rel_disp,
            "_snr": snr_disp,
            "_raw": raw_disp,
        }

    # Bandpass 37–44 Hz
    sos = butter(4, [BP_LOW / nyq, BP_HIGH / nyq], btype="bandpass", output="sos")
    filtered = sosfiltfilt(sos, values)

    # Hilbert envelope → smooth → SNR over full window → trim to display
    envelope = np.abs(hilbert(filtered))
    half_w = max(1, int(rate * SMOOTH_S / 2))
    envelope = boxcar_smooth(envelope, half_w)

    snr_full = snr_normalize(envelope)
    snr_disp = snr_full[disp_mask]

    pre_peak_idx = int(np.argmax(snr_disp[pre_mask]))
    peak_idx = int(np.where(pre_mask)[0][pre_peak_idx])
    return {
        "channel": channel,
        "rate": rate,
        "peak_snr": float(snr_disp[peak_idx]),
        "peak_time_rel": float(t_rel_disp[peak_idx]),
        "method": "bp_envelope",
        "_t": t_rel_disp,
        "_snr": snr_disp,
        "_raw": raw_disp,
    }


def render_batch_plot(
    results: list[dict],
    gps_center: float,
    display_window: float,
    snr_threshold: float,
    batch_offset: int,
    output_path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    plottable = [r for r in results if "_snr" in r]
    n = len(plottable)
    if n == 0:
        log("No plottable results — skipping figure")
        return

    panel_h = 1.6
    fig_h = max(4.0, n * panel_h + 1.2)
    fig, axes = plt.subplots(n, 1, figsize=(14, fig_h), sharex=True)
    if n == 1:
        axes = [axes]

    fig.suptitle(
        f"Fast-channel 37–44 Hz envelope SNR scan\n"
        f"GPS {gps_center:.0f}  ±{display_window:.0f} s  |  batch offset {batch_offset}",
        fontsize=10, y=1.0,
    )

    for ax, r in zip(axes, plottable):
        t = r["_t"]
        snr = r["_snr"]
        above = r["peak_snr"] >= snr_threshold
        color = "tab:red" if above else "tab:gray"

        ax.plot(t, snr, color=color, lw=0.8)
        ax.axhline(snr_threshold, color="black", lw=0.6, ls="--", alpha=0.5)
        ax.axvline(0.0, color="black", lw=0.8, ls="-", alpha=0.4)

        # shade ±window (already is the window, but mark t=0 region clearly)
        ax.axvspan(-display_window, display_window, color="yellow", alpha=0.04, zorder=0)

        if above:
            ax.axvline(r["peak_time_rel"], color="red", lw=0.8, ls=":", alpha=0.7)

        label = r["channel"].replace("V1:", "")
        snr_str = f"  SNR={r['peak_snr']:.1f}"
        method_str = f"  [{r['method']}]"
        ax.text(
            0.01, 0.88, label + snr_str + method_str,
            transform=ax.transAxes, fontsize=7,
            va="top", color=color if above else "dimgray",
            fontweight="bold" if above else "normal",
        )

        ax.set_xlim(-display_window, display_window)
        ax.set_ylabel("SNR", fontsize=7)
        ax.tick_params(labelsize=7)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, integer=True))

        # light grid
        ax.grid(True, which="major", lw=0.3, alpha=0.4)

    axes[-1].set_xlabel("Time relative to glitch [s]", fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved plot: {output_path}  ({output_path.stat().st_size / 1e3:.0f} kB)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--gwf", type=Path, help="Single raw GWF file to scan.")
    g.add_argument("--gwf-dir", type=Path,
                   help="Directory of raw GWF files; all files overlapping the "
                        "requested window are used automatically.")
    p.add_argument("--gps", type=float, required=True, help="Glitch GPS time.")
    p.add_argument("--window", type=float, default=15.0,
                   help="±N seconds around glitch to search for peak (default 15).")
    p.add_argument("--guard", type=float, default=30.0,
                   help="Extra seconds on each side for filter guard band (default 30).")
    p.add_argument("--min-rate", type=int, default=10,
                   help="Minimum channel sample rate Hz to include (default 10).")
    p.add_argument("--n-channels", type=int, default=10,
                   help="Number of channels to process (default 10; 0 = all).")
    p.add_argument("--channel-offset", type=int, default=0,
                   help="Skip the first N channels — for batching.")
    p.add_argument("--snr-threshold", type=float, default=5.0,
                   help="Highlight channels with peak SNR above this (default 5).")
    p.add_argument("--output", type=Path, default=None,
                   help="Optional CSV output path.")
    p.add_argument("--plot-dir", type=Path, default=None,
                   help="Directory for output PNG plots.")
    return p.parse_args()


def main():
    args = parse_args()

    t0_load = args.gps - args.window - args.guard
    t1_load = args.gps + args.window + args.guard

    if args.gwf:
        gwf_files = [args.gwf]
        ref_file = args.gwf
    else:
        gwf_files = find_gwf_files(args.gwf_dir, t0_load, t1_load)
        if not gwf_files:
            raise SystemExit(f"No GWF files found in {args.gwf_dir} for window {t0_load:.0f}–{t1_load:.0f}")
        log(f"Using {len(gwf_files)} GWF files: {[f.name for f in gwf_files]}")
        ref_file = gwf_files[0]

    fast = get_all_channels_with_rate(ref_file, min_rate=args.min_rate)

    start_i = args.channel_offset
    end_i = start_i + args.n_channels if args.n_channels > 0 else len(fast)
    batch = fast[start_i:end_i]

    if not batch:
        log(f"No channels at offset {start_i} (total {len(fast)}). Nothing to do.")
        raise SystemExit(0)

    batch_channels = [ch for ch, _ in batch]
    rate_map = {ch: r for ch, r in batch}

    log(
        f"Batch: channels {start_i}–{start_i + len(batch) - 1} ({len(batch)} total), "
        f"GPS={args.gps:.0f} ±{args.window:.0f}s, guard={args.guard:.0f}s"
    )

    # Read all batch channels in one pass — avoids reopening the GWF N times
    log(f"Reading {len(batch_channels)} channels in one TimeSeriesDict pass ...")
    tsd = read_batch(gwf_files, batch_channels, t0_load, t1_load)
    log(f"Read complete. Processing channels ...")

    results = []
    for i, (channel, rate) in enumerate(batch):
        log(f"  [{i + 1}/{len(batch)}] {channel} ({rate} Hz)")
        entry = tsd.get(channel)
        if entry is None:
            r = {"channel": channel, "rate": rate, "error": "not returned by read"}
        elif isinstance(entry, dict) and "_chan_error" in entry:
            r = {"channel": channel, "rate": rate, "error": entry["_chan_error"]}
        else:
            r = process_channel_from_ts(entry, channel, rate, args.gps, args.window)
        results.append(r)
        if "peak_snr" in r:
            flag = "  ***" if r["peak_snr"] >= args.snr_threshold else ""
            log(
                f"    snr={r['peak_snr']:.2f}  t_rel={r['peak_time_rel']:+.2f}s"
                f"  [{r['method']}]{flag}"
            )
        else:
            log(f"    SKIP: {r.get('error', '?')}")

    # Ranked summary
    ranked = sorted(
        [r for r in results if "peak_snr" in r],
        key=lambda x: x["peak_snr"],
        reverse=True,
    )
    hits = [r for r in ranked if r["peak_snr"] >= args.snr_threshold]

    print("\n=== RANKED RESULTS ===")
    if not hits:
        print(f"No channels exceeded SNR threshold {args.snr_threshold:.1f}")
    else:
        print(f"{'#':<4} {'SNR':>7}  {'t_rel':>8}  {'Rate':>6}  Channel")
        for k, r in enumerate(hits, 1):
            print(
                f"{k:<4} {r['peak_snr']:>7.2f}  {r['peak_time_rel']:>+8.2f}s"
                f"  {r['rate']:>5}Hz  {r['channel']}"
            )

    # Plot
    if args.plot_dir:
        plot_path = args.plot_dir / f"scan_{args.gps:.0f}_offset{start_i:05d}.png"
        render_batch_plot(
            results,
            gps_center=args.gps,
            display_window=args.window,
            snr_threshold=args.snr_threshold,
            batch_offset=start_i,
            output_path=plot_path,
        )

    # CSV
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = ["channel", "rate", "peak_snr", "peak_time_rel", "method", "error"]
        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(results)
        log(f"Wrote {len(results)} rows to {args.output}")


if __name__ == "__main__":
    main()
