#!/usr/bin/env python3
"""
Render glitch-response time-series plots for scan candidate channels.

Reads an aggregated candidates CSV (produced by launch_full_scan.sh),
filters channels above --snr-threshold, reads all of them from raw GWF files
in a single TimeSeriesDict pass, then renders stacked multi-panel plots in
batches of --batch-size channels sorted by peak SNR descending.

Output: one PNG per batch → <plot-dir>/candidates_<gps>_batch<NNN>.png
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
MIN_RATE_FOR_BANDPASS = 100


def log(msg):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{stamp}] {msg}", flush=True)


def find_gwf_files(gwf_dir: Path, t0: float, t1: float) -> list[Path]:
    files = []
    for p in sorted(gwf_dir.glob("V-raw-*.gwf")):
        m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", p.name)
        if not m:
            continue
        fs, fd = float(m.group(1)), float(m.group(2))
        if fs < t1 and fs + fd > t0:
            files.append(p)
    return files


def read_all_channels(gwf_files: list[Path], channels: list[str],
                      t0: float, t1: float) -> dict:
    from gwpy.timeseries import TimeSeriesDict
    src = [str(f) for f in gwf_files]
    log(f"Reading {len(channels)} channels from {len(gwf_files)} GWF files ...")
    tsd = TimeSeriesDict.read(src, channels, start=t0, end=t1)
    log("Read complete.")
    return dict(tsd)


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


def make_envelope(values: np.ndarray, rate: int) -> np.ndarray | None:
    """Return smoothed 37-44 Hz Hilbert envelope, or None if rate too low."""
    from scipy.signal import butter, sosfiltfilt, hilbert
    nyq = rate / 2.0
    if rate < MIN_RATE_FOR_BANDPASS or BP_HIGH >= nyq:
        return None
    sos = butter(4, [BP_LOW / nyq, BP_HIGH / nyq], btype="bandpass", output="sos")
    filtered = sosfiltfilt(sos, values.astype(np.float64))
    envelope = np.abs(hilbert(filtered))
    half_w = max(1, int(rate * SMOOTH_S / 2))
    return boxcar_smooth(envelope, half_w)


def render_batch(batch: list[dict], tsd: dict, gps_center: float,
                 display_window: float, batch_idx: int, out_path: Path,
                 show_envelope: bool = True) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    n = len(batch)
    panel_h = 1.8
    fig, axes = plt.subplots(n, 1, figsize=(14, n * panel_h + 1.4), sharex=True)
    if n == 1:
        axes = [axes]

    fig.suptitle(
        f"Scan candidates — glitch response  |  GPS {gps_center:.0f}"
        f"  ±{display_window:.0f} s  |  batch {batch_idx + 1}",
        fontsize=10, y=1.0,
    )

    for ax, row in zip(axes, batch):
        channel = row["channel"]
        snr     = float(row["peak_snr"])
        t_rel   = float(row["peak_time_rel"])
        method  = row.get("method", "?")
        rate    = int(row.get("rate", 0))

        if channel not in tsd:
            ax.text(0.5, 0.5, f"{channel}\n(not available in GWF)", ha="center",
                    va="center", transform=ax.transAxes, fontsize=8, color="gray")
            ax.set_ylabel("–", fontsize=7)
            ax.set_xlim(-display_window, display_window)
            continue

        times  = np.asarray(tsd[channel].times.value, dtype=np.float64)
        values = np.asarray(tsd[channel].value,       dtype=np.float64)
        disp_mask = (times >= gps_center - display_window) & \
                    (times <= gps_center + display_window)
        t_plot = times[disp_mask] - gps_center
        v_plot = values[disp_mask]

        # colour grade by SNR
        if snr >= 15:
            color = "tab:red"
        elif snr >= 10:
            color = "tab:orange"
        elif snr >= 7:
            color = "goldenrod"
        else:
            color = "tab:blue"

        ax.plot(t_plot, v_plot, color=color, lw=0.7, alpha=0.85, label="raw")

        # overlay 37-44 Hz envelope for fast channels
        if show_envelope and rate >= MIN_RATE_FOR_BANDPASS:
            env = make_envelope(values, rate)
            if env is not None:
                env_disp = env[disp_mask]
                # scale envelope to raw signal amplitude for overlay
                scale = np.percentile(np.abs(v_plot), 95) / (np.percentile(env_disp, 95) + 1e-30)
                ax.plot(t_plot, env_disp * scale, color="black", lw=0.6,
                        ls="--", alpha=0.5, label="37-44Hz env (scaled)")

        ax.axvline(0.0,    color="black", lw=0.8, alpha=0.5, ls="--")
        ax.axvline(t_rel,  color="red",   lw=0.7, alpha=0.6, ls=":")

        label = (f"{channel.replace('V1:', '')}   "
                 f"SNR={snr:.1f}  t={t_rel:+.2f}s  [{method}]  {rate}Hz")
        ax.text(0.01, 0.88, label, transform=ax.transAxes, fontsize=7.5,
                va="top", color=color,
                fontweight="bold" if snr >= 10 else "normal")

        ax.set_xlim(-display_window, display_window)
        ax.tick_params(labelsize=6)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3))
        ax.grid(True, lw=0.3, alpha=0.4)

    axes[-1].set_xlabel("Time relative to glitch [s]", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved {out_path}  ({out_path.stat().st_size / 1e3:.0f} kB)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", type=Path, required=True,
                   help="Candidates CSV from launch_full_scan.sh.")
    p.add_argument("--gwf-dir", type=Path, required=True,
                   help="Directory of raw GWF files.")
    p.add_argument("--gps", type=float, required=True,
                   help="Glitch GPS time.")
    p.add_argument("--window", type=float, default=15.0,
                   help="±N seconds display window (default 15).")
    p.add_argument("--guard", type=float, default=30.0,
                   help="Extra seconds on each side for filter guard band (default 30).")
    p.add_argument("--snr-threshold", type=float, default=5.0,
                   help="Minimum SNR to include (default 5).")
    p.add_argument("--batch-size", type=int, default=10,
                   help="Channels per plot (default 10).")
    p.add_argument("--no-envelope", action="store_true",
                   help="Do not overlay 37-44 Hz envelope.")
    p.add_argument("--plot-dir", type=Path, required=True,
                   help="Output directory for PNG plots.")
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.candidates) as fh:
        all_rows = list(csv.DictReader(fh))

    candidates = [
        r for r in all_rows
        if not r.get("error")
        and float(r.get("peak_snr", 0)) >= args.snr_threshold
    ]
    candidates.sort(key=lambda r: float(r["peak_snr"]), reverse=True)
    log(f"{len(candidates)} candidates at SNR >= {args.snr_threshold:.1f}")

    if not candidates:
        raise SystemExit("No candidates to render.")

    t0 = args.gps - args.window - args.guard
    t1 = args.gps + args.window + args.guard
    gwf_files = find_gwf_files(args.gwf_dir, t0, t1)
    if not gwf_files:
        raise SystemExit(f"No GWF files found in {args.gwf_dir} for window {t0:.0f}–{t1:.0f}")
    log(f"Using {len(gwf_files)} GWF files: {[f.name for f in gwf_files]}")

    # Single read for all candidates
    all_channels = list(dict.fromkeys(r["channel"] for r in candidates))
    tsd = read_all_channels(gwf_files, all_channels, t0, t1)

    n_batches = (len(candidates) + args.batch_size - 1) // args.batch_size
    log(f"Rendering {n_batches} plot(s) of {args.batch_size} channels each ...")
    for i in range(n_batches):
        batch = candidates[i * args.batch_size:(i + 1) * args.batch_size]
        out_path = args.plot_dir / f"candidates_{args.gps:.0f}_batch{i + 1:03d}.png"
        render_batch(
            batch, tsd, args.gps, args.window, i, out_path,
            show_envelope=not args.no_envelope,
        )

    log(f"Done. {n_batches} plot(s) saved to {args.plot_dir}")


if __name__ == "__main__":
    main()
