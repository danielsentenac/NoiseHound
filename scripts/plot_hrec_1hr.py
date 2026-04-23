#!/usr/bin/env python3
"""
Plot a 1-hour amplitude overview of V1:Hrec_hoft_16384Hz around a target GPS.

Reads HoftOnline GWF files, computes 1-second RMS of the bandpassed (or raw)
strain, and marks the target GPS with a red dashed line.

Usage (on CCA):
    /cvmfs/software.igwn.org/conda/envs/igwn/bin/python scripts/plot_hrec_1hr.py \
        --target-gps 1419724815 \
        --gwf-dir frame_cache/hoftonline_1419724818_1419735618 \
        --out usecases/1419724815/hrec_1hr_overview.png
"""
import argparse
import os
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
from glob import glob

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WORKDIR = Path(__file__).parent.parent
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-gps", type=float, required=True,
                        help="Central GPS time to mark with a red line.")
    parser.add_argument("--gwf-dir", type=Path, required=True,
                        help="Directory containing HoftOnline GWF files.")
    parser.add_argument("--channel", default="V1:Hrec_hoft_16384Hz",
                        help="GW strain channel to read.")
    parser.add_argument("--window", type=float, default=1800.0,
                        help="Half-width of the plot window in seconds (default 1800 = 1 hr).")
    parser.add_argument("--seg-s", type=float, default=1.0,
                        help="RMS averaging segment length in seconds (default 1 s).")
    parser.add_argument("--flo", type=float, default=20.0,
                        help="Highpass cutoff Hz (0 to disable).")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output PNG path.")
    return parser.parse_args()


def gps_to_utc(gps: float) -> str:
    return (GPS_EPOCH + pd.to_timedelta(gps, unit="s")).strftime("%Y-%m-%d %H:%M:%S UTC")


def read_rms_trend(gwf_files, channel, t0, t1, seg_s, flo):
    from gwpy.timeseries import TimeSeries, TimeSeriesList

    # Read each file that overlaps the window
    segments = []
    for gwf in sorted(gwf_files):
        # Parse filename timing
        stem = Path(gwf).stem  # e.g. V-HoftOnline-1419724000-2000
        parts = stem.split("-")
        try:
            file_start = int(parts[-2])
            file_dur = int(parts[-1])
            file_end = file_start + file_dur
        except (ValueError, IndexError):
            file_start, file_end = 0, int(1e18)

        if file_end <= t0 or file_start >= t1:
            continue

        seg_t0 = max(t0, file_start)
        seg_t1 = min(t1, file_end)
        try:
            ts = TimeSeries.read(gwf, channel, start=seg_t0, end=seg_t1)
        except Exception as e:
            print(f"  Warning: could not read {gwf}: {e}")
            continue
        segments.append(ts)

    if not segments:
        raise RuntimeError(f"No data could be read for {channel} in [{t0}, {t1}]")

    # Concatenate all segments
    from gwpy.timeseries import TimeSeriesList
    ts = TimeSeriesList(*segments).join(gap="pad")

    # Optional highpass
    if flo > 0:
        ts = ts.highpass(flo)

    # Downsample to 1/seg_s Hz via RMS blocks
    fs = float(ts.sample_rate.value)
    block = int(round(seg_s * fs))
    data = np.asarray(ts.value, dtype=float)
    t_start = float(ts.t0.value)

    n_blocks = len(data) // block
    gps_arr = np.empty(n_blocks)
    rms_arr = np.empty(n_blocks)
    for i in range(n_blocks):
        chunk = data[i * block:(i + 1) * block]
        gps_arr[i] = t_start + (i + 0.5) * seg_s
        rms_arr[i] = np.sqrt(np.mean(chunk**2)) if np.isfinite(chunk).any() else np.nan

    return gps_arr, rms_arr


def main():
    args = parse_args()

    gwf_files = sorted(glob(str(args.gwf_dir / "*.gwf")))
    if not gwf_files:
        raise SystemExit(f"No GWF files found in {args.gwf_dir}")

    t0 = args.target_gps - args.window
    t1 = args.target_gps + args.window

    print(f"Reading {args.channel} over [{t0:.0f}, {t1:.0f}] ...")
    gps_arr, rms_arr = read_rms_trend(gwf_files, args.channel, t0, t1, args.seg_s, args.flo)

    t_rel = gps_arr - args.target_gps

    fig, ax = plt.subplots(figsize=(18, 4))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.87, bottom=0.14)

    ax.plot(t_rel, rms_arr, color="steelblue", lw=0.7, alpha=0.9)
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.9, label=f"Target GPS {args.target_gps:.0f}")
    ax.set_xlabel(f"Time relative to GPS {args.target_gps:.0f}  ({gps_to_utc(args.target_gps)})  [s]",
                  fontsize=9)
    ax.set_ylabel("1-s RMS strain", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(axis="x", ls=":", alpha=0.4)
    ax.tick_params(labelsize=8)

    # x-ticks every 300 s
    step = 300.0
    xticks = np.arange(np.ceil(t_rel[0] / step) * step,
                       np.floor(t_rel[-1] / step) * step + step * 0.01, step)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{v:+.0f}" for v in xticks], fontsize=7)
    ax.set_xlim(t_rel[0], t_rel[-1])

    bp_note = f"highpass {args.flo:.0f} Hz" if args.flo > 0 else "no filter"
    fig.suptitle(
        f"V1:Hrec_hoft_16384Hz — 1-s RMS amplitude overview\n"
        f"GPS {args.target_gps:.0f}  ({gps_to_utc(args.target_gps)})  |  "
        f"Window ±{args.window:.0f} s  |  {bp_note}",
        fontsize=10, fontweight="bold",
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {args.out}")
    print(f"Samples: {len(rms_arr)}, GPS range: [{gps_arr[0]:.0f}, {gps_arr[-1]:.0f}]")


if __name__ == "__main__":
    main()
