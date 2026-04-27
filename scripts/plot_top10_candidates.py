#!/usr/bin/env python3
"""
Plot the 10 most significant pre-glitch candidates in a single stacked figure.
Each panel shows the raw signal (grey) and the 37-44 Hz Hilbert envelope (colour),
with the peak SNR time marked and the full ±window displayed.
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


BP_LOW, BP_HIGH = 37.0, 44.0
SMOOTH_S = 1.0
MIN_RATE_FOR_BANDPASS = 100


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


def read_channels(gwf_files, channels, t0, t1):
    from gwpy.timeseries import TimeSeriesDict
    src = [str(f) for f in gwf_files]
    # Fast path: single file
    if len(src) == 1:
        return dict(TimeSeriesDict.read(src, channels, start=t0, end=t1))
    # Multi-file: read each file only for its own time coverage to avoid
    # gwpy raising on partial-span requests
    parts = []
    for f in gwf_files:
        m = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", f.name)
        if m:
            fs, fd = float(m.group(1)), float(m.group(2))
            ft0, ft1 = max(t0, fs), min(t1, fs + fd)
        else:
            ft0, ft1 = t0, t1
        try:
            parts.append(dict(TimeSeriesDict.read([str(f)], channels, start=ft0, end=ft1)))
        except Exception as e:
            log(f"  Warning: {f.name}: {e}")
            parts.append({})
    result = {}
    for ch in channels:
        segs = [p[ch] for p in parts if ch in p]
        if not segs:
            continue
        ts = segs[0]
        for seg in segs[1:]:
            try:
                ts = ts.append(seg, gap="pad", inplace=False)
            except Exception:
                try:
                    seg.override_unit(ts.unit)
                    ts = ts.append(seg, gap="pad", inplace=False)
                except Exception:
                    pass
        result[ch] = ts
    return result


def envelope(values, rate):
    from scipy.signal import butter, sosfiltfilt, hilbert
    from scipy.ndimage import uniform_filter1d
    nyq = rate / 2.0
    if rate < MIN_RATE_FOR_BANDPASS or BP_HIGH >= nyq:
        return None
    sos = butter(4, [BP_LOW / nyq, BP_HIGH / nyq], btype="bandpass", output="sos")
    filt = sosfiltfilt(sos, values)
    env = np.abs(hilbert(filt))
    hw = max(1, int(rate * SMOOTH_S / 2))
    return uniform_filter1d(env.astype(np.float64), size=2 * hw + 1, mode="reflect")


def snr_norm(arr):
    med = np.median(arr)
    mad = np.median(np.abs(arr - med))
    return (arr - med) / mad if mad > 1e-30 else np.zeros_like(arr)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", type=Path, required=True)
    p.add_argument("--gwf-dir",    type=Path, required=True)
    p.add_argument("--gps",        type=float, required=True)
    p.add_argument("--window",     type=float, default=15.0)
    p.add_argument("--guard",      type=float, default=30.0)
    p.add_argument("--n",          type=int,   default=10)
    p.add_argument("--output",     type=Path,  default=None)
    return p.parse_args()


def main():
    args = parse_args()

    # Load top-N candidates
    rows = []
    with open(args.candidates) as f:
        for row in csv.DictReader(f):
            if row.get("error"):
                continue
            try:
                rows.append((float(row["peak_snr"]), row))
            except (ValueError, KeyError):
                pass
    rows.sort(key=lambda x: x[0], reverse=True)
    top = [r for _, r in rows[:args.n]]
    log(f"Top {len(top)} candidates selected")

    # Load GWF data
    t0 = args.gps - args.window - args.guard
    t1 = args.gps + args.window + args.guard
    gwf_files = find_gwf_files(args.gwf_dir, t0, t1)
    if not gwf_files:
        raise SystemExit(f"No GWF files found in {args.gwf_dir}")
    log(f"Using {len(gwf_files)} GWF files")

    channels = [r["channel"] for r in top]
    log(f"Reading {len(channels)} channels ...")
    tsd = read_channels(gwf_files, channels, t0, t1)
    log("Read complete.")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    n = len(top)
    panel_h = 2.2
    fig, axes = plt.subplots(n, 1, figsize=(15, n * panel_h + 1.5), sharex=True)
    if n == 1:
        axes = [axes]

    fig.suptitle(
        f"Top {n} pre-glitch candidates  —  GPS {args.gps:.0f}  ±{args.window:.0f} s",
        fontsize=11, y=1.002,
    )

    colors = plt.cm.tab10.colors

    for idx, (ax, row) in enumerate(zip(axes, top)):
        ch    = row["channel"]
        rate  = int(row["rate"])
        color = colors[idx % len(colors)]

        if ch not in tsd:
            ax.text(0.5, 0.5, f"{ch}  [no data]", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8)
            continue

        ts     = tsd[ch]
        times  = np.asarray(ts.times.value, dtype=np.float64)
        values = np.asarray(ts.value,       dtype=np.float64)
        disp   = (times >= args.gps - args.window) & (times <= args.gps + args.window)
        t_rel  = times[disp] - args.gps
        raw    = values[disp]

        # Normalise raw to ±1 for display
        rng = np.percentile(np.abs(raw), 99)
        raw_norm = raw / rng if rng > 0 else raw

        # Compute SNR from this glitch's data (full window for baseline, trim to display)
        env_full = envelope(values, rate)
        if env_full is not None:
            snr_full = snr_norm(env_full)
            snr_disp = snr_full[disp]
            method = "bp_env"
            ax.plot(t_rel, snr_disp, color=color, lw=1.0)
            ax.set_ylabel("SNR", fontsize=7)
            ax2 = ax.twinx()
            ax2.plot(t_rel, raw_norm, color="silver", lw=0.5, alpha=0.7, zorder=0)
            ax2.set_ylabel("raw (norm)", fontsize=6, color="silver")
            ax2.tick_params(labelsize=6, colors="silver")
            ax2.set_ylim(-4, 4)
        else:
            snr_disp = snr_norm(np.abs(raw - np.median(raw)))
            method = "absdev"
            ax.plot(t_rel, snr_disp, color=color, lw=1.0)
            ax.set_ylabel("SNR", fontsize=7)

        # Recompute peak time and SNR from t < 0 region of this glitch
        pre_mask = t_rel < 0.0
        if pre_mask.sum() >= 2:
            pre_peak_idx = int(np.argmax(snr_disp[pre_mask]))
            peak_idx = int(np.where(pre_mask)[0][pre_peak_idx])
            t_pk  = float(t_rel[peak_idx])
            snr_v = float(snr_disp[peak_idx])
        else:
            t_pk  = float(row["peak_time_rel"])
            snr_v = float(row["peak_snr"])

        ax.axvline(0,    color="black", lw=0.8, alpha=0.5)
        ax.axvline(t_pk, color=color,  lw=1.0, ls="--", alpha=0.8)
        ax.axhline(5,    color="black", lw=0.5, ls=":", alpha=0.4)

        label = f"#{idx+1}  {ch.replace('V1:','')}   SNR={snr_v:.1f}  t={t_pk:+.2f}s  {rate}Hz  [{method}]"
        ax.text(0.01, 0.90, label, transform=ax.transAxes, fontsize=7.5,
                va="top", color=color, fontweight="bold")

        ax.set_xlim(-args.window, args.window)
        ax.tick_params(labelsize=7)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=3, integer=True))
        ax.grid(True, lw=0.3, alpha=0.4)

    axes[-1].set_xlabel("Time relative to glitch [s]", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 1.0])

    out = args.output or Path(f"top{n}_candidates_{args.gps:.0f}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"Saved: {out}  ({out.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()
