"""
Overlay plot: Hrec spectrogram + MAD-z-score + one panel per non-DARM aux channel.
Usage: python make_overlay_plot.py <trend_gwf> <output_png>
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

from gwpy.timeseries import TimeSeries, TimeSeriesDict

# ── configuration ────────────────────────────────────────────────────────────

GPS_START    = 1419724818
GPS_END      = 1419735618
HREC_FRAMES  = sorted(Path("/pbs/home/s/sentenac/NOISEHOUND/frame_cache/hoftonline_1419724818_1419735618").glob("*.gwf"))
HREC_CH      = "V1:Hrec_hoft_16384Hz"
TRIGGERS_CSV = "/pbs/home/s/sentenac/NOISEHOUND/outputs/detect_1419724818_1419735618/triggers_p95_auto.csv"

FMIN, FMAX  = 20.0, 210.0
STRIDE      = 4.0
FFTLENGTH   = 4.0
OVERLAP     = 2.0
BG_WINDOW_S = 1800.0
THRESHOLD   = 6.0

# Non-DARM channels: (channel_name, short_label, color, rank_score, hit_frac)
AUX_CHANNELS = [
    ("V1:Sc_WI_FF50HZ_P_ERR_mean",                "Sc_WI FF50Hz P_ERR",   "tab:red",    22.8, 1.000),
    ("V1:Sc_WI_FF50HZ_PHASE_mean",                 "Sc_WI FF50Hz PHASE",   "tab:orange", 14.7, 0.875),
    ("V1:Sc_WI_FF50HZ_G_ERR_mean",                 "Sc_WI FF50Hz G_ERR",   "tab:pink",    6.0, 0.625),
    ("V1:ASC_SR_TY_DCP_mag_B1_mag_10Hz_FS_mean",   "ASC SR_TY B1 mag",     "tab:green",   8.3, 0.625),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_mag_10Hz_FS_mean","ASC DIFFp_TY B1 mag",  "tab:cyan",    4.8, 0.625),
    ("V1:SQB1_Cam1_FitWaistY_mean",                "SQB1 Cam1 WaistY",     "tab:purple",  8.0, 0.875),
    ("V1:SDB1_LC_TX_largeError_mean",               "SDB1 LC_TX largeErr",  "tab:brown",   1.6, 0.750),
]

# ── helpers ───────────────────────────────────────────────────────────────────

def rolling_median(values, window_bins):
    s = pd.Series(values)
    return s.rolling(window_bins, center=True, min_periods=max(5, window_bins // 4)).median().bfill().ffill().to_numpy()

def rolling_robust_zscore(values, window_bins, floor=1e-12):
    baseline = rolling_median(values, window_bins)
    residual = np.abs(values - baseline)
    scale = (
        pd.Series(residual)
        .rolling(window_bins, center=True, min_periods=max(5, window_bins // 4))
        .median().bfill().ffill().to_numpy()
    )
    return (values - baseline) / np.maximum(1.4826 * scale, floor)

def draw_triggers(ax, trigger_gps, alpha=0.25, lw=0.9):
    for gps in trigger_gps:
        ax.axvline(gps, color="tab:red", alpha=alpha, linewidth=lw)

# ── 1. Hrec spectrogram ───────────────────────────────────────────────────────

print("Reading Hrec frames …")
frame_paths = [str(f) for f in HREC_FRAMES]
hrec = TimeSeries.read(frame_paths, HREC_CH, start=GPS_START, end=GPS_END)
print(f"  loaded {hrec.duration} s of Hrec")

print("Computing spectrogram …")
spec = hrec.spectrogram(stride=STRIDE, fftlength=FFTLENGTH, overlap=OVERLAP, method="median").crop_frequencies(FMIN, FMAX)
matrix = np.asarray(spec.value, dtype=float)
times  = np.asarray(spec.times.value, dtype=float)
freqs  = np.asarray(spec.frequencies.value, dtype=float)

reference       = np.nanmedian(matrix, axis=0)
floor_val       = max(float(np.nanmin(matrix[matrix > 0])) * 0.1, np.finfo(float).tiny)
relative_excess = np.log10(np.maximum(matrix, floor_val) / np.maximum(reference, floor_val))
broadband_score = np.nanpercentile(relative_excess, 95, axis=1)

bg_bins  = max(5, int(round(BG_WINDOW_S / STRIDE)))
zscore   = rolling_robust_zscore(broadband_score, bg_bins)
dist_bins = max(1, int(round(600 / STRIDE)))
find_peaks(zscore, height=THRESHOLD, distance=dist_bins)

triggers    = pd.read_csv(TRIGGERS_CSV)
trigger_gps = triggers["gps_peak"].values

# ── 2. Trend aux channels ─────────────────────────────────────────────────────

trend_gwf  = sys.argv[1]
output_png = sys.argv[2]

print(f"Reading trend file {trend_gwf} …")
ch_names = [ch for ch, *_ in AUX_CHANNELS]
tdict = TimeSeriesDict.read(trend_gwf, ch_names, start=GPS_START, end=GPS_END)
print(f"  loaded {len(tdict)} channels")

aux_data = {}
for ch, label, color, rank_score, hit_frac in AUX_CHANNELS:
    if ch not in tdict:
        print(f"  WARNING: {ch} not found")
        continue
    ts  = tdict[ch]
    t   = np.asarray(ts.times.value, dtype=float)
    v   = np.asarray(ts.value, dtype=float)
    mask = (t >= GPS_START) & (t <= GPS_END)
    t, v = t[mask], v[mask]
    sr  = float(ts.sample_rate.value)
    win = max(5, int(round(BG_WINDOW_S * sr)))
    z   = rolling_robust_zscore(v, win)
    aux_data[ch] = (t, z, label, color, rank_score, hit_frac)

# ── 3. Plot ───────────────────────────────────────────────────────────────────

n_aux  = len(aux_data)
n_rows = 2 + n_aux          # spectrogram + Hrec z + one row per channel
height_ratios = [3, 1] + [1] * n_aux

fig = plt.figure(figsize=(16, 3 + 2 * n_aux))
grid = fig.add_gridspec(
    n_rows, 2,
    height_ratios=height_ratios,
    width_ratios=[40, 1],
    hspace=0.06,
    wspace=0.04,
)

ax_spec = fig.add_subplot(grid[0, 0])
ax_hrec = fig.add_subplot(grid[1, 0], sharex=ax_spec)
ax_cbar = fig.add_subplot(grid[0, 1])

# -- spectrogram
im = ax_spec.imshow(
    relative_excess.T,
    aspect="auto", origin="lower",
    extent=[times[0], times[-1], freqs[0], freqs[-1]],
    cmap="jet", vmin=-1.5, vmax=2.0,
)
ax_spec.set_ylabel("Frequency [Hz]")
ax_spec.set_title("V1:Hrec_hoft_16384Hz — broadband glitches & auxiliary channel coincidences")
ax_spec.grid(False)
ax_spec.tick_params(axis="x", labelbottom=False)
fig.colorbar(im, cax=ax_cbar, label="log10(power / median)")
draw_triggers(ax_spec, trigger_gps, alpha=0.18)

# -- Hrec z-score
ax_hrec.plot(times, zscore, color="tab:blue", linewidth=1.0)
ax_hrec.axhline(THRESHOLD, color="tab:red", linestyle="--", linewidth=0.8)
draw_triggers(ax_hrec, trigger_gps)
ax_hrec.set_ylabel("Hrec\nMAD-z", fontsize=8)
ax_hrec.tick_params(axis="x", labelbottom=False)
ax_hrec.grid(alpha=0.2)

# -- one panel per aux channel
ax_aux_list = []
for row_idx, (ch, (t, z, label, color, rank_score, hit_frac)) in enumerate(aux_data.items()):
    is_last = (row_idx == n_aux - 1)
    ax = fig.add_subplot(grid[2 + row_idx, 0], sharex=ax_spec)
    ax_aux_list.append(ax)

    ax.plot(t, z, color=color, linewidth=1.1)
    ax.axhline(THRESHOLD, color="gray", linestyle="--", linewidth=0.7, alpha=0.7)
    ax.axhline(0, color="gray", linestyle="-", linewidth=0.4, alpha=0.4)
    draw_triggers(ax, trigger_gps)

    # auto y-limit: clip extreme outliers (±5× 99th percentile)
    p99 = np.nanpercentile(np.abs(z), 99)
    ylim = max(THRESHOLD * 1.5, min(p99 * 1.5, 50))
    ax.set_ylim(-ylim, ylim)

    # right-side annotation
    ax.set_ylabel(f"MAD-z", fontsize=7)
    title = f"{label}   [rank={rank_score:.1f}, hit={hit_frac:.0%}]"
    ax.text(0.01, 0.88, title, transform=ax.transAxes,
            fontsize=7.5, color=color, fontweight="bold",
            va="top", ha="left",
            bbox=dict(fc="white", ec="none", alpha=0.6, pad=1))

    ax.grid(alpha=0.2)
    if not is_last:
        ax.tick_params(axis="x", labelbottom=False)
    else:
        ax.set_xlabel("GPS time [s]")

Path(output_png).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_png, dpi=150, bbox_inches="tight")
print(f"Saved {output_png}")
plt.close(fig)
