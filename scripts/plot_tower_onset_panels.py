#!/usr/bin/env python3
"""
Tower onset panels — 9 panels, one per tower.

Each panel shows all available ERR and CORR channels for that tower
(bandpass 20–150 Hz, normalised to pre-glitch RMS), with V1:Hrec_hoft_16384Hz
superimposed on a twin y-axis (right side, same σ units).

Goal: identify which tower (mirror) reacts first relative to h(t).

Input:
    outputs/tower_onset_probe_{GLITCH_GPS}.csv
    (produced by slurm/nh_tower_onset_probe.slurm)

Output (usecases/25-minute-glitch/):
    tower_onset_panels_{GLITCH_GPS}_{suffix}.png   — three zoom levels
"""
import warnings; warnings.filterwarnings("ignore")
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from scipy.signal import hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ──────────────────────────────────────────────────────────────
WORKDIR          = Path(__file__).parent.parent
GLITCH_GPS_EXACT = 1446106739.098   # 2025-11-02 08:18:59 UTC, SNR=499, f=40.4 Hz
GPS_EPOCH        = pd.Timestamp("1980-01-06", tz="UTC")
OUT_DIR          = WORKDIR / "usecases" / "25-minute-glitch"

CSV = WORKDIR / "outputs" / f"tower_onset_probe_{int(GLITCH_GPS_EXACT)}.csv"

FLO, FHI = 20, 150    # bandpass Hz

# Quiet baseline well before the glitch (within the ±60 s extraction window)
BASELINE_T0 = -55.0
BASELINE_T1 = -45.0

# ── Tower channel definitions ──────────────────────────────────────────────────
# Format: (tower_label, [(label, channel, color, linestyle), ...])
# ERR channels → orange/red palette  |  CORR channels → green/teal palette
TOWERS = [
    ("NI", [
        # ASC_NI_TX/TY_CORR = zero in O4c; actual drive is via Sc_ suspension at 10 kHz
        ("NI TX corr", "V1:Sc_NI_MAR_TX_CORR",  "#2ca02c", "-"),
        ("NI TY corr", "V1:Sc_NI_MAR_TY_CORR",  "#98df8a", "--"),
    ]),
    ("WI", [
        ("WI TX corr", "V1:Sc_WI_MAR_TX_CORR",  "#2ca02c", "-"),
        ("WI TY corr", "V1:Sc_WI_MAR_TY_CORR",  "#98df8a", "--"),
    ]),
    ("NE", [
        ("NE TX corr", "V1:ASC_NE_TX_CORR",  "#2ca02c", "-"),
        ("NE TY corr", "V1:ASC_NE_TY_CORR",  "#98df8a", "--"),
    ]),
    ("WE", [
        ("WE TX corr", "V1:ASC_WE_TX_CORR",  "#2ca02c", "-"),
        ("WE TY corr", "V1:ASC_WE_TY_CORR",  "#98df8a", "--"),
    ]),
    ("BS", [
        ("BS TX err",  "V1:ASC_BS_TX_ERR",   "#d62728", "-"),
        ("BS TY err",  "V1:ASC_BS_TY_ERR",   "#ff9896", "--"),
        ("BS TX corr", "V1:ASC_BS_TX_CORR",  "#2ca02c", "-"),
        ("BS TY corr", "V1:ASC_BS_TY_CORR",  "#98df8a", "--"),
    ]),
    ("PR", [
        ("PR TX err",  "V1:ASC_PR_TX_ERR",   "#d62728", "-"),
        ("PR TY err",  "V1:ASC_PR_TY_ERR",   "#ff9896", "--"),
        ("PR TX corr", "V1:ASC_PR_TX_CORR",  "#2ca02c", "-"),
        ("PR TY corr", "V1:ASC_PR_TY_CORR",  "#98df8a", "--"),
    ]),
    ("SR", [
        # B1p_DCP = raw B1 dark-port demodulated signal, 2000 Hz.
        # This is the optical sensor input to the SR servo BEFORE loop filtering —
        # the most direct measurement of SR angular motion available in raw frames.
        # Step-2 ranked channel ASC_SR_TY_DCP_mag_B1_mag is the 10 Hz envelope of this.
        ("SR TX B1p",  "V1:ASC_SR_TX_B1p_DCP",  "#d62728", "-"),
        ("SR TY B1p",  "V1:ASC_SR_TY_B1p_DCP",  "#ff9896", "--"),
        # ERR = servo error (after input matrix), OUT = actuator correction
        ("SR TX err",  "V1:ASC_SR_TX_ERR",       "#8c564b", "-"),
        ("SR TY err",  "V1:ASC_SR_TY_ERR",       "#c49c94", "--"),
        ("SR TX out",  "V1:ASC_SR_TX_OUT",        "#17becf", "-"),
        ("SR TY out",  "V1:ASC_SR_TY_OUT",        "#9edae5", "--"),
    ]),
    ("IB", [
        ("IB TX corr", "V1:Sc_IB_MAR_TX_CORR",  "#2ca02c", "-"),
        ("IB TY corr", "V1:Sc_IB_MAR_TY_CORR",  "#98df8a", "--"),
    ]),
    ("DET/DIFFp", [
        # DET (OB = Output Bench): no Sc_DET_MAR_* in raw frames.
        # Sa_OB_F0 channels are the suspended bench F0-stage corrections at 500 Hz.
        # DIFFp = differential angular DOF (NE−WE combination); ranked #19 in Step 2.
        # Both shown here: OB as DET proxy, DIFFp as DARM-coupled angular DOF.
        ("OB TY corr",    "V1:Sa_OB_F0_TY_CORR_500Hz",  "#bcbd22", "-"),
        ("DIFFp TX err",  "V1:ASC_DIFFp_TX_ERR",         "#d62728", "-"),
        ("DIFFp TY err",  "V1:ASC_DIFFp_TY_ERR",         "#ff9896", "--"),
        ("DIFFp TX corr", "V1:ASC_DIFFp_TX_CORR",        "#2ca02c", "-"),
        ("DIFFp TY corr", "V1:ASC_DIFFp_TY_CORR",        "#98df8a", "--"),
    ]),
]

HREF_CH   = "V1:Hrec_hoft_16384Hz"
DARM_CORR = "V1:LSC_DARM_CORR"

# Step-2 top channels to overlay in the 40ms onset plot (10th panel)
STEP2_CHANNELS = [
    ("DARM IMC line",   "V1:LSC_DARM_IMC_LINE_mag_100Hz",  "#1f77b4", "-"),   # rank 1
    ("DARM PSTAB0",     "V1:LSC_DARM_PSTAB0_COUPLING_100Hz","#ff7f0e", "-"),  # rank 2
    ("DARM SIB1 line",  "V1:LSC_DARM_SIB1_LINE_mag_100Hz", "#2ca02c", "-"),   # rank 6
    ("WI FF50 P err",   "V1:Sc_WI_FF50HZ_P_ERR",           "#9467bd", "-"),   # rank 7
    ("DARM corr",       "V1:LSC_DARM_CORR",                 "dimgray",  "--"), # already extracted
]

# ── Load data ─────────────────────────────────────────────────────────────────
if not CSV.exists():
    print(f"ERROR: {CSV} not found. Run slurm/nh_tower_onset_probe.slurm first.")
    sys.exit(1)

print(f"Loading {CSV} ...", flush=True)
df = pd.read_csv(CSV).sort_values("gps").set_index("gps")
t_gps = df.index.values
t_rel = t_gps - GLITCH_GPS_EXACT
print(f"Window: t_rel = {t_rel.min():.1f} to {t_rel.max():.1f} s  "
      f"({len(df)} samples, {len(df.columns)} channels)", flush=True)

glitch_utc = GPS_EPOCH + pd.to_timedelta(GLITCH_GPS_EXACT, unit="s")

# ── Signal processing helpers ──────────────────────────────────────────────────
def get_fs(col):
    idx = df[col].dropna().index.values
    if len(idx) < 10:
        return 4096
    return round(1.0 / np.median(np.diff(idx[:200])))

def bandpass_norm(col):
    """Bandpass filter + normalise by pre-glitch baseline RMS.
    Returns (t_ch, normalised_trace) or None if channel absent.

    Works on the non-NaN rows only so the filter sees the true sample rate,
    regardless of how many other channels are merged into the DataFrame.
    For demodulated amplitude channels (fs <= 2*FHI), skip the bandpass and
    just normalise the raw signal."""
    if col not in df.columns:
        return None
    s   = df[col].dropna()
    t_ch = s.index.values - GLITCH_GPS_EXACT   # per-channel relative times
    v    = s.values.astype(float)
    fs   = round(1.0 / np.median(np.diff(t_ch[:200]))) if len(t_ch) >= 10 else 4096
    if fs > 2 * FHI:
        sos = butter(4, [FLO, FHI], btype="band", fs=fs, output="sos")
        bp  = sosfiltfilt(sos, v)
    else:
        bp = v   # demodulated amplitude channel — show as-is
    mask_bl = (t_ch >= BASELINE_T0) & (t_ch <= BASELINE_T1)
    rms = np.std(bp[mask_bl]) if mask_bl.sum() > 0 else np.std(bp)
    norm = bp / rms if rms > 0 else bp
    return t_ch, norm

# Preprocess all channels once
cache = {}
all_cols = (
    [ch for _, chlist in TOWERS for _, ch, _, _ in chlist]
    + [HREF_CH]
    + [ch for _, ch, _, _ in STEP2_CHANNELS]
)
for col in set(all_cols):
    r = bandpass_norm(col)
    if r is not None:
        cache[col] = r
        print(f"  {col:45s}  fs={get_fs(col)} Hz  "
              f"baseline_rms={np.std(cache[col][1][(cache[col][0]>=BASELINE_T0)&(cache[col][0]<=BASELINE_T1)]):.2f}σ",
              flush=True)
    else:
        print(f"  {col:45s}  NOT FOUND", flush=True)

# ── Glitch onset detection ────────────────────────────────────────────────────
# Find the precise glitch onset in h(t): locate peak within ±5 s, then walk
# backwards to find when the Hilbert envelope first exceeds ONSET_THRESH × σ.
ONSET_THRESH = 3.0    # σ threshold for onset detection
ONSET_SEARCH = 2.0    # s: search window half-width around catalog GPS

# ── Hilbert envelope peak timing ─────────────────────────────────────────────
# Use the PEAK of the Hilbert envelope (unambiguous strong feature),
# not a threshold crossing. Search within ±PEAK_SEARCH s of catalog GPS.
PEAK_SEARCH = 5.0   # s

def hilbert_peak_time(t, norm, win_lo=-PEAK_SEARCH, win_hi=PEAK_SEARCH):
    """Time of maximum Hilbert envelope within [win_lo, win_hi]."""
    mask = (t >= win_lo) & (t <= win_hi)
    if mask.sum() < 4:
        return np.nan, np.nan
    env = np.abs(hilbert(norm[mask]))
    idx = np.argmax(env)
    return t[mask][idx], env[idx]   # (peak_time, peak_amplitude_in_σ)

# Find h(t) peak first — all mirror timings are reported relative to this
T_HPEAK = np.nan
AMP_HPEAK = np.nan
if HREF_CH in cache:
    t_arr, norm = cache[HREF_CH]
    T_HPEAK, AMP_HPEAK = hilbert_peak_time(t_arr, norm)

print(f"\nh(t) Hilbert peak: t_rel = {T_HPEAK*1e3:+.1f} ms  "
      f"amplitude = {AMP_HPEAK:.1f} σ  (relative to catalog GPS)")
print(f"\nChannel peak times — Δt = t_peak(channel) − t_peak(h(t)):")
print(f"  {'Channel':45s}  {'peak [ms rel GPS]':>18}  {'Δt vs h(t) [ms]':>17}  {'peak amp [σ]':>13}")
print(f"  {'-'*45}  {'-'*18}  {'-'*17}  {'-'*13}")

peak_table = {}    # col → (t_peak_abs, delta_vs_h, amplitude)
for tower, chlist in TOWERS:
    for label, col, _, _ in chlist:
        if col not in cache:
            continue
        t_arr, norm = cache[col]
        tp, amp = hilbert_peak_time(t_arr, norm)
        dt = (tp - T_HPEAK) * 1e3 if not (np.isnan(tp) or np.isnan(T_HPEAK)) else np.nan
        peak_table[col] = (tp, dt, amp)
        print(f"  {col:45s}  {tp*1e3:+18.1f}  {dt:+17.1f}  {amp:13.1f}")

peak_table[HREF_CH] = (T_HPEAK, 0.0, AMP_HPEAK)
print(f"  {HREF_CH:45s}  {T_HPEAK*1e3:+18.1f}  {'0.0 (ref)':>17}  {AMP_HPEAK:13.1f}")

# Centre tight zoom plots on h(t) peak
T_CENTER = T_HPEAK if not np.isnan(T_HPEAK) else 0.0

# ── Plotting ──────────────────────────────────────────────────────────────────
OUT_DIR.mkdir(parents=True, exist_ok=True)

def make_figure(win_lo, win_hi, suffix, xtick_step):
    """
    3×3 grid of subplots, one per tower.
    Left y-axis: tower ERR+CORR (normalised σ).
    Right y-axis (twin): h(t) (normalised σ), black thick line.
    """
    fig, axes = plt.subplots(3, 3, figsize=(18, 14),
                             sharex=True, sharey=False)
    fig.subplots_adjust(hspace=0.35, wspace=0.35,
                        top=0.90, bottom=0.07, left=0.07, right=0.95)

    ax2_list = []   # collect right axes to enforce uniform h(t) y-scale

    for idx, (tower, chlist) in enumerate(TOWERS):
        row, col = divmod(idx, 3)
        ax  = axes[row, col]
        ax2 = ax.twinx()   # h(t) on right axis
        ax2_list.append(ax2)

        ax.set_facecolor("#f9f9f9")

        # ── Plot tower channels on left axis ──────────────────────────────
        any_data = False
        for label, ch, color, ls in chlist:
            if ch not in cache:
                continue
            t_arr, norm = cache[ch]
            mask = (t_arr >= win_lo) & (t_arr <= win_hi)
            if mask.sum() == 0:
                continue
            # Full channel name + Δt vs h(t) in legend
            tp, dt, amp = peak_table.get(ch, (np.nan, np.nan, np.nan))
            dt_str = f"  Δt={dt:+.0f}ms" if not np.isnan(dt) else ""
            leg_label = f"{ch}{dt_str}"
            ax.plot(t_arr[mask], norm[mask],
                    color=color, lw=1.2, ls=ls, alpha=0.85, label=leg_label)
            # Mark peak with a thin vertical line (only if peak is in window)
            if not np.isnan(tp) and win_lo <= tp <= win_hi:
                ax.axvline(tp, color=color, lw=0.8, ls=":", alpha=0.7, zorder=8)
            any_data = True

        if not any_data:
            ax.text(0.5, 0.5, "no channels\nin raw frame",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color="gray")

        # ── Plot h(t) on right twin axis ──────────────────────────────────
        if HREF_CH in cache:
            t_arr, norm = cache[HREF_CH]
            mask = (t_arr >= win_lo) & (t_arr <= win_hi)
            if mask.sum() > 0:
                ax2.plot(t_arr[mask], norm[mask],
                         color="black", lw=1.8, alpha=0.6,
                         label="h(t)", zorder=5)
                ax2.tick_params(axis="y", labelsize=6, colors="gray")
                ax2.set_ylabel("h(t) [σ]", fontsize=6, color="gray",
                               labelpad=2)

        # ── Reference lines ───────────────────────────────────────────────
        # Catalog GPS marker
        if win_lo <= 0 <= win_hi:
            ax.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8, zorder=10)
        # h(t) peak marker
        if not np.isnan(T_HPEAK) and win_lo <= T_HPEAK <= win_hi:
            ax.axvline(T_HPEAK, color="black", lw=1.2, ls="--", alpha=0.7, zorder=10)
        # Baseline window marker (only if overlaps current window)
        bl_lo = max(BASELINE_T0, win_lo)
        bl_hi = min(BASELINE_T1, win_hi)
        if bl_lo < bl_hi:
            ax.axvspan(bl_lo, bl_hi, color="lightblue", alpha=0.20, zorder=0)

        # ── Cosmetics ──────────────────────────────────────────────────────
        ax.set_title(tower, fontsize=11, fontweight="bold", pad=3)
        ax.set_ylabel("tower [σ]", fontsize=7, labelpad=2)
        ax.tick_params(axis="both", labelsize=7)
        ax.grid(axis="x", ls=":", alpha=0.4)
        ax.axhline(0, color="gray", lw=0.4, alpha=0.4)

        if any_data:
            ax.legend(fontsize=6, loc="upper left",
                      framealpha=0.6, ncol=2, handlelength=1.5)

        if row == 2:
            use_ms = (win_hi - win_lo) <= 2.0
            xlabel = ("Time relative to catalog GPS (ms)" if use_ms
                      else "Time relative to catalog GPS (s)")
            ax.set_xlabel(xlabel, fontsize=8)

    # Uniform h(t) y-scale across all right axes
    if HREF_CH in cache:
        t_arr, norm = cache[HREF_CH]
        mask = (t_arr >= win_lo) & (t_arr <= win_hi)
        if mask.sum() > 0:
            h_max = np.nanmax(np.abs(norm[mask])) * 1.15
            for a2 in ax2_list:
                a2.set_ylim(-h_max, h_max)

    # X-ticks — ms labels for sub-2 s windows
    use_ms = (win_hi - win_lo) <= 2.0
    xt = np.arange(np.ceil(win_lo / xtick_step) * xtick_step,
                   win_hi + xtick_step * 0.01, xtick_step)
    for ax_row in axes:
        for ax in ax_row:
            ax.set_xlim(win_lo, win_hi)   # enforce limits — axvline(0) can blow them out
            ax.set_xticks(xt)
            if use_ms:
                ax.set_xticklabels([f"{v*1e3:+.0f}" for v in xt], fontsize=7)
            else:
                ax.set_xticklabels([f"{v:+.0f}" if v != 0 else "0" for v in xt],
                                    fontsize=7)

    hpeak_str = f"  |  h(t) peak {T_HPEAK*1e3:+.1f} ms rel GPS  ({AMP_HPEAK:.0f}σ)" if not np.isnan(T_HPEAK) else ""
    fig.suptitle(
        f"Tower onset panels — ERR/CORR per mirror + h(t) overlay\n"
        f"GPS {GLITCH_GPS_EXACT:.3f}  "
        f"({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})  "
        f"SNR≈499  f≈40 Hz\n"
        f"Bandpass {FLO}–{FHI} Hz | normalised to baseline [{BASELINE_T0:.0f}, {BASELINE_T1:.0f}] s | "
        f"crimson = catalog GPS | black dashed = h(t) peak | Δt in legend = channel peak − h(t) peak{hpeak_str}",
        fontsize=9, fontweight="bold"
    )

    out = OUT_DIR / f"tower_onset_panels_{int(GLITCH_GPS_EXACT)}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

# Broad context
make_figure(-60.0, 60.0, "120s",  xtick_step=10)
make_figure( -5.0,  5.0, "10s",   xtick_step=1)

# Tight ms-scale zooms centred on h(t) onset — few cycles of 40 Hz signal
# T_CENTER is in seconds relative to catalog GPS
# ±50 ms = 4 cycles at 40 Hz  |  ±20 ms = 1.6 cycles — individual waveforms legible
make_figure(T_CENTER - 0.20, T_CENTER + 0.20, "onset_400ms", xtick_step=0.050)
make_figure(T_CENTER - 0.06, T_CENTER + 0.06, "onset_120ms", xtick_step=0.010)
make_figure(T_CENTER - 0.05, T_CENTER + 0.05, "onset_100ms", xtick_step=0.010)
make_figure(T_CENTER - 0.02, T_CENTER + 0.02, "onset_40ms",  xtick_step=0.005)


def make_figure_step2(win_lo, win_hi, suffix, xtick_step):
    """
    10-panel variant: 9 tower panels (3×3) + 1 Step-2 panel (bottom-left of a 4th row).
    The Step-2 panel overlays top ranked channels from NOISEHOUND Step 2 + h(t).
    Only generated for the tight ms-scale onset plots.
    """
    towers_ext = TOWERS + [("Step 2 / LSC", STEP2_CHANNELS)]
    n_panels   = len(towers_ext)          # 10
    ncols, nrows = 3, 4                   # 12 cells, last 2 empty

    fig, axes_grid = plt.subplots(nrows, ncols, figsize=(18, 18),
                                  sharex=True, sharey=False)
    fig.subplots_adjust(hspace=0.35, wspace=0.35,
                        top=0.90, bottom=0.06, left=0.07, right=0.95)

    ax2_list = []

    use_ms = (win_hi - win_lo) <= 2.0

    for idx in range(nrows * ncols):
        row, col = divmod(idx, ncols)
        ax = axes_grid[row, col]

        if idx >= n_panels:
            ax.set_visible(False)
            continue

        tower, chlist = towers_ext[idx]
        ax2_list.append(ax.twinx())
        ax2 = ax2_list[-1]
        ax.set_facecolor("#f9f9f9")

        any_data = False
        for label, ch, color, ls in chlist:
            if ch not in cache:
                continue
            t_arr, norm = cache[ch]
            mask = (t_arr >= win_lo) & (t_arr <= win_hi)
            if mask.sum() == 0:
                continue
            tp, dt, amp = peak_table.get(ch, (np.nan, np.nan, np.nan))
            dt_str = f"  Δt={dt:+.0f}ms" if not np.isnan(dt) else ""
            ax.plot(t_arr[mask], norm[mask],
                    color=color, lw=1.2, ls=ls, alpha=0.85, label=f"{ch}{dt_str}")
            if not np.isnan(tp) and win_lo <= tp <= win_hi:
                ax.axvline(tp, color=color, lw=0.8, ls=":", alpha=0.7, zorder=8)
            any_data = True

        if not any_data:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8, color="gray")

        if HREF_CH in cache:
            t_arr, norm = cache[HREF_CH]
            mask = (t_arr >= win_lo) & (t_arr <= win_hi)
            if mask.sum() > 0:
                ax2.plot(t_arr[mask], norm[mask],
                         color="black", lw=1.8, alpha=0.6, label="h(t)", zorder=5)
                ax2.tick_params(axis="y", labelsize=6, colors="gray")
                ax2.set_ylabel("h(t) [σ]", fontsize=6, color="gray", labelpad=2)

        if win_lo <= 0 <= win_hi:
            ax.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8, zorder=10)
        if not np.isnan(T_HPEAK) and win_lo <= T_HPEAK <= win_hi:
            ax.axvline(T_HPEAK, color="black", lw=1.2, ls="--", alpha=0.7, zorder=10)

        ax.set_title(tower, fontsize=11, fontweight="bold", pad=3)
        ax.set_ylabel("channel [σ]", fontsize=7, labelpad=2)
        ax.tick_params(axis="both", labelsize=7)
        ax.grid(axis="x", ls=":", alpha=0.4)
        ax.axhline(0, color="gray", lw=0.4, alpha=0.4)
        ax.set_xlim(win_lo, win_hi)

        if any_data:
            ax.legend(fontsize=6, loc="upper left",
                      framealpha=0.6, ncol=2, handlelength=1.5)

        if row == nrows - 1 or idx == n_panels - 1:
            xlabel = "Time relative to catalog GPS (ms)" if use_ms else "Time relative to catalog GPS (s)"
            ax.set_xlabel(xlabel, fontsize=8)

    # Uniform h(t) y-scale across all right axes
    if HREF_CH in cache:
        t_arr, norm = cache[HREF_CH]
        mask = (t_arr >= win_lo) & (t_arr <= win_hi)
        if mask.sum() > 0:
            h_max = np.nanmax(np.abs(norm[mask])) * 1.15
            for a2 in ax2_list:
                a2.set_ylim(-h_max, h_max)

    xt = np.arange(np.ceil(win_lo / xtick_step) * xtick_step,
                   win_hi + xtick_step * 0.01, xtick_step)
    for ax_row in axes_grid:
        for ax in ax_row:
            if not ax.get_visible():
                continue
            ax.set_xlim(win_lo, win_hi)
            ax.set_xticks(xt)
            if use_ms:
                ax.set_xticklabels([f"{v*1e3:+.0f}" for v in xt], fontsize=7)
            else:
                ax.set_xticklabels([f"{v:+.0f}" if v != 0 else "0" for v in xt],
                                    fontsize=7)

    hpeak_str = f"  |  h(t) peak {T_HPEAK*1e3:+.1f} ms rel GPS  ({AMP_HPEAK:.0f}σ)" if not np.isnan(T_HPEAK) else ""
    fig.suptitle(
        f"Tower onset panels + Step-2 top channels — ERR/CORR per mirror + h(t) overlay\n"
        f"GPS {GLITCH_GPS_EXACT:.3f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})  SNR≈499  f≈40 Hz\n"
        f"Bandpass {FLO}–{FHI} Hz | normalised to baseline [{BASELINE_T0:.0f}, {BASELINE_T1:.0f}] s | "
        f"black dashed = h(t) peak | Δt in legend = channel peak − h(t) peak{hpeak_str}",
        fontsize=9, fontweight="bold"
    )

    out = OUT_DIR / f"tower_onset_step2_{int(GLITCH_GPS_EXACT)}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


make_figure_step2(T_CENTER - 0.20, T_CENTER + 0.20, "onset_400ms", xtick_step=0.050)
make_figure_step2(T_CENTER - 0.02, T_CENTER + 0.02, "onset_40ms",  xtick_step=0.005)
