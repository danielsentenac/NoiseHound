#!/usr/bin/env python3
"""
Glitch response plot: did the ASC corrections react to the DARM glitch?

Panel 1 — DARM_CORR bandpassed 30–55 Hz waveform (see the burst)
Panel 2 — DARM_CORR Hilbert amplitude envelope (glitch energy vs time)
Panel 3 — ASC CORR Hilbert envelopes normalised to pre-glitch baseline
           Flat = 1.0 → no response.  Spike → correction reacted.

All panels share the same time axis.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR          = Path(__file__).parent.parent
GLITCH_GPS       = 1415578745
GLITCH_GPS_EXACT = 1415578745.894531   # Omicron peak GPS
CSV              = WORKDIR / "outputs" / f"asc_glitch_probe_{GLITCH_GPS}.csv"
OUT_DIR          = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH        = pd.Timestamp("1980-01-06", tz="UTC")

# Windows to produce (seconds before/after glitch)
WINDOWS = [4.0, 1.0, 0.20]   # 8 s overview, 2 s zoom, 400 ms zoom

# Hilbert envelope smoothing
ENV_SMOOTH_MS = 10   # ms

# Baseline window for normalising envelopes (relative to glitch)
# Use t = -5 s to -3 s — well before glitch activity starts
BASELINE_START = -5.0   # s
BASELINE_END   = -3.0   # s

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV).sort_values("gps").set_index("gps")
print(f"Loaded {df.shape[0]} samples × {df.shape[1]} channels")

# ── Helpers ───────────────────────────────────────────────────────────────────
def bandpass(v, fs, flo=30, fhi=55, order=4):
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)

def envelope(v, fs, smooth_ms=ENV_SMOOTH_MS):
    env = np.abs(hilbert(v))
    n = max(1, int(round(smooth_ms * 1e-3 * fs)))
    if n > 1:
        env = np.convolve(env, np.ones(n) / n, mode="same")
    return env

def get_fs(col):
    idx = df[col].dropna().index
    diffs = np.diff(idx[:200])
    return round(1.0 / np.median(diffs))

# ── Extract and filter all needed channels ─────────────────────────────────
DARM_COL = "V1:LSC_DARM_CORR"

ASC_TOWERS = [
    ("NI", [
        ("V1:ASC_NI_TX_CORR", "tab:blue"),
        ("V1:ASC_NI_TY_CORR", "steelblue"),
    ]),
    ("WI", [
        ("V1:ASC_WI_TX_CORR", "tab:orange"),
        ("V1:ASC_WI_TY_CORR", "goldenrod"),
    ]),
    ("NE", [
        ("V1:ASC_NE_TX_CORR", "tab:green"),
        ("V1:ASC_NE_TY_CORR", "limegreen"),
    ]),
    ("WE", [
        ("V1:ASC_WE_TX_CORR", "tab:red"),
        ("V1:ASC_WE_TY_CORR", "salmon"),
    ]),
]

t_all = df.index.values
glitch_t = GLITCH_GPS_EXACT
t_rel_all = t_all - glitch_t

# ── DARM preparation ─────────────────────────────────────────────────────────
if DARM_COL not in df.columns:
    raise RuntimeError(f"Column {DARM_COL} not in CSV")

darm_raw  = df[DARM_COL].fillna(0).values.astype(float)
darm_fs   = get_fs(DARM_COL)
darm_bp   = bandpass(darm_raw, darm_fs)
darm_env  = envelope(darm_bp, darm_fs)

# Normalise DARM env to its own pre-glitch baseline
mask_bl = (t_rel_all >= BASELINE_START) & (t_rel_all <= BASELINE_END)
darm_baseline = darm_env[mask_bl].mean() if mask_bl.any() else 1.0
darm_env_norm = darm_env / darm_baseline
# Also normalise the waveform for display (divide by pre-glitch RMS)
darm_rms_bl = np.std(darm_bp[mask_bl]) if mask_bl.any() else 1.0
darm_bp_norm = darm_bp / darm_rms_bl

# ── ASC CORR envelopes (normalised to baseline) ───────────────────────────────
asc_towers_traces = []   # list of (tower_label, [(col, color, t_arr, env_norm), ...])
for tower_label, chans in ASC_TOWERS:
    tower_traces = []
    for col, color in chans:
        if col not in df.columns:
            print(f"  Missing: {col}")
            continue
        v = df[col].fillna(0).values.astype(float)
        fs = get_fs(col)
        bp = bandpass(v, fs)
        env = envelope(bp, fs)
        bl_mask = (t_rel_all >= BASELINE_START) & (t_rel_all <= BASELINE_END)
        bl_mean = env[bl_mask].mean() if bl_mask.any() else 1.0
        env_norm = env / bl_mean if bl_mean > 0 else env
        tower_traces.append((col, color, t_rel_all, env_norm))
    asc_towers_traces.append((tower_label, tower_traces))

# ── Plotting ──────────────────────────────────────────────────────────────────
for win in WINDOWS:
    use_ms = win <= 0.5
    scale  = 1e3 if use_ms else 1.0
    unit   = "ms" if use_ms else "s"
    suffix = f"{int(win*2*1000)}ms" if use_ms else f"{int(win*2)}s"

    mask_win = (t_rel_all >= -win) & (t_rel_all <= win)
    t_plot   = t_rel_all[mask_win] * scale

    n_rows = 2 + len(asc_towers_traces)   # DARM waveform + DARM envelope + one per tower
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 3 * n_rows), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.13, right=0.97)

    # ── Panel 1: DARM waveform ────────────────────────────────────────────────
    ax1 = axes[0]
    ax1.plot(t_plot, darm_bp_norm[mask_win],
             color="black", lw=0.8, alpha=0.9, label="V1:LSC_DARM_CORR (30–55 Hz)")
    ax1.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax1.set_ylabel("Amplitude\n(pre-glitch RMS)", fontsize=9)
    ax1.set_title("V1:LSC_DARM_CORR — bandpassed 30–55 Hz waveform", fontsize=9, loc="left")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(axis="x", ls=":", alpha=0.4)

    # ── Panel 2: DARM envelope ────────────────────────────────────────────────
    ax2 = axes[1]
    ax2.plot(t_plot, darm_env_norm[mask_win],
             color="black", lw=1.5, alpha=0.9, label="V1:LSC_DARM_CORR envelope")
    ax2.axhline(1.0, color="black", lw=0.8, ls=":", alpha=0.5, label="baseline = 1")
    ax2.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8, label="catalog GPS")
    ax2.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
    ax2.set_title("V1:LSC_DARM_CORR — Hilbert envelope  (1.0 = pre-glitch level)", fontsize=9, loc="left")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(axis="x", ls=":", alpha=0.4)

    # ── One panel per tower ───────────────────────────────────────────────────
    for i, (tower_label, tower_traces) in enumerate(asc_towers_traces):
        ax = axes[2 + i]
        for col, color, t_arr, env_norm in tower_traces:
            asc_mask = (t_arr >= -win) & (t_arr <= win)
            ax.plot(t_arr[asc_mask] * scale, env_norm[asc_mask],
                    color=color, lw=1.2, alpha=0.9, label=col)
        ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
        ax.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
        ax.set_title(f"{tower_label} — ASC CORR Hilbert envelope", fontsize=9, loc="left")
        ax.legend(fontsize=8, loc="upper right", handlelength=1.5)
        ax.grid(axis="x", ls=":", alpha=0.4)

    # ── Shared x-axis ─────────────────────────────────────────────────────────
    xtick_step = {4.0: 1000, 1.0: 200, 0.20: 50}.get(win, 100)
    xt = np.arange(-win * scale, win * scale + xtick_step * 0.01, xtick_step)
    for ax in axes:
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
    axes[-1].set_xlabel(
        f"Time relative to catalog GPS {glitch_t:.4f}  ({unit})", fontsize=9)

    glitch_utc = GPS_EPOCH + pd.to_timedelta(glitch_t, unit="s")
    fig.suptitle(
        f"25-min glitch — did DARM glitch drive ASC correction responses?\n"
        f"GPS {glitch_t:.4f}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})  "
        f"Window: ±{win*1e3:.0f} ms   Bandpass 30–55 Hz   Envelope smoothed {ENV_SMOOTH_MS} ms",
        fontsize=10, fontweight="bold")

    out = OUT_DIR / f"glitch_response_{GLITCH_GPS}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")
