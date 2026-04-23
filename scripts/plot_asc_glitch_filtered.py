#!/usr/bin/env python3
"""
Bandpass-filtered view of the 25-min glitch in ASC/LSC channels.
Filter: 30–55 Hz bandpass to isolate the ~40 Hz glitch frequency.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, hilbert
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

WORKDIR    = Path(__file__).parent.parent
GLITCH_GPS = 1415578745
GLITCH_GPS_EXACT = 1415578745.894531  # catalog GPS
CSV        = WORKDIR / "outputs" / f"asc_glitch_probe_{GLITCH_GPS}.csv"
OUT_DIR    = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH  = pd.Timestamp("1980-01-06", tz="UTC")

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(CSV).sort_values("gps").set_index("gps")
print(f"Loaded {df.shape[0]} samples × {df.shape[1]} channels")

def bandpass(series, fs, flo=30, fhi=55, order=4):
    """Bandpass filter a pandas Series."""
    v = series.fillna(0).values
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)

def get_fs(col, df):
    """Estimate sample rate from index spacing."""
    idx = df[col].dropna().index
    if len(idx) < 2:
        return 1.0
    diffs = np.diff(idx[:100])
    return round(1.0 / np.median(diffs))

def normalise(v):
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826
    return (v - med) / (mad if mad > 0 else 1.0)

# ── Key channels with their expected sample rates ─────────────────────────────
CHANNELS = [
    # (label, col, color, panel)
    ("DARM corr",  "V1:LSC_DARM_CORR",  "black",       "LSC"),
    ("DARM",       "V1:LSC_DARM",       "dimgray",     "LSC"),
    ("MICH err",   "V1:LSC_MICH_ERR",   "tab:blue",    "LSC"),
    ("PRCL err",   "V1:LSC_PRCL_ERR",   "tab:green",   "LSC"),
    ("SRCL err",   "V1:LSC_SRCL_ERR",   "tab:red",     "LSC"),
    ("CARM err",   "V1:LSC_CARM_ERR",   "tab:purple",  "LSC"),
    ("BS TX err",  "V1:ASC_BS_TX_ERR",  "tab:blue",    "ASC_ERR"),
    ("BS TY err",  "V1:ASC_BS_TY_ERR",  "tab:cyan",    "ASC_ERR"),
    ("PR TX err",  "V1:ASC_PR_TX_ERR",  "tab:orange",  "ASC_ERR"),
    ("PR TY err",  "V1:ASC_PR_TY_ERR",  "gold",        "ASC_ERR"),
    ("SR TX err",  "V1:ASC_SR_TX_ERR",  "tab:red",     "ASC_ERR"),
    ("SR TY err",  "V1:ASC_SR_TY_ERR",  "salmon",      "ASC_ERR"),
    ("NI TX corr", "V1:ASC_NI_TX_CORR", "tab:blue",    "ASC_CORR"),
    ("NI TY corr", "V1:ASC_NI_TY_CORR", "steelblue",   "ASC_CORR"),
    ("WI TX corr", "V1:ASC_WI_TX_CORR", "tab:orange",  "ASC_CORR"),
    ("WI TY corr", "V1:ASC_WI_TY_CORR", "goldenrod",   "ASC_CORR"),
    ("NE TX corr", "V1:ASC_NE_TX_CORR", "tab:green",   "ASC_CORR"),
    ("NE TY corr", "V1:ASC_NE_TY_CORR", "limegreen",   "ASC_CORR"),
    ("WE TX corr", "V1:ASC_WE_TX_CORR", "tab:red",     "ASC_CORR"),
    ("WE TY corr", "V1:ASC_WE_TY_CORR", "salmon",      "ASC_CORR"),
]

PANEL_TITLES = {
    "LSC":      "LSC length signals (DARM_CORR + ERR) — bandpassed 30–55 Hz",
    "ASC_ERR":  "ASC angular error signals (BS, PR, SR) — bandpassed 30–55 Hz",
    "ASC_CORR": "ASC angular correction output (NI, WI, NE, WE) — bandpassed 30–55 Hz",
}

# ── Build filtered traces ─────────────────────────────────────────────────────
# Use full 10s window for filtering, then trim to plot window
traces = {}
for label, col, color, panel in CHANNELS:
    if col not in df.columns:
        continue
    s = df[col].dropna()
    if len(s) < 100:
        continue
    fs = get_fs(col, df)
    filt = bandpass(s, fs)
    norm = normalise(filt)
    traces[col] = (s.index.values, norm, fs, panel, label, color)

print(f"Filtered {len(traces)} channels  (fs values: {set(int(v[2]) for v in traces.values())})")

# ── Use catalog GPS directly; also report DARM peak for info ─────────────────
glitch_t = GLITCH_GPS_EXACT
darm_col = "V1:LSC_DARM_CORR"
if darm_col in traces:
    t, v, fs, *_ = traces[darm_col]
    # Search only within ±2s of catalog GPS to avoid edge artifacts
    mask_search = (t >= glitch_t - 2.0) & (t <= glitch_t + 2.0)
    if mask_search.any():
        idx_local = np.argmax(np.abs(v[mask_search]))
        peak_gps = t[mask_search][idx_local]
        peak_val = v[mask_search][idx_local]
        print(f"DARM peak (±2s search): GPS {peak_gps:.4f}  (offset from catalog: {peak_gps - glitch_t:.3f} s)  value={peak_val:.1f} sigma")
    else:
        print(f"No DARM data within ±2s of catalog GPS — using catalog time")
print(f"Reference time: GPS {glitch_t:.4f} (catalog)")

# ── Plot for each window ──────────────────────────────────────────────────────
glitch_utc = GPS_EPOCH + pd.to_timedelta(glitch_t, unit="s")

PANELS_ORDER = ["LSC", "ASC_ERR", "ASC_CORR"]

# x-axis tick spacing per window: (window_s, tick_step_ms)
XTICK_STEP = {0.5: 100, 0.10: 10, 0.025: 5}

for win, suffix in [(0.5, "1s"), (0.10, "200ms"), (0.025, "50ms")]:
    t0 = glitch_t - win
    t1 = glitch_t + win
    use_ms = win <= 0.10   # plot in ms for zoomed windows

    fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.16, right=0.97)

    for ax, panel_key in zip(axes, PANELS_ORDER):
        offset = 0
        yticks, ylabels = [], []
        for label, col, color, panel in CHANNELS:
            if panel != panel_key or col not in traces:
                continue
            t_arr, v_arr, fs, *_ = traces[col]
            mask = (t_arr >= t0) & (t_arr <= t1)
            t_rel = t_arr[mask] - glitch_t
            t_plot = t_rel * 1e3 if use_ms else t_rel
            v_plot = v_arr[mask]
            if len(t_plot) < 2:
                continue
            ax.plot(t_plot, v_plot + offset, lw=0.9, color=color, alpha=0.9, label=label)
            yticks.append(offset)
            ylabels.append(label)
            offset += max(8, np.percentile(np.abs(v_plot), 99) * 2 + 2)

        # explicit x ticks at round ms (or s) values
        step = XTICK_STEP[win]
        win_plot = win * 1e3 if use_ms else win
        xt = np.arange(-win_plot, win_plot + step * 0.01, step)
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)

        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=7.5)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.9, label="catalog GPS")
        ax.legend(fontsize=6.5, loc="upper right", framealpha=0.7,
                  ncol=2, handlelength=1.2, handletextpad=0.4, borderpad=0.4)
        ax.set_title(PANEL_TITLES[panel_key], fontsize=9, loc="left")
        ax.grid(axis="x", ls=":", alpha=0.4)
        ax.set_ylabel("norm. units + offset", fontsize=7)

    unit = "ms" if use_ms else "s"
    axes[-1].set_xlabel(
        f"Time relative to catalog GPS {glitch_t:.4f}  ({unit})", fontsize=9)
    fig.suptitle(
        f"25-min glitch — BANDPASSED 30–55 Hz — GPS {GLITCH_GPS}  "
        f"({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})\n"
        f"Window ±{win*1e3:.0f} ms",
        fontsize=10, fontweight="bold", y=0.98)

    out = OUT_DIR / f"asc_glitch_filtered_{GLITCH_GPS}_{suffix}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")

# ── Hilbert envelope plot — ±1 s around glitch ────────────────────────────────
# Shows amplitude envelope of each channel to reveal WHICH mirror's motion
# spikes before DARM, without being confused by the 40 Hz carrier phase.
#
# Envelope channels selected for clarity: one per mirror / DOF group
ENV_CHANNELS = [
    # LSC
    ("DARM corr",  "V1:LSC_DARM_CORR",  "black",      "LSC"),
    ("MICH err",   "V1:LSC_MICH_ERR",   "tab:blue",   "LSC"),
    ("PRCL err",   "V1:LSC_PRCL_ERR",   "tab:green",  "LSC"),
    ("SRCL err",   "V1:LSC_SRCL_ERR",   "tab:red",    "LSC"),
    ("CARM err",   "V1:LSC_CARM_ERR",   "tab:purple", "LSC"),
    # ASC error (mirrors with direct sensing)
    ("SR TX err",  "V1:ASC_SR_TX_ERR",  "tab:red",    "ASC_ERR"),
    ("SR TY err",  "V1:ASC_SR_TY_ERR",  "salmon",     "ASC_ERR"),
    ("PR TX err",  "V1:ASC_PR_TX_ERR",  "tab:orange", "ASC_ERR"),
    ("PR TY err",  "V1:ASC_PR_TY_ERR",  "gold",       "ASC_ERR"),
    ("BS TX err",  "V1:ASC_BS_TX_ERR",  "tab:blue",   "ASC_ERR"),
    ("BS TY err",  "V1:ASC_BS_TY_ERR",  "tab:cyan",   "ASC_ERR"),
    # ASC correction (test masses, no direct ERR)
    ("WE TY corr", "V1:ASC_WE_TY_CORR", "tab:red",    "ASC_CORR"),
    ("WE TX corr", "V1:ASC_WE_TX_CORR", "salmon",     "ASC_CORR"),
    ("NE TY corr", "V1:ASC_NE_TY_CORR", "tab:green",  "ASC_CORR"),
    ("NE TX corr", "V1:ASC_NE_TX_CORR", "limegreen",  "ASC_CORR"),
    ("WI TY corr", "V1:ASC_WI_TY_CORR", "goldenrod",  "ASC_CORR"),
    ("WI TX corr", "V1:ASC_WI_TX_CORR", "tab:orange", "ASC_CORR"),
    ("NI TY corr", "V1:ASC_NI_TY_CORR", "steelblue",  "ASC_CORR"),
    ("NI TX corr", "V1:ASC_NI_TX_CORR", "tab:blue",   "ASC_CORR"),
]

ENV_WIN   = 1.0   # ±1 s window for envelope plot
ENV_SMOOTH_MS = 5  # smooth envelope with this many milliseconds rolling mean

fig_env, axes_env = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
fig_env.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.16, right=0.97)

env_t0 = glitch_t - ENV_WIN
env_t1 = glitch_t + ENV_WIN

for ax, panel_key in zip(axes_env, PANELS_ORDER):
    offset = 0
    yticks, ylabels = [], []
    for label, col, color, panel in ENV_CHANNELS:
        if panel != panel_key or col not in traces:
            continue
        t_arr, v_arr, fs, *_ = traces[col]
        mask = (t_arr >= env_t0) & (t_arr <= env_t1)
        t_seg = t_arr[mask] - glitch_t
        v_seg = v_arr[mask]
        if len(v_seg) < 4:
            continue
        # Hilbert envelope
        env = np.abs(hilbert(v_seg))
        # Smooth with rolling mean (convert ms to samples)
        n_smooth = max(1, int(round(ENV_SMOOTH_MS * 1e-3 * fs)))
        if n_smooth > 1:
            env = np.convolve(env, np.ones(n_smooth) / n_smooth, mode="same")
        # Normalise envelope so peak = 1, then offset
        env_peak = env.max()
        if env_peak > 0:
            env_norm = env / env_peak
        else:
            env_norm = env
        # Find envelope peak time (for annotation)
        peak_t = t_seg[np.argmax(env_norm)]
        t_ms = t_seg * 1e3
        ax.plot(t_ms, env_norm + offset, lw=1.2, color=color, alpha=0.9,
                label=f"{label}  ({peak_t*1e3:+.0f} ms)")
        ax.annotate(f"{peak_t*1e3:+.1f} ms", xy=(peak_t*1e3, env_norm.max() + offset),
                    xytext=(3, 2), textcoords="offset points",
                    fontsize=6, color=color, clip_on=True)
        yticks.append(offset + 0.5)
        ylabels.append(label)
        offset += 1.6

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7.5)
    ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.9, label="catalog GPS")
    ax.legend(fontsize=6.5, loc="upper right", framealpha=0.7,
              ncol=2, handlelength=1.2, handletextpad=0.4, borderpad=0.4)
    ax.set_title(PANEL_TITLES[panel_key].replace("bandpassed", "envelope |Hilbert|, smoothed")
                 .replace("30–55 Hz", f"30–55 Hz, {ENV_SMOOTH_MS} ms smooth"),
                 fontsize=9, loc="left")
    xt_env = np.arange(-ENV_WIN * 1e3, ENV_WIN * 1e3 + 100.01, 100)
    ax.set_xticks(xt_env)
    ax.set_xticklabels([f"{v:+.0f}" for v in xt_env], fontsize=8)
    ax.grid(axis="x", ls=":", alpha=0.4)
    ax.set_ylabel("norm. envelope + offset", fontsize=7)
    ax.set_ylim(bottom=0)

axes_env[-1].set_xlabel(f"Time relative to catalog GPS {glitch_t:.4f}  (ms)", fontsize=9)
fig_env.suptitle(
    f"25-min glitch — AMPLITUDE ENVELOPE (Hilbert, 30–55 Hz BP) — GPS {GLITCH_GPS}  "
    f"({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})\n"
    f"Window ±{ENV_WIN}s   |   Annotations = time of peak envelope relative to catalog GPS",
    fontsize=10, fontweight="bold", y=0.98)

out_env = OUT_DIR / f"asc_glitch_envelope_{GLITCH_GPS}.png"
fig_env.savefig(out_env, dpi=150, bbox_inches="tight")
plt.close(fig_env)
print(f"Saved → {out_env}")

# ── Overlay envelope plot — key channels on one panel ─────────────────────────
# Highlight DARM_CORR vs ASC angular channels, all normalised to peak=1,
# so timing relationships are immediately visible.
OVERLAY_CHANNELS = [
    ("DARM corr",  "V1:LSC_DARM_CORR",  "black",       2.5),
    ("SR TX err",  "V1:ASC_SR_TX_ERR",  "tab:red",     1.5),
    ("SR TY err",  "V1:ASC_SR_TY_ERR",  "salmon",      1.5),
    ("PR TX err",  "V1:ASC_PR_TX_ERR",  "tab:orange",  1.0),
    ("PR TY err",  "V1:ASC_PR_TY_ERR",  "gold",        1.0),
    ("BS TX err",  "V1:ASC_BS_TX_ERR",  "tab:blue",    1.0),
    ("WE TY corr", "V1:ASC_WE_TY_CORR", "tab:green",   1.0),
    ("NE TY corr", "V1:ASC_NE_TY_CORR", "limegreen",   1.0),
    ("NI TX corr", "V1:ASC_NI_TX_CORR", "steelblue",   1.0),
]

OVL_WIN = 0.5   # ±500 ms
ovl_t0  = glitch_t - OVL_WIN
ovl_t1  = glitch_t + OVL_WIN

fig_ov, ax_ov = plt.subplots(figsize=(16, 6))
fig_ov.subplots_adjust(top=0.88, bottom=0.10, left=0.07, right=0.97)

for label, col, color, lw in OVERLAY_CHANNELS:
    if col not in traces:
        continue
    t_arr, v_arr, fs, *_ = traces[col]
    mask = (t_arr >= ovl_t0) & (t_arr <= ovl_t1)
    t_seg = t_arr[mask] - glitch_t
    v_seg = v_arr[mask]
    if len(v_seg) < 4:
        continue
    env = np.abs(hilbert(v_seg))
    n_smooth = max(1, int(round(ENV_SMOOTH_MS * 1e-3 * fs)))
    if n_smooth > 1:
        env = np.convolve(env, np.ones(n_smooth) / n_smooth, mode="same")
    env_peak = env.max()
    env_norm = env / env_peak if env_peak > 0 else env
    peak_t_ms = t_seg[np.argmax(env_norm)] * 1e3
    ax_ov.plot(t_seg * 1e3, env_norm, lw=lw, color=color, alpha=0.85,
               label=f"{label}  (peak {peak_t_ms:+.0f} ms)")

ax_ov.axvline(0, color="crimson", lw=2, ls="--", alpha=0.9, label="catalog GPS")
ax_ov.set_xlabel("Time relative to catalog GPS (ms)", fontsize=10)
ax_ov.set_ylabel("Normalised amplitude envelope  (peak = 1)", fontsize=9)
ax_ov.set_xlim(-OVL_WIN * 1e3, OVL_WIN * 1e3)
xt_ov = np.arange(-OVL_WIN * 1e3, OVL_WIN * 1e3 + 50.01, 50)
ax_ov.set_xticks(xt_ov)
ax_ov.set_xticklabels([f"{v:+.0f}" for v in xt_ov], fontsize=9)
ax_ov.set_ylim(0, 1.15)
ax_ov.grid(ls=":", alpha=0.4)
ax_ov.legend(fontsize=8, ncol=2, loc="upper left")
fig_ov.suptitle(
    f"25-min glitch — ENVELOPE OVERLAY (Hilbert, 30–55 Hz BP, {ENV_SMOOTH_MS} ms smooth)\n"
    f"GPS {GLITCH_GPS}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})   Window ±{int(OVL_WIN*1e3)} ms",
    fontsize=10, fontweight="bold")

out_ov = OUT_DIR / f"asc_glitch_envelope_overlay_{GLITCH_GPS}.png"
fig_ov.savefig(out_ov, dpi=150, bbox_inches="tight")
plt.close(fig_ov)
print(f"Saved → {out_ov}")

# ── Peak-time timeline (lollipop chart) ───────────────────────────────────────
# For every channel in OVERLAY_CHANNELS, compute the Hilbert envelope peak
# time (ms) within ±500 ms of catalog GPS and plot as horizontal lollipop.
# Channels are sorted by peak time so lead/lag pattern is immediate.
TIMELINE_WIN = 0.5   # ±500 ms search window

peak_records = []
for label, col, color, lw in OVERLAY_CHANNELS:
    if col not in traces:
        continue
    t_arr, v_arr, fs, *_ = traces[col]
    mask = (t_arr >= glitch_t - TIMELINE_WIN) & (t_arr <= glitch_t + TIMELINE_WIN)
    t_seg = t_arr[mask] - glitch_t
    v_seg = v_arr[mask]
    if len(v_seg) < 4:
        continue
    env = np.abs(hilbert(v_seg))
    n_smooth = max(1, int(round(ENV_SMOOTH_MS * 1e-3 * fs)))
    if n_smooth > 1:
        env = np.convolve(env, np.ones(n_smooth) / n_smooth, mode="same")
    peak_t_ms = t_seg[np.argmax(env)] * 1e3
    # classify sharpness: DARM_CORR is "sharp", all ASC are "broad"
    sharp = (col == "V1:LSC_DARM_CORR")
    subsystem = ("LSC" if col.startswith("V1:LSC") else
                 "ASC ERR" if col.endswith("_ERR") else "ASC CORR")
    peak_records.append((peak_t_ms, label, col, color, sharp, subsystem))

# sort by peak time
peak_records.sort(key=lambda x: x[0])

fig_tl, ax_tl = plt.subplots(figsize=(11, max(5, len(peak_records) * 0.55 + 1.5)))
fig_tl.subplots_adjust(left=0.28, right=0.88, top=0.88, bottom=0.12)

# subsystem marker shapes
markers = {"LSC": "D", "ASC ERR": "o", "ASC CORR": "s"}
subsys_labels_added = set()

for i, (peak_ms, label, col, color, sharp, subsys) in enumerate(peak_records):
    mkr = markers[subsys]
    ms  = 12 if sharp else 9
    # stem line from t=0 to peak
    ax_tl.plot([0, peak_ms], [i, i], color=color, lw=1.5, alpha=0.7)
    # marker at peak
    leg_label = subsys if subsys not in subsys_labels_added else "_nolegend_"
    subsys_labels_added.add(subsys)
    ax_tl.scatter(peak_ms, i, color=color, s=ms**2, marker=mkr,
                  zorder=5, label=leg_label,
                  edgecolors="black" if sharp else "none", linewidths=1.2)
    # numeric annotation
    ax_tl.text(peak_ms + (8 if peak_ms >= 0 else -8), i,
               f"{peak_ms:+.0f} ms",
               va="center", ha="left" if peak_ms >= 0 else "right",
               fontsize=8, color=color)

# channel labels on y-axis
ax_tl.set_yticks(range(len(peak_records)))
ax_tl.set_yticklabels([r[1] for r in peak_records], fontsize=9)

# color y-tick labels to match trace color
for ytl, rec in zip(ax_tl.get_yticklabels(), peak_records):
    ytl.set_color(rec[3])

ax_tl.axvline(0, color="crimson", lw=2, ls="--", alpha=0.9, label="catalog GPS (t = 0)")
ax_tl.axvspan(-TIMELINE_WIN * 1e3, 0, alpha=0.04, color="steelblue")
ax_tl.axvspan(0, TIMELINE_WIN * 1e3, alpha=0.04, color="tomato")
ax_tl.text(-TIMELINE_WIN * 1e3 * 0.97, len(peak_records) - 0.3,
           "← leads DARM", fontsize=8, color="steelblue", ha="left")
ax_tl.text(TIMELINE_WIN * 1e3 * 0.97, len(peak_records) - 0.3,
           "lags DARM →", fontsize=8, color="tomato", ha="right")

xt_tl = np.arange(-500, 501, 100)
ax_tl.set_xticks(xt_tl)
ax_tl.set_xticklabels([f"{v:+.0f}" for v in xt_tl], fontsize=9)
ax_tl.set_xlim(-TIMELINE_WIN * 1e3 * 1.15, TIMELINE_WIN * 1e3 * 1.15)
ax_tl.set_xlabel("Hilbert envelope peak time relative to catalog GPS  (ms)", fontsize=10)
ax_tl.grid(axis="x", ls=":", alpha=0.4)
ax_tl.legend(fontsize=8, loc="lower right", framealpha=0.8,
             title="subsystem / shape", title_fontsize=7)

fig_tl.suptitle(
    f"25-min glitch — ENVELOPE PEAK TIMES per mirror/DOF\n"
    f"GPS {GLITCH_GPS}  ({glitch_utc.strftime('%Y-%m-%d %H:%M:%S UTC')})   "
    f"Hilbert |BP 30–55 Hz|, {ENV_SMOOTH_MS} ms smooth, ±{int(TIMELINE_WIN*1e3)} ms window\n"
    f"Diamond = LSC (sharp), circle = ASC ERR, square = ASC CORR  |  "
    f"Black edge = impulsive (DARM only)",
    fontsize=9, fontweight="bold")

out_tl = OUT_DIR / f"asc_glitch_peak_times_{GLITCH_GPS}.png"
fig_tl.savefig(out_tl, dpi=150, bbox_inches="tight")
plt.close(fig_tl)
print(f"Saved → {out_tl}")
