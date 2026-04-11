"""
plot_sr_ty_input_timeseries.py — ASC_SR_TY_INPUT (SR tilt proxy) + glitch properties, full O4.

ASC_SR_TY_INPUT is the DCP-based SR TY tilt monitor: the ITF_LOCK automation
servos SR TY to hold this channel at 190 Hz (LN3, tilted) or 0 Hz (LN3_ALIGNED).
It is the ground-truth SR angular position in the optical basis.

3 panels sharing the same x-axis (full O4, 2023-2026):
  Panel 1:  ASC_SR_TY_INPUT  (SR tilt proxy, [Hz DCP])
  Panel 2:  Glitch SNR per event
  Panel 3:  Glitch Q-factor per event

Usage:
    python scripts/plot_sr_ty_input_timeseries.py \\
        --sr-csv    outputs/sr_ty_input_merged.csv \\
        --glitches  data/full_25min_glitches_ER16-O4b.csv \\
        --output    usecases/25-minute-glitch/sr_ty_input_timeseries.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")


def gps_to_dt(gps: np.ndarray) -> np.ndarray:
    return np.array(
        [GPS_EPOCH + pd.to_timedelta(float(g), unit="s") for g in gps],
        dtype="datetime64[ms]",
    )


def detect_steps(tx: np.ndarray, gps: np.ndarray,
                 threshold: float = 20.0, min_gap_h: int = 48) -> list[dict]:
    """Detect persistent level changes in ASC_SR_TY_INPUT (sliding-window median comparison)."""
    steps = []
    half = 6
    last_step_idx = -min_gap_h
    for i in range(half, len(tx) - half):
        before = np.nanmedian(tx[max(0, i - half):i])
        after  = np.nanmedian(tx[i:i + half])
        delta  = after - before
        if abs(delta) > threshold and (i - last_step_idx) > min_gap_h:
            steps.append({"gps": int(gps[i]), "delta_tx": delta})
            last_step_idx = i
    return steps




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sr-csv",   default="outputs/sr_ty_input_merged.csv")
    ap.add_argument("--glitches", default="data/full_25min_glitches_ER16-O4b.csv")
    ap.add_argument("--output",   default="usecases/25-minute-glitch/sr_ty_input_timeseries.png")
    args = ap.parse_args()

    # ── Load ────────────────────────────────────────────────────────────────
    sr = pd.read_csv(args.sr_csv).sort_values("gps_bin").reset_index(drop=True)
    sr["asc_sr_ty_input"] = sr["asc_sr_ty_input"].ffill().bfill()
    gl = pd.read_csv(args.glitches).sort_values("time").reset_index(drop=True)

    sr_dt = gps_to_dt(sr["gps_bin"].values)
    tx    = sr["asc_sr_ty_input"].values

    # Per-glitch time, SNR, Q
    gl_dt = gps_to_dt(gl["time"].values)
    gl_snr = gl["amplitude"].values
    gl_q   = gl["q"].values

    last_gl_gps = gl["time"].max()
    last_gl_dt  = np.datetime64(
        GPS_EPOCH + pd.to_timedelta(float(last_gl_gps), unit="s"), "ms")

    # Detect MAR TX step events
    steps = detect_steps(tx, sr["gps_bin"].values)
    step_dts   = gps_to_dt([s["gps"] for s in steps])
    step_deltas = np.array([s["delta_tx"] for s in steps])
    print(f"Detected {len(steps)} TX steps")

    xfmt = mdates.DateFormatter("%Y-%m")

    # ── LN3 locked spans from lock_frac + TX gate ───────────────────────────
    # Full O4 lock index (merged from quarterly extraction + step7).
    lock_csv = Path("outputs/itf_lock_full_o4.csv")
    if lock_csv.exists():
        lk = pd.read_csv(lock_csv).sort_values("gps_bin")
        lk_dt = gps_to_dt(lk["gps_bin"].values)
        lm = lk["lock_mean"].values
        # ITF_LOCK state machine (from metatron_nodes.json):
        #   1–38   : pre-arm lock (DRMI, PRMI, MICH…)
        #   39–103 : arms locked, CARM/DARM acquisition
        #   104–120: CARM_NULL, OMC, DC readout
        #   124–125: LOW_NOISE_1
        #   129–130: LOW_NOISE_2
        #   134–135: LOW_NOISE_3
        #   144–145: LOW_NOISE_3_ALIGNED
        #   149–150: LOW_NOISE_3_SQZ
        ln3_mask         = (lm >= 134) & (lm <= 135)
        ln3_aligned_mask = (lm >= 144) & (lm <= 145)
        ln3_sqz_mask     = (lm >= 149) & (lm <= 150)
        ln2_mask         = (lm >= 129) & (lm <= 130)
        ln1_mask         = (lm >= 124) & (lm <= 125)
        dc_lock_mask     = (lm >= 104) & (lm < 124)   # DC readout / OMC
        arms_lock_mask   = (lm >= 39)  & (lm < 104)   # arms locked
        drmi_mask        = (lm >  1)   & (lm < 39)    # DRMI/PRMI/pre-arms
    else:
        lk_dt = np.array([])
        ln3_mask = ln3_aligned_mask = ln3_sqz_mask = np.array([])
        ln2_mask = ln1_mask = dc_lock_mask = arms_lock_mask = drmi_mask = np.array([])

    # 8 lock levels cold→hot using turbo (max perceptual separation, blue→red)
    _cmap = plt.get_cmap("turbo")
    c = [_cmap(v) for v in np.linspace(0.08, 0.95, 8)]

    def add_lock_spans(ax):
        """Draw lock-state background spans."""
        if len(lk_dt) == 0:
            return
        for mask, color in [
            (drmi_mask,        c[0]),   # DRMI/PRMI — cold blue
            (arms_lock_mask,   c[1]),   # arms locked
            (dc_lock_mask,     c[2]),   # DC readout / OMC
            (ln1_mask,         c[3]),   # LN1
            (ln2_mask,         c[4]),   # LN2
            (ln3_mask,         c[5]),   # LN3
            (ln3_sqz_mask,     c[6]),   # LN3+SQZ
            (ln3_aligned_mask, c[7]),   # LN3_ALIGNED — hot red
        ]:
            # group consecutive True bins into spans
            in_span = False
            for i, val in enumerate(mask):
                if val and not in_span:
                    span_start = lk_dt[i]
                    in_span = True
                elif not val and in_span:
                    ax.axvspan(span_start, lk_dt[i], color=color,
                               alpha=0.45, zorder=0, lw=0)
                    in_span = False
            if in_span:
                ax.axvspan(span_start, lk_dt[-1], color=color,
                           alpha=0.45, zorder=0, lw=0)

    # ── Figure: 3 panels, shared x ──────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True,
                             gridspec_kw={"hspace": 0.08, "height_ratios": [1.4, 1, 1]})
    fig.suptitle("ASC_SR_TY_INPUT (SR tilt proxy) + glitch properties — full O4 (2023–2026)",
                 fontsize=11, y=0.995)

    # ── Panel 0: MAR TX ─────────────────────────────────────────────────────
    ax0 = axes[0]
    ax0.plot(sr_dt, tx, color="steelblue", lw=0.8, label="ASC_SR_TY_INPUT")
    ax0.axhline(190, color="darkorange", lw=0.9, ls="--", alpha=0.7,
                label="LN3 tilt setpoint (190 Hz)")
    ax0.axhline(0, color="darkgreen", lw=0.9, ls="--", alpha=0.7,
                label="LN3_ALIGNED setpoint (0 Hz)")
    ax0.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
                label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax0)
    ax0.set_ylabel("ASC_SR_TY_INPUT [Hz DCP]", fontsize=9)
    ax0.tick_params(labelsize=8)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    step_handles = [
        Patch(fc=c[0], alpha=0.7, label="DRMI/PRMI"),
        Patch(fc=c[1], alpha=0.7, label="Arms locked"),
        Patch(fc=c[2], alpha=0.7, label="DC readout / OMC"),
        Patch(fc=c[3], alpha=0.7, label="LOW_NOISE_1"),
        Patch(fc=c[4], alpha=0.7, label="LOW_NOISE_2"),
        Patch(fc=c[5], alpha=0.7, label="LOW_NOISE_3"),
        Patch(fc=c[6], alpha=0.7, label="LOW_NOISE_3 + SQZ"),
        Patch(fc=c[7], alpha=0.7, label="LOW_NOISE_3_ALIGNED"),
    ]
    handles0, labels0 = ax0.get_legend_handles_labels()
    ax0.legend(handles=handles0 + step_handles, fontsize=7, loc="lower left")
    ax0.grid(alpha=0.25)

    # ── Panel 1: glitch SNR (one point per glitch) ───────────────────────────
    ax1 = axes[1]
    ax1.plot(gl_dt, gl_snr, ',', color="tab:red", alpha=0.6,
             rasterized=True, label="amplitude per glitch")
    ax1.set_yscale("log")
    ax1.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
                label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax1)
    ax1.set_ylabel("glitch amplitude\n(strain)", fontsize=9)
    ax1.tick_params(labelsize=8)
    ax1.legend(fontsize=7, loc="upper left")
    ax1.grid(alpha=0.25)

    # ── Panel 2: glitch Q-factor (one point per glitch) ──────────────────────
    ax2 = axes[2]
    ax2.plot(gl_dt, gl_q, ',', color="tab:green", alpha=0.4,
             rasterized=True, label="Q per glitch")
    ax2.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
                label="Last glitch")
    add_lock_spans(ax2)
    ax2.set_ylabel("Q-factor", fontsize=9)
    ax2.tick_params(axis="x", rotation=20, labelsize=8)
    ax2.grid(alpha=0.25)
    ax2.xaxis.set_major_formatter(xfmt)
    ax2.legend(fontsize=7, loc="upper left")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
