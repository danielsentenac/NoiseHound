"""
plot_ni_thermal_timeseries.py — NI tower thermal channels + ASC_SR_TY_INPUT + glitches, full O4.

7 panels sharing the same x-axis:
  Panel 0: ASC_SR_TY_INPUT  (SR tilt proxy, [Hz DCP])
  Panel 1: NI temperatures  — ni_bottom_te1, ni_mir_coil_ul_te, ni_rh_te
  Panel 2: Ring heater      — ni_rh_out (clipped ±800), ni_rh_set
  Panel 3: CO2 laser        — ni_co2_pwrlas, ni_co2_pwrin
  Panel 4: NI TX angular    — ni_tx_mag
  Panel 5: Glitch SNR       (one point per glitch)
  Panel 6: Glitch Q-factor  (one point per glitch)

Usage:
    python scripts/plot_ni_thermal_timeseries.py \\
        --ni-csv   outputs/ni_thermal_merged.csv \\
        --sr-csv   outputs/sr_ty_input_merged.csv \\
        --glitches data/full_25min_glitches_ER16-O4b.csv \\
        --output   usecases/25-minute-glitch/ni_thermal_timeseries.png
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ni-csv",   default="outputs/ni_thermal_merged.csv")
    ap.add_argument("--sr-csv",   default="outputs/sr_ty_input_merged.csv")
    ap.add_argument("--glitches", default="data/full_25min_glitches_ER16-O4b.csv")
    ap.add_argument("--output",   default="usecases/25-minute-glitch/ni_thermal_timeseries.png")
    args = ap.parse_args()

    # ── Load ────────────────────────────────────────────────────────────────
    ni = pd.read_csv(args.ni_csv).sort_values("gps_bin").reset_index(drop=True)
    sr = pd.read_csv(args.sr_csv).sort_values("gps_bin").reset_index(drop=True)
    sr["asc_sr_ty_input"] = sr["asc_sr_ty_input"].ffill().bfill()
    gl = pd.read_csv(args.glitches).sort_values("time").reset_index(drop=True)

    ni_dt  = gps_to_dt(ni["gps_bin"].values)
    sr_dt  = gps_to_dt(sr["gps_bin"].values)
    tx     = sr["asc_sr_ty_input"].values
    gl_dt  = gps_to_dt(gl["time"].values)
    gl_snr = gl["amplitude"].values
    gl_q   = gl["q"].values

    last_gl_gps = gl["time"].max()
    last_gl_dt  = np.datetime64(
        GPS_EPOCH + pd.to_timedelta(float(last_gl_gps), unit="s"), "ms")

    # Clip ring heater output to remove extreme outliers
    rh_out = ni["ni_rh_out"].clip(-800, 800).values

    xfmt = mdates.DateFormatter("%Y-%m")

    # ── Lock spans ──────────────────────────────────────────────────────────
    lock_csv = Path("outputs/itf_lock_full_o4.csv")
    if lock_csv.exists():
        lk = pd.read_csv(lock_csv).sort_values("gps_bin")
        lk_dt = gps_to_dt(lk["gps_bin"].values)
        lm = lk["lock_mean"].values
        ln3_mask         = (lm >= 134) & (lm <= 135)
        ln3_aligned_mask = (lm >= 144) & (lm <= 145)
        ln3_sqz_mask     = (lm >= 149) & (lm <= 150)
        ln2_mask         = (lm >= 129) & (lm <= 130)
        ln1_mask         = (lm >= 124) & (lm <= 125)
        dc_lock_mask     = (lm >= 104) & (lm < 124)
        arms_lock_mask   = (lm >= 39)  & (lm < 104)
        drmi_mask        = (lm >  1)   & (lm < 39)
    else:
        lk_dt = np.array([])
        ln3_mask = ln3_aligned_mask = ln3_sqz_mask = np.array([])
        ln2_mask = ln1_mask = dc_lock_mask = arms_lock_mask = drmi_mask = np.array([])

    _cmap = plt.get_cmap("turbo")
    c = [_cmap(v) for v in np.linspace(0.08, 0.95, 8)]

    def add_lock_spans(ax):
        if len(lk_dt) == 0:
            return
        for mask, color in [
            (drmi_mask,        c[0]),
            (arms_lock_mask,   c[1]),
            (dc_lock_mask,     c[2]),
            (ln1_mask,         c[3]),
            (ln2_mask,         c[4]),
            (ln3_mask,         c[5]),
            (ln3_sqz_mask,     c[6]),
            (ln3_aligned_mask, c[7]),
        ]:
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

    # ── Figure ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(7, 1, figsize=(18, 18), sharex=True,
                             gridspec_kw={"hspace": 0.07,
                                          "height_ratios": [1.2, 1.4, 1.4, 1, 1, 1, 1]})
    fig.suptitle("NI tower thermal channels vs ASC_SR_TY_INPUT — full O4 (2023–2026)",
                 fontsize=11, y=0.995)

    # ── Panel 0: ASC_SR_TY_INPUT ─────────────────────────────────────────
    ax0 = axes[0]
    ax0.plot(sr_dt, tx, color="steelblue", lw=0.7, label="ASC_SR_TY_INPUT")
    ax0.axhline(190, color="darkorange", lw=0.9, ls="--", alpha=0.7,
                label="LN3 setpoint (190 Hz)")
    ax0.axhline(0, color="darkgreen", lw=0.9, ls="--", alpha=0.7,
                label="LN3_ALIGNED (0 Hz)")
    ax0.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
                label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax0)
    ax0.set_ylabel("SR TY INPUT\n[Hz DCP]", fontsize=8)
    ax0.tick_params(labelsize=8)
    ax0.legend(fontsize=7, loc="lower left", ncol=3)
    ax0.grid(alpha=0.25)

    # ── Panel 1: NI temperatures ─────────────────────────────────────────
    ax1 = axes[1]
    ax1.plot(ni_dt, ni["ni_bottom_te1"].values,      color="tab:red",    lw=0.7,
             label="NI bottom TE1 (tower base)")
    ax1.plot(ni_dt, ni["ni_mir_coil_ul_te"].values,  color="tab:orange", lw=0.7,
             label="NI mirror coil UL TE")
    ax1.plot(ni_dt, ni["ni_rh_te"].values,           color="tab:purple", lw=0.7,
             label="NI RH TE")
    ax1.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax1)
    ax1.set_ylabel("Temperature [°C]", fontsize=8)
    ax1.tick_params(labelsize=8)
    ax1.legend(fontsize=7, loc="upper left", ncol=3)
    ax1.grid(alpha=0.25)

    # ── Panel 2: Ring heater ──────────────────────────────────────────────
    ax2 = axes[2]
    ax2.plot(ni_dt, rh_out,                    color="tab:blue",  lw=0.5,
             label="NI RH OUT (clipped ±800)")
    ax2.plot(ni_dt, ni["ni_rh_set"].values,    color="tab:green", lw=0.8,
             label="NI RH SET")
    ax2.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax2)
    ax2.set_ylabel("Ring heater [arb]", fontsize=8)
    ax2.tick_params(labelsize=8)
    ax2.legend(fontsize=7, loc="upper left", ncol=2)
    ax2.grid(alpha=0.25)

    # ── Panel 3: CO2 laser ────────────────────────────────────────────────
    ax3 = axes[3]
    ax3.plot(ni_dt, ni["ni_co2_pwrlas"].values, color="tab:red",    lw=0.7,
             label="CO2 laser power (PWRLAS)")
    ax3.plot(ni_dt, ni["ni_co2_pwrin"].values,  color="tab:orange", lw=0.7, alpha=0.7,
             label="CO2 input power (PWRIN)")
    ax3.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax3)
    ax3.set_ylabel("CO2 power [W]", fontsize=8)
    ax3.tick_params(labelsize=8)
    ax3.legend(fontsize=7, loc="upper left", ncol=2)
    ax3.grid(alpha=0.25)

    # ── Panel 4: Floor environment temperatures ───────────────────────────
    ax4 = axes[4]
    ax4.plot(ni_dt, ni["ni_f0_te1"].values, color="tab:blue",   lw=0.6, label="F0 TE1")
    ax4.plot(ni_dt, ni["ni_f0_te2"].values, color="tab:cyan",   lw=0.6, label="F0 TE2")
    ax4.plot(ni_dt, ni["ni_f4_te1"].values, color="tab:green",  lw=0.6, label="F4 TE1")
    ax4.plot(ni_dt, ni["ni_f4_te2"].values, color="limegreen",  lw=0.6, label="F4 TE2")
    ax4.plot(ni_dt, ni["ni_f7_te1"].values, color="tab:brown",  lw=0.6, label="F7 TE1")
    ax4.plot(ni_dt, ni["ni_f7_te2"].values, color="goldenrod",  lw=0.6, label="F7 TE2")
    ax4.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax4)
    ax4.set_ylabel("Floor env temp [°C]", fontsize=8)
    ax4.tick_params(labelsize=8)
    ax4.legend(fontsize=7, loc="upper left", ncol=3)
    ax4.grid(alpha=0.25)

    # ── Panel 5: Glitch SNR ───────────────────────────────────────────────
    ax5 = axes[5]
    ax5.plot(gl_dt, gl_snr, ",", color="tab:red", alpha=0.6,
             rasterized=True, label="amplitude per glitch")
    ax5.set_yscale("log")
    ax5.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
                label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax5)
    ax5.set_ylabel("glitch amplitude\n(strain)", fontsize=8)
    ax5.tick_params(labelsize=8)
    ax5.legend(fontsize=7, loc="upper left")
    ax5.grid(alpha=0.25)

    # ── Panel 6: Glitch Q-factor ──────────────────────────────────────────
    ax6 = axes[6]
    ax6.plot(gl_dt, gl_q, ",", color="tab:green", alpha=0.4,
             rasterized=True, label="Q per glitch")
    ax6.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax6)
    ax6.set_ylabel("Q-factor", fontsize=8)
    ax6.tick_params(axis="x", rotation=20, labelsize=8)
    ax6.tick_params(axis="y", labelsize=8)
    ax6.xaxis.set_major_formatter(xfmt)
    ax6.legend(fontsize=7, loc="upper left")
    ax6.grid(alpha=0.25)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
