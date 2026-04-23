"""
plot_tower_thermal_timeseries.py — per-tower thermal channels + ASC_SR_TY_INPUT + glitches.

6 panels sharing the same x-axis:
  Panel 0: ASC_SR_TY_INPUT  (SR tilt proxy)
  Panel 1: Core temperatures (mir_coil, rh_te, bottom_te where available)
  Panel 2: TCS actuators     (rh_out clipped, rh_set, co2_pwrlas where available)
  Panel 3: Floor env temps   (F0, F4, F7 TE1/2)
  Panel 4: Glitch SNR
  Panel 5: Glitch Q-factor

Usage:
    python scripts/plot_tower_thermal_timeseries.py --tower WI
    python scripts/plot_tower_thermal_timeseries.py --tower SR --output usecases/...
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

# ── Per-tower channel config ─────────────────────────────────────────────────
# Each entry: list of (column_key, label, color)
TOWER_CONFIG = {
    "WI": {
        "core": [
            ("wi_bottom_te1",      "WI bottom TE1",       "tab:red"),
            ("wi_mir_coil_dr_te",  "WI mirror coil DR",   "tab:orange"),
            ("wi_rh_te",           "WI RH TE",            "tab:purple"),
        ],
        "tcs": [
            ("wi_rh_out",          "WI RH OUT (clip±800)", "tab:blue"),
            ("wi_rh_set",          "WI RH SET",            "tab:green"),
            ("wi_co2_pwrlas",      "WI CO2 PWRLAS",        "tab:red"),
        ],
        "floors": [
            ("wi_f0_te1", "F0 TE1", "tab:blue"),   ("wi_f0_te2", "F0 TE2", "tab:cyan"),
            ("wi_f4_te1", "F4 TE1", "tab:green"),  ("wi_f4_te2", "F4 TE2", "limegreen"),
            ("wi_f7_te1", "F7 TE1", "tab:brown"),  ("wi_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": True,
    },
    "WE": {
        "core": [
            ("we_mir_coil_ul_te",  "WE mirror coil UL",   "tab:orange"),
            ("we_rh_te",           "WE RH TE",            "tab:purple"),
        ],
        "tcs": [
            ("we_rh_set",          "WE RH SET",           "tab:green"),
        ],
        "floors": [
            ("we_f0_te1", "F0 TE1", "tab:blue"),   ("we_f0_te2", "F0 TE2", "tab:cyan"),
            ("we_f4_te1", "F4 TE1", "tab:green"),  ("we_f4_te2", "F4 TE2", "limegreen"),
            ("we_f7_te1", "F7 TE1", "tab:brown"),  ("we_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "NE": {
        "core": [
            ("ne_mir_coil_ul_te",  "NE mirror coil UL",   "tab:orange"),
            ("ne_rh_te",           "NE RH TE",            "tab:purple"),
        ],
        "tcs": [
            ("ne_rh_set",          "NE RH SET",           "tab:green"),
        ],
        "floors": [
            ("ne_f0_te1", "F0 TE1", "tab:blue"),   ("ne_f0_te2", "F0 TE2", "tab:cyan"),
            ("ne_f4_te1", "F4 TE1", "tab:green"),  ("ne_f4_te2", "F4 TE2", "limegreen"),
            ("ne_f7_te1", "F7 TE1", "tab:brown"),  ("ne_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "BS": {
        "core": [],
        "tcs":  [],
        "floors": [
            ("bs_f0_te1", "F0 TE1", "tab:blue"),   ("bs_f0_te2", "F0 TE2", "tab:cyan"),
            ("bs_f4_te1", "F4 TE1", "tab:green"),  ("bs_f4_te2", "F4 TE2", "limegreen"),
            ("bs_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "PR": {
        "core": [
            ("pr_mir_coil_ul_te",  "PR mirror coil UL",   "tab:orange"),
            ("pr_rh_te",           "PR RH TE",            "tab:purple"),
        ],
        "tcs":  [],
        "floors": [
            ("pr_f0_te1", "F0 TE1", "tab:blue"),   ("pr_f0_te2", "F0 TE2", "tab:cyan"),
            ("pr_f4_te1", "F4 TE1", "tab:green"),  ("pr_f4_te2", "F4 TE2", "limegreen"),
            ("pr_f7_te1", "F7 TE1", "tab:brown"),  ("pr_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "SR": {
        "core": [
            ("sr_mir_coil_ul_te",  "SR mirror coil UL",   "tab:orange"),
            ("sr_rh_te",           "SR RH TE",            "tab:purple"),
        ],
        "tcs": [
            ("sr_rh_set",          "SR RH SET",           "tab:green"),
        ],
        "floors": [
            ("sr_f0_te1", "F0 TE1", "tab:blue"),   ("sr_f0_te2", "F0 TE2", "tab:cyan"),
            ("sr_f4_te1", "F4 TE1", "tab:green"),  ("sr_f4_te2", "F4 TE2", "limegreen"),
            ("sr_f7_te1", "F7 TE1", "tab:brown"),  ("sr_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "IB": {
        "core": [
            ("ib_bench_te1", "IB bench TE1", "tab:red"),
            ("ib_bench_te2", "IB bench TE2", "tab:orange"),
            ("ib_box_te1",   "IB box TE1",   "tab:purple"),
            ("ib_box_te2",   "IB box TE2",   "tab:pink"),
        ],
        "tcs":  [],
        "floors": [
            ("ib_f0_te1", "F0 TE1", "tab:blue"),   ("ib_f0_te2", "F0 TE2", "tab:cyan"),
            ("ib_f4_te1", "F4 TE1", "tab:green"),  ("ib_f4_te2", "F4 TE2", "limegreen"),
            ("ib_f7_te1", "F7 TE1", "tab:brown"),  ("ib_f7_te2", "F7 TE2", "goldenrod"),
        ],
        "rh_clip": False,
    },
    "DET": {
        "core": [
            ("det_dlab_te",    "DET lab TE",     "tab:red"),
            ("det_sas_te_in",  "DET SAS TE in",  "tab:orange"),
            ("det_sas_te_out", "DET SAS TE out", "tab:purple"),
        ],
        "tcs": [
            ("det_te_in",  "DET TE in",  "tab:blue"),
            ("det_te_out", "DET TE out", "tab:cyan"),
            ("det_te_set", "DET TE set", "tab:green"),
        ],
        "floors": [],
        "rh_clip": False,
    },
}


def gps_to_dt(gps: np.ndarray) -> np.ndarray:
    return np.array(
        [GPS_EPOCH + pd.to_timedelta(float(g), unit="s") for g in gps],
        dtype="datetime64[ms]",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tower",    required=True, choices=list(TOWER_CONFIG))
    ap.add_argument("--tw-csv",   default="outputs/tower_thermal_merged.csv")
    ap.add_argument("--sr-csv",   default="outputs/sr_ty_input_merged.csv")
    ap.add_argument("--glitches", default="data/full_25min_glitches_ER16-O4b.csv")
    ap.add_argument("--output",   default=None)
    args = ap.parse_args()

    tower = args.tower.upper()
    cfg   = TOWER_CONFIG[tower]
    out_path = args.output or f"usecases/25-minute-glitch/{tower.lower()}_thermal_timeseries.png"

    # ── Load ────────────────────────────────────────────────────────────────
    tw = pd.read_csv(args.tw_csv).sort_values("gps_bin").reset_index(drop=True)
    sr = pd.read_csv(args.sr_csv).sort_values("gps_bin").reset_index(drop=True)
    sr["asc_sr_ty_input"] = sr["asc_sr_ty_input"].ffill().bfill()
    gl = pd.read_csv(args.glitches).sort_values("time").reset_index(drop=True)

    tw_dt  = gps_to_dt(tw["gps_bin"].values)
    sr_dt  = gps_to_dt(sr["gps_bin"].values)
    tx     = sr["asc_sr_ty_input"].values
    gl_dt  = gps_to_dt(gl["time"].values)
    gl_snr = gl["amplitude"].values
    gl_q   = gl["q"].values

    last_gl_gps = gl["time"].max()
    last_gl_dt  = np.datetime64(
        GPS_EPOCH + pd.to_timedelta(float(last_gl_gps), unit="s"), "ms")

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

    # Decide which panels to draw
    has_core   = len(cfg["core"]) > 0
    has_tcs    = len(cfg["tcs"]) > 0
    has_floors = len(cfg["floors"]) > 0

    n_panels = 1 + has_core + has_tcs + has_floors + 2  # SR_TY + core + tcs + floors + SNR + Q
    ratios = ([1.2]
              + ([1.2] if has_core   else [])
              + ([1.2] if has_tcs    else [])
              + ([1.2] if has_floors else [])
              + [1.0, 1.0])

    fig, axes = plt.subplots(n_panels, 1, figsize=(18, 3 * n_panels), sharex=True,
                             gridspec_kw={"hspace": 0.07, "height_ratios": ratios})
    fig.suptitle(f"{tower} tower thermal channels vs ASC_SR_TY_INPUT — full O4 (2023–2026)",
                 fontsize=11, y=0.995)

    panel = iter(axes)

    # ── Panel 0: ASC_SR_TY_INPUT ─────────────────────────────────────────
    ax = next(panel)
    ax.plot(sr_dt, tx, color="steelblue", lw=0.7, label="ASC_SR_TY_INPUT")
    ax.axhline(190, color="darkorange", lw=0.9, ls="--", alpha=0.7, label="LN3 (190 Hz)")
    ax.axhline(0,   color="darkgreen",  lw=0.9, ls="--", alpha=0.7, label="LN3_ALIGNED (0 Hz)")
    ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
               label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax)
    ax.set_ylabel("SR TY INPUT\n[Hz DCP]", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, loc="lower left", ncol=4)
    ax.grid(alpha=0.25)

    # ── Panel: core temperatures ─────────────────────────────────────────
    if has_core:
        ax = next(panel)
        for col, lbl, clr in cfg["core"]:
            if col in tw.columns:
                ax.plot(tw_dt, tw[col].values, color=clr, lw=0.7, label=lbl)
        ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
        add_lock_spans(ax)
        ax.set_ylabel("Temperature [°C]", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper left", ncol=3)
        ax.grid(alpha=0.25)

    # ── Panel: TCS actuators ─────────────────────────────────────────────
    if has_tcs:
        ax = next(panel)
        for col, lbl, clr in cfg["tcs"]:
            if col not in tw.columns:
                continue
            vals = tw[col].values
            if cfg["rh_clip"] and "rh_out" in col:
                vals = np.clip(vals, -800, 800)
                lbl  = lbl + " (clip±800)"
            ax.plot(tw_dt, vals, color=clr, lw=0.6, label=lbl)
        ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
        add_lock_spans(ax)
        ylabel = "TCS actuators [arb/W]"
        if tower == "DET":
            ylabel = "Temperature [°C]"
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper left", ncol=3)
        ax.grid(alpha=0.25)

    # ── Panel: floor env temperatures ───────────────────────────────────
    if has_floors:
        ax = next(panel)
        for col, lbl, clr in cfg["floors"]:
            if col in tw.columns:
                ax.plot(tw_dt, tw[col].values, color=clr, lw=0.6, label=lbl)
        ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
        add_lock_spans(ax)
        ax.set_ylabel("Floor env temp [°C]", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc="upper left", ncol=3)
        ax.grid(alpha=0.25)

    # ── Panel: glitch SNR ────────────────────────────────────────────────
    ax = next(panel)
    ax.plot(gl_dt, gl_snr, ",", color="tab:red", alpha=0.6,
            rasterized=True, label="amplitude per glitch")
    ax.set_yscale("log")
    ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--",
               label=f"Last glitch {str(last_gl_dt)[:10]}")
    add_lock_spans(ax)
    ax.set_ylabel("glitch amplitude\n(strain)", fontsize=8)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.25)

    # ── Panel: glitch Q ──────────────────────────────────────────────────
    ax = next(panel)
    ax.plot(gl_dt, gl_q, ",", color="tab:green", alpha=0.4,
            rasterized=True, label="Q per glitch")
    ax.axvline(last_gl_dt, color="crimson", lw=1.2, ls="--")
    add_lock_spans(ax)
    ax.set_ylabel("Q-factor", fontsize=8)
    ax.tick_params(axis="x", rotation=20, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.xaxis.set_major_formatter(xfmt)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.25)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
