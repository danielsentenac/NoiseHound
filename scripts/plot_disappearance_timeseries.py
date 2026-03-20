"""
plot_disappearance_timeseries.py — Step 7 multiplot.

Plots the glitch rate and key thermal / electrical channels over the
disappearance epoch (Oct 2025 – Apr 2026), with:
  - Yellow band: ITF shutdown (Christmas 2025-12-12 – 2026-01-08)
  - Orange band: general power outage (2026-02-28, ~10-day recovery)
  - ITF lock fraction panel (META_ITF_LOCK_index) replacing recurrence

Key GPS markers:
  - LAST_TRIGGER_GPS  = 1449582668  (2025-12-12 13:51 UTC) — yellow band start
  - ITF_DOWN_END      = 1451908330  (2026-01-08)           — yellow band end
  - POWER_OUTAGE_GPS  = 1456288198  (2026-02-28 04:30 UTC) — orange band start
  - POWER_RECOVERY    = POWER_OUTAGE_GPS + 10 days         — orange band end

Usage (run on CCA after jobs complete):
    python scripts/plot_disappearance_timeseries.py \\
        --binned      outputs/rate_correlation/binned_summary.csv \\
        --step4-dir   outputs/rate_correlation_step4 \\
        --triggers    data/full_25min_glitches_ER16-O4b.csv \\
        --itf-lock    outputs/itf_lock_step7.csv \\
        --output      usecases/25-minute-glitch/disappearance_timeseries.png
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

EPOCH_START_GPS  = 1443657600   # 2025-10-01
EPOCH_END_GPS    = 1459468800   # 2026-04-01

LAST_TRIGGER_GPS = 1449582668   # 2025-12-12 13:51 UTC — yellow zone starts HERE
ITF_DOWN_END     = 1451908330   # 2026-01-08

POWER_OUTAGE_GPS = 1456288198   # 2026-02-28 04:30 UTC
POWER_RECOVERY   = POWER_OUTAGE_GPS + 10 * 86400   # ~2026-03-10


def gps_to_dt(gps):
    return GPS_EPOCH + pd.to_timedelta(np.asarray(gps, float), unit="s")


def load_data(binned_csv: str, step4_dir: str) -> pd.DataFrame:
    """Merge step4 partials (full coverage, 4 channels) with binned_summary
    (all channels, partial coverage) for the Oct 2025 – Apr 2026 epoch."""
    # Step4: full coverage, 4 channels
    files = sorted(glob.glob(f"{step4_dir}/partial_*.csv"))
    df4 = (pd.concat([pd.read_csv(f) for f in files])
           .sort_values("gps_bin").drop_duplicates("gps_bin"))
    df4 = df4[(df4.gps_bin >= EPOCH_START_GPS) & (df4.gps_bin < EPOCH_END_GPS)].reset_index(drop=True)

    # binned_summary: all channels, partial coverage — merge extra columns
    dfs = pd.read_csv(binned_csv)
    dfs = dfs[(dfs.gps_bin >= EPOCH_START_GPS) & (dfs.gps_bin < EPOCH_END_GPS)]
    extra = [c for c in dfs.columns if c not in df4.columns]
    if extra:
        df4 = df4.merge(dfs[["gps_bin"] + extra], on="gps_bin", how="left")

    print(f"Bins: {len(df4)}, with triggers: {(df4.n_triggers > 0).sum()}")
    return df4


def load_lock(itf_lock_csv: str) -> pd.DataFrame | None:
    if not itf_lock_csv or not Path(itf_lock_csv).exists():
        print(f"  ITF lock CSV not found: {itf_lock_csv}")
        return None
    df = pd.read_csv(itf_lock_csv)
    df = df[(df.gps_bin >= EPOCH_START_GPS) & (df.gps_bin < EPOCH_END_GPS)]
    print(f"Lock bins: {len(df)}")
    return df


def plot(binned_csv: str, step4_dir: str, triggers_csv: str,
         itf_lock_csv: str, output: str) -> None:

    df   = load_data(binned_csv, step4_dir)
    lock = load_lock(itf_lock_csv)
    dt   = gps_to_dt(df["gps_bin"].values)

    # Recurrence (kept for fallback if lock data unavailable)
    trig  = pd.read_csv(triggers_csv)
    t_all = np.sort(trig["time"].values)
    keep  = np.concatenate([[True], np.diff(t_all) > 300])
    t     = t_all[keep]
    t_ep  = t[(t >= EPOCH_START_GPS) & (t < EPOCH_END_GPS)]
    intervals = np.diff(t_ep)
    t_mid     = t_ep[:-1] + intervals / 2
    smooth    = (pd.Series(intervals)
                 .rolling(20, center=True, min_periods=5)
                 .median().bfill().ffill().values)

    last_trig_dt    = gps_to_dt(LAST_TRIGGER_GPS)
    itf_end_dt      = gps_to_dt(ITF_DOWN_END)
    outage_dt       = gps_to_dt(POWER_OUTAGE_GPS)
    recovery_dt     = gps_to_dt(POWER_RECOVERY)

    xfmt = DateFormatter("%b %Y")

    def decorate(ax, legend=True):
        """Add both event bands + red dashed line to an axis."""
        ax.axvspan(last_trig_dt, itf_end_dt,
                   color="gold", alpha=0.35, zorder=0,
                   label="SR tower opening (Christmas)")
        ax.axvspan(outage_dt, recovery_dt,
                   color="tomato", alpha=0.25, zorder=0,
                   label="Power outage (~10-day recovery)")
        ax.axvline(last_trig_dt, color="crimson", lw=1.2, ls="--")
        ax.xaxis.set_major_formatter(xfmt)
        ax.grid(alpha=0.25)
        if legend:
            ax.legend(fontsize=7, loc="upper right")

    def draw_timeline(ax):
        """Dedicated thin panel: colored bands + clear date labels at edges."""
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.axvspan(last_trig_dt, itf_end_dt, color="gold",   alpha=0.6, zorder=0)
        ax.axvspan(outage_dt,    recovery_dt, color="tomato", alpha=0.5, zorder=0)
        ax.axvline(last_trig_dt, color="crimson", lw=1.2, ls="--")
        ax.xaxis.set_major_formatter(xfmt)
        ax.grid(alpha=0.15)
        kw = dict(fontsize=8, va="center", ha="left", rotation=90,
                  fontweight="bold", clip_on=False, zorder=7)
        # Yellow band edges
        ax.text(last_trig_dt + pd.Timedelta(hours=4), 0.5,
                last_trig_dt.strftime("%Y-%m-%d %H:%M UTC"),
                color="darkgoldenrod", **kw)
        ax.text(itf_end_dt   + pd.Timedelta(hours=4), 0.5,
                itf_end_dt.strftime("%Y-%m-%d"),
                color="darkgoldenrod", **kw)
        # Orange band edges
        ax.text(outage_dt    + pd.Timedelta(hours=4), 0.5,
                outage_dt.strftime("%Y-%m-%d %H:%M UTC"),
                color="darkred", **kw)
        ax.text(recovery_dt  + pd.Timedelta(hours=4), 0.5,
                recovery_dt.strftime("%Y-%m-%d (+10 d)"),
                color="darkred", **kw)
        # Legend patches
        from matplotlib.patches import Patch
        ax.legend(handles=[
            Patch(fc="gold",   alpha=0.7, label="SR tower opening (Christmas)"),
            Patch(fc="tomato", alpha=0.6, label="General power outage (~10-day recovery)"),
        ], fontsize=8, loc="upper left", framealpha=0.85)

    # ── figure: 9 panels (timeline + 8 data) ─────────────────────────────────
    fig, axes = plt.subplots(9, 1, figsize=(16, 26), sharex=True,
                             gridspec_kw={"height_ratios": [0.4]+[1]*8})
    fig.suptitle("Glitch rate and thermal channels: Oct 2025 – Apr 2026",
                 fontsize=13)

    # ── Panel 0: event timeline ───────────────────────────────────────────────
    draw_timeline(axes[0])

    # ── Panel 1: glitch rate ──────────────────────────────────────────────────
    ax = axes[1]
    rate = df["n_triggers"].values.astype(float)
    ax.bar(dt.values, rate, width=np.timedelta64(3600, "s"),
           color="steelblue", alpha=0.7)
    ax.annotate(f"Last glitch\n{last_trig_dt.strftime('%Y-%m-%d %H:%M UTC')}",
                xy=(last_trig_dt, rate.max() * 0.9),
                xytext=(last_trig_dt + pd.Timedelta(days=3), rate.max() * 0.75),
                color="crimson", fontsize=7,
                arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8))
    decorate(ax, legend=False)
    ax.set_ylabel("Rate [/h]", fontsize=9)
    ax.set_title("25-min glitch rate (1-h bins)", fontsize=9)

    # ── Panel 2: ITF lock fraction ────────────────────────────────────────────
    ax = axes[2]
    if lock is not None:
        lock_dt = gps_to_dt(lock["gps_bin"].values)
        ax.plot(lock_dt.values, lock["lock_frac"].values,
                color="tab:green", lw=0.8, label="ITF lock fraction")
        ax.set_ylabel("Lock fraction", fontsize=9)
        ax.set_ylim(-0.05, 1.05)
    else:
        ax.scatter(gps_to_dt(t_mid).values, smooth / 60,
                   s=2, alpha=0.4, color="tab:blue",
                   label="Recurrence [min]")
        ax.set_ylabel("Recurrence [min]", fontsize=9)
    decorate(ax)

    # ── Panel 3: tower bottom temperatures ───────────────────────────────────
    ax = axes[3]
    for col, lab, color in [("ni_bottom_te1", "NI_BOTTOM_TE1", "tab:orange"),
                             ("wi_bottom_te1", "WI_BOTTOM_TE1", "tab:red")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    decorate(ax)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 4: NI/WI CO2 bench ambient temperatures ─────────────────────────
    ax = axes[4]
    for col, lab, color in [("ni_co2_env_te", "NI CO2 ambient", "tab:blue"),
                             ("wi_co2_env_te", "WI CO2 ambient", "tab:cyan")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    decorate(ax)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 5: mirror coil temperatures ────────────────────────────────────
    ax = axes[5]
    for col, lab, color in [("ni_mir_coil_te", "NI mirror coil", "tab:blue"),
                             ("wi_mir_coil_te", "WI mirror coil", "tab:purple")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    decorate(ax)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 6: ring heater setpoints ───────────────────────────────────────
    ax = axes[6]
    for col, lab, color in [("ni_rh_set", "NI RH setpoint", "tab:green"),
                             ("wi_rh_set", "WI RH setpoint", "tab:olive")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    decorate(ax)
    ax.set_ylabel("[W]", fontsize=9)

    # ── Panel 7: CO2 laser body temperatures ─────────────────────────────────
    ax = axes[7]
    for col, lab, color in [("ni_co2_tc", "NI CO2 laser body", "tab:green"),
                             ("wi_co2_tc", "WI CO2 laser body", "tab:olive")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    decorate(ax)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 8: CEB UPS current ──────────────────────────────────────────────
    ax = axes[8]
    if "ceb_ups_curr_r" in df.columns:
        ax.plot(dt.values, df["ceb_ups_curr_r"].values,
                color="tab:red", lw=0.6, label="CEB UPS curr R")
    decorate(ax)
    ax.set_ylabel("[A]", fontsize=9)
    ax.set_xlabel("Date", fontsize=9)

    fig.autofmt_xdate(rotation=20)
    fig.tight_layout()
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binned",    required=True)
    p.add_argument("--step4-dir", required=True)
    p.add_argument("--triggers",  required=True)
    p.add_argument("--itf-lock",  default="")
    p.add_argument("--output",
                   default="usecases/25-minute-glitch/disappearance_timeseries.png")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot(args.binned, args.step4_dir, args.triggers, args.itf_lock, args.output)
