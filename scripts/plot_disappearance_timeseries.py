"""
plot_disappearance_timeseries.py — Step 7 multiplot.

Plots the glitch rate and key thermal / electrical channels over the
disappearance epoch (Oct 2025 – Apr 2026), with a yellow band marking the
ITF shutdown period (Christmas 2025, Dec 12 – Jan 8).

The ITF down period is derived from the trigger catalog: the last trigger
before the 26.9-day gap is GPS 1449582668 (2025-12-12) and the first trigger
after the gap is GPS 1451908330 (2026-01-08).

Usage (run on CCA after Step 4 merge):
    python scripts/plot_disappearance_timeseries.py \\
        --binned   outputs/rate_correlation/binned_summary.csv \\
        --triggers data/full_25min_glitches_ER16-O4b.csv \\
        --output   usecases/25-minute-glitch/disappearance_timeseries.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.dates import DateFormatter

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

# Step-7 epoch
EPOCH_START_GPS = 1443657600   # 2025-10-01
EPOCH_END_GPS   = 1459468800   # 2026-04-01

# ITF down period derived from the 26.9-day gap in the trigger catalog.
# Last trigger before gap: GPS 1449582668 (2025-12-12)
# First trigger after gap: GPS 1451908330 (2026-01-08)
LAST_TRIGGER_GPS = 1449582668
ITF_DOWN_START   = 1449532800   # 2025-12-12 00:00 UTC (day boundary)
ITF_DOWN_END     = 1451952000   # 2026-01-08 12:00 UTC (conservative end)


def gps_to_dt(gps):
    return GPS_EPOCH + pd.to_timedelta(np.asarray(gps, float), unit="s")


def plot(binned_csv: str, triggers_csv: str, output: str) -> None:
    # Load and filter binned data to epoch
    df = pd.read_csv(binned_csv)
    df = df[(df["gps_bin"] >= EPOCH_START_GPS) & (df["gps_bin"] < EPOCH_END_GPS)].copy()
    df = df.sort_values("gps_bin").reset_index(drop=True)
    dt = gps_to_dt(df["gps_bin"].values)
    print(f"Bins: {len(df)}")

    # Load triggers for recurrence period panel
    trig = pd.read_csv(triggers_csv)
    t_all = np.sort(trig["time"].values)
    keep = np.concatenate([[True], np.diff(t_all) > 300])
    t = t_all[keep]
    t_ep = t[(t >= EPOCH_START_GPS) & (t < EPOCH_END_GPS)]
    intervals = np.diff(t_ep)
    t_mid = t_ep[:-1] + intervals / 2
    smooth = (pd.Series(intervals)
              .rolling(20, center=True, min_periods=5)
              .median().bfill().ffill().values)
    print(f"Triggers in epoch: {len(t_ep)}")

    last_trig_dt = gps_to_dt(LAST_TRIGGER_GPS)
    itf_start_dt = gps_to_dt(ITF_DOWN_START)
    itf_end_dt   = gps_to_dt(ITF_DOWN_END)

    xfmt = DateFormatter("%b %Y")

    def band_and_line(ax):
        ax.axvspan(itf_start_dt, itf_end_dt, color="gold", alpha=0.35, zorder=0)
        ax.axvline(last_trig_dt, color="crimson", lw=1.2, ls="--")
        ax.xaxis.set_major_formatter(xfmt)
        ax.grid(alpha=0.25)

    fig, axes = plt.subplots(8, 1, figsize=(16, 24), sharex=True)
    fig.suptitle("Glitch rate and thermal channels: Oct 2025 – Apr 2026",
                 fontsize=13)

    # ── Panel 1: glitch rate ───────────────────────────────────────────────────
    ax = axes[0]
    rate = df["n_triggers"].values.astype(float)
    ax.bar(dt.values, rate, width=np.timedelta64(3600, "s"),
           color="steelblue", alpha=0.7)
    band_and_line(ax)
    ylim = ax.get_ylim()
    ax.text(itf_start_dt + (itf_end_dt - itf_start_dt) * 0.1,
            max(ylim[1] * 0.7, 0.5),
            "ITF down\n(Christmas)", fontsize=7, va="top", color="saddlebrown")
    ax.set_ylabel("Rate [/h]", fontsize=9)
    ax.set_title("25-min glitch rate (1-h bins)", fontsize=9)

    # ── Panel 2: recurrence period ─────────────────────────────────────────────
    ax = axes[1]
    ax.scatter(gps_to_dt(t_mid).values, smooth / 60, s=2, alpha=0.4,
               color="tab:blue")
    band_and_line(ax)
    ax.set_ylabel("Recurrence [min]", fontsize=9)

    # ── Panel 3: tower bottom temperatures ────────────────────────────────────
    ax = axes[2]
    for col, lab, color in [("ni_bottom_te1", "NI_BOTTOM_TE1", "tab:orange"),
                             ("wi_bottom_te1", "WI_BOTTOM_TE1", "tab:red")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 4: mirror coil temperatures ─────────────────────────────────────
    ax = axes[3]
    for col, lab, color in [("ni_mir_coil_te", "NI mirror coil", "tab:blue"),
                             ("wi_mir_coil_te", "WI mirror coil", "tab:cyan")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 5: ring heater setpoints ────────────────────────────────────────
    ax = axes[4]
    for col, lab, color in [("ni_rh_set", "NI RH setpoint", "tab:green"),
                             ("wi_rh_set", "WI RH setpoint", "tab:olive")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[W]", fontsize=9)

    # ── Panel 6: NI CO2 laser power + NI ring heater input ────────────────────
    ax = axes[5]
    for col, lab, color in [("ni_co2_pwr", "NI CO2 laser power [W]", "tab:purple"),
                             ("ni_rh_in",  "NI RH input [W]",         "tab:brown")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[W]", fontsize=9)

    # ── Panel 7: CO2 laser body temperatures ──────────────────────────────────
    ax = axes[6]
    for col, lab, color in [("ni_co2_tc", "NI CO2 laser body", "tab:green"),
                             ("wi_co2_tc", "WI CO2 laser body", "tab:olive")]:
        if col in df.columns:
            ax.plot(dt.values, df[col].values, color=color, lw=0.8, label=lab)
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[°C]", fontsize=9)

    # ── Panel 8: CEB UPS current ───────────────────────────────────────────────
    ax = axes[7]
    if "ceb_ups_curr_r" in df.columns:
        ax.plot(dt.values, df["ceb_ups_curr_r"].values,
                color="tab:red", lw=0.6, label="CEB UPS curr R")
    band_and_line(ax)
    ax.legend(fontsize=7)
    ax.set_ylabel("[A]", fontsize=9)
    ax.set_xlabel("Date", fontsize=9)

    fig.autofmt_xdate(rotation=20)
    fig.tight_layout()
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {out}")
    print(f"\nITF down band: {itf_start_dt.date()} – {itf_end_dt.date()}")
    print(f"Last trigger:  {last_trig_dt.date()}  GPS {LAST_TRIGGER_GPS}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binned",   required=True,
                   help="binned_summary.csv from rate_correlation_direct.py")
    p.add_argument("--triggers", required=True,
                   help="full_25min_glitches_ER16-O4b.csv")
    p.add_argument("--output",
                   default="usecases/25-minute-glitch/disappearance_timeseries.png")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot(args.binned, args.triggers, args.output)
