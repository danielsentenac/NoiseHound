#!/usr/bin/env python3
"""
plot_observations_summary.py — compact summary figure for the 25-minute glitch study.

The figure combines:
  - SR angular operating-point change across four key windows
  - Exact LSC/sideband medians from narrow CCA pulls
  - NI-bottom threshold test showing the old thermal range persists post-Jan-14

The exact medians below were measured from the CCA narrow extraction jobs
34759608, 34759609, 34759610, 34759611 and then folded into the README.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
LAST_TRIGGER_GPS = 1449582668
POST_BAFFLE_GPS = 1452456000

WINDOWS = [
    {
        "key": "oct12",
        "label": "Oct 12\n2025",
        "start": "2025-10-12",
        "end": "2025-10-13",
        "states": ["LN3"],
        "color": "#4C78A8",
    },
    {
        "key": "dec11",
        "label": "Dec 11\n2025",
        "start": "2025-12-11",
        "end": "2025-12-12",
        "states": ["LN2", "LN3"],
        "color": "#F58518",
    },
    {
        "key": "jan8",
        "label": "Jan 8\n2026",
        "start": "2026-01-08",
        "end": "2026-01-09",
        "states": ["LN2", "LN3"],
        "color": "#B279A2",
    },
    {
        "key": "jan14_29",
        "label": "Jan 14-29\n2026",
        "start": "2026-01-14",
        "end": "2026-01-30",
        "states": ["LN2", "LN3", "LN3_ALIGNED"],
        "color": "#54A24B",
    },
]

# Exact medians from CCA narrow pulls.
EXACT = {
    "oct12": {
        "b4_112mhz": 0.026303,
        "mich_set_tot": 13.336378,
        "srcl_set_tot": 53.820241,
    },
    "dec11": {
        "b4_112mhz": 0.021563,
        "mich_set_tot": 51.183179,
        "srcl_set_tot": 19.196669,
    },
    "jan8": {
        "b4_112mhz": 0.024611,
        "mich_set_tot": 8.733583,
        "srcl_set_tot": 25.071546,
    },
    "jan14_29": {
        "b4_112mhz": 0.023942,
        "mich_set_tot": 22.372792,
        "srcl_set_tot": 3.761551,
    },
}


def state_name(v: float) -> str:
    if pd.isna(v):
        return "NA"
    if 134 <= v <= 135:
        return "LN3"
    if 144 <= v <= 145:
        return "LN3_ALIGNED"
    if 149 <= v <= 150:
        return "LN3_SQZ"
    if 129 <= v <= 130:
        return "LN2"
    if 124 <= v <= 125:
        return "LN1"
    if 104 <= v < 124:
        return "DC_READOUT"
    if 39 <= v < 104:
        return "ARMS"
    if 1 < v < 39:
        return "DRMI"
    return "OTHER"


def dedupe_glitches(gl: pd.DataFrame, min_sep_s: float = 300.0) -> pd.DataFrame:
    gl = gl.sort_values("time").reset_index(drop=True)
    keep = np.concatenate([[True], np.diff(gl["time"].values) > min_sep_s])
    return gl.loc[keep].reset_index(drop=True)


def build_base(lock_csv: str, sr_csv: str, srty_csv: str, ni_csv: str, glitch_csv: str) -> pd.DataFrame:
    lock = pd.read_csv(lock_csv)
    sr = pd.read_csv(sr_csv)
    srty = pd.read_csv(srty_csv)
    ni = pd.read_csv(ni_csv)
    gl = dedupe_glitches(pd.read_csv(glitch_csv))

    lock["state"] = lock["lock_mean"].map(state_name)
    base = (lock[["gps_bin", "lock_mean", "state"]]
            .merge(sr[["gps_bin", "sr_mar_tx_set"]], on="gps_bin", how="left")
            .merge(srty[["gps_bin", "asc_sr_ty_input"]], on="gps_bin", how="left")
            .merge(ni[["gps_bin", "ni_bottom_te1"]], on="gps_bin", how="left"))

    for col in ["sr_mar_tx_set", "asc_sr_ty_input", "ni_bottom_te1"]:
        base[col] = base[col].ffill().bfill()

    gl["gps_bin"] = (gl["time"] // 3600 * 3600).astype(int)
    rate = gl.groupby("gps_bin").size().rename("n_triggers").reset_index()
    base = base.merge(rate, on="gps_bin", how="left").fillna({"n_triggers": 0})
    base["time"] = GPS_EPOCH + pd.to_timedelta(base["gps_bin"], unit="s")
    return base


def compute_window_table(base: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for win in WINDOWS:
        sub = base[
            (base["time"] >= pd.Timestamp(win["start"], tz="UTC"))
            & (base["time"] < pd.Timestamp(win["end"], tz="UTC"))
            & (base["state"].isin(win["states"]))
        ]
        rows.append({
            "key": win["key"],
            "label": win["label"],
            "hours": len(sub),
            "glitch_rate": float(sub["n_triggers"].mean()),
            "tx": float(sub["sr_mar_tx_set"].median()),
            "ty_input": float(sub["asc_sr_ty_input"].median()),
            "ni_bottom": float(sub["ni_bottom_te1"].median()),
            "b4_112mhz": EXACT[win["key"]]["b4_112mhz"],
            "mich_set_tot": EXACT[win["key"]]["mich_set_tot"],
            "srcl_set_tot": EXACT[win["key"]]["srcl_set_tot"],
            "color": win["color"],
        })
    return pd.DataFrame(rows)


def compute_threshold_table(base: pd.DataFrame) -> pd.DataFrame:
    pre = base[(base["gps_bin"] < LAST_TRIGGER_GPS) & (base["state"] == "LN3")]
    post = base[(base["gps_bin"] >= POST_BAFFLE_GPS)
                & (base["state"].isin(["LN2", "LN3", "LN3_ALIGNED"]))]

    rows = []
    for thr in [24, 25, 26, 27]:
        sub_pre = pre[pre["ni_bottom_te1"] >= thr]
        sub_post = post[post["ni_bottom_te1"] >= thr]
        rows.append({
            "thr": thr,
            "pre_rate": float(sub_pre["n_triggers"].mean()),
            "post_rate": float(sub_post["n_triggers"].mean()) if len(sub_post) else np.nan,
            "pre_hours": int(len(sub_pre)),
            "post_hours": int(len(sub_post)),
        })
    return pd.DataFrame(rows)


def annotate_points(ax, xs, ys, fmt: str, dy_frac: float = 0.04) -> None:
    vals = [y for y in ys if np.isfinite(y)]
    if not vals:
        return
    ymin, ymax = min(vals), max(vals)
    span = ymax - ymin
    if span == 0:
        span = max(abs(ymax), 1.0)
    dy = span * dy_frac
    for x, y in zip(xs, ys):
        if np.isfinite(y):
            ax.text(x, y + dy, fmt.format(y), ha="center", va="bottom",
                    fontsize=8, color="black")


def add_metric_panel(ax, df: pd.DataFrame, column: str, title: str,
                     yfmt: str, yscale: str = "linear",
                     hlines: list[tuple[float, str, str]] | None = None) -> None:
    x = np.arange(len(df))
    colors = df["color"].tolist()
    y = df[column].values.astype(float)

    ax.plot(x, y, color="#666666", lw=1.0, alpha=0.8, zorder=1)
    ax.scatter(x, y, c=colors, s=70, edgecolors="black", linewidths=0.6, zorder=2)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_xticks(x, df["label"].tolist(), fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(alpha=0.25, axis="y")
    if yscale != "linear":
        ax.set_yscale(yscale)
    if hlines:
        for yv, ls, lbl in hlines:
            ax.axhline(yv, color="#444444", ls=ls, lw=0.9, alpha=0.6, label=lbl)
    annotate_points(ax, x, y, yfmt)


def plot_summary(base: pd.DataFrame, output: str) -> Path:
    window_df = compute_window_table(base)
    thr_df = compute_threshold_table(base)

    fig = plt.figure(figsize=(16, 8.5))
    fig.subplots_adjust(top=0.86)
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.15], wspace=0.35, hspace=0.38)

    ax_tx = fig.add_subplot(gs[0, 0])
    ax_ty = fig.add_subplot(gs[0, 1])
    ax_rate = fig.add_subplot(gs[0, 2])
    ax_b4 = fig.add_subplot(gs[1, 0])
    ax_mich = fig.add_subplot(gs[1, 1])
    ax_srcl = fig.add_subplot(gs[1, 2])
    ax_thr = fig.add_subplot(gs[:, 3])

    add_metric_panel(ax_tx, window_df, "tx", "SR TX Operating Point", "{:.0f}")
    add_metric_panel(
        ax_ty, window_df, "ty_input", "SR TY Input", "{:.0f}",
        yscale="symlog",
        hlines=[(190.0, "--", "LN3 setpoint"), (0.0, ":", "LN3_ALIGNED setpoint")],
    )
    add_metric_panel(ax_rate, window_df, "glitch_rate", "25-Minute Glitch Rate", "{:.2f}/h")
    add_metric_panel(ax_b4, window_df, "b4_112mhz", "B4 112 MHz Level", "{:.5f}")
    add_metric_panel(ax_mich, window_df, "mich_set_tot", "MICH_SET_TOT", "{:.1f}")
    add_metric_panel(ax_srcl, window_df, "srcl_set_tot", "SRCL_SET_TOT", "{:.1f}")

    ax_ty.legend(fontsize=7, loc="upper right", framealpha=0.8)

    x = np.arange(len(thr_df))
    w = 0.36
    ax_thr.bar(x - w / 2, thr_df["pre_rate"], width=w, color="#E45756", label="Pre-Christmas LN3")
    ax_thr.bar(x + w / 2, thr_df["post_rate"], width=w, color="#54A24B", label="Post-Jan-14 science-like")
    ax_thr.set_xticks(x, [f"{t}C" for t in thr_df["thr"]], fontsize=9)
    ax_thr.set_ylabel("glitch rate [per hour]", fontsize=9)
    ax_thr.set_title("NI Bottom Threshold Test", fontsize=10, loc="left")
    ax_thr.grid(alpha=0.25, axis="y")
    ax_thr.tick_params(axis="y", labelsize=8)
    ax_thr.legend(fontsize=8, loc="upper right", framealpha=0.8)

    for i, row in thr_df.iterrows():
        ax_thr.text(i - w / 2, row["pre_rate"] + 0.05, f"{row['pre_hours']}h",
                    ha="center", va="bottom", fontsize=7, color="#A33D3C", rotation=90)
        ax_thr.text(i + w / 2, max(row["post_rate"], 0.0) + 0.05, f"{row['post_hours']}h",
                    ha="center", va="bottom", fontsize=7, color="#2F6F35", rotation=90)

    ax_thr.text(
        0.02, 0.02,
        "Same NI thermal range persists after Jan 14,\n"
        "but the 25-minute family stays at 0/h.",
        transform=ax_thr.transAxes, ha="left", va="bottom", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#cccccc"),
    )

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=w["color"],
               markeredgecolor="black", label=w["label"], markersize=8)
        for w in WINDOWS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=9, frameon=False,
               bbox_to_anchor=(0.47, 0.955))

    fig.suptitle("25-Minute Glitch Disappearance: Summary of Main Observations",
                 fontsize=13, y=0.985)
    fig.text(
        0.5, 0.925,
        "SR angular state jumps at Christmas 2025; Jan 8 shows an intermediate lock point; "
        "stable post-Jan-14 operation is glitch-free despite NI temperatures returning to the old range.",
        ha="center", va="top", fontsize=9
    )

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--itf-lock", default="outputs/itf_lock_full_o4.csv")
    ap.add_argument("--sr", default="outputs/sr_o4_merged.csv")
    ap.add_argument("--srty", default="outputs/sr_ty_input_merged.csv")
    ap.add_argument("--ni", default="outputs/ni_thermal_merged.csv")
    ap.add_argument("--glitches", default="data/full_25min_glitches_ER16-O4b.csv")
    ap.add_argument("--output", default="usecases/25-minute-glitch/observations_summary.png")
    args = ap.parse_args()

    base = build_base(args.itf_lock, args.sr, args.srty, args.ni, args.glitches)
    out = plot_summary(base, args.output)
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
