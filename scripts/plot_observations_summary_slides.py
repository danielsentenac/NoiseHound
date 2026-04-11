#!/usr/bin/env python3
"""
plot_observations_summary_slides.py — slide-style narrative summary.

This version is intentionally more presentation-oriented than
plot_observations_summary.py.  It focuses on the key narrative:
  1. Christmas shutdown / Jan 8 relock chronology
  2. SR angular-state jump with explicit Jan 8 LN3 -> LN2 split
  3. Exact B4 / MICH / SRCL lock-point comparison from CCA narrow pulls
  4. NI-bottom threshold test showing the driver-side thermal condition survives

All values are the medians already documented in README section 4.8 and in
outputs/lockpoint_exact_comparison.txt.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import pandas as pd


COLORS = {
    "oct12": "#4C78A8",
    "dec11": "#F58518",
    "jan8": "#B279A2",
    "post": "#54A24B",
}

FULL_WINDOWS = [
    {"label": "Oct 12", "date": "2025-10-12", "color": COLORS["oct12"],
     "tx": -22.23, "ty": 190.248069, "rate": 2.0833,
     "b4": 0.026303, "mich": 13.336378, "srcl": 53.820241},
    {"label": "Dec 11", "date": "2025-12-11", "color": COLORS["dec11"],
     "tx": 2.5216, "ty": 189.761261, "rate": 1.7692,
     "b4": 0.021563, "mich": 51.183179, "srcl": 19.196669},
    {"label": "Jan 8", "date": "2026-01-08", "color": COLORS["jan8"],
     "tx": -747.0, "ty": -0.000794, "rate": 0.25,
     "b4": 0.024611, "mich": 8.733583, "srcl": 25.071546},
    {"label": "Jan 14-29", "date": "2026-01-22", "color": COLORS["post"],
     "tx": -758.7, "ty": 0.000229, "rate": 0.0,
     "b4": 0.023942, "mich": 22.372792, "srcl": 3.761551},
]

JAN8_SPLIT = [
    {"label": "Jan 8 LN3", "date": "2026-01-08 12:30", "color": COLORS["jan8"],
     "tx": -747.0, "ty": 158.949241, "rate": 0.5,
     "b4": 0.024090, "mich": 24.781122, "srcl": 27.966505},
    {"label": "Jan 8 LN2", "date": "2026-01-08 20:30", "color": COLORS["jan8"],
     "tx": -747.0, "ty": -0.003176, "rate": 0.1667,
     "b4": 0.024694, "mich": 8.247099, "srcl": 24.121096},
]

THRESHOLDS = [
    {"thr": 24, "pre_rate": 1.9705, "post_rate": 0.0, "pre_hours": 5865, "post_hours": 223},
    {"thr": 25, "pre_rate": 1.9843, "post_rate": 0.0, "pre_hours": 4912, "post_hours": 185},
    {"thr": 26, "pre_rate": 1.9995, "post_rate": 0.0, "pre_hours": 3869, "post_hours": 134},
    {"thr": 27, "pre_rate": 2.0013, "post_rate": 0.0, "pre_hours": 3029, "post_hours": 91},
]

CHRISTMAS_START = pd.Timestamp("2025-12-12 13:51", tz="UTC")
JAN8_RELOCK = pd.Timestamp("2026-01-08 12:00", tz="UTC")
STABLE_START = pd.Timestamp("2026-01-14 00:00", tz="UTC")
TIMELINE_END = pd.Timestamp("2026-01-30 00:00", tz="UTC")


def ts(s: str) -> pd.Timestamp:
    return pd.Timestamp(s, tz="UTC")


def make_plot(output: str) -> Path:
    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[0.55, 1.2, 1.0],
                  width_ratios=[1.08, 1.0], wspace=0.28, hspace=0.32)

    ax_time = fig.add_subplot(gs[0, :])
    ax_state = fig.add_subplot(gs[1, 0])
    ax_text = fig.add_subplot(gs[2, 0])
    ax_thr = fig.add_subplot(gs[2, 1])
    metric_gs = gs[1, 1].subgridspec(3, 1, hspace=0.18)
    ax_b4 = fig.add_subplot(metric_gs[0, 0])
    ax_mich = fig.add_subplot(metric_gs[1, 0])
    ax_srcl = fig.add_subplot(metric_gs[2, 0])

    # Timeline
    ax_time.set_title("Chronology", fontsize=11, loc="left")
    ax_time.axhspan(-0.2, 0.2, color="#f3f3f3", zorder=0)
    ax_time.axvspan(CHRISTMAS_START, JAN8_RELOCK, color="gold", alpha=0.35, zorder=0)
    ax_time.axvspan(STABLE_START, TIMELINE_END, color="#D9F2D9", alpha=0.7, zorder=0)
    ax_time.axhline(0, color="#444444", lw=1.0)
    for row in FULL_WINDOWS[:2]:
        t = ts(row["date"])
        ax_time.scatter(t, 0, s=90, color=row["color"], edgecolors="black", zorder=3)
        ax_time.text(t, 0.28, row["label"], ha="center", va="bottom", fontsize=9)
    for row in JAN8_SPLIT:
        t = ts(row["date"])
        ax_time.scatter(t, 0, s=85, marker="D", color=row["color"], edgecolors="black", zorder=3)
        ax_time.text(t, -0.34 if "LN2" in row["label"] else 0.28,
                     row["label"], ha="center",
                     va="top" if "LN2" in row["label"] else "bottom",
                     fontsize=9)
    post = FULL_WINDOWS[3]
    ax_time.scatter(ts(post["date"]), 0, s=90, color=post["color"], edgecolors="black", zorder=3)
    ax_time.text(ts(post["date"]), 0.28, post["label"], ha="center", va="bottom", fontsize=9)
    ax_time.text(CHRISTMAS_START + (JAN8_RELOCK - CHRISTMAS_START) / 2, -0.34,
                 "Christmas SR shutdown", ha="center", va="top", fontsize=9)
    ax_time.text(STABLE_START + (TIMELINE_END - STABLE_START) / 2, -0.34,
                 "stable glitch-free regime", ha="center", va="top", fontsize=9, color="#2F6F35")
    ax_time.set_ylim(-0.55, 0.55)
    ax_time.set_yticks([])
    ax_time.set_xlim(ts("2025-10-05"), ts("2026-01-30"))
    ax_time.xaxis.set_major_locator(mdates.MonthLocator())
    ax_time.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax_time.tick_params(axis="x", labelsize=9)
    ax_time.spines[["left", "right", "top"]].set_visible(False)

    # SR state map
    ax_state.set_title("SR Angular Working-Point Jump", fontsize=11, loc="left")
    state_points = [
        (FULL_WINDOWS[0]["tx"], FULL_WINDOWS[0]["ty"], FULL_WINDOWS[0]["label"], FULL_WINDOWS[0]["color"], "o"),
        (FULL_WINDOWS[1]["tx"], FULL_WINDOWS[1]["ty"], FULL_WINDOWS[1]["label"], FULL_WINDOWS[1]["color"], "o"),
        (JAN8_SPLIT[0]["tx"], JAN8_SPLIT[0]["ty"], JAN8_SPLIT[0]["label"], JAN8_SPLIT[0]["color"], "D"),
        (JAN8_SPLIT[1]["tx"], JAN8_SPLIT[1]["ty"], JAN8_SPLIT[1]["label"], JAN8_SPLIT[1]["color"], "D"),
        (FULL_WINDOWS[3]["tx"], FULL_WINDOWS[3]["ty"], FULL_WINDOWS[3]["label"], FULL_WINDOWS[3]["color"], "o"),
    ]
    for i in range(len(state_points) - 1):
        x0, y0, *_ = state_points[i]
        x1, y1, *_ = state_points[i + 1]
        ax_state.annotate("", xy=(x1, y1), xytext=(x0, y0),
                          arrowprops=dict(arrowstyle="->", color="#777777", lw=1.4))
    for x, y, label, color, marker in state_points:
        ax_state.scatter(x, y, s=120, marker=marker, color=color, edgecolors="black", zorder=3)
        ax_state.text(x, y, f" {label}", fontsize=9, ha="left", va="bottom")
    ax_state.set_xlabel("SAT_SR_MAR_TX_SET", fontsize=9)
    ax_state.set_ylabel("ASC_SR_TY_INPUT", fontsize=9)
    ax_state.set_yscale("symlog", linthresh=1.0)
    ax_state.axhline(190, color="#666666", ls="--", lw=0.9, alpha=0.7)
    ax_state.axhline(0, color="#666666", ls=":", lw=0.9, alpha=0.7)
    ax_state.grid(alpha=0.25)
    ax_state.tick_params(labelsize=8)
    ax_state.text(0.02, 0.03,
                  "Jan 8 is transitional:\nLN3 still has large TY,\nLN2 is already near TY ~ 0.",
                  transform=ax_state.transAxes, ha="left", va="bottom", fontsize=8,
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#cccccc"))

    # Metric triptych
    x = [0, 1, 2, 3]
    labels = ["Oct 12", "Dec 11", "Jan 8", "Jan 14-29"]
    split_x = [1.88, 2.12]
    metrics = [
        ("B4 112 MHz", "b4", ax_b4, "{:.5f}"),
        ("MICH_SET_TOT", "mich", ax_mich, "{:.1f}"),
        ("SRCL_SET_TOT", "srcl", ax_srcl, "{:.1f}"),
    ]
    for title, key, ax, fmt in metrics:
        y = [w[key] for w in FULL_WINDOWS]
        ax.plot(x, y, color="#777777", lw=1.2, zorder=1)
        for xi, w in zip(x, FULL_WINDOWS):
            ax.scatter(xi, w[key], s=80, color=w["color"], edgecolors="black", zorder=2)
            ax.text(xi, w[key], fmt.format(w[key]), fontsize=8, ha="center", va="bottom")
        for xi, row, marker in zip(split_x, JAN8_SPLIT, ["^", "v"]):
            ax.scatter(xi, row[key], s=75, marker=marker, facecolors="white",
                       edgecolors=row["color"], linewidths=1.8, zorder=3)
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(alpha=0.25, axis="y")
        ax.tick_params(labelsize=8)
        ax.set_xticks(x, labels if ax is ax_srcl else [""] * 4, fontsize=8)
    metric_legend = [
        Line2D([0], [0], marker="o", color="#777777", markerfacecolor="white",
               markeredgecolor="black", label="full-window median"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="white",
               markeredgecolor=COLORS["jan8"], label="Jan 8 LN3"),
        Line2D([0], [0], marker="v", color="none", markerfacecolor="white",
               markeredgecolor=COLORS["jan8"], label="Jan 8 LN2"),
    ]
    ax_b4.legend(handles=metric_legend, fontsize=7, loc="upper right", framealpha=0.8)

    # Text panel
    ax_text.axis("off")
    ax_text.set_title("What The Figure Says", fontsize=11, loc="left")
    callout = (
        "1. Christmas 2025 moves Virgo into a new SR family.\n\n"
        "2. Jan 8 is not the final post-baffle state:\n"
        "   the same day contains an SR-tilted LN3 and an almost aligned LN2.\n\n"
        "3. The exact CCA pulls confirm the commissioning picture:\n"
        "   sidebands sit between Oct and Dec on Jan 8,\n"
        "   MICH drops strongly, and SRCL then collapses further\n"
        "   in the stable glitch-free regime.\n\n"
        "4. NI-bottom still revisits 24–27 C after Jan 14,\n"
        "   but the 25-minute family stays at 0/h."
    )
    ax_text.text(0.02, 0.98, callout, ha="left", va="top", fontsize=10,
                 bbox=dict(boxstyle="round,pad=0.45", facecolor="#F8F8F8", edgecolor="#DDDDDD"))

    # Threshold panel
    ax_thr.set_title("NI Bottom Threshold Test", fontsize=11, loc="left")
    thr_x = range(len(THRESHOLDS))
    w = 0.34
    ax_thr.bar([i - w / 2 for i in thr_x], [r["pre_rate"] for r in THRESHOLDS],
               width=w, color="#E45756", label="Pre-Christmas LN3")
    ax_thr.bar([i + w / 2 for i in thr_x], [r["post_rate"] for r in THRESHOLDS],
               width=w, color=COLORS["post"], label="Post-Jan-14 science-like")
    ax_thr.set_xticks(list(thr_x), [f'{r["thr"]}C' for r in THRESHOLDS], fontsize=9)
    ax_thr.set_ylabel("glitch rate [per hour]", fontsize=9)
    ax_thr.grid(alpha=0.25, axis="y")
    ax_thr.tick_params(labelsize=8)
    for i, r in enumerate(THRESHOLDS):
        ax_thr.text(i - w / 2, r["pre_rate"] + 0.05, f'{r["pre_hours"]}h',
                    ha="center", va="bottom", fontsize=7, rotation=90, color="#A33D3C")
        ax_thr.text(i + w / 2, r["post_rate"] + 0.05, f'{r["post_hours"]}h',
                    ha="center", va="bottom", fontsize=7, rotation=90, color="#2F6F35")
    ax_thr.legend(fontsize=8, loc="upper right", framealpha=0.85)
    ax_thr.text(0.03, 0.05,
                "Old NI thermal range survives.\nGlitch family does not.",
                transform=ax_thr.transAxes, ha="left", va="bottom", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.9, edgecolor="#cccccc"))

    fig.suptitle("25-Minute Glitch: Narrative Summary", fontsize=15, y=0.985)
    fig.text(0.5, 0.955,
             "The disappearance is best explained by a new SR/LSC working point after the Christmas intervention, "
             "not by the disappearance of the slow NI-side thermal condition.",
             ha="center", va="top", fontsize=10)

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    out = make_plot("usecases/25-minute-glitch/observations_summary_slides.png")
    print(f"Saved -> {out}")


if __name__ == "__main__":
    main()
