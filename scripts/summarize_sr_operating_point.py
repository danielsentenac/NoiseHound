"""
summarize_sr_operating_point.py — compact report on the Christmas 2025 SR change.

Builds a 1-hour merged table from the local extraction CSVs and summarizes:
  - the SR operating-point jump across January 2026,
  - the post-Christmas disappearance of the 25-minute family,
  - whether NI-bottom temperatures still reach the old glitch-friendly range.

Usage:
    python scripts/summarize_sr_operating_point.py

Optional:
    python scripts/summarize_sr_operating_point.py --output outputs/sr_operating_point_summary.txt
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

LAST_TRIGGER_GPS = 1449582668     # 2025-12-12 13:51 UTC
POST_BAFFLE_GPS = 1452456000      # 2026-01-14 00:00 UTC

PRE_BAFFLE = (1448323200, 1448928000)   # 2025-11-24 -> 2025-12-01
POST_WIDE = (1452456000, 1458000000)    # 2026-01-14 -> 2026-03-20


def gps_to_dt(gps: float | int) -> pd.Timestamp:
    return GPS_EPOCH + pd.to_timedelta(float(gps), unit="s")


def dedupe_glitches(gl: pd.DataFrame, min_sep_s: float = 300.0) -> pd.DataFrame:
    gl = gl.sort_values("time").reset_index(drop=True)
    keep = np.concatenate([[True], np.diff(gl["time"].values) > min_sep_s])
    return gl.loc[keep].reset_index(drop=True)


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


def detect_steps(tx: np.ndarray, gps: np.ndarray,
                 threshold: float = 15.0,
                 min_gap_h: int = 48) -> pd.DataFrame:
    rows = []
    half = 6
    last_idx = -min_gap_h
    for i in range(half, len(tx) - half):
        before = np.nanmedian(tx[max(0, i - half):i])
        after = np.nanmedian(tx[i:i + half])
        delta = after - before
        if abs(delta) > threshold and (i - last_idx) > min_gap_h:
            rows.append({
                "gps": int(gps[i]),
                "tx_before": before,
                "tx_after": after,
                "delta_tx": delta,
            })
            last_idx = i
    return pd.DataFrame(rows)


def build_base(gl: pd.DataFrame,
               lock: pd.DataFrame,
               ni: pd.DataFrame,
               sr: pd.DataFrame,
               srty: pd.DataFrame) -> pd.DataFrame:
    base_bins = sorted(
        set(lock["gps_bin"])
        .union(ni["gps_bin"])
        .union(sr["gps_bin"])
        .union(srty["gps_bin"])
    )
    df = pd.DataFrame({"gps_bin": base_bins})

    rate = (pd.DataFrame({"gps_bin": (gl["time"] // 3600 * 3600).astype(int),
                          "n_triggers": 1})
            .groupby("gps_bin", as_index=False)["n_triggers"].sum())

    df = (df.merge(rate, on="gps_bin", how="left")
            .fillna({"n_triggers": 0})
            .merge(lock, on="gps_bin", how="left")
            .merge(ni, on="gps_bin", how="left")
            .merge(sr[["gps_bin", "sr_mar_tx_set"]], on="gps_bin", how="left")
            .merge(srty, on="gps_bin", how="left"))

    for col in ["sr_mar_tx_set", "asc_sr_ty_input", "ni_bottom_te1"]:
        df[col] = df[col].ffill().bfill()
    df["state"] = df["lock_mean"].map(state_name)
    return df


def epoch_summary(df: pd.DataFrame) -> pd.DataFrame:
    return (df.groupby("state")
              .agg(hours=("gps_bin", "size"),
                   glitches=("n_triggers", "sum"),
                   rate_per_h=("n_triggers", "mean"),
                   tx_med=("sr_mar_tx_set", "median"),
                   ty_med=("asc_sr_ty_input", "median"),
                   ni_bottom_med=("ni_bottom_te1", "median"))
              .sort_values("glitches", ascending=False))


def canonical_post_last(gl: pd.DataFrame,
                        sr: pd.DataFrame,
                        srty: pd.DataFrame,
                        lock: pd.DataFrame) -> pd.DataFrame:
    sub = gl[(gl["time"] >= LAST_TRIGGER_GPS)
             & (gl["q"] < 10)
             & (gl["frequency"].between(35, 55))].copy()
    base = (sr.merge(srty, on="gps_bin", how="outer")
              .merge(lock, on="gps_bin", how="outer")
              .sort_values("gps_bin")
              .reset_index(drop=True))
    base["sr_mar_tx_set"] = base["sr_mar_tx_set"].ffill().bfill()
    base["asc_sr_ty_input"] = base["asc_sr_ty_input"].ffill().bfill()
    idx = np.searchsorted(base["gps_bin"].values, sub["time"].values, side="right") - 1
    idx = np.clip(idx, 0, len(base) - 1)
    sub["gps_bin"] = base["gps_bin"].values[idx]
    sub["lock_mean"] = base["lock_mean"].values[idx]
    sub["state"] = [state_name(v) for v in sub["lock_mean"].values]
    sub["sr_mar_tx_set"] = base["sr_mar_tx_set"].values[idx]
    sub["asc_sr_ty_input"] = base["asc_sr_ty_input"].values[idx]
    sub["dt_utc"] = [gps_to_dt(t).strftime("%Y-%m-%d %H:%M:%S") for t in sub["time"].values]
    cols = ["dt_utc", "snr", "frequency", "q",
            "state", "lock_mean", "sr_mar_tx_set", "asc_sr_ty_input"]
    return sub[cols].reset_index(drop=True)


def threshold_summary(df: pd.DataFrame, thr: float) -> str:
    pre = df[(df["gps_bin"] < LAST_TRIGGER_GPS) & (df["state"] == "LN3")]
    post = df[(df["gps_bin"] >= POST_BAFFLE_GPS)
              & (df["state"].isin(["LN2", "LN3", "LN3_ALIGNED"]))]

    sub_pre = pre[pre["ni_bottom_te1"] >= thr]
    sub_post = post[post["ni_bottom_te1"] >= thr]

    lines = [
        f"NI_BOTTOM >= {thr:.0f} C",
        (f"  pre-Christmas LN3      : hours={len(sub_pre):4d}  "
         f"glitches={int(sub_pre['n_triggers'].sum()):5d}  "
         f"rate={sub_pre['n_triggers'].mean():.3f}/h  "
         f"TX_med={sub_pre['sr_mar_tx_set'].median():.1f}  "
         f"TY_med={sub_pre['asc_sr_ty_input'].median():.3f}"),
        (f"  post-14-Jan science-ish: hours={len(sub_post):4d}  "
         f"glitches={int(sub_post['n_triggers'].sum()):5d}  "
         f"rate={sub_post['n_triggers'].mean() if len(sub_post) else np.nan:.3f}/h  "
         f"TX_med={sub_post['sr_mar_tx_set'].median() if len(sub_post) else np.nan:.1f}  "
         f"TY_med={sub_post['asc_sr_ty_input'].median() if len(sub_post) else np.nan:.3f}"),
    ]
    return "\n".join(lines)


def render(args) -> str:
    raw_gl = pd.read_csv(args.glitches)
    gl = dedupe_glitches(raw_gl)
    lock = pd.read_csv(args.itf_lock).sort_values("gps_bin")
    ni = pd.read_csv(args.ni).sort_values("gps_bin")
    sr = pd.read_csv(args.sr).sort_values("gps_bin")
    srty = pd.read_csv(args.srty).sort_values("gps_bin")

    base = build_base(gl, lock, ni, sr, srty)

    pre = base[(base["gps_bin"] >= PRE_BAFFLE[0]) & (base["gps_bin"] < PRE_BAFFLE[1])]
    post = base[(base["gps_bin"] >= POST_WIDE[0]) & (base["gps_bin"] < POST_WIDE[1])]

    steps = detect_steps(sr["sr_mar_tx_set"].ffill().bfill().values,
                         sr["gps_bin"].values)
    steps["dt_utc"] = [gps_to_dt(g).strftime("%Y-%m-%d %H:%M:%S") for g in steps["gps"].values]
    steps["abs_delta"] = steps["delta_tx"].abs()
    major_steps = steps.sort_values("abs_delta", ascending=False).head(5).copy()

    for side in ["before", "after"]:
        col = f"n_{side}_72h"
        major_steps[col] = [
            len(gl[(gl["time"] >= gps - 72 * 3600) & (gl["time"] < gps)]) if side == "before"
            else len(gl[(gl["time"] > gps) & (gl["time"] <= gps + 72 * 3600)])
            for gps in major_steps["gps"].values
        ]

    lines = []
    lines.append("SR Operating-Point Summary")
    lines.append("")
    lines.append(f"Trigger file      : {args.glitches}")
    lines.append(f"Deduped triggers  : {len(gl)} from {len(raw_gl)} raw rows")
    lines.append(f"Last pre-gap GPS  : {LAST_TRIGGER_GPS} ({gps_to_dt(LAST_TRIGGER_GPS)})")
    lines.append("")
    lines.append("Window summaries")
    lines.append("")
    lines.append("Pre-baffle window (2025-11-24 -> 2025-12-01)")
    lines.append(epoch_summary(pre).to_string())
    lines.append("")
    lines.append("Post-baffle window (2026-01-14 -> 2026-03-20)")
    lines.append(epoch_summary(post).to_string())
    lines.append("")
    lines.append("NI-bottom threshold test")
    for thr in [24, 25, 26, 27]:
        lines.append(threshold_summary(base, thr))
        lines.append("")
    lines.append("Canonical post-gap residual triggers")
    lines.append(canonical_post_last(gl, sr, srty, lock).to_string(index=False))
    lines.append("")
    lines.append("Largest SR TX steps")
    lines.append(major_steps[["dt_utc", "tx_before", "tx_after", "delta_tx",
                              "n_before_72h", "n_after_72h"]].to_string(index=False))
    lines.append("")
    lines.append("Interpretation")
    lines.append("- The SR operating point does not return to the pre-Christmas regime after 2026-01-08.")
    lines.append("- Post-2026-01-14 science-like states sit near TX ~ -759 to -806 and TY_INPUT ~ 0, with zero 25-minute glitches.")
    lines.append("- NI_BOTTOM still revisits the old 26-27+ C range post-Christmas, so NI temperature alone is not sufficient once the SR configuration changes.")
    lines.append("- Only a few canonical residual events remain on 2026-01-08 during relock/reacquisition; none survive in the stable post-2026-01-14 configuration.")
    return "\n".join(lines) + "\n"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--glitches", default="data/full_25min_glitches_ER16-O4b.csv")
    p.add_argument("--itf-lock", default="outputs/itf_lock_full_o4.csv")
    p.add_argument("--ni",       default="outputs/ni_thermal_merged.csv")
    p.add_argument("--sr",       default="outputs/sr_o4_merged.csv")
    p.add_argument("--srty",     default="outputs/sr_ty_input_merged.csv")
    p.add_argument("--output",   default="")
    return p.parse_args()


def main():
    args = parse_args()
    text = render(args)
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
