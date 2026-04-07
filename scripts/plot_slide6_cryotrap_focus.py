#!/usr/bin/env python3
"""
Focused Jan-Feb 2026 view around the slide-6 cryogenic-trap-leak period.

Builds a multi-panel figure around the aligned SR operation used in the
presentation:
  - DET_B1p_DC (B1p carrier proxy)
  - Hrec_Range_BNS sensitivity trend
  - NI thermal channels tied to the glitch investigation
  - Additional NI tower thermal channels
  - NI laser body temperatures near the etalon thermal actuators
  - Merged NI etalon raw-value channels

The script accepts two alignment masks:
  - proxy:  lock metric >= threshold
  - strict: lock metric in the chosen aligned band

The lock metric defaults to the historical hourly `lock_mean`, but can also be
an hourly summary derived from 1 Hz `V1:META_ITF_LOCK_index`, such as
`lock_mode`, `frac_145`, or `frac_144_145`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

CHANNEL_LABELS_1HZ = {
    "b1p_dc": "V1:LSC_B1p_DC_mean",
    "det_b1p_dc": "V1:DET_B1p_DC_mean",
    "hrec_range_bns": "V1:Hrec_Range_BNS",
    "ni_bottom_te1": "V1:INF_NI_BOTTOM_TE1",
    "ni_mir_coil_ul_te": "V1:INF_NI_MIR_COIL_UL_TE",
    "ni_rh_te": "V1:INF_TCS_NI_RH_TE",
    "env_co2_ni_te": "V1:ENV_TCS_CO2_NI_TE",
    "ni_f0_te1": "V1:ENV_NI_F0_TE1",
    "ni_f0_te2": "V1:ENV_NI_F0_TE2",
    "ni_f4_te1": "V1:ENV_NI_F4_TE1",
    "ni_f4_te2": "V1:ENV_NI_F4_TE2",
    "ni_f7_te1": "V1:ENV_NI_F7_TE1",
    "ni_f7_te2": "V1:ENV_NI_F7_TE2",
    "tcs_ni_co2laser_te": "V1:TCS_NI_TE_CO2Laser",
    "tcs_ni_auxlaser_te": "V1:TCS_NI_TE_AUXLaser",
    "ni_rh_out": "V1:LSC_Etalon_NI_RH_OUT_mean",
    "ni_rh_set": "V1:LSC_Etalon_NI_RH_SET_mean",
    "rh_ni_in": "V1:LSC_Etalon_NI_RH_IN_mean",
    "rh_ni_err": "V1:LSC_Etalon_NI_RH_ERR_mean",
}

CHANNEL_LABELS_50HZ = {
    **CHANNEL_LABELS_1HZ,
    "b1p_dc": "V1:LSC_B1p_DC_50Hz",
    "det_b1p_dc": "V1:DET_B1p_DC_50Hz",
    "ni_rh_out": "V1:LSC_Etalon_NI_RH_OUT",
    "ni_rh_set": "V1:LSC_Etalon_NI_RH_SET",
    "rh_ni_in": "V1:LSC_Etalon_NI_RH_IN",
    "rh_ni_err": "V1:LSC_Etalon_NI_RH_ERR",
}


def gps_to_dt(gps: pd.Series | np.ndarray) -> pd.Series:
    return GPS_EPOCH + pd.to_timedelta(pd.Series(gps, dtype=float), unit="s")


def robust_z(series: pd.Series) -> pd.Series:
    vals = pd.Series(series, dtype=float).copy()
    med = np.nanmedian(vals)
    mad = np.nanmedian(np.abs(vals - med))
    if not np.isfinite(mad) or mad == 0:
        return vals - med
    return (vals - med) / (1.4826 * mad)


def break_large_gaps(
    df: pd.DataFrame, cols: list[str], max_gap_hours: float = 8.0
) -> pd.DataFrame:
    out = df[["dt", *cols]].dropna(subset=["dt"]).sort_values("dt").copy()
    if out.empty:
        return out
    gap_h = out["dt"].diff().dt.total_seconds().div(3600.0)
    break_idx = np.flatnonzero(gap_h.fillna(0).to_numpy() > max_gap_hours)
    if len(break_idx) == 0:
        return out

    pieces: list[pd.DataFrame] = []
    start = 0
    for idx in break_idx:
        pieces.append(out.iloc[start:idx])
        nan_row = {c: np.nan for c in cols}
        nan_row["dt"] = out.iloc[idx - 1]["dt"] + pd.Timedelta(seconds=1)
        pieces.append(pd.DataFrame([nan_row]))
        start = idx
    pieces.append(out.iloc[start:])
    return pd.concat(pieces, ignore_index=True)


def load_csv(path: str | Path, cols: list[str] | None = None) -> pd.DataFrame:
    df = pd.read_csv(path)
    if cols is not None:
        keep = [c for c in cols if c in df.columns]
        df = df[keep]
    return df.sort_values("gps_bin").reset_index(drop=True)


def load_b1p_series(patterns: list[str]) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for pattern in patterns:
        for path in sorted(Path(".").glob(pattern)):
            df = pd.read_csv(path)
            cols = [c for c in ["gps_bin", "b1p_dc", "det_b1p_dc"] if c in df.columns]
            if "gps_bin" in cols and len(cols) > 1:
                parts.append(df[cols])
    if not parts:
        raise FileNotFoundError(
            "No B1p CSV matched the provided patterns with a b1p_dc or det_b1p_dc column."
        )
    out = (
        pd.concat(parts, ignore_index=True)
        .drop_duplicates(subset=["gps_bin"])
        .sort_values("gps_bin")
        .reset_index(drop=True)
    )
    return out


def build_mask(
    df: pd.DataFrame,
    mode: str,
    lock_field: str,
    lock_min: float,
    ty_abs_max: float,
    strict_low: float,
    strict_high: float,
) -> pd.Series:
    if lock_field not in df.columns:
        raise KeyError(f"Lock field {lock_field!r} not found in dataframe columns.")
    lock = df[lock_field]
    if mode == "strict":
        return lock.between(strict_low, strict_high, inclusive="both")
    if mode == "proxy":
        return lock >= lock_min
    raise ValueError(f"Unsupported alignment mode: {mode}")


def build_piecewise_mask(
    df: pd.DataFrame,
    cutoff: pd.Timestamp,
    lock_field: str,
    pre_lock_min: float,
    mode: str,
    post_lock_min: float,
    ty_abs_max: float,
    strict_low: float,
    strict_high: float,
) -> pd.Series:
    if lock_field not in df.columns:
        raise KeyError(f"Lock field {lock_field!r} not found in dataframe columns.")
    pre = (df["dt"] < cutoff) & (df[lock_field] >= pre_lock_min)
    post = build_mask(
        df,
        mode=mode,
        lock_field=lock_field,
        lock_min=post_lock_min,
        ty_abs_max=ty_abs_max,
        strict_low=strict_low,
        strict_high=strict_high,
    ) & (df["dt"] >= cutoff)
    return pre | post


def add_marker(ax: plt.Axes, when: pd.Timestamp, label: str) -> None:
    ax.axvline(when, color="black", lw=1.0, ls="--", alpha=0.65)
    ax.text(
        when,
        0.98,
        label,
        transform=ax.get_xaxis_transform(),
        rotation=90,
        ha="right",
        va="top",
        fontsize=8,
        color="black",
        backgroundcolor="white",
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lock-csv", default="outputs/itf_lock_full_o4.csv")
    ap.add_argument(
        "--lock-mask-field",
        default="lock_mean",
        help="Column from --lock-csv used to build the alignment mask.",
    )
    ap.add_argument("--sr-csv", default="outputs/sr_ty_input_merged.csv")
    ap.add_argument(
        "--hrec-csv",
        default="outputs/hrec_range_bns_1447200000_1455580800_fromcca.csv",
    )
    ap.add_argument("--ni-csv", default="outputs/ni_thermal_merged.csv")
    ap.add_argument(
        "--ni-laser-csv",
        default="outputs/ni_laser_temps_1443312000_1456358400_fromcca.csv",
    )
    ap.add_argument("--logbook-csv", default="outputs/logbook_channels_full_o4.csv")
    ap.add_argument(
        "--archive50-csv",
        default=None,
        help="Optional CSV derived from /data/archive_50Hz to override DET/etalon mean channels.",
    )
    ap.add_argument(
        "--channel-label-mode",
        choices=["1hz", "50hz"],
        default="1hz",
    )
    ap.add_argument(
        "--b1p-pattern",
        action="append",
        default=None,
        help="Glob pattern(s) for local B1p trend CSVs; may be passed more than once.",
    )
    ap.add_argument("--start", default="2025-10-01")
    ap.add_argument("--end", default="2026-02-28 23:59:59")
    ap.add_argument("--alignment-mode", choices=["proxy", "strict"], default="proxy")
    ap.add_argument(
        "--reference-alignment-mode", choices=["proxy", "strict"], default="proxy"
    )
    ap.add_argument("--aligned-lock-min", type=float, default=100.0)
    ap.add_argument("--aligned-ty-abs-max", type=float, default=10.0)
    ap.add_argument("--pre-christmas-lock-min", type=float, default=134.0)
    ap.add_argument("--christmas-cutoff", default="2025-12-24")
    ap.add_argument("--strict-lock-low", type=float, default=143.0)
    ap.add_argument("--strict-lock-high", type=float, default=146.0)
    ap.add_argument(
        "--uniform-mask",
        action="store_true",
        help="Apply the same alignment mask over the full period instead of the historical pre/post Christmas split.",
    )
    ap.add_argument("--reference-start", default="2025-11-15")
    ap.add_argument("--reference-end", default="2026-01-14")
    ap.add_argument("--fixed-b1p-70", type=float, default=None)
    ap.add_argument("--fixed-b1p-50", type=float, default=None)
    ap.add_argument(
        "--output",
        default="usecases/25-minute-glitch/slide6_cryotrap_focus_proxy.png",
    )
    args = ap.parse_args()

    # Keep the full lock CSV so custom hourly state summaries can be used
    # directly as --lock-mask-field without having to update a fixed allowlist.
    lock = load_csv(args.lock_csv)
    sr = load_csv(args.sr_csv, ["gps_bin", "asc_sr_ty_input"])
    hrec = load_csv(args.hrec_csv, ["gps_bin", "hrec_range_bns"])
    ni = load_csv(
        args.ni_csv,
        [
            "gps_bin",
            "ni_bottom_te1",
            "ni_mir_coil_ul_te",
            "ni_rh_te",
            "ni_f0_te1",
            "ni_f0_te2",
            "ni_f4_te1",
            "ni_f4_te2",
            "ni_f7_te1",
            "ni_f7_te2",
            "ni_rh_set",
            "ni_rh_out",
        ],
    )
    ni_laser_path = Path(args.ni_laser_csv)
    if ni_laser_path.exists():
        ni_laser = load_csv(
            ni_laser_path,
            ["gps_bin", "tcs_ni_co2laser_te", "tcs_ni_auxlaser_te"],
        )
        ni_laser = ni_laser.rename(
            columns={"tcs_ni_co2laser_te": "tcs_ni_co2laser_te_extra"}
        )
    else:
        ni_laser = pd.DataFrame(
            columns=["gps_bin", "tcs_ni_co2laser_te_extra", "tcs_ni_auxlaser_te"]
        )
    logbook = load_csv(
        args.logbook_csv,
        [
            "gps_bin",
            "rh_ni_in",
            "rh_ni_err",
            "env_co2_ni_te",
            "tcs_ni_co2laser_te",
        ],
    )
    b1p_patterns = args.b1p_pattern or ["outputs/lockpoint_compare_*_fromcca.csv"]
    b1p = load_b1p_series(b1p_patterns)
    archive50 = None
    if args.archive50_csv:
        archive50_path = Path(args.archive50_csv)
        if archive50_path.exists():
            archive50 = load_csv(
                archive50_path,
                [
                    "gps_bin",
                    "b1p_dc",
                    "det_b1p_dc",
                    "ni_rh_set",
                    "ni_rh_out",
                    "rh_ni_in",
                    "rh_ni_err",
                ],
            )

    df = lock.merge(sr, on="gps_bin", how="outer")
    df = df.merge(hrec, on="gps_bin", how="outer")
    df = df.merge(ni, on="gps_bin", how="outer")
    df = df.merge(ni_laser, on="gps_bin", how="outer")
    df = df.merge(logbook, on="gps_bin", how="outer")
    df = df.merge(b1p, on="gps_bin", how="outer")
    if archive50 is not None:
        df = df.merge(
            archive50,
            on="gps_bin",
            how="outer",
            suffixes=("", "_archive50"),
        )
    df = df.sort_values("gps_bin").drop_duplicates(subset=["gps_bin"]).reset_index(drop=True)
    if archive50 is not None:
        for col in ["b1p_dc", "det_b1p_dc", "ni_rh_set", "ni_rh_out", "rh_ni_in", "rh_ni_err"]:
            extra = f"{col}_archive50"
            if extra in df.columns:
                if col in df.columns:
                    df[col] = df[extra].combine_first(df[col])
                else:
                    df[col] = df[extra]
                df = df.drop(columns=[extra])
    if "tcs_ni_co2laser_te_extra" in df.columns:
        extra_co2 = df["tcs_ni_co2laser_te_extra"]
        if extra_co2.notna().any() and "tcs_ni_co2laser_te" in df.columns:
            df["tcs_ni_co2laser_te"] = df["tcs_ni_co2laser_te_extra"].combine_first(
                df["tcs_ni_co2laser_te"]
            )
        elif extra_co2.notna().any():
            df["tcs_ni_co2laser_te"] = df["tcs_ni_co2laser_te_extra"]
        df = df.drop(columns=["tcs_ni_co2laser_te_extra"])
    df["dt"] = gps_to_dt(df["gps_bin"])

    start = pd.Timestamp(args.start, tz="UTC")
    end = pd.Timestamp(args.end, tz="UTC")
    christmas_cutoff = pd.Timestamp(args.christmas_cutoff, tz="UTC")
    ref_start = pd.Timestamp(args.reference_start, tz="UTC")
    ref_end = pd.Timestamp(args.reference_end, tz="UTC")

    df = df[(df["dt"] >= start) & (df["dt"] <= end)].copy()
    if args.uniform_mask:
        plot_mask = build_mask(
            df,
            mode=args.alignment_mode,
            lock_field=args.lock_mask_field,
            lock_min=args.aligned_lock_min,
            ty_abs_max=args.aligned_ty_abs_max,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )
        ref_mask = build_mask(
            df,
            mode=args.reference_alignment_mode,
            lock_field=args.lock_mask_field,
            lock_min=args.aligned_lock_min,
            ty_abs_max=args.aligned_ty_abs_max,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )
    else:
        plot_mask = build_piecewise_mask(
            df,
            cutoff=christmas_cutoff,
            lock_field=args.lock_mask_field,
            pre_lock_min=args.pre_christmas_lock_min,
            mode=args.alignment_mode,
            post_lock_min=args.aligned_lock_min,
            ty_abs_max=args.aligned_ty_abs_max,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )
        ref_mask = build_piecewise_mask(
            df,
            cutoff=christmas_cutoff,
            lock_field=args.lock_mask_field,
            pre_lock_min=args.pre_christmas_lock_min,
            mode=args.reference_alignment_mode,
            post_lock_min=args.aligned_lock_min,
            ty_abs_max=args.aligned_ty_abs_max,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )

    plot_df = df.loc[plot_mask].copy()
    ref_df = df.loc[ref_mask & df["dt"].between(ref_start, ref_end, inclusive="left")].copy()
    ref_series = ref_df["det_b1p_dc"].dropna()
    if ref_series.empty:
        fallback = df.loc[ref_mask & df["det_b1p_dc"].notna(), ["dt", "det_b1p_dc"]].head(24)
        ref_series = fallback["det_b1p_dc"]
        ref_msg = (
            "reference fallback used: first 24 reference-mask B1p points in the plotted window"
        )
    else:
        ref_msg = "reference window used as requested"

    if ref_series.empty:
        raise RuntimeError("Could not derive any B1p reference points for the first-panel reference.")

    b1p_ref = float(ref_series.median())
    b1p_70 = 0.70 * b1p_ref
    b1p_50 = 0.50 * b1p_ref
    using_fixed_b1p_lines = (
        args.fixed_b1p_70 is not None and args.fixed_b1p_50 is not None
    )
    if using_fixed_b1p_lines:
        b1p_70 = float(args.fixed_b1p_70)
        b1p_50 = float(args.fixed_b1p_50)

    first_diaphragm = pd.Timestamp("2025-12-25", tz="UTC")
    second_diaphragm = pd.Timestamp("2026-02-07", tz="UTC")
    cryotrap_leak_start = pd.Timestamp("2026-01-10", tz="UTC")
    cryotrap_leak_end = second_diaphragm
    channel_labels = CHANNEL_LABELS_50HZ if args.channel_label_mode == "50hz" else CHANNEL_LABELS_1HZ

    fig, axes = plt.subplots(
        6,
        1,
        figsize=(16, 16),
        sharex=True,
        gridspec_kw={"height_ratios": [1.5, 1.0, 1.2, 1.0, 0.9, 1.1], "hspace": 0.08},
    )

    # Panel 0: B1p carrier (same mask as the lower panels)
    ax = axes[0]
    b1p_sel = plot_df[plot_df["det_b1p_dc"].notna() & (plot_df["det_b1p_dc"] > 50)]
    b1p_pre = b1p_sel[b1p_sel["dt"] < first_diaphragm]["det_b1p_dc"]
    b1p_between = b1p_sel[
        (b1p_sel["dt"] >= first_diaphragm) & (b1p_sel["dt"] < second_diaphragm)
    ]["det_b1p_dc"]
    b1p_after = b1p_sel[b1p_sel["dt"] >= second_diaphragm]["det_b1p_dc"]
    b1p_daily_pre = (
        b1p_sel[b1p_sel["dt"] < first_diaphragm]
        .assign(day=lambda x: x["dt"].dt.floor("D"))
        .groupby("day", as_index=False)["det_b1p_dc"]
        .mean()
    )
    b1p_daily = (
        b1p_sel[
            (b1p_sel["dt"] >= cryotrap_leak_start) & (b1p_sel["dt"] < cryotrap_leak_end)
        ]
        .assign(day=lambda x: x["dt"].dt.floor("D"))
        .groupby("day", as_index=False)["det_b1p_dc"]
        .mean()
    )
    b1p_daily_post = (
        b1p_sel[b1p_sel["dt"] >= second_diaphragm]
        .assign(day=lambda x: x["dt"].dt.floor("D"))
        .groupby("day", as_index=False)["det_b1p_dc"]
        .mean()
    )
    ax.plot(
        b1p_sel["dt"],
        b1p_sel["det_b1p_dc"],
        color="#1f77b4",
        ls="None",
        marker="o",
        ms=2.6,
        label=channel_labels["det_b1p_dc"],
    )
    if not b1p_daily_pre.empty:
        ax.plot(
            b1p_daily_pre["day"],
            b1p_daily_pre["det_b1p_dc"],
            color="indianred",
            lw=1.4,
            marker="o",
            ms=2.8,
            ls=":",
            label="Daily mean pre-Christmas",
            zorder=4,
        )
        pre_factor = float(
            b1p_daily_pre["det_b1p_dc"].max() / b1p_daily_pre["det_b1p_dc"].min()
        )
        pre_mid_time = start + (first_diaphragm - start) / 2
        ax.text(
            pre_mid_time,
            float(b1p_daily_pre["det_b1p_dc"].max()) + 14,
            f"x{pre_factor:.1f} modulation",
            color="indianred",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    if not b1p_daily.empty:
        ax.plot(
            b1p_daily["day"],
            b1p_daily["det_b1p_dc"],
            color="crimson",
            lw=1.6,
            marker="o",
            ms=3.0,
            label="Daily mean in leak window",
            zorder=4,
        )
        leak_factor = float(b1p_daily["det_b1p_dc"].max() / b1p_daily["det_b1p_dc"].min())
        leak_mid_time = cryotrap_leak_start + (cryotrap_leak_end - cryotrap_leak_start) / 2
        ax.text(
            leak_mid_time,
            float(b1p_daily["det_b1p_dc"].max()) + 14,
            f"x{leak_factor:.1f} modulation",
            color="crimson",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    if not b1p_daily_post.empty:
        ax.plot(
            b1p_daily_post["day"],
            b1p_daily_post["det_b1p_dc"],
            color="firebrick",
            lw=1.6,
            marker="o",
            ms=3.0,
            ls="--",
            label="Daily mean after 2nd baffle",
            zorder=4,
        )
        post_factor = float(
            b1p_daily_post["det_b1p_dc"].max() / b1p_daily_post["det_b1p_dc"].min()
        )
        post_mid_time = second_diaphragm + (end - second_diaphragm) / 2
        ax.text(
            post_mid_time,
            float(b1p_daily_post["det_b1p_dc"].max()) + 14,
            f"x{post_factor:.1f} modulation",
            color="firebrick",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    if not b1p_pre.empty:
        mean_pre = float(b1p_pre.mean())
        ax.plot(
            [start, first_diaphragm],
            [mean_pre, mean_pre],
            color="crimson",
            lw=1.8,
            zorder=3,
        )
        ax.text(
            start + (first_diaphragm - start) / 2,
            mean_pre + 8,
            f"{mean_pre:.1f} mW",
            color="crimson",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    if not b1p_between.empty:
        mean_between = float(b1p_between.mean())
        ax.plot(
            [first_diaphragm, second_diaphragm],
            [mean_between, mean_between],
            color="crimson",
            lw=1.8,
            zorder=3,
        )
        ax.text(
            first_diaphragm + (second_diaphragm - first_diaphragm) / 2,
            mean_between + 8,
            f"{mean_between:.1f} mW",
            color="crimson",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    if not b1p_after.empty:
        mean_after = float(b1p_after.mean())
        ax.plot(
            [second_diaphragm, end],
            [mean_after, mean_after],
            color="crimson",
            lw=1.8,
            zorder=3,
        )
        ax.text(
            second_diaphragm + 0.84 * (end - second_diaphragm),
            mean_after + 8,
            f"{mean_after:.1f} mW",
            color="crimson",
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )
    top_label_y = 1.05
    if args.uniform_mask:
        ax.text(
            start + (end - start) / 2,
            top_label_y,
            "SR ALIGNED",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=12,
            color="black",
            clip_on=False,
        )
    else:
        ax.text(
            start + (first_diaphragm - start) / 2,
            top_label_y,
            "SR NOT ALIGNED",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=11,
            color="black",
            clip_on=False,
        )
        ax.text(
            first_diaphragm + (second_diaphragm - first_diaphragm) / 2,
            top_label_y,
            "SR ALIGNED",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=11,
            color="black",
            clip_on=False,
        )
        ax.text(
            second_diaphragm + (end - second_diaphragm) / 2,
            top_label_y,
            "SR ALIGNED",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=11,
            color="black",
            clip_on=False,
        )
    add_marker(ax, first_diaphragm, "First diaphragm\n2025-12-25")
    add_marker(ax, second_diaphragm, "Second diaphragm\n2026-02-07")
    ax.set_ylabel("B1p carrier\n[mW]")
    ax.set_ylim(0, 500)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=8, ncol=2, framealpha=0.9)

    # Panel 1: Hrec range BNS
    ax = axes[1]
    hrec_sel = plot_df.dropna(subset=["hrec_range_bns"])
    ax.plot(
        hrec_sel["dt"],
        hrec_sel["hrec_range_bns"],
        color="#2c7fb8",
        ls="None",
        marker="o",
        ms=2.4,
        label=channel_labels["hrec_range_bns"],
    )
    add_marker(ax, first_diaphragm, "First diaphragm")
    add_marker(ax, second_diaphragm, "Second diaphragm")
    ax.set_ylabel("Hrec range BNS\n[Mpc]")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.9)

    # Panel 2: NI temperatures
    ax = axes[2]
    temp_specs = [
        ("ni_bottom_te1", channel_labels["ni_bottom_te1"], "tab:red", "-"),
        ("ni_mir_coil_ul_te", channel_labels["ni_mir_coil_ul_te"], "tab:orange", "-"),
        ("ni_rh_te", channel_labels["ni_rh_te"], "tab:purple", "-"),
        ("env_co2_ni_te", channel_labels["env_co2_ni_te"], "darkgreen", "--"),
    ]
    for col, label, color, ls in temp_specs:
        s = plot_df[["dt", col]].dropna()
        if not s.empty:
            ax.plot(
                s["dt"],
                s[col],
                color=color,
                ls="None",
                marker="o",
                ms=2.4,
                label=label,
            )
    add_marker(ax, first_diaphragm, "First diaphragm")
    add_marker(ax, second_diaphragm, "Second diaphragm")
    ax.set_ylabel("Temperature\n[°C]")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left", fontsize=8, ncol=2, framealpha=0.9)

    # Panel 3: Additional NI thermal channels
    ax = axes[3]
    extra_temp_specs = [
        ("ni_f0_te1", channel_labels["ni_f0_te1"], "tab:brown"),
        ("ni_f0_te2", channel_labels["ni_f0_te2"], "peru"),
        ("ni_f4_te1", channel_labels["ni_f4_te1"], "tab:green"),
        ("ni_f4_te2", channel_labels["ni_f4_te2"], "limegreen"),
        ("ni_f7_te1", channel_labels["ni_f7_te1"], "tab:cyan"),
        ("ni_f7_te2", channel_labels["ni_f7_te2"], "deepskyblue"),
    ]
    for col, label, color in extra_temp_specs:
        s = plot_df[["dt", col]].dropna()
        if not s.empty:
            ax.plot(
                s["dt"],
                s[col],
                color=color,
                ls="None",
                marker="o",
                ms=2.0,
                label=label,
            )
    add_marker(ax, first_diaphragm, "First diaphragm")
    add_marker(ax, second_diaphragm, "Second diaphragm")
    ax.set_ylabel("Additional NI thermal\n[degC]")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left", fontsize=8, ncol=3, framealpha=0.9)

    # Panel 4: NI laser body temperatures near the etalon thermal actuators
    ax = axes[4]
    for col, label, color in [
        ("tcs_ni_co2laser_te", channel_labels["tcs_ni_co2laser_te"], "tab:blue"),
        ("tcs_ni_auxlaser_te", channel_labels["tcs_ni_auxlaser_te"], "tab:olive"),
    ]:
        s = plot_df[["dt", col]].dropna()
        if not s.empty:
            ax.plot(
                s["dt"],
                s[col],
                color=color,
                ls="None",
                marker="o",
                ms=2.2,
                label=label,
            )
    add_marker(ax, first_diaphragm, "First diaphragm")
    add_marker(ax, second_diaphragm, "Second diaphragm")
    ax.set_ylabel("NI laser TE\n[degC]")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", fontsize=8, framealpha=0.9)

    # Panel 5: merged NI etalon raw-value channels
    ax = axes[5]
    for col, label, color in [
        ("ni_rh_out", channel_labels["ni_rh_out"], "tab:blue"),
        ("ni_rh_set", channel_labels["ni_rh_set"], "tab:green"),
        ("rh_ni_in", channel_labels["rh_ni_in"], "goldenrod"),
        ("rh_ni_err", channel_labels["rh_ni_err"], "crimson"),
    ]:
        s = plot_df[["dt", col]].dropna()
        if not s.empty:
            ax.plot(
                s["dt"],
                s[col],
                color=color,
                ls="None",
                marker="o",
                ms=2.2,
                label=label,
            )
    add_marker(ax, first_diaphragm, "First diaphragm")
    add_marker(ax, second_diaphragm, "Second diaphragm")
    ax.set_ylabel("NI etalon\nraw value")
    ax.set_ylim(0, 3000)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=8, ncol=2, framealpha=0.9)

    for ax in axes:
        ax.set_xlim(start, end)

    axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[-1].tick_params(axis="x", rotation=25, labelsize=9)
    for ax in axes[:-1]:
        ax.tick_params(labelsize=9)

    post_selection_label = {
        "proxy": f"{args.lock_mask_field} >= {args.aligned_lock_min:.3g}",
        "strict": (
            f"{args.strict_lock_low:.3g} <= {args.lock_mask_field} <= "
            f"{args.strict_lock_high:.3g}"
        ),
    }[args.alignment_mode]
    if args.uniform_mask:
        selection_label = f"uniform full-period mask: {post_selection_label}"
    else:
        selection_label = (
            f"pre-Christmas {args.lock_mask_field} >= {args.pre_christmas_lock_min:.3g}; "
            f"post-Christmas {post_selection_label}"
        )
    ref_label = (
        f"reference: same piecewise logic with {args.reference_alignment_mode} post-Christmas mode, "
        f"{ref_start.date()} to {ref_end.date()}, median={b1p_ref:.1f} mW"
    )
    fig.suptitle(
        "Slide-6 cryogenic-trap-leak focus around SR aligned operation\n"
        f"{selection_label}; {ref_label}",
        fontsize=12,
        y=0.995,
    )
    fig.text(
        0.995,
        0.005,
        ref_msg,
        ha="right",
        va="bottom",
        fontsize=8,
        color="0.35",
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved -> {out}")
    print(f"Aligned plot mode: {args.alignment_mode}")
    print(f"Reference mode: {args.reference_alignment_mode}")
    print(f"Reference points: {len(ref_series)}")
    print(f"Reference B1p median: {b1p_ref:.3f} mW")
    print(f"Plotted aligned rows: {len(plot_df)}")
    print(ref_msg)


if __name__ == "__main__":
    main()
