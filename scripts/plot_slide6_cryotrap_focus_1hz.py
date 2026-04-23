#!/usr/bin/env python3
"""
Build a slide-6-style cryotrap-focus plot directly from 1 Hz Virgo trend files.

Unlike plot_slide6_cryotrap_focus.py, this script does not read hourly CSVs.
It stages daily V-trend GWF files, reads the relevant channels at 1 Hz, applies
the same piecewise lock-state mask, and plots one point per second.

Intended to run on CCA where rfcp to HPSS is available.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import astropy.units.quantity as _aq
    from astropy.units import UnitConversionError as _UCE
    from astropy.units import dimensionless_unscaled as _dless

    _orig_cau = _aq.converters_and_unit

    def _safe_cau(function, method, *inputs):
        try:
            return _orig_cau(function, method, *inputs)
        except (_UCE, ValueError):
            return [None] * len(inputs), _dless

    _aq.converters_and_unit = _safe_cau
except Exception:
    pass

from gwpy.timeseries import TimeSeries, TimeSeriesDict

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"

CHANNELS = {
    "lock_state": "V1:META_ITF_LOCK_index",
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


def gps_to_dt(gps: np.ndarray | pd.Series) -> pd.Series:
    return GPS_EPOCH + pd.to_timedelta(pd.Series(gps, dtype=np.int64), unit="s")


def dt_to_gps(ts: pd.Timestamp) -> int:
    return int((ts - GPS_EPOCH) / pd.Timedelta(seconds=1))


def stage_day(gps_day: int, stage_dir: Path) -> Path | None:
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local = stage_dir / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"rfcp failed for {gps_day}: {r.stderr[:160]}", flush=True)
            return None
    return local


def _series_to_aligned_values(ts: TimeSeries, gps_grid: np.ndarray) -> np.ndarray:
    times = np.rint(np.asarray(ts.times.value, float)).astype(np.int64)
    vals = np.asarray(ts.value, dtype=float)
    if len(times) == 0:
        return np.full(len(gps_grid), np.nan, dtype=float)
    s = pd.Series(vals, index=times)
    s = s[~s.index.duplicated(keep="first")]
    return s.reindex(gps_grid).to_numpy(dtype=float)


def read_day_1hz(gwf: Path, t0: int, t1: int) -> pd.DataFrame:
    gps_grid = np.arange(t0, t1, dtype=np.int64)
    out = pd.DataFrame({"gps": gps_grid})
    requested = list(CHANNELS.values())
    ts_map: dict[str, TimeSeries] = {}

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tdict = TimeSeriesDict.read(str(gwf), requested, start=t0, end=t1)
        for key, ch in CHANNELS.items():
            if ch in tdict:
                ts_map[key] = tdict[ch]
    except Exception as exc:
        print(f"bulk read failed for {gwf.name}: {exc}", flush=True)

    for key, ch in CHANNELS.items():
        if key in ts_map:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts_map[key] = TimeSeries.read(str(gwf), ch, start=t0, end=t1)
        except Exception as exc:
            print(f"{gwf.name}: missing {ch} ({exc})", flush=True)

    for key, ts in ts_map.items():
        out[key] = _series_to_aligned_values(ts, gps_grid)
    out["dt"] = gps_to_dt(out["gps"].to_numpy())
    return out


def build_piecewise_mask(
    df: pd.DataFrame,
    cutoff: pd.Timestamp,
    pre_lock_min: float,
    mode: str,
    post_lock_min: float,
    strict_low: float,
    strict_high: float,
) -> pd.Series:
    lock = df["lock_state"]
    pre = (df["dt"] < cutoff) & lock.notna() & (lock >= pre_lock_min)
    if mode == "strict":
        post = (
            (df["dt"] >= cutoff)
            & lock.notna()
            & lock.between(strict_low, strict_high, inclusive="both")
        )
    elif mode == "proxy":
        post = (df["dt"] >= cutoff) & lock.notna() & (lock >= post_lock_min)
    else:
        raise ValueError(f"Unsupported mode: {mode}")
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


def plot_points(
    ax: plt.Axes,
    dt: pd.Series,
    vals: pd.Series,
    color: str,
    label: str | None,
    ms: float = 0.65,
) -> None:
    if vals.empty:
        return
    ax.plot(
        dt.to_numpy(),
        vals.to_numpy(),
        ls="None",
        marker=",",
        ms=ms,
        color=color,
        alpha=0.65,
        label=label,
        rasterized=True,
    )


def finalize_legend(ax: plt.Axes, **kwargs):
    leg = ax.legend(**kwargs)
    handles = getattr(leg, "legend_handles", None)
    if handles is None:
        handles = getattr(leg, "legendHandles", [])
    for handle in handles:
        if not hasattr(handle, "get_marker"):
            continue
        marker = handle.get_marker()
        if marker in (None, "", " ", "None", "none"):
            continue
        handle.set_alpha(1.0)
        handle.set_marker("o")
        handle.set_markersize(6.5)
        if hasattr(handle, "set_markeredgewidth"):
            handle.set_markeredgewidth(0.0)
    return leg


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2025-10-01 00:00:00")
    ap.add_argument("--end", default="2026-02-28 23:59:59")
    ap.add_argument("--alignment-mode", choices=["proxy", "strict"], default="proxy")
    ap.add_argument(
        "--reference-alignment-mode", choices=["proxy", "strict"], default="proxy"
    )
    ap.add_argument("--aligned-lock-min", type=float, default=100.0)
    ap.add_argument("--pre-christmas-lock-min", type=float, default=134.0)
    ap.add_argument("--christmas-cutoff", default="2025-12-24 00:00:00")
    ap.add_argument("--strict-lock-low", type=float, default=143.0)
    ap.add_argument("--strict-lock-high", type=float, default=146.0)
    ap.add_argument("--reference-start", default="2025-11-15 00:00:00")
    ap.add_argument("--reference-end", default="2026-01-14 00:00:00")
    ap.add_argument("--fixed-b1p-70", type=float, default=145.0)
    ap.add_argument("--fixed-b1p-50", type=float, default=75.0)
    ap.add_argument("--stage-dir", default=os.environ.get("STAGE_DIR", "/tmp/nh_slide6_1hz"))
    ap.add_argument(
        "--output",
        default="usecases/25-minute-glitch/slide6_cryotrap_focus_proxy_1hz.png",
    )
    args = ap.parse_args()

    start = pd.Timestamp(args.start, tz="UTC")
    end = pd.Timestamp(args.end, tz="UTC")
    cutoff = pd.Timestamp(args.christmas_cutoff, tz="UTC")
    ref_start = pd.Timestamp(args.reference_start, tz="UTC")
    ref_end = pd.Timestamp(args.reference_end, tz="UTC")
    first_diaphragm = pd.Timestamp("2025-12-25 00:00:00", tz="UTC")
    second_diaphragm = pd.Timestamp("2026-02-07 00:00:00", tz="UTC")
    cryotrap_leak_start = pd.Timestamp("2026-01-10 00:00:00", tz="UTC")
    cryotrap_leak_end = second_diaphragm

    gps_start = dt_to_gps(start)
    gps_end_inclusive = dt_to_gps(end)
    gps_end_exclusive = gps_end_inclusive + 1

    stage_dir = Path(args.stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        6,
        1,
        figsize=(16, 16),
        sharex=True,
        gridspec_kw={"height_ratios": [1.5, 1.0, 1.2, 1.0, 0.9, 1.1], "hspace": 0.08},
    )

    legend_seen: set[tuple[int, str]] = set()

    b1p_ref_arrays: list[np.ndarray] = []
    b1p_pre_sum = 0.0
    b1p_pre_count = 0
    b1p_between_sum = 0.0
    b1p_between_count = 0
    b1p_after_sum = 0.0
    b1p_after_count = 0
    daily_pre: list[tuple[pd.Timestamp, float]] = []
    daily_leak: list[tuple[pd.Timestamp, float]] = []
    daily_post: list[tuple[pd.Timestamp, float]] = []

    day_start = (gps_start // 86400) * 86400
    day_end = (gps_end_inclusive // 86400) * 86400

    for gps_day in range(day_start, day_end + 86400, 86400):
        gwf = stage_day(gps_day, stage_dir)
        if gwf is None:
            continue
        t0 = max(gps_start, gps_day)
        t1 = min(gps_end_exclusive, gps_day + 86400)
        if t1 <= t0:
            gwf.unlink(missing_ok=True)
            continue

        print(f"Processing day {gps_day} [{t0}, {t1})", flush=True)
        day = read_day_1hz(gwf, t0, t1)
        gwf.unlink(missing_ok=True)

        plot_mask = build_piecewise_mask(
            day,
            cutoff=cutoff,
            pre_lock_min=args.pre_christmas_lock_min,
            mode=args.alignment_mode,
            post_lock_min=args.aligned_lock_min,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )
        ref_mask = build_piecewise_mask(
            day,
            cutoff=cutoff,
            pre_lock_min=args.pre_christmas_lock_min,
            mode=args.reference_alignment_mode,
            post_lock_min=args.aligned_lock_min,
            strict_low=args.strict_lock_low,
            strict_high=args.strict_lock_high,
        )

        plot_day = day.loc[plot_mask].copy()
        ref_day = day.loc[
            ref_mask & day["dt"].between(ref_start, ref_end, inclusive="left")
        ].copy()

        if not ref_day.empty and "det_b1p_dc" in ref_day:
            ref_vals = ref_day["det_b1p_dc"].dropna().to_numpy(dtype=float)
            if len(ref_vals):
                b1p_ref_arrays.append(ref_vals)

        # Panel 1 raw 1 Hz B1p.
        b1p_sel = plot_day.loc[plot_day["det_b1p_dc"].notna() & (plot_day["det_b1p_dc"] > 50)]
        if not b1p_sel.empty:
            key = (0, CHANNELS["det_b1p_dc"])
            plot_points(
                axes[0],
                b1p_sel["dt"],
                b1p_sel["det_b1p_dc"],
                color="#1f77b4",
                label=CHANNELS["det_b1p_dc"] if key not in legend_seen else None,
                ms=0.6,
            )
            legend_seen.add(key)

            pre_mask = b1p_sel["dt"] < first_diaphragm
            between_mask = (b1p_sel["dt"] >= first_diaphragm) & (b1p_sel["dt"] < second_diaphragm)
            after_mask = b1p_sel["dt"] >= second_diaphragm
            if pre_mask.any():
                vals = b1p_sel.loc[pre_mask, "det_b1p_dc"].to_numpy(dtype=float)
                b1p_pre_sum += float(np.nansum(vals))
                b1p_pre_count += int(np.isfinite(vals).sum())
                daily_pre.append((b1p_sel["dt"].iloc[0].floor("D"), float(np.nanmean(vals))))
            if between_mask.any():
                vals = b1p_sel.loc[between_mask, "det_b1p_dc"].to_numpy(dtype=float)
                b1p_between_sum += float(np.nansum(vals))
                b1p_between_count += int(np.isfinite(vals).sum())
                d0 = b1p_sel.loc[between_mask, "dt"].iloc[0].floor("D")
                daily_leak.append((d0, float(np.nanmean(vals))))
            if after_mask.any():
                vals = b1p_sel.loc[after_mask, "det_b1p_dc"].to_numpy(dtype=float)
                b1p_after_sum += float(np.nansum(vals))
                b1p_after_count += int(np.isfinite(vals).sum())
                d0 = b1p_sel.loc[after_mask, "dt"].iloc[0].floor("D")
                daily_post.append((d0, float(np.nanmean(vals))))

        # Panel 2 Hrec.
        s = plot_day[["dt", "hrec_range_bns"]].dropna()
        if not s.empty:
            key = (1, CHANNELS["hrec_range_bns"])
            plot_points(
                axes[1],
                s["dt"],
                s["hrec_range_bns"],
                color="#2c7fb8",
                label=CHANNELS["hrec_range_bns"] if key not in legend_seen else None,
            )
            legend_seen.add(key)

        # Panel 3 NI temperatures.
        for col, color in [
            ("ni_bottom_te1", "tab:red"),
            ("ni_mir_coil_ul_te", "tab:orange"),
            ("ni_rh_te", "tab:purple"),
            ("env_co2_ni_te", "darkgreen"),
        ]:
            s = plot_day[["dt", col]].dropna()
            if s.empty:
                continue
            key = (2, CHANNELS[col])
            plot_points(
                axes[2],
                s["dt"],
                s[col],
                color=color,
                label=CHANNELS[col] if key not in legend_seen else None,
            )
            legend_seen.add(key)

        # Panel 4 additional NI thermals.
        for col, color in [
            ("ni_f0_te1", "tab:brown"),
            ("ni_f0_te2", "peru"),
            ("ni_f4_te1", "tab:green"),
            ("ni_f4_te2", "limegreen"),
            ("ni_f7_te1", "tab:cyan"),
            ("ni_f7_te2", "deepskyblue"),
        ]:
            s = plot_day[["dt", col]].dropna()
            if s.empty:
                continue
            key = (3, CHANNELS[col])
            plot_points(
                axes[3],
                s["dt"],
                s[col],
                color=color,
                label=CHANNELS[col] if key not in legend_seen else None,
            )
            legend_seen.add(key)

        # Panel 5 NI laser temperatures.
        for col, color in [
            ("tcs_ni_co2laser_te", "tab:blue"),
            ("tcs_ni_auxlaser_te", "tab:olive"),
        ]:
            s = plot_day[["dt", col]].dropna()
            if s.empty:
                continue
            key = (4, CHANNELS[col])
            plot_points(
                axes[4],
                s["dt"],
                s[col],
                color=color,
                label=CHANNELS[col] if key not in legend_seen else None,
            )
            legend_seen.add(key)

        # Panel 6 etalon.
        for col, color in [
            ("ni_rh_out", "tab:blue"),
            ("ni_rh_set", "tab:green"),
            ("rh_ni_in", "goldenrod"),
            ("rh_ni_err", "crimson"),
        ]:
            s = plot_day[["dt", col]].dropna()
            if s.empty:
                continue
            key = (5, CHANNELS[col])
            plot_points(
                axes[5],
                s["dt"],
                s[col],
                color=color,
                label=CHANNELS[col] if key not in legend_seen else None,
            )
            legend_seen.add(key)

    if not b1p_ref_arrays:
        raise RuntimeError("No per-second B1p reference samples were collected.")

    b1p_ref = float(np.nanmedian(np.concatenate(b1p_ref_arrays)))
    b1p_70 = float(args.fixed_b1p_70) if args.fixed_b1p_70 is not None else 0.70 * b1p_ref
    b1p_50 = float(args.fixed_b1p_50) if args.fixed_b1p_50 is not None else 0.50 * b1p_ref

    # First-panel overlays.
    ax = axes[0]

    if b1p_pre_count > 0:
        mean_pre = b1p_pre_sum / b1p_pre_count
        ax.plot([start, first_diaphragm], [mean_pre, mean_pre], color="crimson", lw=1.8, zorder=3)
        ax.text(start + (first_diaphragm - start) / 2, mean_pre + 8, f"{mean_pre:.1f} mW", color="crimson", fontsize=10, ha="center", va="bottom", backgroundcolor="white")
    if b1p_between_count > 0:
        mean_between = b1p_between_sum / b1p_between_count
        ax.plot([first_diaphragm, second_diaphragm], [mean_between, mean_between], color="crimson", lw=1.8, zorder=3)
        ax.text(first_diaphragm + (second_diaphragm - first_diaphragm) / 2, mean_between + 8, f"{mean_between:.1f} mW", color="crimson", fontsize=10, ha="center", va="bottom", backgroundcolor="white")
    if b1p_after_count > 0:
        mean_after = b1p_after_sum / b1p_after_count
        ax.plot([second_diaphragm, end], [mean_after, mean_after], color="crimson", lw=1.8, zorder=3)
        ax.text(second_diaphragm + 0.84 * (end - second_diaphragm), mean_after + 8, f"{mean_after:.1f} mW", color="crimson", fontsize=10, ha="center", va="bottom", backgroundcolor="white")

    def add_mod_label(day_vals: list[tuple[pd.Timestamp, float]], x0: pd.Timestamp, x1: pd.Timestamp, color: str) -> None:
        if len(day_vals) < 2:
            return
        vals = np.asarray([v for _, v in day_vals], dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        if len(vals) < 2:
            return
        factor = float(vals.max() / vals.min())
        ax.text(
            x0 + (x1 - x0) / 2,
            float(vals.max()) + 14,
            f"x{factor:.1f} modulation",
            color=color,
            fontsize=10,
            ha="center",
            va="bottom",
            backgroundcolor="white",
        )

    add_mod_label(daily_pre, start, first_diaphragm, "indianred")
    add_mod_label(daily_leak, cryotrap_leak_start, cryotrap_leak_end, "crimson")
    add_mod_label(daily_post, second_diaphragm, end, "firebrick")

    top_label_y = 1.05
    ax.text(start + (first_diaphragm - start) / 2, top_label_y, "SR NOT ALIGNED", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11, color="black", clip_on=False)
    ax.text(first_diaphragm + (second_diaphragm - first_diaphragm) / 2, top_label_y, "SR ALIGNED", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11, color="black", clip_on=False)
    ax.text(second_diaphragm + (end - second_diaphragm) / 2, top_label_y, "SR ALIGNED", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=11, color="black", clip_on=False)

    for idx, axis in enumerate(axes):
        add_marker(axis, first_diaphragm, "First diaphragm" if idx else "First diaphragm\n2025-12-25")
        add_marker(axis, second_diaphragm, "Second diaphragm" if idx else "Second diaphragm\n2026-02-07")
        axis.set_xlim(start, end)
        axis.grid(alpha=0.25)

    axes[0].set_ylabel("B1p carrier\n[mW]")
    axes[0].set_ylim(0, 500)
    finalize_legend(axes[0], loc="upper left", fontsize=8, ncol=2, framealpha=0.9)

    axes[1].set_ylabel("Hrec range BNS\n[Mpc]")
    finalize_legend(axes[1], loc="lower left", fontsize=8, framealpha=0.9)

    axes[2].set_ylabel("Temperature\n[°C]")
    finalize_legend(axes[2], loc="lower left", fontsize=8, ncol=2, framealpha=0.9)

    axes[3].set_ylabel("Additional NI thermal\n[degC]")
    finalize_legend(axes[3], loc="lower left", fontsize=8, ncol=3, framealpha=0.9)

    axes[4].set_ylabel("NI laser TE\n[degC]")
    finalize_legend(axes[4], loc="center left", fontsize=8, framealpha=0.9)

    axes[5].set_ylabel("NI etalon\nraw value")
    axes[5].set_ylim(0, 3000)
    finalize_legend(axes[5], loc="upper left", fontsize=8, ncol=2, framealpha=0.9)

    axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[-1].tick_params(axis="x", rotation=25, labelsize=9)
    for axis in axes[:-1]:
        axis.tick_params(labelsize=9)

    post_selection_label = {
        "proxy": f"post-Christmas proxy mask: lock_state >= {args.aligned_lock_min:.0f}",
        "strict": (
            f"post-Christmas strict mask: {args.strict_lock_low:.0f} <= lock_state <= "
            f"{args.strict_lock_high:.0f}"
        ),
    }[args.alignment_mode]
    fig.suptitle(
        "Slide-6 cryogenic-trap-leak focus from 1 Hz trend points\n"
        f"pre-Christmas lock_state >= {args.pre_christmas_lock_min:.0f}; {post_selection_label}",
        fontsize=12,
        y=0.995,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved 1 Hz plot -> {out}", flush=True)
    print(f"Reference B1p median: {b1p_ref:.3f} mW", flush=True)


if __name__ == "__main__":
    main()
