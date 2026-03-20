"""
rate_correlation_direct.py — stage daily trend GWF files and compute
rate-vs-temperature correlations in a single pass. No intermediate CSV.

Each SLURM array task covers a 3-month block. Tasks write a partial
summary CSV (one row per time bin). The final task (or a manual step)
merges all partial results and produces plots + correlation table.

Usage (single block):
    python rate_correlation_direct.py \
        --triggers   /path/to/full_25min_glitches_ER16-O4b.csv \
        --gps-start  1388534400 \
        --gps-end    1396310400 \
        --output     /path/to/outputs/rate_correlation/partial_3.csv \
        --hpss-base  cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend \
        --bin-size-s 3600

Merge + plot (after all tasks done):
    python rate_correlation_direct.py --merge-only \
        --triggers   /path/to/full_25min_glitches_ER16-O4b.csv \
        --output-dir /path/to/outputs/rate_correlation
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import correlate, correlation_lags

# gwpy/astropy: patch converters_and_unit so channels with unit 'degC'
# (an UnrecognizedUnit) can be read without crashing when gwpy applies
# the GWF bias term (new += bias, plain float + Quantity).
try:
    import astropy.units.quantity as _aq
    from astropy.units import dimensionless_unscaled as _dless
    from astropy.units import UnitConversionError as _UCE
    _orig_cau = _aq.converters_and_unit
    def _safe_cau(function, method, *inputs):
        try:
            return _orig_cau(function, method, *inputs)
        except (_UCE, ValueError):
            return [None] * len(inputs), _dless
    _aq.converters_and_unit = _safe_cau
except Exception:
    pass

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

CHANNELS = [
    # INF infrastructure temperatures (no _mean suffix)
    ("V1:INF_NI_BOTTOM_TE1",           "ni_bottom_te1"),  # ★ logbook r=-0.72
    ("V1:INF_WI_BOTTOM_TE1",           "wi_bottom_te1"),
    ("V1:INF_NI_MIR_COIL_UL_TE",      "ni_mir_coil_te"),
    ("V1:INF_WI_MIR_COIL_DR_TE",      "wi_mir_coil_te"),
    # TCS optics temperatures
    ("V1:TCS_HWS_NI_TE1_mean",         "ni_hws_te1"),
    ("V1:TCS_HWS_NI_TE2_mean",         "ni_hws_te2"),
    ("V1:ENV_TCS_CO2_NI_TE",            "ni_co2_env_te"),   # NI CO2 bench ambient temp
    ("V1:TCS_NI_TE_CO2Laser",          "ni_co2_laser_te"), # NI CO2 laser body temp
    ("V1:TCS_HWS_WI_TE1_mean",         "wi_hws_te1"),
    ("V1:TCS_HWS_WI_TE2_mean",         "wi_hws_te2"),
    ("V1:ENV_TCS_CO2_WI_TE",           "wi_co2_env_te"),   # WI CO2 bench ambient temp
    ("V1:TCS_WI_TE_CO2Laser",          "wi_co2_laser_te"), # WI CO2 laser body temp
    ("V1:TCS_HWS_NE_TE1_mean",         "ne_hws_te1"),     # control
    ("V1:TCS_NI_CO2_PWRLAS_mean",      "ni_co2_pwr"),
    # Ring heater power
    ("V1:LSC_Etalon_NI_RH_SET_mean",   "ni_rh_set"),
    ("V1:LSC_Etalon_NI_RH_OUT_mean",   "ni_rh_out"),
    ("V1:LSC_Etalon_NI_RH_IN_mean",    "ni_rh_in"),
    ("V1:LSC_Etalon_NI_RH_ERR_mean",   "ni_rh_err"),
    ("V1:LSC_Etalon_WI_RH_SET_mean",   "wi_rh_set"),
    ("V1:LSC_Etalon_WI_RH_OUT_mean",   "wi_rh_out"),
    ("V1:LSC_Etalon_WI_RH_ERR_mean",   "wi_rh_err"),
    # Electrical / mains
    ("V1:ENV_NEB_UPS_VOLT_R_mean",     "neb_ups_volt_r"),
    ("V1:ENV_CEB_UPS_VOLT_R_mean",     "ceb_ups_volt_r"),
    ("V1:ENV_WEB_UPS_VOLT_R_mean",     "web_ups_volt_r"),
    ("V1:ENV_MCB_IPS_CURR_T_mean",     "mcb_ips_curr_t"),
    ("V1:ENV_CEB_UPS_CURR_R_mean",     "ceb_ups_curr_r"),  # Excavator rank 15 per-event
    # Suspended end bench geophones (Excavator: SNEB W-E rank 5, N-S rank 20)
    ("V1:SBE_SNEB_GEO_GRWE_raw_mean", "sneb_geo_we"),     # North End Bench W-E
    ("V1:SBE_SNEB_GEO_GRNS_raw_mean", "sneb_geo_ns"),     # North End Bench N-S
    ("V1:SBE_SWEB_GEO_GRWE_raw_mean", "sweb_geo_we"),     # West End Bench W-E (control)
    ("V1:SBE_SWEB_GEO_GRNS_raw_mean", "sweb_geo_ns"),     # West End Bench N-S (control)
    # Ring heater thermistor temperatures (logbook-backed, GWF-confirmed)
    ("V1:INF_TCS_NI_RH_TE",          "ni_rh_te"),        # NI ring heater thermistor °C
    ("V1:INF_TCS_WI_RH_TE",          "wi_rh_te"),        # WI ring heater thermistor °C
    ("V1:ENV_CEB_N_TE",              "ceb_n_te"),         # CEB north ambient temp °C
]

CH_LABELS = {
    "ni_bottom_te1":  "NI tower bottom TE1 [°C]  ★",
    "wi_bottom_te1":  "WI tower bottom TE1 [°C]",
    "ni_mir_coil_te": "NI mirror coil TE [°C]",
    "wi_mir_coil_te": "WI mirror coil TE [°C]",
    "ni_hws_te1":     "NI HWS TE1 [°C]",
    "ni_hws_te2":     "NI HWS TE2 [°C]",
    "ni_co2_env_te":   "NI CO2 bench ambient [°C]",
    "ni_co2_laser_te": "NI CO2 laser body [°C]",
    "wi_hws_te1":     "WI HWS TE1 [°C]",
    "wi_hws_te2":     "WI HWS TE2 [°C]",
    "wi_co2_env_te":   "WI CO2 bench ambient [°C]",
    "wi_co2_laser_te": "WI CO2 laser body [°C]",
    "ne_hws_te1":     "NE HWS TE1 [°C] (control)",
    "ni_co2_pwr":     "NI CO2 laser power [W]",
    "ni_rh_set":      "NI ring heater setpoint [W]",
    "ni_rh_out":      "NI ring heater output [W]",
    "ni_rh_in":       "NI ring heater input [W]",
    "ni_rh_err":      "NI ring heater error",
    "wi_rh_set":      "WI ring heater setpoint [W]",
    "wi_rh_out":      "WI ring heater output [W]",
    "wi_rh_err":      "WI ring heater error",
    "neb_ups_volt_r": "NEB UPS voltage R [V]",
    "ceb_ups_volt_r": "CEB UPS voltage R [V]",
    "web_ups_volt_r": "WEB UPS voltage R [V]",
    "mcb_ips_curr_t": "MCB IPS current T [A]",
    "ceb_ups_curr_r": "CEB UPS current R [A]  (Excavator #15)",
    "sneb_geo_we":    "SNEB geophone W-E [counts]  (Excavator #5)",
    "sneb_geo_ns":    "SNEB geophone N-S [counts]  (Excavator #20)",
    "sweb_geo_we":    "SWEB geophone W-E [counts] (control)",
    "sweb_geo_ns":    "SWEB geophone N-S [counts] (control)",
    "ni_rh_te":      "NI ring heater thermistor [°C]",
    "wi_rh_te":      "WI ring heater thermistor [°C]",
    "ceb_n_te":      "CEB north ambient temp [°C]",
}


# ── helpers ────────────────────────────────────────────────────────────────────

def gps_to_dt(gps):
    return GPS_EPOCH + pd.to_timedelta(np.asarray(gps, float), unit="s")


def gps_to_year(gps: float) -> int:
    return (GPS_EPOCH + pd.to_timedelta(gps, unit="s")).year


def stage(hpss_base: str, gps_day: int, dest: Path) -> Path | None:
    year = gps_to_year(gps_day)
    fname = f"V-trend-{gps_day}-86400.gwf"
    dst = dest / fname
    if dst.exists():
        return dst
    ret = subprocess.run(
        ["rfcp", f"{hpss_base}/{year}/{fname}", str(dst)],
        capture_output=True
    )
    if ret.returncode != 0:
        print(f"  rfcp failed: {ret.stderr.decode()[:120]}", flush=True)
        return None
    return dst


def bin_day(gwf: Path, gps_day: int, bin_size: int,
            trigger_times: np.ndarray) -> pd.DataFrame | None:
    """Read GWF channel by channel, compute per-bin medians + trigger count."""
    try:
        from gwpy.timeseries import TimeSeries
    except ImportError:
        sys.exit("gwpy not found")

    # Read each channel individually to avoid unit-mixing errors
    ch_data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ch_name, label in CHANNELS:
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts = TimeSeries.read(
                    str(gwf), ch_name,
                    start=gps_day, end=gps_day + 86400,
                )
            ch_data[ch_name] = (
                np.asarray(ts.times.value, float),
                np.asarray(ts.value, float),
            )
        except Exception:
            ch_data[ch_name] = (np.array([]), np.array([]))

    n_ok = sum(1 for t, v in ch_data.values() if len(t) > 0)
    print(f"  {n_ok}/{len(CHANNELS)} channels read", flush=True)
    if n_ok == 0:
        return None

    rows = []
    for bin_start in range(gps_day, gps_day + 86400, bin_size):
        bin_end = bin_start + bin_size
        row = {"gps_bin": bin_start}
        # Trigger count in this bin
        row["n_triggers"] = int(np.sum(
            (trigger_times >= bin_start) & (trigger_times < bin_end)
        ))
        # Channel medians
        for ch_name, label in CHANNELS:
            t, v = ch_data[ch_name]
            if len(t) == 0:
                row[label] = np.nan
                continue
            mask = (t >= bin_start) & (t < bin_end)
            row[label] = float(np.nanmedian(v[mask])) if mask.any() else np.nan
        rows.append(row)

    return pd.DataFrame(rows)


# ── correlation + plots ────────────────────────────────────────────────────────

def run_analysis(df: pd.DataFrame, triggers_csv: str, out_dir: Path,
                 bin_size: int) -> None:
    from matplotlib import pyplot as plt
    from matplotlib.dates import DateFormatter

    bin_h = bin_size / 3600.0
    df = df[df["n_triggers"] >= 0].copy()
    df["rate"] = df["n_triggers"] / bin_h
    df = df[df["rate"] > 0]

    print(f"\nBins with triggers: {len(df)}")
    labels = [l for _, l in CHANNELS]
    available = [l for l in labels if l in df.columns]

    # ── Correlation table ──────────────────────────────────────────────────────
    rows = []
    for col in available:
        x = df[col].values
        y = df["rate"].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 10:
            continue
        pr, pp = stats.pearsonr(x[mask], y[mask])
        sr, sp = stats.spearmanr(x[mask], y[mask])
        rows.append(dict(channel=col, label=CH_LABELS.get(col, col),
                         pearson_r=pr, pearson_p=pp,
                         spearman_r=sr, spearman_p=sp, n=mask.sum()))

    corr_df = pd.DataFrame(rows).sort_values("pearson_r", key=abs, ascending=False)
    corr_df.to_csv(out_dir / "rate_correlation_table.csv", index=False, float_format="%.4f")
    print("\nCorrelation table:")
    print(corr_df[["label", "pearson_r", "spearman_r", "n"]].to_string(index=False))

    # ── Lag scan ───────────────────────────────────────────────────────────────
    max_lag = int(round(7 * 24 * 3600 / bin_size))
    lag_rows = []
    for col in available:
        x = df[col].values
        y = df["rate"].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 10:
            continue
        xn = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
        yn = (y - np.nanmean(y)) / (np.nanstd(y) + 1e-12)
        xn[~mask] = 0; yn[~mask] = 0
        cc = correlate(yn, xn, mode="full") / mask.sum()
        lags = correlation_lags(len(yn), len(xn), mode="full")
        m = np.abs(lags) <= max_lag
        best = np.argmax(np.abs(cc[m]))
        lag_rows.append(dict(channel=col,
                             best_lag_h=lags[m][best] * bin_h,
                             best_r=cc[m][best]))

    lag_df = pd.DataFrame(lag_rows).sort_values("best_r", key=abs, ascending=False)
    lag_df.to_csv(out_dir / "lag_scan.csv", index=False, float_format="%.3f")
    print("\nLag scan (positive = sensor leads rate):")
    print(lag_df.to_string(index=False))

    # ── Parabolic lag refinement ────────────────────────────────────────────────
    # For channels with |integer lag| < 10 h, fit a parabola to the three-point
    # neighbourhood of the cross-correlation peak to refine the lag to sub-bin
    # precision.  Formula: fractional offset = (y_{-1} - y_{+1}) /
    # (2*(y_{-1} - 2*y_0 + y_{+1})) in units of one bin step.
    refine_rows = []
    for col in available:
        x = df[col].values
        y = df["rate"].values
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 10:
            continue
        xn = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
        yn = (y - np.nanmean(y)) / (np.nanstd(y) + 1e-12)
        xn[~mask] = 0; yn[~mask] = 0
        cc = correlate(yn, xn, mode="full") / mask.sum()
        lags = correlation_lags(len(yn), len(xn), mode="full")
        m = np.abs(lags) <= max_lag
        cc_m = cc[m]; lags_m = lags[m]
        best = np.argmax(np.abs(cc_m))
        int_lag_bins = int(lags_m[best])
        int_lag_h = int_lag_bins * bin_h
        best_r = float(cc_m[best])
        # Only refine physically motivated lags (< 10 h)
        if abs(int_lag_h) < 10 and 0 < best < len(cc_m) - 1:
            y0, ym1, yp1 = cc_m[best], cc_m[best - 1], cc_m[best + 1]
            denom = ym1 - 2 * y0 + yp1
            frac = (ym1 - yp1) / (2 * denom) if abs(denom) > 1e-12 else 0.0
            refined_lag_h = (int_lag_bins + frac) * bin_h
        else:
            refined_lag_h = int_lag_h
        refine_rows.append(dict(channel=col,
                                int_lag_h=int_lag_h,
                                refined_lag_h=refined_lag_h,
                                best_r=best_r))

    refine_df = (pd.DataFrame(refine_rows)
                 .sort_values("best_r", key=abs, ascending=False))
    refine_df.to_csv(out_dir / "lag_scan_refined.csv", index=False,
                     float_format="%.4f")
    print("\nParabolic lag refinement (|int lag| < 10 h only):")
    phys = refine_df[refine_df["int_lag_h"].abs() < 10]
    print(phys.to_string(index=False))

    # ── Plot 1: time series overview ───────────────────────────────────────────
    dt = gps_to_dt(df["gps_bin"].values)

    # Load full trigger list for recurrence period
    trig = pd.read_csv(triggers_csv)
    t_all = np.sort(trig["time"].values)
    keep = np.concatenate([[True], np.diff(t_all) > 300])
    t = t_all[keep]
    intervals = np.diff(t)
    t_mid = t[:-1] + intervals / 2
    smooth = (pd.Series(intervals).rolling(20, center=True, min_periods=5)
              .median().bfill().ffill().values)

    key_cols = [c for c in ["ni_bottom_te1", "ni_rh_set", "ni_co2_tc"] if c in df.columns]
    n = 2 + len(key_cols)
    fig, axes = plt.subplots(n, 1, figsize=(16, 3 * n), sharex=False)

    axes[0].scatter(gps_to_dt(t_mid), smooth / 60, s=1.5, alpha=0.3, color="tab:blue")
    axes[0].set_ylabel("Recurrence [min]")
    axes[0].set_title("25-min glitch recurrence period (20-trigger rolling median)")
    axes[0].xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    axes[0].grid(alpha=0.3)

    axes[1].plot(dt, df["rate"].values, color="tab:green", lw=0.7)
    axes[1].set_ylabel(f"Rate [/h]")
    axes[1].set_title(f"25-min glitch rate ({bin_size//3600:.0f}-h bins)")
    axes[1].xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    axes[1].grid(alpha=0.3)

    colors = ["tab:red", "tab:orange", "tab:brown"]
    for i, col in enumerate(key_cols):
        ax = axes[2 + i]
        ax.plot(dt, df[col].values, color=colors[i % len(colors)], lw=0.8)
        ax.set_ylabel(CH_LABELS.get(col, col), fontsize=8)
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m"))
        ax.grid(alpha=0.3)

    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out_dir / "timeseries_overview.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── Plot 2: scatter matrix ─────────────────────────────────────────────────
    top = corr_df.head(9)["channel"].tolist()
    nc = min(3, len(top)); nr = int(np.ceil(len(top) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(5 * nc, 4 * nr))
    axes = np.array(axes).flatten()
    for i, col in enumerate(top):
        ax = axes[i]
        x = df[col].values; y = df["rate"].values
        mask = np.isfinite(x) & np.isfinite(y)
        r = corr_df.loc[corr_df["channel"] == col, "pearson_r"].values[0]
        ax.scatter(x[mask], y[mask], s=3, alpha=0.3, color="tab:blue")
        if mask.sum() > 5:
            m_, b_ = np.polyfit(x[mask], y[mask], 1)
            xl = np.array([np.nanmin(x), np.nanmax(x)])
            ax.plot(xl, m_ * xl + b_, "r-", lw=1.5)
        ax.set_xlabel(CH_LABELS.get(col, col), fontsize=8)
        ax.set_ylabel("Rate [/h]", fontsize=8)
        ax.set_title(f"r = {r:.3f}", fontsize=9)
        ax.grid(alpha=0.2)
    for j in range(len(top), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f"Glitch rate vs. slow sensors  ({bin_size//3600:.0f}-h bins)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "scatter_rate_vs_sensors.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── Plot 3: lag scan for top 4 channels ────────────────────────────────────
    top4 = [r for r in lag_df.head(4)["channel"].tolist() if r in available]
    if top4:
        fig, axes = plt.subplots(2, 2, figsize=(14, 8))
        axes = axes.flatten()
        for i, col in enumerate(top4):
            ax = axes[i]
            x = df[col].values; y = df["rate"].values
            mask = np.isfinite(x) & np.isfinite(y)
            xn = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
            yn = (y - np.nanmean(y)) / (np.nanstd(y) + 1e-12)
            xn[~mask] = 0; yn[~mask] = 0
            cc = correlate(yn, xn, mode="full") / mask.sum()
            lags = correlation_lags(len(yn), len(xn), mode="full")
            m = np.abs(lags) <= max_lag
            lag_h = lags[m] * bin_h
            ax.plot(lag_h, cc[m], color="tab:blue", lw=1.2)
            ax.axvline(0, color="gray", lw=0.8, ls="--")
            bi = np.argmax(np.abs(cc[m]))
            ax.axvline(lag_h[bi], color="tab:red", lw=1.0, ls=":",
                       label=f"lag={lag_h[bi]:.1f}h  r={cc[m][bi]:.3f}")
            ax.legend(fontsize=8)
            ax.set_xlabel("Lag [h]  (+ = sensor leads rate)")
            ax.set_ylabel("Cross-correlation r")
            ax.set_title(CH_LABELS.get(col, col), fontsize=9)
            ax.grid(alpha=0.2)
        fig.suptitle("Lag scan: slow sensors vs. glitch rate", fontsize=11)
        fig.tight_layout()
        fig.savefig(out_dir / "lag_scan.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"\nAll outputs saved to {out_dir}")


# ── main ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--triggers",    required=True)
    p.add_argument("--gps-start",   type=float)
    p.add_argument("--gps-end",     type=float)
    p.add_argument("--output",      help="Partial output CSV (extraction mode)")
    p.add_argument("--output-dir",  help="Output directory (merge+plot mode)")
    p.add_argument("--hpss-base",   default="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend")
    p.add_argument("--stage-dir",   default=None)
    p.add_argument("--bin-size-s",  type=int, default=3600)
    p.add_argument("--merge-only",  action="store_true",
                   help="Skip extraction, just merge partial CSVs and plot")
    p.add_argument("--channels",    nargs="+", metavar="LABEL",
                   help="Restrict to these channel short-labels (default: all)")
    p.add_argument("--source-dir",  help="Read binned_summary.csv or partial_*.csv from here "
                   "(merge-only mode); write analysis outputs to --output-dir")
    p.add_argument("--epoch-start", type=float,
                   help="GPS start of analysis epoch (filter binned data, merge-only mode)")
    p.add_argument("--epoch-end",   type=float,
                   help="GPS end of analysis epoch (filter binned data, merge-only mode)")
    return p.parse_args()


def main():
    args = parse_args()

    # ── Optional channel subset ────────────────────────────────────────────────
    if args.channels:
        global CHANNELS, CH_LABELS
        CHANNELS = [(ch, lbl) for ch, lbl in CHANNELS if lbl in args.channels]
        if not CHANNELS:
            sys.exit(f"No matching channels for labels: {args.channels}")
        print(f"Channel filter active: {[lbl for _, lbl in CHANNELS]}", flush=True)

    # ── Merge + plot mode ──────────────────────────────────────────────────────
    if args.merge_only:
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        src_dir = Path(args.source_dir) if args.source_dir else out_dir
        # Prefer an already-merged binned_summary.csv to avoid re-merging
        summary_csv = src_dir / "binned_summary.csv"
        if summary_csv.exists():
            df = pd.read_csv(summary_csv)
            print(f"Loaded {len(df)} bins from {summary_csv}")
        else:
            partials = sorted(src_dir.glob("partial_*.csv"))
            if not partials:
                sys.exit(f"No partial_*.csv or binned_summary.csv found in {src_dir}")
            df = pd.concat([pd.read_csv(f) for f in partials], ignore_index=True)
            df = df.sort_values("gps_bin").drop_duplicates("gps_bin")
            print(f"Merged {len(partials)} partials → {len(df)} bins")
        # Optional epoch filter
        if args.epoch_start:
            df = df[df["gps_bin"] >= args.epoch_start]
        if args.epoch_end:
            df = df[df["gps_bin"] < args.epoch_end]
        if args.epoch_start or args.epoch_end:
            print(f"Epoch filter: {len(df)} bins remaining")
        merged_csv = out_dir / "binned_summary.csv"
        df.to_csv(merged_csv, index=False, float_format="%.6g")
        print(f"Saved {len(df)} bins → {merged_csv}")
        run_analysis(df, args.triggers, out_dir, args.bin_size_s)
        return

    # ── Extraction mode ────────────────────────────────────────────────────────
    assert args.gps_start and args.gps_end and args.output, \
        "Need --gps-start, --gps-end, --output for extraction mode"

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Load triggers once
    trig = pd.read_csv(args.triggers)
    t_all = np.sort(trig["time"].values)
    keep = np.concatenate([[True], np.diff(t_all) > 300])
    trigger_times = t_all[keep]
    print(f"Triggers loaded: {len(trigger_times)}", flush=True)

    # Resume: skip days already in output
    done_days: set[int] = set()
    if out_path.exists():
        existing = pd.read_csv(out_path)
        done_days = set(int(g // 86400) * 86400 for g in existing["gps_bin"])
        print(f"Resuming: {len(done_days)} days already done", flush=True)
    else:
        existing = pd.DataFrame()

    stage_dir = Path(args.stage_dir) if args.stage_dir \
        else Path(tempfile.mkdtemp(prefix="nh_direct_"))
    stage_dir.mkdir(parents=True, exist_ok=True)

    day_start = int(args.gps_start // 86400) * 86400
    day_end   = int(args.gps_end   // 86400) * 86400
    gps_days  = list(range(day_start, day_end + 86400, 86400))
    print(f"Days to process: {len(gps_days)}", flush=True)

    all_frames = [existing] if not existing.empty else []

    for gps_day in gps_days:
        if gps_day in done_days:
            continue
        dt = GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")
        print(f"\n[{dt.strftime('%Y-%m-%d')}] GPS {gps_day}", flush=True)

        gwf = stage(args.hpss_base, gps_day, stage_dir)
        if gwf is None:
            continue

        df_day = bin_day(gwf, gps_day, args.bin_size_s, trigger_times)
        gwf.unlink(missing_ok=True)

        if df_day is not None:
            all_frames.append(df_day)
            pd.concat(all_frames, ignore_index=True)\
              .sort_values("gps_bin")\
              .to_csv(out_path, index=False, float_format="%.6g")
            print(f"  {len(df_day)} bins written → {out_path}", flush=True)

    print(f"\nDone. Total bins: {sum(len(f) for f in all_frames)}", flush=True)


if __name__ == "__main__":
    main()
