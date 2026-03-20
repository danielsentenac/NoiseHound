"""
zoom_asc_sr_ty.py — High-resolution zoom + periodogram for V1:ASC_SR_TY_ERR_mean.

Extracts the channel at 1-second resolution for two reference weeks:
  BEFORE: 2025-11-24 – 2025-12-01  (active glitch period, 2 weeks before shutdown)
  AFTER:  2026-01-12 – 2026-01-19  (post-relock, 4 days after ITF back)

For each window:
  - Time series at 5-minute averages
  - Lomb-Scargle periodogram (periods 5 min – 3 h)
  - Highlights the ~25-min glitch recurrence period

Usage (on CCA):
    python scripts/zoom_asc_sr_ty.py \
        --output outputs/zoom_asc_sr_ty.png \
        --stage-dir /tmp/nh_zoom_sr
"""
from __future__ import annotations

import argparse
import subprocess
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

# astropy degC patch (same as other scripts)
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

GPS_EPOCH  = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE  = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
CHANNEL    = "V1:ASC_SR_TY_ERR_mean"   # default; overridden by --channel arg
GLITCH_PERIOD_MIN = 25.0   # known ~25-min glitch recurrence

# Two reference windows (GPS)
WINDOWS = {
    "Before shutdown\n(24 Nov – 01 Dec 2025)": (1448323200, 1448928000),  # 7 days
    "After relock — quiet period\n(01 Feb – 08 Feb 2026)": (1454284800, 1454889600),  # 7 days
}
BIN_MIN = 5   # 5-minute averages for time series


def gps_to_dt(gps):
    return GPS_EPOCH + pd.to_timedelta(np.asarray(gps, float), unit="s")


def stage(gps_day: int, stage_dir: Path) -> Optional[Path]:
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = stage_dir / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  rfcp failed: {r.stderr[:120]}", flush=True)
            return None
    return local


def extract_window(gps_start: int, gps_end: int,
                   stage_dir: Path,
                   channel: str = CHANNEL) -> Optional[pd.DataFrame]:
    """Extract channel at 1-s resolution for [gps_start, gps_end), binned to BIN_MIN."""
    from gwpy.timeseries import TimeSeries

    day_start = (gps_start // 86400) * 86400
    day_end   = (gps_end   // 86400) * 86400

    times, vals = [], []
    for gps_day in range(day_start, day_end + 86400, 86400):
        gwf = stage(gps_day, stage_dir)
        if gwf is None:
            continue
        t0 = max(gps_start, gps_day)
        t1 = min(gps_end,   gps_day + 86400)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts = TimeSeries.read(str(gwf), channel, start=t0, end=t1)
            times.append(np.asarray(ts.times.value, float))
            vals.append(np.asarray(ts.value, float))
        except Exception as e:
            print(f"  read failed {gps_day}: {e}", flush=True)

    if not times:
        return None

    t = np.concatenate(times)
    v = np.concatenate(vals)
    df = pd.DataFrame({"gps": t, "val": v})

    # Bin to BIN_MIN-minute averages
    bin_s = BIN_MIN * 60
    df["bin"] = (df["gps"] // bin_s) * bin_s
    binned = df.groupby("bin")["val"].mean().reset_index()
    binned.columns = ["gps", "val"]
    return binned


def lombscargle(gps: np.ndarray, val: np.ndarray,
                period_min_range=(5.0, 180.0)):
    """Return (periods_min, power) Lomb-Scargle normalised to unit variance."""
    from astropy.timeseries import LombScargle
    t_s   = gps - gps[0]
    y     = val - np.nanmean(val)
    mask  = np.isfinite(y)
    t_s, y = t_s[mask], y[mask]
    freq_min = 1.0 / (period_min_range[1] * 60)
    freq_max = 1.0 / (period_min_range[0] * 60)
    freq = np.linspace(freq_min, freq_max, 5000)
    ls   = LombScargle(t_s, y)
    power = ls.power(freq)
    periods = 1.0 / freq / 60.0   # in minutes
    return periods, power


def plot(output: str, stage_dir: str,
         channel: str = CHANNEL, ylabel: Optional[str] = None) -> None:
    sd = Path(stage_dir)
    sd.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(f"{channel} — zoom & periodogram\n"
                 f"Before vs after SR baffle installation (Christmas 2025 shutdown)",
                 fontsize=11)
    gs = GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

    colors = ["tab:blue", "tab:orange"]
    data = {}

    for idx, (label, (gps_start, gps_end)) in enumerate(WINDOWS.items()):
        print(f"\n[{label.splitlines()[0]}]", flush=True)
        df = extract_window(gps_start, gps_end, sd, channel=channel)
        if df is None or len(df) < 10:
            print("  No data.", flush=True)
            continue
        data[label] = df
        color = colors[idx]
        dt = gps_to_dt(df["gps"].values)

        # Time series
        ax_ts = fig.add_subplot(gs[idx, 0])
        ax_ts.plot(dt.values, df["val"].values, lw=0.7, color=color)
        ax_ts.set_title(label, fontsize=9)
        ax_ts.set_ylabel(ylabel or channel.split(":")[-1], fontsize=8)
        ax_ts.xaxis.set_major_formatter(
            plt.matplotlib.dates.DateFormatter("%d %b\n%H:%M"))
        ax_ts.tick_params(axis="x", labelsize=7)
        ax_ts.grid(alpha=0.25)

        # Periodogram
        ax_pg = fig.add_subplot(gs[idx, 1])
        periods, power = lombscargle(df["gps"].values, df["val"].values)
        ax_pg.plot(periods, power, lw=0.8, color=color)
        ax_pg.axvline(GLITCH_PERIOD_MIN, color="crimson", lw=1.2, ls="--",
                      label=f"~{GLITCH_PERIOD_MIN:.0f} min")
        ax_pg.set_xlabel("Period [min]", fontsize=8)
        ax_pg.set_ylabel("LS power", fontsize=8)
        ax_pg.set_title(f"Periodogram — {label.splitlines()[0]}", fontsize=9)
        ax_pg.legend(fontsize=7)
        ax_pg.grid(alpha=0.25)
        ax_pg.set_xlim(5, 180)

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved → {out}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",    default="outputs/zoom_asc_sr_ty.png")
    p.add_argument("--stage-dir", default="/tmp/nh_zoom_sr")
    p.add_argument("--channel",   default=None,
                   help="Override channel (default: V1:ASC_SR_TY_ERR_mean)")
    p.add_argument("--ylabel",    default=None,
                   help="Y-axis label for time-series panels")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot(args.output, args.stage_dir,
         channel=args.channel or CHANNEL,
         ylabel=args.ylabel)
