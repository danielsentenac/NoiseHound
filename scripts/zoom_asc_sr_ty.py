"""
zoom_asc_sr_ty.py — Periodogram zoom for any ITF ASC channel.

Produces a 3-row plot comparing:
  Row 0 — Pre-baffle  LOW_NOISE_3         (SR misaligned, glitches present)
  Row 1 — Post-baffle LOW_NOISE_3         (SR misaligned, same config, glitches gone?)
  Row 2 — Post-baffle LOW_NOISE_3_ALIGNED (SR aligned, probe glitch/mystery-noise coupling)

Lock states (V1:META_ITF_LOCK_index):
  134 = ACQUIRE_LOW_NOISE_3       135 = LOW_NOISE_3
  144 = ACQUIRE_LOW_NOISE_3_ALIGNED  145 = LOW_NOISE_3_ALIGNED

Usage (on CCA):
    python scripts/zoom_asc_sr_ty.py \\
        --channel   V1:ASC_SR_TY_ERR_mean \\
        --ylabel    "SR TY error [rad]" \\
        --output    outputs/zoom_asc_sr_ty.png \\
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

# astropy degC patch
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
CHANNEL      = "V1:ASC_SR_TY_ERR_mean"
LOCK_CHANNEL = "V1:META_ITF_LOCK_index"
GLITCH_PERIOD_MIN = 25.0

# GPS windows
GPS_BEFORE = (1448323200, 1448928000)   # 24 Nov – 01 Dec 2025 (pre-baffle)
GPS_AFTER  = (1452456000, 1458000000)   # 14 Jan – 20 Mar 2026 (post-baffle, wide)

# Lock index ranges per state
LN3         = (133.0, 136.0)   # LOW_NOISE_3 (SR misaligned)
LN3_ALIGNED = (143.0, 146.0)   # LOW_NOISE_3_ALIGNED (SR aligned)

BIN_MIN = 5   # 5-minute averages


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
                   channel: str = CHANNEL,
                   lock_range: Optional[tuple] = None) -> Optional[pd.DataFrame]:
    """Extract channel for [gps_start, gps_end), binned to BIN_MIN minutes.

    If lock_range=(lo, hi) is given, also reads V1:META_ITF_LOCK_index and
    keeps only 5-min bins where the mean lock index is in [lo, hi].
    """
    from gwpy.timeseries import TimeSeries

    day_start = (gps_start // 86400) * 86400
    day_end   = (gps_end   // 86400) * 86400

    times, vals, locks = [], [], []
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
            t_arr = np.asarray(ts.times.value, float)
            v_arr = np.asarray(ts.value, float)
            times.append(t_arr)
            vals.append(v_arr)
            if lock_range is not None:
                try:
                    lk = TimeSeries.read(str(gwf), LOCK_CHANNEL, start=t0, end=t1)
                    locks.append(np.asarray(lk.value, float))
                except Exception:
                    locks.append(np.full(len(t_arr), np.nan))
        except Exception as e:
            print(f"  read failed {gps_day}: {e}", flush=True)

    if not times:
        return None

    t = np.concatenate(times)
    v = np.concatenate(vals)
    df = pd.DataFrame({"gps": t, "val": v})
    if lock_range is not None and locks:
        df["lock"] = np.concatenate(locks)

    # Bin to BIN_MIN-minute averages
    bin_s = BIN_MIN * 60
    df["bin"] = (df["gps"] // bin_s) * bin_s

    if lock_range is not None and "lock" in df.columns:
        agg = df.groupby("bin").agg(val=("val", "mean"), lock=("lock", "mean")).reset_index()
        agg.columns = ["gps", "val", "lock"]
        lo, hi = lock_range
        n_total = len(agg)
        agg = agg[(agg["lock"] >= lo) & (agg["lock"] <= hi)].drop(columns="lock")
        n_kept = len(agg)
        print(f"  lock filter [{lo:.0f},{hi:.0f}]: {n_kept}/{n_total} 5-min bins "
              f"({100*n_kept/max(n_total,1):.1f}%)", flush=True)
        if n_kept < 20:
            print("  WARNING: very few bins in this lock state", flush=True)
            return None
        return agg.rename(columns={"bin": "gps"}) if "bin" in agg.columns else agg
    else:
        binned = df.groupby("bin")["val"].mean().reset_index()
        binned.columns = ["gps", "val"]
        return binned


def lombscargle(gps: np.ndarray, val: np.ndarray,
                period_min_range=(5.0, 180.0)):
    """Return (periods_min, power) Lomb-Scargle normalised to unit variance."""
    from astropy.timeseries import LombScargle
    t_s = gps - gps[0]
    y   = val - np.nanmean(val)
    mask = np.isfinite(y)
    t_s, y = t_s[mask], y[mask]
    freq_min = 1.0 / (period_min_range[1] * 60)
    freq_max = 1.0 / (period_min_range[0] * 60)
    freq  = np.linspace(freq_min, freq_max, 5000)
    power = LombScargle(t_s, y).power(freq)
    return 1.0 / freq / 60.0, power   # periods in minutes


def plot(output: str, stage_dir: str,
         channel: str = CHANNEL, ylabel: Optional[str] = None,
         after_gps: Optional[tuple] = None) -> None:
    sd = Path(stage_dir)
    sd.mkdir(parents=True, exist_ok=True)

    after_start, after_end = after_gps if after_gps else GPS_AFTER

    from astropy.time import Time as AstroTime
    def fmt(gps): return AstroTime(gps, format="gps").strftime("%d %b %Y")

    windows = [
        # (label, gps_start, gps_end, lock_range, color)
        ("Pre-baffle — LOW_NOISE_3 (SR misaligned)\n"
         f"({fmt(GPS_BEFORE[0])} – {fmt(GPS_BEFORE[1])})",
         GPS_BEFORE[0], GPS_BEFORE[1], LN3, "tab:blue"),

        ("Post-baffle — LOW_NOISE_3 (SR misaligned)\n"
         f"({fmt(after_start)} – {fmt(after_end)}, locked only)",
         after_start, after_end, LN3, "tab:orange"),

        ("Post-baffle — LOW_NOISE_3_ALIGNED (SR aligned)\n"
         f"({fmt(after_start)} – {fmt(after_end)}, locked only)",
         after_start, after_end, LN3_ALIGNED, "tab:green"),
    ]

    fig = plt.figure(figsize=(14, 15))
    fig.suptitle(f"{channel} — zoom & periodogram\n"
                 f"Before vs after SR baffle (Christmas 2025) — by lock state",
                 fontsize=11)
    gs = GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.35)

    for idx, (label, gps_start, gps_end, lock_range, color) in enumerate(windows):
        print(f"\n[{label.splitlines()[0]}]", flush=True)
        df = extract_window(gps_start, gps_end, sd, channel=channel,
                            lock_range=lock_range)
        if df is None or len(df) < 10:
            print("  No data.", flush=True)
            continue
        dt = gps_to_dt(df["gps"].values)

        ax_ts = fig.add_subplot(gs[idx, 0])
        ax_ts.plot(dt.values, df["val"].values, lw=0.7, color=color)
        ax_ts.set_title(label, fontsize=8)
        ax_ts.set_ylabel(ylabel or channel.split(":")[-1], fontsize=8)
        ax_ts.xaxis.set_major_formatter(
            plt.matplotlib.dates.DateFormatter("%d %b\n%H:%M"))
        ax_ts.tick_params(axis="x", labelsize=7)
        ax_ts.grid(alpha=0.25)

        ax_pg = fig.add_subplot(gs[idx, 1])
        periods, power = lombscargle(df["gps"].values, df["val"].values)
        ax_pg.plot(periods, power, lw=0.8, color=color)
        ax_pg.axvline(GLITCH_PERIOD_MIN, color="crimson", lw=1.2, ls="--",
                      label=f"~{GLITCH_PERIOD_MIN:.0f} min")
        ax_pg.set_xlabel("Period [min]", fontsize=8)
        ax_pg.set_ylabel("LS power", fontsize=8)
        ax_pg.set_title(f"Periodogram — {label.splitlines()[0]}", fontsize=8)
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
    p.add_argument("--channel",   default=None)
    p.add_argument("--ylabel",    default=None)
    p.add_argument("--after-gps-start", type=int, default=None,
                   help="GPS start of post-baffle search window (default: 14 Jan 2026)")
    p.add_argument("--after-gps-end",   type=int, default=None,
                   help="GPS end of post-baffle search window (default: 01 Mar 2026)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    after_gps = None
    if args.after_gps_start and args.after_gps_end:
        after_gps = (args.after_gps_start, args.after_gps_end)
    elif args.after_gps_start or args.after_gps_end:
        raise SystemExit("ERROR: provide both --after-gps-start and --after-gps-end")
    plot(args.output, args.stage_dir,
         channel=args.channel or CHANNEL,
         ylabel=args.ylabel,
         after_gps=after_gps)
