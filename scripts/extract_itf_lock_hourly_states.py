#!/usr/bin/env python3
"""
Extract hourly summaries of V1:META_ITF_LOCK_index from Virgo trend files.

For each hour, compute:
- lock_mean: arithmetic mean of the 1 Hz lock-state index
- lock_mode: most frequent integer state in the hour
- mode_frac: fraction of 1 Hz samples equal to lock_mode
- frac_145: fraction of samples in LOW_NOISE_3_ALIGNED
- frac_144_145: fraction of samples in ACQUIRE/LOW_NOISE_3_ALIGNED
- frac_143_146: fraction of samples in the old broad strict interval
"""

from __future__ import annotations

import argparse
import subprocess
import warnings
from pathlib import Path

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

from gwpy.timeseries import TimeSeries

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
CHANNEL = "V1:META_ITF_LOCK_index"


def stage_day(gps_day: int, stage_dir: Path) -> Path | None:
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local = stage_dir / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"rfcp failed {gps_day}: {r.stderr[:120]}", flush=True)
            return None
    return local


def summarize_hour(bin_start: int, vals: np.ndarray) -> dict[str, float | int]:
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return {}
    states = np.rint(vals).astype(int)
    uniq, counts = np.unique(states, return_counts=True)
    mode = int(uniq[np.argmax(counts)])
    mode_frac = float(counts.max() / len(states))
    return {
        "gps_bin": int(bin_start),
        "lock_mean": float(np.mean(vals)),
        "lock_mode": mode,
        "mode_frac": mode_frac,
        "frac_145": float(np.mean(states == 145)),
        "frac_144_145": float(np.mean((states >= 144) & (states <= 145))),
        "frac_143_146": float(np.mean((states >= 143) & (states <= 146))),
        "n_samples": int(len(states)),
    }


def extract_day(gwf: Path, t0: int, t1: int) -> list[dict[str, float | int]]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ts = TimeSeries.read(str(gwf), CHANNEL, start=t0, end=t1)
    t_arr = np.asarray(ts.times.value, float)
    v_arr = np.asarray(ts.value, float)
    out: list[dict[str, float | int]] = []
    for bin_start in range((t0 // 3600) * 3600, t1, 3600):
        mask = (t_arr >= bin_start) & (t_arr < bin_start + 3600)
        if not np.any(mask):
            continue
        row = summarize_hour(bin_start, v_arr[mask])
        if row:
            out.append(row)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gps-start", type=int, required=True)
    ap.add_argument("--gps-end", type=int, required=True)
    ap.add_argument("--stage-dir", default="/tmp/nh_itf_lock_hourly_states")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    stage_dir = Path(args.stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, float | int]] = []
    day_start = (args.gps_start // 86400) * 86400
    day_end = (args.gps_end // 86400) * 86400
    for gps_day in range(day_start, day_end + 86400, 86400):
        gwf = stage_day(gps_day, stage_dir)
        if gwf is None:
            continue
        t0 = max(args.gps_start, gps_day)
        t1 = min(args.gps_end, gps_day + 86400)
        if t1 <= t0:
            gwf.unlink(missing_ok=True)
            continue
        try:
            rows.extend(extract_day(gwf, t0, t1))
            print(f"day {gps_day} OK", flush=True)
        except Exception as exc:
            print(f"{CHANNEL} @ {gps_day}: failed ({exc})", flush=True)
        gwf.unlink(missing_ok=True)

    if not rows:
        raise SystemExit("No hourly lock-state summaries extracted.")

    df = (
        pd.DataFrame(rows)
        .sort_values("gps_bin")
        .drop_duplicates(subset=["gps_bin"])
        .reset_index(drop=True)
    )
    df.to_csv(output, index=False)
    print(f"Saved {len(df)} hourly rows -> {output}", flush=True)


if __name__ == "__main__":
    main()
