"""
extract_itf_lock.py — Extract V1:META_ITF_LOCK_index for Oct 2025 – Apr 2026.

Produces a CSV with hourly bins and the mean lock fraction (fraction of
1-second samples where META_ITF_LOCK_index > 0, i.e. ITF in locked state).

Usage:
    python scripts/extract_itf_lock.py \
        --output outputs/itf_lock_step7.csv \
        --gps-start 1443657600 --gps-end 1459468800
"""
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
CHANNEL   = "V1:META_ITF_LOCK_index"


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


def extract_day(gwf: Path, gps_day: int) -> list[dict]:
    from gwpy.timeseries import TimeSeries
    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = TimeSeries.read(str(gwf), CHANNEL,
                                 start=gps_day, end=gps_day + 86400)
    except Exception as e:
        print(f"  read failed: {e}", flush=True)
        return []
    t = np.asarray(ts.times.value, float)
    v = np.asarray(ts.value, float)
    rows = []
    for bin_start in range(gps_day, gps_day + 86400, 3600):
        mask = (t >= bin_start) & (t < bin_start + 3600)
        if not mask.any():
            continue
        vals = v[mask]
        rows.append(dict(
            gps_bin   = bin_start,
            lock_mean = float(np.nanmean(vals)),
            lock_frac = float(np.mean(vals > 0)),   # fraction of samples locked
        ))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",    default="outputs/itf_lock_step7.csv")
    p.add_argument("--gps-start", type=float, default=1443657600)  # 2025-10-01
    p.add_argument("--gps-end",   type=float, default=1459468800)  # 2026-04-01
    p.add_argument("--stage-dir", default="/tmp/nh_itf_lock")
    args = p.parse_args()

    stage_dir = Path(args.stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Resume support
    done_days: set[int] = set()
    existing_rows: list[dict] = []
    if out.exists():
        df_ex = pd.read_csv(out)
        existing_rows = df_ex.to_dict("records")
        done_days = {int(g // 86400) * 86400 for g in df_ex["gps_bin"]}
        print(f"Resuming: {len(done_days)} days already done", flush=True)

    day_start = int(args.gps_start // 86400) * 86400
    day_end   = int(args.gps_end   // 86400) * 86400
    all_rows  = list(existing_rows)

    for gps_day in range(day_start, day_end, 86400):
        if gps_day in done_days:
            continue
        dt = GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")
        print(f"[{dt.date()}] GPS {gps_day}", flush=True)

        gwf = stage(gps_day, stage_dir)
        if gwf is None:
            continue

        rows = extract_day(gwf, gps_day)
        gwf.unlink(missing_ok=True)

        if rows:
            all_rows.extend(rows)
            df = pd.DataFrame(all_rows).sort_values("gps_bin").drop_duplicates("gps_bin")
            df.to_csv(out, index=False, float_format="%.4f")
            print(f"  {len(rows)} bins written → {out}", flush=True)

    print(f"\nDone. Total bins: {len(all_rows)}", flush=True)


if __name__ == "__main__":
    main()
