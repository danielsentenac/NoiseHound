#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeries

warnings.filterwarnings("ignore")

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
BIN_S = 3600
CHANNELS = [
    ("V1:TCS_NI_TE_CO2Laser", "tcs_ni_co2laser_te"),
    ("V1:TCS_NI_TE_AUXLaser", "tcs_ni_auxlaser_te"),
]


def stage_gwf(stage_dir: Path, gps_day: int) -> Path | None:
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local = stage_dir / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        proc = subprocess.run(
            ["rfcp", remote, str(local)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            msg = proc.stderr.strip() or proc.stdout.strip()
            print(f"rfcp failed for {gps_day}: {msg}", file=sys.stderr, flush=True)
            return None
    return local


def read_hourly(gwf: Path, channel: str, start: int, end: int) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ts = TimeSeries.read(str(gwf), channel, start=start, end=end)
    t_arr = np.asarray(ts.times.value, float)
    v_arr = np.asarray(ts.value, float)
    bins = (t_arr // BIN_S).astype(int) * BIN_S
    return (
        pd.DataFrame({"gps_bin": bins, "val": v_arr})
        .groupby("gps_bin", as_index=False)["val"]
        .mean()
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gps-start", type=int, required=True)
    ap.add_argument("--gps-end", type=int, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    stage_dir = Path(
        tempfile.mkdtemp(prefix="nh_ni_laser_", dir=os.environ.get("TMPDIR", "/tmp"))
    )
    try:
        day_start = (args.gps_start // 86400) * 86400
        day_end = (args.gps_end // 86400) * 86400
        records: list[dict[str, float]] = []
        selected: dict[str, str] = {}

        for gps_day in range(day_start, day_end + 86400, 86400):
            gwf = stage_gwf(stage_dir, gps_day)
            if gwf is None:
                continue
            t0 = max(args.gps_start, gps_day)
            t1 = min(args.gps_end, gps_day + 86400)
            any_ok = False
            for channel, key in CHANNELS:
                try:
                    df_ch = read_hourly(gwf, channel, t0, t1)
                    for gps_bin, val in df_ch.itertuples(index=False):
                        records.append({"gps_bin": int(gps_bin), key: val})
                    selected[key] = channel
                    any_ok = True
                except Exception as exc:
                    print(f"{channel} @ {gps_day}: {exc}", file=sys.stderr, flush=True)
            if any_ok:
                print(f"day {gps_day} OK", flush=True)
            gwf.unlink(missing_ok=True)

        if not records:
            print("No data extracted", file=sys.stderr)
            return 1

        df = (
            pd.DataFrame(records)
            .groupby("gps_bin", as_index=False)
            .mean()
            .sort_values("gps_bin")
        )
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        print(f"Saved {len(df)} bins -> {out}", flush=True)
        for key, channel in selected.items():
            print(f"{key}: {channel}", flush=True)
        return 0
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
