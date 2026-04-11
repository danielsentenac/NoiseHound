#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
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
    ("V1:TCS_HWS_NI_TE1_mean", "tcs_hws_ni_te1"),
    ("V1:TCS_HWS_NI_TE2_mean", "tcs_hws_ni_te2"),
    ("V1:TCS_NI_TE_CO2Laser", "tcs_ni_co2laser_te"),
]


def stage(stage_dir: Path, gps_day: int) -> Path | None:
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local = stage_dir / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)], capture_output=True, text=True)
        if r.returncode != 0:
            return None
    return local


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gps-start", type=int, required=True)
    ap.add_argument("--gps-end", type=int, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    stage_dir = Path(os.environ.get("STAGE_DIR", "/tmp/nh_ni_tcs_extra"))
    stage_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, float]] = []
    day_start = (args.gps_start // 86400) * 86400
    day_end = (args.gps_end // 86400) * 86400

    for gps_day in range(day_start, day_end + 86400, 86400):
        gwf = stage(stage_dir, gps_day)
        if gwf is None:
            continue
        t0 = max(args.gps_start, gps_day)
        t1 = min(args.gps_end, gps_day + 86400)
        for ch, key in CHANNELS:
            try:
                ts = TimeSeries.read(str(gwf), ch, start=t0, end=t1)
            except Exception:
                continue
            t_arr = np.asarray(ts.times.value, float)
            v_arr = np.asarray(ts.value, float)
            bins = (t_arr // BIN_S).astype(int) * BIN_S
            df_ch = (
                pd.DataFrame({"gps_bin": bins, "val": v_arr})
                .groupby("gps_bin", as_index=False)["val"]
                .mean()
            )
            for row in df_ch.itertuples(index=False):
                records.append({"gps_bin": int(row.gps_bin), key: float(row.val)})
        gwf.unlink(missing_ok=True)

    out = (
        pd.DataFrame(records)
        .groupby("gps_bin", as_index=False)
        .mean()
        .sort_values("gps_bin")
        .reset_index(drop=True)
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(out_path)
    print(list(out.columns))
    print(len(out))


if __name__ == "__main__":
    main()
