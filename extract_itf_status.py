"""Extract V1:META_ITF_STATUS_index for Oct 2025 – Apr 2026 (hourly bins)."""
import subprocess, os, sys
import numpy as np
import pandas as pd
from pathlib import Path
from gwpy.timeseries import TimeSeries

GPS_EPOCH  = pd.Timestamp("1980-01-06", tz="UTC")
HPSS_BASE  = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
CHANNEL    = "V1:META_ITF_STATUS_index"
GPS_START  = 1443657600   # 2025-10-01
GPS_END    = 1459468800   # 2026-04-01
OUTPUT     = "/home/sentenac/NOISEHOUND/outputs/itf_status_step7.csv"
STAGE_DIR  = Path("/tmp/nh_itf_status")
STAGE_DIR.mkdir(exist_ok=True)

def stage(gps_day):
    year = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            return None
    return local

day_start = int(GPS_START // 86400) * 86400
day_end   = int(GPS_END   // 86400) * 86400
rows = []

for gps_day in range(day_start, day_end, 86400):
    dt = GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")
    print(f"{dt.date()} ...", flush=True)
    gwf = stage(gps_day)
    if gwf is None:
        print("  rfcp failed", flush=True)
        continue
    try:
        ts = TimeSeries.read(str(gwf), CHANNEL)
    except Exception as e:
        print(f"  read failed: {e}", flush=True)
        continue
    t = ts.times.value
    v = ts.value
    for bin_start in range(gps_day, gps_day + 86400, 3600):
        mask = (t >= bin_start) & (t < bin_start + 3600)
        if not mask.any():
            continue
        vals = v[mask]
        rows.append(dict(
            gps_bin        = bin_start,
            itf_status_mean= float(np.nanmean(vals)),
            science_frac   = float(np.mean(vals == 130)),
        ))

df = pd.DataFrame(rows)
df.to_csv(OUTPUT, index=False, float_format="%.4f")
print(f"\nSaved {len(df)} rows → {OUTPUT}")
