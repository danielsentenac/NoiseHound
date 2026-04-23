#!/usr/bin/env python3
import sys, os
from pathlib import Path
import numpy as np
from astropy.timeseries import LombScargle
sys.path.insert(0, os.path.expanduser("~/NOISEHOUND"))
from scripts.zoom_asc_sr_ty import extract_window

GPS_BEFORE = (1448323200, 1448928000)
STAGE_DIR  = Path("/tmp/peak_period_stage")
STAGE_DIR.mkdir(exist_ok=True)

print("Extracting NI_BOTTOM_TE1 pre-baffle LN3...", flush=True)
t, y = extract_window(GPS_BEFORE[0], GPS_BEFORE[1], STAGE_DIR,
                      "V1:INF_NI_BOTTOM_TE1", lock_range=(133.0, 136.0))
print(f"  {len(t)} bins, {(t[-1]-t[0])/3600:.1f} h span", flush=True)

freqs = np.linspace(1/40, 1/15, 5000)  # cycles/min
ls = LombScargle(t / 60, y)
power = ls.power(freqs)
periods = 1.0 / freqs
idx_sort = np.argsort(power)[::-1]
print("\nTop-10 LS peaks in 15-40 min range:")
for i in idx_sort[:10]:
    print(f"  period = {periods[i]:.2f} min   power = {power[i]:.5f}")
