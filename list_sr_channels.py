"""
List all channels containing 'SR' from a sample GWF trend file.
Stages one day from HPSS, dumps channel names, then cleans up.
"""
import subprocess, sys
from pathlib import Path

GPS_EPOCH_UNIX = 315964800  # 1980-01-06 UTC in Unix seconds

# Pick a day we know works: 2025-10-01 = GPS 1443657600
GPS_DAY = 1443657600
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
YEAR = 2025
STAGE = Path("/tmp/sr_channel_probe.gwf")

remote = f"{HPSS_BASE}/{YEAR}/V-trend-{GPS_DAY}-86400.gwf"
if not STAGE.exists():
    print(f"Staging {remote} ...", flush=True)
    r = subprocess.run(["rfcp", remote, str(STAGE)], capture_output=True, text=True)
    if r.returncode != 0:
        print("rfcp failed:", r.stderr[:200])
        sys.exit(1)
    print("Staged.", flush=True)

from gwpy.io.gwf import iter_channel_names
names = sorted(iter_channel_names(str(STAGE)))
sr = [n for n in names if "SR" in n.upper()]
print(f"Total channels in file: {len(names)}")
print(f"Channels containing 'SR': {len(sr)}")
for n in sr:
    print(" ", n)
