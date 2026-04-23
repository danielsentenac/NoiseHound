import subprocess, sys
from pathlib import Path
from gwpy.io.gwf import iter_channel_names

GPS_DAY  = 1443657600
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
STAGE    = Path("/tmp/asc_probe.gwf")

if not STAGE.exists():
    r = subprocess.run(["rfcp",
        f"{HPSS_BASE}/2025/V-trend-{GPS_DAY}-86400.gwf", str(STAGE)],
        capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr)

names = sorted(iter_channel_names(str(STAGE)))

mirrors = ["NI","WI","NE","WE","PR","BS","SR"]
dofs    = ["TX","TY","TZ"]

for m in mirrors:
    # Try _ERR_mean / _OUT_mean / _IN_mean / _mean with TX/TY/TZ
    found = [n for n in names
             if f":ASC_{m}_" in n
             and any(f"_{d}" in n for d in dofs)
             and n.endswith("_mean")]
    print(f"\n=== {m} ===")
    for n in found:
        print(" ", n)
