#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path

STAGE = Path("/tmp/nh_probe_o4start")
STAGE.mkdir(exist_ok=True)
gwf = STAGE / "probe.gwf"

if not gwf.exists():
    r = subprocess.run(
        ["rfcp",
         "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2023/V-trend-1364774400-86400.gwf",
         str(gwf)],
        capture_output=True, text=True)
    if r.returncode != 0:
        print("rfcp failed:", r.stderr[:200])
        sys.exit(1)

from gwpy.io.gwf import iter_channel_names
chs = list(iter_channel_names(str(gwf)))
target = "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean"
print(f"Total channels: {len(chs)}")
print(f"Target '{target}' present: {target in chs}")
imc = [c for c in chs if "IMC_LINE" in c]
print(f"IMC_LINE channels: {len(imc)}")
for c in imc:
    print(" ", c)
