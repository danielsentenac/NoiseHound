#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path

STAGE = Path("/tmp/nh_probe_lsc_channels")
STAGE.mkdir(exist_ok=True)
gwf = STAGE / "probe.gwf"

if not gwf.exists():
    r = subprocess.run(
        ["rfcp",
         "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2023/V-trend-1369094400-86400.gwf",
         str(gwf)],
        capture_output=True, text=True)
    if r.returncode != 0:
        print("rfcp failed:", r.stderr[:200])
        sys.exit(1)
    print(f"Staged {gwf}")

from gwpy.io.gwf import iter_channel_names
chs = list(iter_channel_names(str(gwf)))
print(f"Total channels in file: {len(chs)}")

lsc  = [c for c in chs if "LSC" in c]
imc  = [c for c in lsc if "IMC" in c]
darm = [c for c in lsc if "DARM" in c]

print(f"\nLSC channels total: {len(lsc)}")
print(f"LSC+IMC channels: {len(imc)}")
for c in imc[:40]:
    print(" ", c)

print(f"\nLSC+DARM channels: {len(darm)}")
for c in darm[:40]:
    print(" ", c)

# Also check for the specific channel
target = "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean"
print(f"\nExact match '{target}': {target in chs}")
# fuzzy
candidates = [c for c in chs if "IMC_LINE" in c or "100Hz" in c]
print(f"Channels with 'IMC_LINE' or '100Hz': {len(candidates)}")
for c in candidates[:20]:
    print(" ", c)
