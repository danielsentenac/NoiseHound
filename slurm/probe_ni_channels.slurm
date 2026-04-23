#!/usr/bin/env bash
#SBATCH -J nh-probe-ni
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-probe-ni-%j.out

IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_probe_ni_${SLURM_JOB_ID}"
mkdir -p "$STAGE_DIR"

"$IGWN_PYTHON" - << 'PYEOF'
import subprocess, os
from pathlib import Path

GPS_DAY = 1448928000   # 2025-12-01
STAGE_DIR = Path(os.environ.get("STAGE_DIR", "/tmp/nh_probe_ni"))
local = STAGE_DIR / f"V-trend-{GPS_DAY}-86400.gwf"

if not local.exists():
    r = subprocess.run(["rfcp",
        f"cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2025/V-trend-{GPS_DAY}-86400.gwf",
        str(local)], capture_output=True, text=True)
    print("rfcp:", r.returncode, r.stderr[:300], flush=True)

from gwpy.io.gwf import get_channel_names
channels = get_channel_names(str(local))

ni = [c for c in channels if "_NI_" in c or c.startswith("V1:NI_")]
print(f"\nTotal NI channels: {len(ni)}")

import re
te_pat  = re.compile(r"TE|TEMP|temp|therm|THERM", re.I)
pwr_pat = re.compile(r"PWR|CURR|VOLT|RH|HEAT|_I$|_I_|_P$|_P_|_W$", re.I)

ni_te  = [c for c in ni if te_pat.search(c)]
ni_pwr = [c for c in ni if pwr_pat.search(c)]

print(f"\n=== ALL NI CHANNELS ({len(ni)}) ===")
for c in sorted(ni): print(c)

print(f"\n=== NI TEMPERATURE CHANNELS ({len(ni_te)}) ===")
for c in sorted(ni_te): print(c)

print(f"\n=== NI POWER/CURRENT CHANNELS ({len(ni_pwr)}) ===")
for c in sorted(ni_pwr): print(c)
PYEOF

rm -rf "$STAGE_DIR"
