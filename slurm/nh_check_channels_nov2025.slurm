#!/usr/bin/env bash
#SBATCH -J nh-ch-check
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-ch-check-%j.out

# Quickly stage one Nov 2025 raw frame and dump its channel list.
# Specifically check for Sc_DET_MAR_TX_CORR / Sc_DET_MAR_TY_CORR.
#
# Raw file: V-raw-1446106700-100.gwf  (covers GPS 1446106700–1446106800)
# HPSS path: cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw/1446/

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
STAGE_DIR="/tmp/nh_ch_check_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

GWF_FILE="V-raw-1446106700-100.gwf"
SUBDIR="1446"
REMOTE="${HPSS_BASE}/${SUBDIR}/${GWF_FILE}"
LOCAL="${STAGE_DIR}/${GWF_FILE}"
OUT_CH="${WORKDIR}/outputs/asc_glitch_channels_1446106739.txt"

echo "=== Channel check started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    Staging ${REMOTE} ..."
rfcp "${REMOTE}" "${LOCAL}"
echo "    Staged ($(du -sh ${LOCAL} | cut -f1))"

"$IGWN_PYTHON" - << PYEOF
from gwpy.io.gwf import iter_channel_names
from pathlib import Path

local = "${LOCAL}"
out   = "${OUT_CH}"

print("Listing channels ...", flush=True)
ch = list(iter_channel_names(local))
print(f"Total channels: {len(ch)}", flush=True)

Path(out).write_text("\n".join(ch))
print(f"Full list saved -> {out}", flush=True)

# --- targeted check ---
targets = [
    "V1:Sc_DET_MAR_TX_CORR",
    "V1:Sc_DET_MAR_TY_CORR",
    "V1:Sc_IB_MAR_TX_CORR",
    "V1:Sc_IB_MAR_TY_CORR",
    "V1:Sa_OB_F0_TY_CORR_500Hz",
    "V1:Sa_OB_F0_X_CORR_500Hz",
    "V1:ASC_SR_TX_OUT",
    "V1:ASC_SR_TY_OUT",
    "V1:ASC_SR_TX_CORR",
    "V1:ASC_SR_TY_CORR",
]
ch_set = set(ch)
print("\n--- Targeted channel check ---")
for t in targets:
    status = "FOUND" if t in ch_set else "absent"
    print(f"  {status:6s}  {t}")

# Also show all Sc_DET and Sa_OB channels found
print("\n--- All Sc_DET_* channels ---")
for c in sorted(c for c in ch if "Sc_DET" in c):
    print(f"  {c}")
print("\n--- All Sa_OB_* channels ---")
for c in sorted(c for c in ch if "Sa_OB" in c):
    print(f"  {c}")
PYEOF

rm -rf "$STAGE_DIR"
echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
