#!/usr/bin/env bash
#SBATCH -J nh-zoom-sr-temps
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-sr-temps-%j.out

# Periodogram zoom for SR tower temperature channels
# Tests whether INF_SR_MIR_COIL_UL_TE and INF_TCS_SR_RH_TE oscillate at ~25 min
# before the Christmas shutdown, matching the glitch recurrence period.

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_sr_temps_$$"
mkdir -p "$STAGE_DIR"

echo "=== SR temps zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/zoom_asc_sr_ty.py" \
    --channel   "V1:INF_SR_MIR_COIL_UL_TE" \
    --ylabel    "SR mirror coil UL TE [°C]" \
    --output    "$WORKDIR/outputs/zoom_sr_mir_coil_te.png" \
    --stage-dir "$STAGE_DIR"

"$IGWN_PYTHON" "$WORKDIR/scripts/zoom_asc_sr_ty.py" \
    --channel   "V1:INF_TCS_SR_RH_TE" \
    --ylabel    "SR RH thermistor [°C]" \
    --output    "$WORKDIR/outputs/zoom_sr_rh_te.png" \
    --stage-dir "$STAGE_DIR"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
