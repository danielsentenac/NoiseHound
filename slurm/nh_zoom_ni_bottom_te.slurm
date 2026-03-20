#!/usr/bin/env bash
#SBATCH -J nh-zoom-ni-te
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-ni-te-%j.out

# Periodogram zoom for V1:INF_NI_BOTTOM_TE1
# Tests whether this channel oscillates at ~25 min before the shutdown
# (same two windows as ASC_SR_TY and SR temperature analysis).

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_ni_te_$$"
mkdir -p "$STAGE_DIR"

echo "=== NI bottom TE zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/zoom_asc_sr_ty.py" \
    --channel   "V1:INF_NI_BOTTOM_TE1" \
    --ylabel    "NI tower bottom TE1 [°C]" \
    --output    "$WORKDIR/outputs/zoom_ni_bottom_te1.png" \
    --stage-dir "$STAGE_DIR"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
