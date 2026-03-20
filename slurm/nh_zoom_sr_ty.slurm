#!/usr/bin/env bash
#SBATCH -J nh-zoom-sr-ty
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-sr-ty-%j.out

# High-resolution zoom + Lomb-Scargle periodogram for V1:ASC_SR_TY_ERR_mean
# Compares 1 week before Christmas shutdown vs 1 week after ITF relock

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_sr_$$"
mkdir -p "$STAGE_DIR"

echo "=== SR TY zoom job started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/zoom_asc_sr_ty.py" \
    --output    "$WORKDIR/outputs/zoom_asc_sr_ty.png" \
    --stage-dir "$STAGE_DIR"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
