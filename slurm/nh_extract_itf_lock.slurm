#!/usr/bin/env bash
#SBATCH -J nh-itf-lock
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-itf-lock-%j.out

# Extract V1:META_ITF_LOCK_index for Oct 2025 – Apr 2026 (hourly lock fraction)

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
PYTHON="$HOME/.venvs/noisehound/bin/python3"
STAGE_DIR="/tmp/nh_itf_lock_$$"
mkdir -p "$STAGE_DIR"

echo "=== ITF lock extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

"$PYTHON" "$WORKDIR/scripts/extract_itf_lock.py" \
    --output    "$WORKDIR/outputs/itf_lock_step7.csv" \
    --gps-start 1443657600 \
    --gps-end   1459468800 \
    --stage-dir "$STAGE_DIR"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
