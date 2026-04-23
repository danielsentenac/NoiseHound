#!/usr/bin/env bash
#SBATCH -J nh-itf-lock-hourly
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-itf-lock-hourly-%j.out

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
GPS_START="${GPS_START:?GPS_START is required}"
GPS_END="${GPS_END:?GPS_END is required}"
OUTPUT="${OUTPUT:?OUTPUT is required}"
STAGE_DIR="/tmp/nh_itf_lock_hourly_${SLURM_JOB_ID:-$$}"

mkdir -p "$STAGE_DIR"

echo "=== ITF lock hourly-state extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "GPS_START=$GPS_START"
echo "GPS_END=$GPS_END"
echo "OUTPUT=$OUTPUT"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/extract_itf_lock_hourly_states.py" \
    --gps-start "$GPS_START" \
    --gps-end "$GPS_END" \
    --output "$OUTPUT" \
    --stage-dir "$STAGE_DIR"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
