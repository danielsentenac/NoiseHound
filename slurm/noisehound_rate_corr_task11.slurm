#!/usr/bin/env bash
#SBATCH -J nh-rate-corr-t11
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o slurm-nh-rate-corr-t11-%j.out
# Task 11 split: covers 2025-10-01 to 2026-04-01 (180 days)
# Run as array: sbatch --array=0-1 noisehound_rate_corr_task11.slurm

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
TEMP_CSV="$WORKDIR/outputs/rate_correlation/temperature_series.csv"
STAGE_DIR="/tmp/nh_trend_stage_t11_${SLURM_ARRAY_TASK_ID:-0}_$$"

mkdir -p "$WORKDIR/outputs/rate_correlation" "$STAGE_DIR"

BLOCK_STARTS=(1443657600 1451606400)   # 2025-10-01  2026-01-01
BLOCK_ENDS=(  1451606400 1459468800)   # 2026-01-01  2026-04-01

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
GPS_START="${BLOCK_STARTS[$TASK_ID]}"
GPS_END="${BLOCK_ENDS[$TASK_ID]}"

echo "=== Task 11.$TASK_ID: GPS $GPS_START – $GPS_END $(date -u '+%H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/extract_temperature_series.py" \
    --gps-start "$GPS_START" \
    --gps-end   "$GPS_END" \
    --output    "$TEMP_CSV" \
    --hpss-base "$HPSS_BASE" \
    --stage-dir "$STAGE_DIR" \
    --bin-size-s 3600

echo "=== Done $(date -u '+%H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
