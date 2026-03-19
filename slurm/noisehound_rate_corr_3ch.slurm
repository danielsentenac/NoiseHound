#!/usr/bin/env bash
#SBATCH -J nh-rate-3ch
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH -o slurm-nh-rate-3ch-%j.out
#
# Follow-up array job for the 3 logbook-confirmed ring heater / CEB channels
# that were missing from the original 30-channel run.
#
# Submit: sbatch --array=0-11 noisehound_rate_corr_3ch.slurm

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
TRIGGERS="$WORKDIR/data/full_25min_glitches_ER16-O4b.csv"
OUTPUT_DIR="$WORKDIR/outputs/rate_correlation_3ch"
STAGE_DIR="/tmp/nh_3ch_${SLURM_ARRAY_TASK_ID:-0}_$$"

mkdir -p "$OUTPUT_DIR" "$STAGE_DIR"

BLOCK_STARTS=(
    1364774400   # 2023-04-01
    1372636800   # 2023-07-01
    1380585600   # 2023-10-01
    1388534400   # 2024-01-01
    1396310400   # 2024-04-01
    1404172800   # 2024-07-01
    1412121600   # 2024-10-01
    1420070400   # 2025-01-01
    1427846400   # 2025-04-01
    1435708800   # 2025-07-01
    1443657600   # 2025-10-01
    1451606400   # 2026-01-01
)
BLOCK_ENDS=(
    1372636800   # 2023-07-01
    1380585600   # 2023-10-01
    1388534400   # 2024-01-01
    1396310400   # 2024-04-01
    1404172800   # 2024-07-01
    1412121600   # 2024-10-01
    1420070400   # 2025-01-01
    1427846400   # 2025-04-01
    1435708800   # 2025-07-01
    1443657600   # 2025-10-01
    1451606400   # 2026-01-01
    1459468800   # 2026-04-01
)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
GPS_START="${BLOCK_STARTS[$TASK_ID]}"
GPS_END="${BLOCK_ENDS[$TASK_ID]}"

echo "=== Task $TASK_ID: GPS $GPS_START – $GPS_END $(date -u '+%H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/rate_correlation_direct.py" \
    --triggers   "$TRIGGERS" \
    --gps-start  "$GPS_START" \
    --gps-end    "$GPS_END" \
    --output     "$OUTPUT_DIR/partial_${TASK_ID}.csv" \
    --hpss-base  "$HPSS_BASE" \
    --stage-dir  "$STAGE_DIR" \
    --bin-size-s 3600 \
    --channels   ni_rh_te wi_rh_te ceb_n_te

echo "=== Task $TASK_ID done $(date -u '+%H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
