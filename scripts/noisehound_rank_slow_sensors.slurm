#!/usr/bin/env bash
#SBATCH -J nh-rank-slow
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o slurm-nh-rank-slow-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_TREND="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2025/V-trend-1419724800-86400.gwf"
DEST_DIR="${DEST_DIR:-/tmp/noisehound_trend_$$}"
GPS_START=1419724818
GPS_END=1419735618
TRIGGERS="$WORKDIR/outputs/detect_1419724818_1419735618/triggers_p95_auto.csv"
OUTPUT_DIR="$WORKDIR/outputs/rank_1419724818_1419735618_slow_sensors"

mkdir -p "$DEST_DIR" "$OUTPUT_DIR"
TREND_FILE="$DEST_DIR/V-trend-1419724800-86400.gwf"

echo "=== staging trend file $(date -u '+%H:%M:%SZ') ==="
rfcp "$HPSS_TREND" "$TREND_FILE"
echo "=== trend staged $(date -u '+%H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

# --- Step 1: extended pre-glitch window ranking (slow sensors) ---
echo "--- pre-glitch window ranking (1800s lookback) ---"
"$IGWN_PYTHON" -m noisehound.cli rank \
    "$TREND_FILE" \
    --triggers "$TRIGGERS" \
    --include "^V1:(INF_|ENV_|HVAC_|VAC_)" \
    --gps-start "$GPS_START" \
    --gps-end   "$GPS_END" \
    --window-before-s  1800 \
    --search-before-s  1800 \
    --search-after-s   30 \
    --control-offset-s 3600 \
    --hit-threshold    3.0 \
    --output "$OUTPUT_DIR/ranking_slow_1800s.csv" \
    --top 30

echo "--- standard window ranking (slow sensors, for comparison) ---"
"$IGWN_PYTHON" -m noisehound.cli rank \
    "$TREND_FILE" \
    --triggers "$TRIGGERS" \
    --include "^V1:(INF_|ENV_|HVAC_|VAC_)" \
    --gps-start "$GPS_START" \
    --gps-end   "$GPS_END" \
    --window-before-s  60 \
    --search-before-s  5 \
    --search-after-s   5 \
    --control-offset-s 600 \
    --hit-threshold    3.0 \
    --output "$OUTPUT_DIR/ranking_slow_5s.csv" \
    --top 30

# --- Step 2: causality analysis ---
echo "--- running causality analysis $(date -u '+%H:%M:%SZ') ---"
"$IGWN_PYTHON" "$WORKDIR/scripts/causality_analysis.py" \
    --trend     "$TREND_FILE" \
    --hrec-dir  "$WORKDIR/frame_cache/hoftonline_1419724818_1419735618" \
    --triggers  "$TRIGGERS" \
    --ranking   "$WORKDIR/outputs/rank_trend_1419724800_full/ranking.csv" \
    --gps-start "$GPS_START" \
    --gps-end   "$GPS_END" \
    --output    "$OUTPUT_DIR"

echo "=== done $(date -u '+%H:%M:%SZ') ==="
rm -rf "$DEST_DIR"
