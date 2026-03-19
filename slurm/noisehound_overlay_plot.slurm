#!/usr/bin/env bash
#SBATCH -J nh-overlay-plot
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o slurm-nh-overlay-plot-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_HOST="${HPSS_HOST:-cchpss0}"
HPSS_TREND="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2025/V-trend-1419724800-86400.gwf"
DEST_DIR="${DEST_DIR:-/tmp/noisehound_trend_$$}"
OUTPUT_PNG="${OUTPUT_PNG:-$WORKDIR/outputs/detect_1419724818_1419735618/overlay_aux_channels.png}"

mkdir -p "$DEST_DIR"
TREND_FILE="$DEST_DIR/V-trend-1419724800-86400.gwf"

echo "=== staging trend file $(date -u '+%H:%M:%SZ') ==="
rfcp "$HPSS_TREND" "$TREND_FILE"
echo "=== trend staged, running plot $(date -u '+%H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"
"$IGWN_PYTHON" "$WORKDIR/scripts/make_overlay_plot.py" "$TREND_FILE" "$OUTPUT_PNG"

echo "=== done $(date -u '+%H:%M:%SZ') ==="
rm -rf "$DEST_DIR"
