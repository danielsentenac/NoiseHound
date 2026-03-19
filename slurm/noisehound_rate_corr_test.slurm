#!/usr/bin/env bash
#SBATCH -J nh-rate-corr-test
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o slurm-nh-rate-corr-test-%j.out
#
# Short validation job: 3 days (Apr 1-4 2023) to confirm per-channel read fix.
# Submit: sbatch noisehound_rate_corr_test.slurm

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
TRIGGERS="$WORKDIR/data/full_25min_glitches_ER16-O4b.csv"
OUTPUT_DIR="$WORKDIR/outputs/rate_correlation"
STAGE_DIR="/tmp/nh_test_$$"

# 3 days: 2023-04-01 00:00 UTC → 2023-04-04 00:00 UTC
GPS_START=1364774400
GPS_END=1365033600

mkdir -p "$OUTPUT_DIR" "$STAGE_DIR"

echo "=== TEST JOB: GPS $GPS_START – $GPS_END (3 days) $(date -u '+%H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/rate_correlation_direct.py" \
    --triggers   "$TRIGGERS" \
    --gps-start  "$GPS_START" \
    --gps-end    "$GPS_END" \
    --output     "$OUTPUT_DIR/partial_test.csv" \
    --hpss-base  "$HPSS_BASE" \
    --stage-dir  "$STAGE_DIR" \
    --bin-size-s 3600

echo "=== TEST done $(date -u '+%H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
