#!/usr/bin/env bash
#SBATCH -J nh-fill-gap
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-fill-gap-%j.out

# Fill the 80-day data hole in binned_summary.csv
# Gap: 2025-10-17 -> 2026-01-05 (GPS 1444694400 -> 1451606400)
# All 33 channels, 1-hour bins

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
STAGE_DIR="/tmp/nh_gap_fill_$$"
mkdir -p "$STAGE_DIR"

echo "=== Fill gap job started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/rate_correlation_direct.py" \
    --triggers   "$WORKDIR/data/full_25min_glitches_ER16-O4b.csv" \
    --gps-start  1444694400 \
    --gps-end    1451606400 \
    --output     "$WORKDIR/outputs/rate_correlation/partial_gap.csv" \
    --hpss-base  "$HPSS_BASE" \
    --stage-dir  "$STAGE_DIR" \
    --bin-size-s 3600

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
