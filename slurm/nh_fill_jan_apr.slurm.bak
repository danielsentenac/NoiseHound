#!/usr/bin/env bash
#SBATCH -J nh-fill-janapr
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 08:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-fill-janapr-%j.out

# Fill missing Jan 17 – Apr 1 2026 (block 11 timed out at Jan 17)
# GPS 1452038400 (2026-01-17) -> 1459468800 (2026-04-01)
# All 33 channels, 1-hour bins

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
STAGE_DIR="/tmp/nh_fill_janapr_$$"
mkdir -p "$STAGE_DIR"

echo "=== Fill Jan-Apr job started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/rate_correlation_direct.py" \
    --triggers   "$WORKDIR/data/full_25min_glitches_ER16-O4b.csv" \
    --gps-start  1452038400 \
    --gps-end    1459468800 \
    --output     "$WORKDIR/outputs/rate_correlation/partial_janapr.csv" \
    --hpss-base  "$HPSS_BASE" \
    --stage-dir  "$STAGE_DIR" \
    --bin-size-s 3600

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
