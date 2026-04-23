#!/usr/bin/env bash
#SBATCH -J nh-slide6-1hz
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-slide6-1hz-%j.out

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
ALIGNMENT_MODE="${ALIGNMENT_MODE:-proxy}"
REFERENCE_ALIGNMENT_MODE="${REFERENCE_ALIGNMENT_MODE:-$ALIGNMENT_MODE}"
export STAGE_DIR="/tmp/nh_slide6_1hz_${SLURM_JOB_ID:-local}_${ALIGNMENT_MODE}"
mkdir -p "$STAGE_DIR"

echo "=== slide6 1Hz plot started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    mode=$ALIGNMENT_MODE ref_mode=$REFERENCE_ALIGNMENT_MODE"

cd "$WORKDIR"
"$IGWN_PYTHON" scripts/plot_slide6_cryotrap_focus_1hz.py \
  --alignment-mode "$ALIGNMENT_MODE" \
  --reference-alignment-mode "$REFERENCE_ALIGNMENT_MODE" \
  --output "usecases/25-minute-glitch/slide6_cryotrap_focus_${ALIGNMENT_MODE}_1hz.png"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
