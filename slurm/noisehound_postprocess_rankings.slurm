#!/usr/bin/env bash
#SBATCH -J nh-post-rank
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 01:00:00
#SBATCH -o slurm-nh-post-rank-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SUMMARY_DIR="${SUMMARY_DIR:-$WORKDIR/outputs/rank_postprocess}"
SUMMARY_PATH="${SUMMARY_PATH:-$SUMMARY_DIR/summary.md}"
MERGED_CSV="${MERGED_CSV:-$SUMMARY_DIR/merged_rankings.csv}"
RESULT_DIRS="${RESULT_DIRS:?set RESULT_DIRS to a space-separated list of ranking output directories}"

cd "$WORKDIR"
mkdir -p "$SUMMARY_DIR"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

echo "start $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "python bin: $IGWN_PYTHON"
echo "summary dir: $SUMMARY_DIR"
echo "result dirs: $RESULT_DIRS"

read -r -a RESULT_DIR_ARRAY <<< "$RESULT_DIRS"

"$IGWN_PYTHON" scripts/summarize_rankings.py \
  --summary-path "$SUMMARY_PATH" \
  --merged-csv "$MERGED_CSV" \
  "${RESULT_DIR_ARRAY[@]}"

echo "done $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
