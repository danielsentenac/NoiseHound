#!/usr/bin/env bash
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o slurm-nh-rank-%x-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
AUX_SOURCE_LIST="${AUX_SOURCE_LIST:?set AUX_SOURCE_LIST}"
TRIGGERS="${TRIGGERS:?set TRIGGERS}"
STAGED_LIST="${STAGED_LIST:-$WORKDIR/frame_lists/aux_raw_${SLURM_JOB_ID}_local.txt}"
DEST_DIR="${DEST_DIR:-}"
HPSS_HOST="${HPSS_HOST:-cchpss0}"
HPSS_BASE="${HPSS_BASE:-/hpss/in2p3.fr/group/virgo}"
FORCE="${FORCE:-0}"
CHANNEL_FILE="${CHANNEL_FILE:-}"
INCLUDE="${INCLUDE:-}"
EXCLUDE="${EXCLUDE:-}"
LIMIT="${LIMIT:-}"
WINDOW_BEFORE_S="${WINDOW_BEFORE_S:-60}"
SEARCH_BEFORE_S="${SEARCH_BEFORE_S:-5}"
SEARCH_AFTER_S="${SEARCH_AFTER_S:-5}"
CONTROL_OFFSET_S="${CONTROL_OFFSET_S:-600}"
HIT_THRESHOLD="${HIT_THRESHOLD:-6}"
MAX_SAMPLE_RATE="${MAX_SAMPLE_RATE:-}"
TOP="${TOP:-20}"
OUTPUT_DIR="${OUTPUT_DIR:-$WORKDIR/outputs/rank_${SLURM_JOB_ID}}"

cd "$WORKDIR"
mkdir -p "$OUTPUT_DIR"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

echo "=== start $(date -u '+%Y-%m-%dT%H:%M:%SZ') job=$SLURM_JOB_ID name=$SLURM_JOB_NAME ==="
echo "python:    $IGWN_PYTHON"
echo "triggers:  $TRIGGERS"
echo "channels:  ${CHANNEL_FILE:-pattern=${INCLUDE:-all}}"
echo "output:    $OUTPUT_DIR"

stage_env=(
  SOURCE_LIST="$AUX_SOURCE_LIST"
  OUTPUT_LIST="$STAGED_LIST"
  HPSS_HOST="$HPSS_HOST"
  HPSS_BASE="$HPSS_BASE"
  FORCE="$FORCE"
)
[[ -n "$DEST_DIR" ]] && stage_env+=(DEST_DIR="$DEST_DIR")

echo "--- staging frames $(date -u '+%H:%M:%SZ') ---"
env "${stage_env[@]}" ./scripts/cc_stage_hpss.sh

mapfile -t FRAMES < "$STAGED_LIST"
echo "--- staged ${#FRAMES[@]} frames, starting rank $(date -u '+%H:%M:%SZ') ---"
[[ "${#FRAMES[@]}" -eq 0 ]] && { echo "no frames staged" >&2; exit 2; }

ARGS=(
  -m noisehound.cli rank
  --triggers "$TRIGGERS"
  --window-before-s "$WINDOW_BEFORE_S"
  --search-before-s "$SEARCH_BEFORE_S"
  --search-after-s "$SEARCH_AFTER_S"
  --control-offset-s "$CONTROL_OFFSET_S"
  --hit-threshold "$HIT_THRESHOLD"
  --top "$TOP"
  --output "$OUTPUT_DIR/ranking.csv"
)
[[ -n "$CHANNEL_FILE" ]]    && ARGS+=(--channel-file "$CHANNEL_FILE")
[[ -n "$INCLUDE" ]]         && ARGS+=(--include "$INCLUDE")
[[ -n "$EXCLUDE" ]]         && ARGS+=(--exclude "$EXCLUDE")
[[ -n "$LIMIT" ]]           && ARGS+=(--limit "$LIMIT")
[[ -n "$MAX_SAMPLE_RATE" ]] && ARGS+=(--max-sample-rate "$MAX_SAMPLE_RATE")

"$IGWN_PYTHON" "${ARGS[@]}" "${FRAMES[@]}"

echo "=== done $(date -u '+%Y-%m-%dT%H:%M:%SZ') ==="
