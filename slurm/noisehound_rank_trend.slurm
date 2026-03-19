#!/usr/bin/env bash
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o slurm-nh-rank-trend-%x-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_TREND="${HPSS_TREND:?set HPSS_TREND to the HPSS path of the trend gwf file}"
TRIGGERS="${TRIGGERS:?set TRIGGERS}"
HPSS_HOST="${HPSS_HOST:-cchpss0}"
DEST_DIR="${DEST_DIR:-/tmp/noisehound_trend}"
CHANNEL_FILE="${CHANNEL_FILE:-}"
INCLUDE="${INCLUDE:-}"
EXCLUDE="${EXCLUDE:-}"
LIMIT="${LIMIT:-}"
WINDOW_BEFORE_S="${WINDOW_BEFORE_S:-60}"
SEARCH_BEFORE_S="${SEARCH_BEFORE_S:-5}"
SEARCH_AFTER_S="${SEARCH_AFTER_S:-5}"
CONTROL_OFFSET_S="${CONTROL_OFFSET_S:-600}"
HIT_THRESHOLD="${HIT_THRESHOLD:-6}"
TOP="${TOP:-30}"
OUTPUT_DIR="${OUTPUT_DIR:-$WORKDIR/outputs/rank_trend_${SLURM_JOB_ID}}"

cd "$WORKDIR"
mkdir -p "$OUTPUT_DIR" "$DEST_DIR"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

echo "=== start $(date -u '+%Y-%m-%dT%H:%M:%SZ') job=$SLURM_JOB_ID name=$SLURM_JOB_NAME ==="
echo "trend src: $HPSS_TREND"
echo "triggers:  $TRIGGERS"
echo "output:    $OUTPUT_DIR"

TREND_FILE="$DEST_DIR/$(basename "$HPSS_TREND")"

if [[ ! -f "$TREND_FILE" ]]; then
    echo "--- staging trend file $(date -u '+%H:%M:%SZ') ---"
    rfcp "${HPSS_HOST}:${HPSS_TREND}" "$TREND_FILE"
    echo "--- staged $(du -sh "$TREND_FILE" | cut -f1) $(date -u '+%H:%M:%SZ') ---"
else
    echo "--- trend file already staged: $TREND_FILE ---"
fi

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

echo "--- starting rank $(date -u '+%H:%M:%SZ') ---"
"$IGWN_PYTHON" "${ARGS[@]}" "$TREND_FILE"

echo "=== done $(date -u '+%Y-%m-%dT%H:%M:%SZ') ==="
