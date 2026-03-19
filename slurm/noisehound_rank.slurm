#!/usr/bin/env bash
#SBATCH -J nh-rank
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o slurm-nh-rank-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
VENV_DIR="${VENV_DIR:-$HOME/.venvs/noisehound}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
PYTHON_BIN="${PYTHON_BIN:-}"
AUX_FRAME_LIST_FILE="${AUX_FRAME_LIST_FILE:?set AUX_FRAME_LIST_FILE to a file containing one auxiliary frame per line}"
TRIGGERS="${TRIGGERS:?set TRIGGERS to a triggers.csv path}"
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
OUTPUT_DIR="${OUTPUT_DIR:-$WORKDIR/outputs/rank}"

resolve_python_bin() {
  if [[ -n "$PYTHON_BIN" ]]; then
    printf '%s\n' "$PYTHON_BIN"
    return
  fi
  if [[ -x "$IGWN_PYTHON" ]]; then
    printf '%s\n' "$IGWN_PYTHON"
    return
  fi
  if [[ -x "$VENV_DIR/bin/python" ]]; then
    printf '%s\n' "$VENV_DIR/bin/python"
    return
  fi
  echo "no usable Python interpreter found; set PYTHON_BIN explicitly" >&2
  exit 2
}

PYTHON_BIN="$(resolve_python_bin)"
if [[ "$PYTHON_BIN" == "$VENV_DIR/bin/python" ]]; then
  source "$VENV_DIR/bin/activate"
fi

cd "$WORKDIR"
mkdir -p "$OUTPUT_DIR"
export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

echo "python bin: $PYTHON_BIN"

mapfile -t FRAMES < "$AUX_FRAME_LIST_FILE"
FIRST_FRAME="${FRAMES[0]:-}"

case "$FIRST_FRAME" in
  osdf://*|https://*)
    echo "AUX_FRAME_LIST_FILE contains remote URLs, but rank expects local .gwf paths." >&2
    echo "Stage the frames first or replace the list with local files accessible on CC." >&2
    exit 2
    ;;
  chpss0:*|cchpss0:*|/hpss/*)
    echo "AUX_FRAME_LIST_FILE contains HPSS paths, but rank expects local .gwf paths." >&2
    echo "Stage them first with scripts/cc_stage_hpss.sh or replace the list with local files." >&2
    exit 2
    ;;
esac

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

if [[ -n "$CHANNEL_FILE" ]]; then
  ARGS+=(--channel-file "$CHANNEL_FILE")
fi
if [[ -n "$INCLUDE" ]]; then
  ARGS+=(--include "$INCLUDE")
fi
if [[ -n "$EXCLUDE" ]]; then
  ARGS+=(--exclude "$EXCLUDE")
fi
if [[ -n "$LIMIT" ]]; then
  ARGS+=(--limit "$LIMIT")
fi
if [[ -n "$MAX_SAMPLE_RATE" ]]; then
  ARGS+=(--max-sample-rate "$MAX_SAMPLE_RATE")
fi

"$PYTHON_BIN" "${ARGS[@]}" "${FRAMES[@]}"
