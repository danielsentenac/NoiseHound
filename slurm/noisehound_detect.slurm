#!/usr/bin/env bash
#SBATCH -J nh-detect
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o slurm-nh-detect-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
VENV_DIR="${VENV_DIR:-$HOME/.venvs/noisehound}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
PYTHON_BIN="${PYTHON_BIN:-}"
FRAME_LIST_FILE="${FRAME_LIST_FILE:?set FRAME_LIST_FILE to a file containing one Hrec frame per line}"
TARGET_CHANNEL="${TARGET_CHANNEL:-V1:Hrec_hoft_16384Hz}"
GPS_START="${GPS_START:?set GPS_START}"
GPS_END="${GPS_END:?set GPS_END}"
FMIN="${FMIN:-20}"
FMAX="${FMAX:-210}"
STRIDE="${STRIDE:-4}"
FFTLENGTH="${FFTLENGTH:-4}"
OVERLAP="${OVERLAP:-2}"
THRESHOLD="${THRESHOLD:-6}"
BACKGROUND_WINDOW_S="${BACKGROUND_WINDOW_S:-1800}"
MIN_SEPARATION_S="${MIN_SEPARATION_S:-600}"
OUTPUT_DIR="${OUTPUT_DIR:-$WORKDIR/outputs/detect_${GPS_START}_${GPS_END}}"

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
echo "target channel: $TARGET_CHANNEL"

mapfile -t FRAMES < "$FRAME_LIST_FILE"
FIRST_FRAME="${FRAMES[0]:-}"

case "$FIRST_FRAME" in
  osdf://*|https://*)
    echo "FRAME_LIST_FILE contains remote URLs, but detect expects local .gwf paths." >&2
    echo "Stage the frames first or replace the list with local files accessible on CC." >&2
    exit 2
    ;;
  chpss0:*|cchpss0:*|/hpss/*)
    echo "FRAME_LIST_FILE contains HPSS paths, but detect expects local .gwf paths." >&2
    echo "Stage them first with scripts/cc_stage_hpss.sh or replace the list with local files." >&2
    exit 2
    ;;
esac

"$PYTHON_BIN" -m noisehound.cli detect \
  --channel "$TARGET_CHANNEL" \
  --gps-start "$GPS_START" \
  --gps-end "$GPS_END" \
  --fmin "$FMIN" \
  --fmax "$FMAX" \
  --stride "$STRIDE" \
  --fftlength "$FFTLENGTH" \
  --overlap "$OVERLAP" \
  --threshold "$THRESHOLD" \
  --background-window-s "$BACKGROUND_WINDOW_S" \
  --min-separation-s "$MIN_SEPARATION_S" \
  --output "$OUTPUT_DIR/triggers.csv" \
  --plot "$OUTPUT_DIR/triggers.png" \
  "${FRAMES[@]}"
