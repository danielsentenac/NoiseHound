#!/usr/bin/env bash
#SBATCH -J nh-stage
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 24:00:00
#SBATCH -o slurm-nh-stage-%j.out

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
SOURCE_LIST="${SOURCE_LIST:?set SOURCE_LIST to a file containing one HPSS path per line}"
DEST_DIR="${DEST_DIR:-}"
OUTPUT_LIST="${OUTPUT_LIST:-$WORKDIR/frame_lists/local_paths.txt}"
HPSS_HOST="${HPSS_HOST:-cchpss0}"
HPSS_BASE="${HPSS_BASE:-/hpss/in2p3.fr/group/virgo}"
FORCE="${FORCE:-0}"
LOG_DIR="${LOG_DIR:-$WORKDIR/logs}"

mkdir -p "$LOG_DIR"
cd "$WORKDIR"

echo "start $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "source list: $SOURCE_LIST"
echo "dest dir: $DEST_DIR"
echo "output list: $OUTPUT_LIST"
echo "hpss host: $HPSS_HOST"
echo "force: $FORCE"
echo "slurm tmpdir: ${SLURM_TMPDIR:-}"
echo "tmpdir: ${TMPDIR:-}"

stage_env=(
  SOURCE_LIST="$SOURCE_LIST"
  OUTPUT_LIST="$OUTPUT_LIST"
  HPSS_HOST="$HPSS_HOST"
  HPSS_BASE="$HPSS_BASE"
  FORCE="$FORCE"
)

if [[ -n "$DEST_DIR" ]]; then
  stage_env+=(DEST_DIR="$DEST_DIR")
fi

env "${stage_env[@]}" ./scripts/cc_stage_hpss.sh

echo "done $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
