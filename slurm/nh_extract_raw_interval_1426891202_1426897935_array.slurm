#!/usr/bin/env bash
#SBATCH -J nh-raw-1426-shard
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-raw-1426-shard-%A_%a.out

# Array job: one task per 100-second raw GWF file for interval 1426891202–1426897935.
#
# Submit:
#   sbatch --array=0-67 slurm/nh_extract_raw_interval_1426891202_1426897935_array.slurm
#
# Output shards:
#   /scratch/sentenac/noisehound_raw/1426/shards_1426891202_1426897935/raw_interval_1426891202_1426897935_<file_start>.pkl.gz

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935}"
HPSS_BASE="${HPSS_BASE:-cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw}"
GPS_START="${GPS_START:-1426891202}"
GPS_END="${GPS_END:-1426897935}"
BASE_FILE_START="${BASE_FILE_START:-1426891200}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
FILE_START="$((BASE_FILE_START + 100 * TASK_ID))"
SUBDIR="${FILE_START:0:4}"
GWF_NAME="V-raw-${FILE_START}-100.gwf"
REMOTE_GWF="$HPSS_BASE/$SUBDIR/$GWF_NAME"
STAGE_DIR="${STAGE_DIR:-/tmp/nh_raw_interval_${SLURM_ARRAY_JOB_ID:-local}_${TASK_ID}}"
GWF_PATH="$STAGE_DIR/$GWF_NAME"
OUTPUT_PATH="$SHARD_DIR/raw_interval_1426891202_1426897935_${FILE_START}.pkl.gz"

echo "=== raw shard task started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    task_id     = $TASK_ID"
echo "    file_start  = $FILE_START"
echo "    remote_gwf  = $REMOTE_GWF"
echo "    stage_dir   = $STAGE_DIR"
echo "    gwf_path    = $GWF_PATH"
echo "    output_path = $OUTPUT_PATH"

mkdir -p "$STAGE_DIR"
mkdir -p "$SHARD_DIR"
cleanup() {
    rm -rf "$STAGE_DIR"
}
trap cleanup EXIT

echo "Staging $REMOTE_GWF ..."
rfcp "$REMOTE_GWF" "$GWF_PATH"

cd "$WORKDIR"
"$IGWN_PYTHON" "$WORKDIR/scripts/extract_raw_interval_shard.py" \
    --gwf "$GWF_PATH" \
    --gps-start "$GPS_START" \
    --gps-end "$GPS_END" \
    --output "$OUTPUT_PATH"

echo "=== raw shard task finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
