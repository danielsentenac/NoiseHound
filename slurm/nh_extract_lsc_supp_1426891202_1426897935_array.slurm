#!/usr/bin/env bash
#SBATCH -J nh-lsc-supp-1426
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-lsc-supp-1426-%A_%a.out

# Supplemental shard extraction: LSC longitudinal correction channels only.
# Fetches the same raw GWF files from HPSS but extracts only the new LSC
# channels (LSC_DARM_ERR, LSC_MICH_CORR, LSC_PRCL_CORR, LSC_SRCL_CORR,
# LSC_CARM_CORR, LSC_BS/PR/SR/NI/WI/NE/WE_CORR) into a separate shard dir.
# The render job merges both shard sets.
#
# Submit:
#   sbatch --array=0-67 slurm/nh_extract_lsc_supp_1426891202_1426897935_array.slurm

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-/scratch/sentenac/noisehound_raw/raw_interval_shards_1426891202_1426897935_lsc_supp}"
GWF_DIR="${GWF_DIR:-/scratch/sentenac/noisehound_raw/1426}"
GPS_START="${GPS_START:-1426891202}"
GPS_END="${GPS_END:-1426897935}"
BASE_FILE_START="${BASE_FILE_START:-1426891200}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
FILE_START="$((BASE_FILE_START + 100 * TASK_ID))"
GWF_PATH="$GWF_DIR/V-raw-${FILE_START}-100.gwf"
OUTPUT_PATH="$SHARD_DIR/lsc_supp_1426891202_1426897935_${FILE_START}.pkl.gz"

echo "=== LSC supp shard task started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    task_id     = $TASK_ID"
echo "    file_start  = $FILE_START"
echo "    gwf_path    = $GWF_PATH"
echo "    output_path = $OUTPUT_PATH"

if [[ ! -f "$GWF_PATH" ]]; then
    echo "ERROR: GWF file not found: $GWF_PATH" >&2; exit 1
fi

mkdir -p "$SHARD_DIR"
cd "$WORKDIR"
"$IGWN_PYTHON" "$WORKDIR/scripts/extract_raw_interval_shard.py" \
    --gwf "$GWF_PATH" \
    --gps-start "$GPS_START" \
    --gps-end "$GPS_END" \
    --channel-set lsc-supp \
    --output "$OUTPUT_PATH"

echo "=== LSC supp shard task finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
