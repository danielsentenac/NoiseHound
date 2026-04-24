#!/usr/bin/env bash
#SBATCH --job-name=nh_scan_4550
#SBATCH --array=0-6
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/sps/virgo/USERS/sentenac/scan_lal_1426894550/slurm_scan_%A_%a.out
#SBATCH --error=/sps/virgo/USERS/sentenac/scan_lal_1426894550/slurm_scan_%A_%a.err

set -euo pipefail
cd ~/NOISEHOUND

export GWF_DIR=/scratch/sentenac/noisehound_raw/1426
export GPS=1426894550
export OUT_DIR=/sps/virgo/USERS/sentenac/scan_lal_${GPS}
export BATCH_SIZE=100
export OFFSET_START=$(( SLURM_ARRAY_TASK_ID * 1000 ))
export OFFSET_END=$(( (SLURM_ARRAY_TASK_ID + 1) * 1000 ))
export AGGREGATE=0

mkdir -p "$OUT_DIR/batches" "$OUT_DIR/scan_plots"

bash scripts/launch_full_scan.sh
