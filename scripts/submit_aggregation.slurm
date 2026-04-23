#!/usr/bin/env bash
#SBATCH --job-name=noisehound_agg
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/sps/virgo/USERS/sentenac/scan_lal_1426892892/slurm_agg_%j.out
#SBATCH --error=/sps/virgo/USERS/sentenac/scan_lal_1426892892/slurm_agg_%j.err

# Run aggregation only — must be submitted with --dependency=afterok:<array_jobid>

set -euo pipefail
cd ~/NOISEHOUND

export GWF_DIR=/sps/virgo/USERS/sentenac/gwf_1426
export GPS=1426892892
export OUT_DIR=/sps/virgo/USERS/sentenac/scan_lal_${GPS}
export AGGREGATE=1
export SNR_THRESHOLD=5

# OFFSET_END=1 forces the scan loop to exit immediately (offset 0 >= 1),
# so only the aggregation block runs.
export OFFSET_START=0
export OFFSET_END=1

bash scripts/launch_full_scan.sh
