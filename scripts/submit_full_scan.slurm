#!/usr/bin/env bash
#SBATCH --job-name=noisehound_scan
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/sps/virgo/USERS/sentenac/scan_full_1426892892/slurm_%j.out
#SBATCH --error=/sps/virgo/USERS/sentenac/scan_full_1426892892/slurm_%j.err

set -euo pipefail
cd ~/NOISEHOUND
export OUT_DIR=/sps/virgo/USERS/sentenac/scan_full_1426892892
export GWF_DIR=/sps/virgo/USERS/sentenac/gwf_1426
bash scripts/launch_full_scan.sh
