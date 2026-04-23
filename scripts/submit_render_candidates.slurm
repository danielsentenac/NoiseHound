#!/usr/bin/env bash
#SBATCH --job-name=noisehound_render
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/sps/virgo/USERS/sentenac/scan_full_1426892892/slurm_render_%j.out
#SBATCH --error=/sps/virgo/USERS/sentenac/scan_full_1426892892/slurm_render_%j.err

set -euo pipefail
/cvmfs/software.igwn.org/conda/envs/igwn/bin/python \
    ~/NOISEHOUND/scripts/render_candidates_from_gwf.py \
    --candidates /sps/virgo/USERS/sentenac/scan_full_1426892892/candidates.csv \
    --gwf-dir    /sps/virgo/USERS/sentenac/gwf_1426 \
    --gps        1426892892 \
    --snr-threshold 5 \
    --batch-size 10 \
    --plot-dir   /sps/virgo/USERS/sentenac/scan_full_1426892892/render_plots
