#!/usr/bin/env bash
#SBATCH --job-name=nh_watchlist
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/sps/virgo/USERS/sentenac/scan_lal_%j.out
#SBATCH --error=/sps/virgo/USERS/sentenac/scan_lal_%j.err

# Usage: sbatch --export=GPS=<gps>,GWF_DIR=<dir>,OUT=<file.png> submit_watchlist_plot.slurm

set -euo pipefail
cd ~/NOISEHOUND

: "${GPS:?set GPS}"
: "${GWF_DIR:?set GWF_DIR}"
: "${OUT:=watchlist_${GPS}.png}"

python scripts/plot_top10_candidates.py \
    --candidates data/watchlist_1426892892.csv \
    --gwf-dir    "$GWF_DIR" \
    --gps        "$GPS" \
    --window     15 \
    --guard      30 \
    --n          23 \
    --output     "$OUT"
