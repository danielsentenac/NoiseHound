#!/usr/bin/env bash
#SBATCH -J nh-sr-channels
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 10:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-sr-channels-%j.out

# Slow-channel rate correlation for 6 new SR channels (baffle hypothesis)
# GPS 1443657600 (2025-10-01) -> 1459468800 (2026-04-01)
# Channels: INF_SR_MIR_COIL_UL_TE, INF_TCS_SR_RH_TE, TCS_SR_RH_SET,
#           INF_TCS_SR_CHROCC_TE_IntHeater, ASC_SR_TY_ERR_mean, ASC_SR_TX_ERR_mean

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
HPSS_BASE="cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
STAGE_DIR="/tmp/nh_sr_channels_$$"
mkdir -p "$STAGE_DIR"

echo "=== SR channels job started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" "$WORKDIR/scripts/rate_correlation_direct.py" \
    --triggers   "$WORKDIR/data/full_25min_glitches_ER16-O4b.csv" \
    --gps-start  1443657600 \
    --gps-end    1459468800 \
    --output     "$WORKDIR/outputs/rate_correlation/partial_sr.csv" \
    --hpss-base  "$HPSS_BASE" \
    --stage-dir  "$STAGE_DIR" \
    --bin-size-s 3600 \
    --channels   V1:INF_SR_MIR_COIL_UL_TE \
                 V1:INF_TCS_SR_RH_TE \
                 V1:TCS_SR_RH_SET \
                 V1:INF_TCS_SR_CHROCC_TE_IntHeater \
                 V1:ASC_SR_TY_ERR_mean \
                 V1:ASC_SR_TX_ERR_mean

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
