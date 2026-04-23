#!/usr/bin/env bash
#SBATCH -J nh-zoom-ni-temps
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 05:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-ni-temps-%j.out

# 3-row periodogram zoom for NI temperature channels
# Pre-baffle LN3 / Post-baffle LN3 / Post-baffle LN3_ALIGNED
# Channels: INF_NI_BOTTOM_TE1 and INF_TCS_NI_RH_TE

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_ni_temps_$$"
mkdir -p "$STAGE_DIR"

echo "=== NI temps zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

run_zoom() {
    local channel="$1"
    local ylabel="$2"
    local outfile="$3"
    echo "--- Processing $channel ---"
    "$IGWN_PYTHON" "$WORKDIR/scripts/zoom_asc_sr_ty.py" \
        --channel          "$channel" \
        --ylabel           "$ylabel" \
        --output           "$WORKDIR/outputs/$outfile" \
        --stage-dir        "$STAGE_DIR" \
        --after-gps-start  1452456000 \
        --after-gps-end    1458000000 \
    && echo "  OK -> outputs/$outfile" \
    || echo "  FAILED: $channel"
}

run_zoom "V1:INF_NI_BOTTOM_TE1"  "NI tower bottom TE1 [°C]"      "zoom_ni_bottom_te1.png"
run_zoom "V1:INF_TCS_NI_RH_TE"   "NI ring heater thermistor [°C]" "zoom_ni_rh_te.png"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
