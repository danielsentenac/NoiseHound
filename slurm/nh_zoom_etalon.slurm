#!/usr/bin/env bash
#SBATCH -J nh-zoom-etalon
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 05:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-etalon-%j.out

# 3-row periodogram zoom for Etalon NI ring heater channels
# Testing whether LSC_Etalon_NI_RH_SET oscillates at 25 min —
# candidate root driver of NI tower thermal oscillation
#
# Pre-baffle LN3 / Post-baffle LN3 / Post-baffle LN3_ALIGNED

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_etalon_$$"
mkdir -p "$STAGE_DIR"

echo "=== Etalon NI zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

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

# Etalon NI ring heater — setpoint, output, input
run_zoom "V1:LSC_Etalon_NI_RH_SET_mean" "Etalon NI RH setpoint [arb]"  "zoom_etalon_ni_rh_set.png"
run_zoom "V1:LSC_Etalon_NI_RH_OUT_mean" "Etalon NI RH output [arb]"    "zoom_etalon_ni_rh_out.png"
run_zoom "V1:LSC_Etalon_NI_RH_IN_mean"  "Etalon NI RH input [arb]"     "zoom_etalon_ni_rh_in.png"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
