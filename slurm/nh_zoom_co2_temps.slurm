#!/usr/bin/env bash
#SBATCH -J nh-zoom-co2-temps
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 05:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-co2-temps-%j.out

# 3-row periodogram zoom for CO2 bench ambient temperature channels
# These are the strongest correlators with glitch rate (|r|~0.14-0.18, lag~-1h)
# Testing whether they oscillate at 25 min pre/post baffle
#
# Pre-baffle LN3 / Post-baffle LN3 / Post-baffle LN3_ALIGNED

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_co2_$$"
mkdir -p "$STAGE_DIR"

echo "=== CO2 bench temps zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

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

# CO2 bench ambient temperatures (strongest correlators with glitch rate)
run_zoom "V1:ENV_TCS_CO2_NI_TE"   "CO2 NI bench ambient TE [°C]"  "zoom_co2_ni_bench_te.png"
run_zoom "V1:TCS_NI_TE_CO2Laser"  "CO2 NI laser body TE [°C]"     "zoom_co2_ni_laser_te.png"
run_zoom "V1:ENV_TCS_CO2_WI_TE"   "CO2 WI bench ambient TE [°C]"  "zoom_co2_wi_bench_te.png"
run_zoom "V1:TCS_WI_TE_CO2Laser"  "CO2 WI laser body TE [°C]"     "zoom_co2_wi_laser_te.png"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
