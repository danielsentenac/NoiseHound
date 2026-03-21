#!/usr/bin/env bash
#SBATCH -J nh-zoom-mirrors
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-zoom-mirrors-%j.out

# Periodogram zoom for all ITF mirror ASC angular channels
# Checks whether the ~25-min oscillation seen in ASC_SR_TY_ERR is SR-specific
# or a global Alignment/Control response to a common temperature perturbation.
#
# Mirrors: NI, WI, NE, WE, PR, BS (SR already done separately)
# DOFs:    TX, TY  (no TZ channels found in trend files)
# Channel strategy:
#   NI/WI: _LF_B7/B8_I_10Hz_mean  (ERR/CORR are zero — use LF beam sensing)
#   NE/WE: _CORR_mean              (no ERR channel in trend; CORR is active)
#   PR/BS: _ERR_mean               (disturbance signal; CORR=ERR numerically but ERR is the drive)
#
# Two windows:
#   Before shutdown: GPS 1448323200-1448928000 (24 Nov – 01 Dec 2025)
#   After baffle:    GPS 1452456000-1456790400 (14 Jan – 01 Mar 2026), locked only
#     Best locked blocks: Jan 14-18 (36+57h), Jan 21-25 (24+46h), Feb 20-23 (63h)
#     Lock filter uses META_ITF_LOCK_index_mean > 10

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
STAGE_DIR="/tmp/nh_zoom_mirrors_$$"
mkdir -p "$STAGE_DIR"

echo "=== All-mirrors ASC zoom started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

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

# NI (North Input) — LF beam sensing (CORR is zero, no ERR in trend)
run_zoom "V1:ASC_NI_TX_LF_B7_I_10Hz_mean" "NI TX LF B7 I [rad]"  "zoom_asc_ni_tx.png"
run_zoom "V1:ASC_NI_TY_LF_B7_I_10Hz_mean" "NI TY LF B7 I [rad]"  "zoom_asc_ni_ty.png"

# WI (West Input) — LF beam sensing (CORR is zero, no ERR in trend)
run_zoom "V1:ASC_WI_TX_LF_B8_I_10Hz_mean" "WI TX LF B8 I [rad]"  "zoom_asc_wi_tx.png"
run_zoom "V1:ASC_WI_TY_LF_B8_I_10Hz_mean" "WI TY LF B8 I [rad]"  "zoom_asc_wi_ty.png"

# NE (North End) — CORR (no ERR in trend; CORR is active)
run_zoom "V1:ASC_NE_TX_CORR_mean" "NE TX correction [rad]"  "zoom_asc_ne_tx.png"
run_zoom "V1:ASC_NE_TY_CORR_mean" "NE TY correction [rad]"  "zoom_asc_ne_ty.png"

# WE (West End) — CORR (no ERR in trend; CORR is active)
run_zoom "V1:ASC_WE_TX_CORR_mean" "WE TX correction [rad]"  "zoom_asc_we_tx.png"
run_zoom "V1:ASC_WE_TY_CORR_mean" "WE TY correction [rad]"  "zoom_asc_we_ty.png"

# PR (Power Recycling) — ERR (disturbance signal)
run_zoom "V1:ASC_PR_TX_ERR_mean" "PR TX error [rad]"  "zoom_asc_pr_tx.png"
run_zoom "V1:ASC_PR_TY_ERR_mean" "PR TY error [rad]"  "zoom_asc_pr_ty.png"

# BS (Beam Splitter) — ERR (disturbance signal)
run_zoom "V1:ASC_BS_TX_ERR_mean" "BS TX error [rad]"  "zoom_asc_bs_tx.png"
run_zoom "V1:ASC_BS_TY_ERR_mean" "BS TY error [rad]"  "zoom_asc_bs_ty.png"

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
