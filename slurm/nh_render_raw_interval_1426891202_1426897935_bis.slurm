#!/usr/bin/env bash
#SBATCH -J nh-raw-1426-bis
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-raw-1426-bis-%j.out

# Render the full raw interval 1426891202–1426897935 with the "bis" settings:
# - uniform-bandpass mode
# - 35–50 Hz band-pass
# - glitch markers from full_25min_glitches_ER16-O4b.csv
#
# Submit:
#   sbatch slurm/nh_render_raw_interval_1426891202_1426897935_bis.slurm
#
# Main outputs:
#   /scratch/sentenac/noisehound_raw/1426/render_1426891202_1426897935.log
#   /pbs/home/s/sentenac/NOISEHOUND/usecases/1426891202_1426897935/raw_interval_1426891202_1426897935_bp35_50.png

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
WRAPPER="${WRAPPER:-$HOME/run_envelope_1426891202_1426897935_bis.sh}"
SCRDIR="${SCRDIR:-/scratch/sentenac/noisehound_raw/1426}"
RENDER_LOG="${RENDER_LOG:-$SCRDIR/render_1426891202_1426897935.log}"

echo "=== raw interval render started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    job_id     = ${SLURM_JOB_ID:-local}"
echo "    workdir    = $WORKDIR"
echo "    wrapper    = $WRAPPER"
echo "    render_log = $RENDER_LOG"

cd "$WORKDIR"

if [[ ! -x "$WRAPPER" ]]; then
    echo "Wrapper missing or not executable: $WRAPPER" >&2
    exit 1
fi

mkdir -p "$SCRDIR"

bash "$WRAPPER"

echo "=== raw interval render finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
if [[ -f "$RENDER_LOG" ]]; then
    echo "=== tail of render log ==="
    tail -n 40 "$RENDER_LOG"
fi
