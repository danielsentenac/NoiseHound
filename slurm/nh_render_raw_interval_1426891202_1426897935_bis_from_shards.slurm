#!/usr/bin/env bash
#SBATCH -J nh-raw-1426-render
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-raw-1426-render-%j.out

# Final render step for the interval 1426891202–1426897935 using per-file shard bundles.
# This keeps the per-file raw extraction distributed across the array tasks.
# The final renderer rebuilds the q-transform directly from the Hrec shard data,
# so it does not need the raw GWF files again.
#
# Submit after successful shard extraction:
#   sbatch slurm/nh_render_raw_interval_1426891202_1426897935_bis_from_shards.slurm
#
# Main outputs:
#   /scratch/sentenac/noisehound_raw/1426/render_1426891202_1426897935_from_shards.log
#   /pbs/home/s/sentenac/NOISEHOUND/usecases/1426891202_1426897935/raw_interval_1426891202_1426897935_bp35_50_from_shards.png

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935}"
OUTDIR="${OUTDIR:-$WORKDIR/usecases/1426891202_1426897935}"
LOG="${LOG:-$WORKDIR/outputs/render_1426891202_1426897935_from_shards.log}"
OUTFILE="${OUTFILE:-$OUTDIR/raw_interval_1426891202_1426897935_allenv_bp37_44_from_shards.png}"
GLITCH_CSV="${GLITCH_CSV:-$WORKDIR/data/full_25min_glitches_ER16-O4b.csv}"

mkdir -p "$(dirname "$LOG")"
mkdir -p "$OUTDIR"

shopt -s nullglob
SHARD_FILES=( "$SHARD_DIR"/raw_interval_1426891202_1426897935_*.pkl.gz )
if [[ ${#SHARD_FILES[@]} -ne 68 ]]; then
    echo "Expected 68 shard files in $SHARD_DIR, found ${#SHARD_FILES[@]}" >&2
    exit 1
fi

{
    echo "=== shard render started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
    echo "    shard_dir   = $SHARD_DIR"
    echo "    render_log  = $LOG"
    echo "    outdir      = $OUTDIR"
    echo "    outfile     = $OUTFILE"
    echo "    shard_count = ${#SHARD_FILES[@]}"

    cd "$WORKDIR"
    "$IGWN_PYTHON" "$WORKDIR/scripts/render_raw_glitch_response_from_shards.py" \
        --gps-start 1426891202 \
        --gps-end 1426897935 \
        --panel-mode all-envelope \
        --aux-bandpass-low 37 \
        --aux-bandpass-high 44 \
        --aux-envelope-smooth-s 1 \
        --glitch-times "$GLITCH_CSV" \
        --q-fmin 10 \
        --q-fmax 500 \
        --q-qmin 4 \
        --q-qmax 64 \
        --q-max-rate 4096 \
        --max-plot-points 8000 \
        --shard "${SHARD_FILES[@]}" \
        --out "$OUTFILE"

    echo "=== shard render finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
} 2>&1 | tee "$LOG"
