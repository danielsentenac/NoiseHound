#!/usr/bin/env bash
#SBATCH -J nh-lsc-render-1426
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-lsc-render-1426-%j.out

# Render glitch 1426892892 with LSC longitudinal channels only (all-envelope mode).
# Merges original shards + LSC supplemental shards.
#
# Submit after lsc_supp extraction completes:
#   sbatch slurm/nh_render_lsc_glitch_1426892892.slurm
#
# Output:
#   usecases/1426891202_1426897935/raw_glitch_1426892892_lsc_allenv_bp37_44_60s.png

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935}"
SUPP_DIR="${SUPP_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935_lsc_supp}"
OUTDIR="${OUTDIR:-$WORKDIR/usecases/1426891202_1426897935}"
GLITCH_CSV="${GLITCH_CSV:-$WORKDIR/data/full_25min_glitches_ER16-O4b.csv}"
OUTFILE="${OUTFILE:-$OUTDIR/raw_glitch_1426892892_lsc_allenv_bp37_44_60s.png}"
LOG="${LOG:-$WORKDIR/outputs/render_lsc_glitch_1426892892.log}"

shopt -s nullglob
SHARD_FILES=( "$SHARD_DIR"/raw_interval_1426891202_1426897935_*.pkl.gz )
SUPP_FILES=( "$SUPP_DIR"/lsc_supp_1426891202_1426897935_*.pkl.gz )

if [[ ${#SHARD_FILES[@]} -eq 0 ]]; then
    echo "No original shard files found in $SHARD_DIR" >&2; exit 1
fi
if [[ ${#SUPP_FILES[@]} -eq 0 ]]; then
    echo "No LSC supplemental shard files found in $SUPP_DIR" >&2; exit 1
fi

echo "Original shards: ${#SHARD_FILES[@]}, LSC supp shards: ${#SUPP_FILES[@]}"
mkdir -p "$OUTDIR" "$(dirname "$LOG")"
cd "$WORKDIR"

{
    echo "=== started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
    "$IGWN_PYTHON" "$WORKDIR/scripts/render_raw_glitch_response_from_shards.py" \
        --glitch-gps 1426892892.417968 \
        --window 60 \
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
        --shard "${SHARD_FILES[@]}" "${SUPP_FILES[@]}" \
        --out "$OUTFILE"
    echo "=== finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
} 2>&1 | tee "$LOG"
