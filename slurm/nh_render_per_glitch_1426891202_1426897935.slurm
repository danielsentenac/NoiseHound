#!/usr/bin/env bash
#SBATCH -J nh-raw-1426-per-glitch
#SBATCH -p htc
#SBATCH -c 2
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-raw-1426-per-glitch-%j.out

# Per-glitch renders for interval 1426891202–1426897935.
# For each of the 5 glitches: 60 s window, 35-50 Hz bandpass, 1 s smoothing.
# Reuses the shard files already extracted by nh_extract_raw_interval_1426891202_1426897935_array.slurm.
#
# Submit:
#   sbatch slurm/nh_render_per_glitch_1426891202_1426897935.slurm
#
# Outputs:
#   /pbs/home/s/sentenac/NOISEHOUND/usecases/1426891202_1426897935/raw_glitch_<GPS>_bp35_50_60s_smooth1s.png

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935}"
OUTDIR="${OUTDIR:-$WORKDIR/usecases/1426891202_1426897935}"
GLITCH_CSV="${GLITCH_CSV:-$WORKDIR/data/full_25min_glitches_ER16-O4b.csv}"

shopt -s nullglob
SHARD_FILES=( "$SHARD_DIR"/raw_interval_1426891202_1426897935_*.pkl.gz )
if [[ ${#SHARD_FILES[@]} -ne 68 ]]; then
    echo "Expected 68 shard files in $SHARD_DIR, found ${#SHARD_FILES[@]}" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
cd "$WORKDIR"

GLITCH_TIMES=(
    1426892892.417968
)

for GPS in "${GLITCH_TIMES[@]}"; do
    GPS_INT=$(printf '%.0f' "$GPS")
    OUTFILE="$OUTDIR/raw_glitch_${GPS_INT}_allenv_bp37_44_60s.png"
    LOG="$WORKDIR/outputs/render_glitch_${GPS_INT}_allenv_from_shards.log"
    echo "=== rendering GPS ${GPS} -> $OUTFILE ==="
    {
        echo "=== started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
        "$IGWN_PYTHON" "$WORKDIR/scripts/render_raw_glitch_response_from_shards.py" \
            --glitch-gps "$GPS" \
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
            --shard "${SHARD_FILES[@]}" \
            --out "$OUTFILE"
        echo "=== finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
    } 2>&1 | tee "$LOG"
done

echo "All per-glitch renders done."
