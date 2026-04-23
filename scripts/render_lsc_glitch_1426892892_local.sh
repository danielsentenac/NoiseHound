#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935}"
SUPP_DIR="${SUPP_DIR:-$WORKDIR/outputs/raw_interval_shards_1426891202_1426897935_lsc_supp}"
OUTDIR="$WORKDIR/usecases/1426891202_1426897935"
LOG="$WORKDIR/outputs/render_lsc_glitch_1426892892.log"

mkdir -p "$OUTDIR" "$(dirname "$LOG")"
cd "$WORKDIR"

echo "=== render lsc glitch 1426892892 started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

"$IGWN_PYTHON" "$WORKDIR/scripts/render_raw_glitch_response_from_shards.py" \
    --glitch-gps 1426892892.417968 \
    --window 60 \
    --panel-mode all-envelope \
    --aux-bandpass-low 37 \
    --aux-bandpass-high 44 \
    --aux-envelope-smooth-s 1 \
    --glitch-times "$WORKDIR/data/full_25min_glitches_ER16-O4b.csv" \
    --q-fmin 10 --q-fmax 500 --q-qmin 4 --q-qmax 64 \
    --q-max-rate 4096 --max-plot-points 8000 \
    --shard \
        "$SHARD_DIR/raw_interval_1426891202_1426897935_1426892800.pkl.gz" \
        "$SHARD_DIR/raw_interval_1426891202_1426897935_1426892900.pkl.gz" \
        "$SUPP_DIR/lsc_supp_1426891202_1426897935_1426892800.pkl.gz" \
        "$SUPP_DIR/lsc_supp_1426891202_1426897935_1426892900.pkl.gz" \
    --out "$OUTDIR/raw_glitch_1426892892_lsc_allenv_bp37_44_60s.png" \
    2>&1 | tee "$LOG"

echo "=== render lsc glitch 1426892892 finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
