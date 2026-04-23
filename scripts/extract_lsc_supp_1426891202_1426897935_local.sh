#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
SHARD_DIR="${SHARD_DIR:-/scratch/sentenac/noisehound_raw/raw_interval_shards_1426891202_1426897935_lsc_supp}"
GWF_DIR="${GWF_DIR:-/scratch/sentenac/noisehound_raw/1426}"
LOG_DIR="${LOG_DIR:-/scratch/sentenac/noisehound_raw/lsc_supp_logs}"

mkdir -p "$SHARD_DIR" "$LOG_DIR"
cd "$WORKDIR"

echo "=== LSC supp extraction started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    shard_dir = $SHARD_DIR"
echo "    gwf_dir   = $GWF_DIR"

for i in $(seq 0 67); do
    FILE_START=$((1426891200 + 100 * i))
    GWF="$GWF_DIR/V-raw-${FILE_START}-100.gwf"
    OUT="$SHARD_DIR/lsc_supp_1426891202_1426897935_${FILE_START}.pkl.gz"
    LOG="$LOG_DIR/lsc_supp_${FILE_START}.log"
    if [[ -f "$GWF" ]]; then
        "$IGWN_PYTHON" scripts/extract_raw_interval_shard.py \
            --gwf "$GWF" \
            --gps-start 1426891202 \
            --gps-end 1426897935 \
            --channel-set lsc-supp \
            --output "$OUT" \
            > "$LOG" 2>&1 &
    else
        echo "SKIP: $GWF not found"
    fi
done

wait
echo "=== LSC supp extraction finished $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
