#!/usr/bin/env bash
# Run a partial or full fast-channel scan in sequential batches of BATCH_SIZE.
# Designed to be called directly or via a SLURM array job.
#
# Usage:
#   bash scripts/launch_full_scan.sh
#
# Environment overrides (all optional):
#   IGWN_PYTHON       path to Python with gwpy
#   GWF_DIR           directory of raw GWF files
#   GPS               glitch GPS time
#   WINDOW            ±seconds display window (default 15)
#   GUARD             filter guard seconds (default 30)
#   MIN_RATE          minimum channel sample rate Hz (default 10)
#   BATCH_SIZE        channels per scan batch (default 100)
#   SNR_THRESHOLD     minimum SNR to flag as candidate (default 5)
#   OUT_DIR           output directory
#   OFFSET_START      first channel offset to process (default 0)
#   OFFSET_END        stop before this offset; 0 = run to end (default 0)
#   AGGREGATE         run aggregation step after scan? 1=yes 0=no (default 0)
#                     Set to 1 only in the dedicated aggregation job.

set -euo pipefail

IGWN_PYTHON=${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}
WORKDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
GWF_DIR=${GWF_DIR:-/sps/virgo/USERS/sentenac/gwf_1426}
GPS=${GPS:-1426892892}
WINDOW=${WINDOW:-15}
GUARD=${GUARD:-30}
MIN_RATE=${MIN_RATE:-10}
BATCH_SIZE=${BATCH_SIZE:-100}
SNR_THRESHOLD=${SNR_THRESHOLD:-5}
OUT_DIR=${OUT_DIR:-/sps/virgo/USERS/sentenac/scan_full_${GPS}}
OFFSET_START=${OFFSET_START:-0}
OFFSET_END=${OFFSET_END:-0}       # 0 means run to end
AGGREGATE=${AGGREGATE:-0}

# Per-task log so parallel jobs don't interleave into one file
TASK_ID=${SLURM_ARRAY_TASK_ID:-solo}
LOG="$OUT_DIR/scan_task_${TASK_ID}.log"

mkdir -p "$OUT_DIR/batches" "$OUT_DIR/scan_plots"
cd "$WORKDIR"

tee_log() { echo "$*" | tee -a "$LOG"; }
tee_log "=== Scan task $TASK_ID started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
tee_log "    GPS=$GPS  window=+-${WINDOW}s  guard=${GUARD}s"
tee_log "    min_rate=${MIN_RATE}Hz  batch_size=$BATCH_SIZE  snr_threshold=$SNR_THRESHOLD"
tee_log "    gwf_dir=$GWF_DIR  out_dir=$OUT_DIR"
tee_log "    offset_start=$OFFSET_START  offset_end=${OFFSET_END:-end}"

OFFSET=$OFFSET_START
while true; do
    # Stop if we've reached OFFSET_END (when non-zero)
    if [ "$OFFSET_END" -gt 0 ] && [ "$OFFSET" -ge "$OFFSET_END" ]; then
        tee_log "[$(date -u '+%H:%M:%SZ')] Reached offset_end=$OFFSET_END — task done."
        break
    fi

    CSV="$OUT_DIR/batches/scan_offset$(printf '%05d' $OFFSET).csv"

    if [ -f "$CSV" ] && [ -s "$CSV" ] && [ "$(tail -n +2 "$CSV" | wc -l)" -gt 0 ]; then
        NROWS=$(tail -n +2 "$CSV" | wc -l)
        tee_log "[$(date -u '+%H:%M:%SZ')] Skip offset=$OFFSET (already done, $NROWS rows)"
    else
        tee_log "[$(date -u '+%H:%M:%SZ')] Scanning offset=$OFFSET ..."
        set +e
        $IGWN_PYTHON scripts/scan_fast_channels.py \
            --gwf-dir "$GWF_DIR" \
            --gps "$GPS" \
            --window "$WINDOW" \
            --guard "$GUARD" \
            --min-rate "$MIN_RATE" \
            --n-channels "$BATCH_SIZE" \
            --channel-offset "$OFFSET" \
            --snr-threshold "$SNR_THRESHOLD" \
            --plot-dir "$OUT_DIR/scan_plots" \
            --output "$CSV" 2>&1 | tee -a "$LOG"
        EXIT_CODE=${PIPESTATUS[0]}
        set -e

        if [ $EXIT_CODE -eq 0 ] && [ -f "$CSV" ] && [ -s "$CSV" ]; then
            NROWS=$(tail -n +2 "$CSV" | wc -l)
            tee_log "[$(date -u '+%H:%M:%SZ')] Offset=$OFFSET done: $NROWS rows."
        else
            tee_log "[$(date -u '+%H:%M:%SZ')] Offset=$OFFSET: scan exited (code=$EXIT_CODE) — end of channel list."
            break
        fi
    fi

    NROWS=$(tail -n +2 "$CSV" | wc -l)
    if [ "$NROWS" -lt "$BATCH_SIZE" ]; then
        tee_log "[$(date -u '+%H:%M:%SZ')] Last batch ($NROWS < $BATCH_SIZE). Task complete."
        break
    fi
    OFFSET=$((OFFSET + BATCH_SIZE))
done

tee_log "=== Scan task $TASK_ID done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

# ── Aggregation (only when AGGREGATE=1) ───────────────────────────────────────
if [ "$AGGREGATE" = "1" ]; then
    tee_log "[$(date -u '+%H:%M:%SZ')] Aggregating candidates (SNR >= $SNR_THRESHOLD) ..."

    OUT_DIR_AGG="$OUT_DIR" SNR_THRESH_AGG="$SNR_THRESHOLD" \
    $IGWN_PYTHON - << 'PYEOF'
import csv, glob, os, sys

out_dir   = os.environ["OUT_DIR_AGG"]
threshold = float(os.environ["SNR_THRESH_AGG"])

candidates = []
for f in sorted(glob.glob(os.path.join(out_dir, "batches", "scan_offset*.csv"))):
    for row in csv.DictReader(open(f)):
        if row.get("error"):
            continue
        try:
            snr = float(row["peak_snr"])
        except (ValueError, KeyError):
            continue
        if snr >= threshold:
            candidates.append(row)

candidates.sort(key=lambda r: float(r["peak_snr"]), reverse=True)

out_path = os.path.join(out_dir, "candidates.csv")
if not candidates:
    print(f"No candidates above SNR={threshold:.1f}", file=sys.stderr)
    sys.exit(0)

fieldnames = list(candidates[0].keys())
with open(out_path, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(candidates)

print(f"Wrote {len(candidates)} candidates → {out_path}")
print(f"{'Rank':<5} {'SNR':>7}  {'t_rel':>8}  {'Rate':>6}  Channel")
for k, r in enumerate(candidates[:30], 1):
    print(f"{k:<5} {float(r['peak_snr']):>7.2f}  {float(r['peak_time_rel']):>+8.2f}s"
          f"  {r['rate']:>5}Hz  {r['channel']}")
if len(candidates) > 30:
    print(f"  … and {len(candidates) - 30} more")
PYEOF

    tee_log "=== Aggregation done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
fi
