#!/usr/bin/env bash

set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
ARRAY_SPEC="${ARRAY_SPEC:-0-67}"
EXTRACT_JOB="$WORKDIR/slurm/nh_extract_raw_interval_1426891202_1426897935_array.slurm"
RENDER_JOB="$WORKDIR/slurm/nh_render_raw_interval_1426891202_1426897935_bis_from_shards.slurm"

cd "$WORKDIR"

extract_id="$(sbatch --parsable --array="$ARRAY_SPEC" "$EXTRACT_JOB")"
render_id="$(sbatch --parsable --dependency="afterok:${extract_id}" "$RENDER_JOB")"

printf 'Submitted shard array job: %s\n' "$extract_id"
printf 'Submitted dependent render job: %s\n' "$render_id"
printf 'Shard logs: %s\n' "$WORKDIR/slurm-nh-raw-1426-shard-<arrayjob>_<task>.out"
printf 'Render log: %s\n' "$WORKDIR/outputs/render_1426891202_1426897935_from_shards.log"
printf 'Final PNG: %s\n' "$WORKDIR/usecases/1426891202_1426897935/raw_interval_1426891202_1426897935_bp35_50_from_shards.png"
