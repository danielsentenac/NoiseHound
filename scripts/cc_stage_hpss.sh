#!/usr/bin/env bash
set -euo pipefail

SOURCE_LIST="${SOURCE_LIST:?set SOURCE_LIST to a file containing one HPSS file per line}"
DEST_DIR="${DEST_DIR:-}"
OUTPUT_LIST="${OUTPUT_LIST:-}"
HPSS_HOST="${HPSS_HOST:-cchpss0}"
HPSS_BASE="${HPSS_BASE:-/hpss/in2p3.fr/group/virgo}"
FORCE="${FORCE:-0}"
FALLBACK_BASE="${FALLBACK_BASE:-$HOME/NOISEHOUND/staging}"

resolve_dest_dir() {
  if [[ -n "$DEST_DIR" ]]; then
    mkdir -p "$DEST_DIR"
    printf '%s\n' "$DEST_DIR"
    return
  fi

  local candidates=()
  [[ -n "${SLURM_TMPDIR:-}" ]] && candidates+=("${SLURM_TMPDIR%/}/noisehound_frames")
  [[ -n "${TMPDIR:-}" ]] && candidates+=("${TMPDIR%/}/noisehound_frames")
  candidates+=("/tmp/$USER/noisehound_frames")
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    candidates+=("${FALLBACK_BASE%/}/${SLURM_JOB_ID}")
  fi
  candidates+=("${FALLBACK_BASE%/}/manual")

  local candidate
  for candidate in "${candidates[@]}"; do
    if mkdir -p "$candidate" 2>/dev/null; then
      printf '%s\n' "$candidate"
      return
    fi
  done

  echo "unable to create a writable staging directory" >&2
  echo "tried: ${candidates[*]}" >&2
  exit 1
}

DEST_DIR="$(resolve_dest_dir)"
if [[ -z "${OUTPUT_LIST:-}" ]]; then
  OUTPUT_LIST="$DEST_DIR/local_paths.txt"
fi

mkdir -p "$(dirname "$OUTPUT_LIST")"
: > "$OUTPUT_LIST"

printf 'resolved dest dir: %s\n' "$DEST_DIR"
printf 'output list: %s\n' "$OUTPUT_LIST"
printf 'slurm tmpdir: %s\n' "${SLURM_TMPDIR:-}"
printf 'tmpdir: %s\n' "${TMPDIR:-}"

normalize_remote_path() {
  local entry="$1"
  if [[ "$entry" == *:* ]]; then
    printf '%s\n' "${entry#*:}"
  elif [[ "$entry" == /* ]]; then
    printf '%s\n' "$entry"
  else
    printf '%s/%s\n' "${HPSS_BASE%/}" "${entry#/}"
  fi
}

build_remote_spec() {
  local entry="$1"
  local path
  path="$(normalize_remote_path "$entry")"
  printf '%s:%s\n' "$HPSS_HOST" "$path"
}

build_local_path() {
  local remote_path="$1"
  local relative_path
  if [[ "$remote_path" == "${HPSS_BASE%/}/"* ]]; then
    relative_path="${remote_path#${HPSS_BASE%/}/}"
  else
    relative_path="${remote_path#/}"
  fi
  printf '%s/%s\n' "${DEST_DIR%/}" "$relative_path"
}

while IFS= read -r entry || [[ -n "$entry" ]]; do
  [[ -z "$entry" ]] && continue
  [[ "$entry" =~ ^[[:space:]]*# ]] && continue

  remote_path="$(normalize_remote_path "$entry")"
  remote_spec="$(build_remote_spec "$entry")"
  local_path="$(build_local_path "$remote_path")"

  mkdir -p "$(dirname "$local_path")"

  if [[ -f "$local_path" && "$FORCE" != "1" ]]; then
    printf 'reuse %s\n' "$local_path"
  else
    printf 'rfcp %s -> %s\n' "$remote_spec" "$local_path"
    rfcp "$remote_spec" "$local_path"
  fi

  printf '%s\n' "$local_path" >> "$OUTPUT_LIST"
done < "$SOURCE_LIST"

printf 'Wrote %s\n' "$OUTPUT_LIST"
