#!/usr/bin/env bash
set -euo pipefail

REMOTE_HOST="${1:-sentenac@cca020.in2p3.fr}"
REMOTE_DIR="${2:-/pbs/home/s/sentenac/NOISEHOUND}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ssh "$REMOTE_HOST" "mkdir -p '$REMOTE_DIR'"
rsync -avz \
  --exclude "__pycache__" \
  --exclude ".pytest_cache" \
  --exclude "*.pyc" \
  --exclude ".git" \
  "$ROOT_DIR"/ "$REMOTE_HOST":"$REMOTE_DIR"/

printf 'Synced %s to %s:%s\n' "$ROOT_DIR" "$REMOTE_HOST" "$REMOTE_DIR"
