#!/usr/bin/env bash
set -euo pipefail

WORKDIR="${WORKDIR:-$HOME/NOISEHOUND}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.venvs/noisehound/bin/python}"
SERVER="${SERVER:-${GWDATAFIND_SERVER:-}}"
OBSERVATORY="${OBSERVATORY:-V}"
TOKEN_FILE="${TOKEN_FILE:-${BEARER_TOKEN_FILE:-/run/user/$(id -u)/bt_u$(id -u)}}"
URL_TYPE="${URL_TYPE:-osdf}"
HREC_FRAME_TYPE="${HREC_FRAME_TYPE:-HoftOnline}"
AUX_FRAME_TYPE="${AUX_FRAME_TYPE:-raw}"
GPS_START="${GPS_START:?set GPS_START}"
GPS_END="${GPS_END:?set GPS_END}"
OUTPUT_DIR="${OUTPUT_DIR:-$WORKDIR/frame_lists}"

mkdir -p "$OUTPUT_DIR"

"$PYTHON_BIN" -m noisehound.cli urls \
  --observatory "$OBSERVATORY" \
  --frame-type "$HREC_FRAME_TYPE" \
  --gps-start "$GPS_START" \
  --gps-end "$GPS_END" \
  --url-type "$URL_TYPE" \
  --token-file "$TOKEN_FILE" \
  ${SERVER:+--server "$SERVER"} \
  > "$OUTPUT_DIR/hrec_${HREC_FRAME_TYPE}_${GPS_START}_${GPS_END}_${URL_TYPE}.txt"

"$PYTHON_BIN" -m noisehound.cli urls \
  --observatory "$OBSERVATORY" \
  --frame-type "$AUX_FRAME_TYPE" \
  --gps-start "$GPS_START" \
  --gps-end "$GPS_END" \
  --url-type "$URL_TYPE" \
  --token-file "$TOKEN_FILE" \
  ${SERVER:+--server "$SERVER"} \
  > "$OUTPUT_DIR/aux_${AUX_FRAME_TYPE}_${GPS_START}_${GPS_END}_${URL_TYPE}.txt"

printf 'Wrote %s\n' "$OUTPUT_DIR/hrec_${HREC_FRAME_TYPE}_${GPS_START}_${GPS_END}_${URL_TYPE}.txt"
printf 'Wrote %s\n' "$OUTPUT_DIR/aux_${AUX_FRAME_TYPE}_${GPS_START}_${GPS_END}_${URL_TYPE}.txt"
