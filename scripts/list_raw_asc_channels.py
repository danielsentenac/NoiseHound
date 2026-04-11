#!/usr/bin/env python3
"""
List raw ASC ERR/CORR channels for the main ITF towers from one Virgo raw frame.

The script stages a single raw GWF file on `cca.in2p3.fr` with `rfcp`, lists
its channels with `FrChannels`, then filters for ASC alignment channels ending
in `ERR` or `CORR`.

Examples
--------
List channels for a known pre-Christmas 2025 glitch frame:
    python scripts/list_raw_asc_channels.py --gps 1446106739 --run O4c

Try O4c then O4b automatically:
    python scripts/list_raw_asc_channels.py --gps 1446106739 --run auto

Inspect an older O4b frame and save the full raw channel inventory:
    python scripts/list_raw_asc_channels.py \
        --raw-start 1415578700 \
        --run O4b \
        --save-all-channels
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import subprocess
import sys
from pathlib import Path


WORKDIR = Path(__file__).resolve().parents[1]
OUTPUTS_DIR = WORKDIR / "outputs"
DEFAULT_HOST = "cca.in2p3.fr"
DEFAULT_HPSS_ROOT = "cchpss0:/hpss/in2p3.fr/group/virgo/Run"
DEFAULT_TOWERS = ["NI", "WI", "NE", "WE", "BS", "PR", "SR"]
DEFAULT_TYPES = ["ERR", "CORR"]
ALL_ASC_KIND_RE = re.compile(
    r"^(V1:ASC_(?P<scope>[A-Za-z0-9]+)_(?P<rest>.+)_(?P<kind>ERR|CORR))$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="List raw ASC ERR/CORR channels for one Virgo raw frame."
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument(
        "--gps",
        type=float,
        help="Glitch/event GPS time. The raw frame start is floor(gps/100)*100.",
    )
    src.add_argument(
        "--raw-start",
        type=int,
        help="Exact raw frame start GPS for V-raw-<start>-100.gwf.",
    )
    src.add_argument(
        "--remote-path",
        help="Full HPSS path to one raw frame, e.g. "
        "cchpss0:/.../Run/O4c/raw/1446/V-raw-1446106700-100.gwf",
    )
    parser.add_argument(
        "--run",
        default="auto",
        choices=["auto", "O4b", "O4c"],
        help="Run to search under when building the HPSS path.",
    )
    parser.add_argument(
        "--host",
        default=DEFAULT_HOST,
        help=f"SSH host that provides rfdir/rfcp/FrChannels (default: {DEFAULT_HOST}).",
    )
    parser.add_argument(
        "--hpss-root",
        default=DEFAULT_HPSS_ROOT,
        help=f"HPSS Run root (default: {DEFAULT_HPSS_ROOT}).",
    )
    parser.add_argument(
        "--towers",
        nargs="+",
        default=DEFAULT_TOWERS,
        help="Tower list to keep, default: NI WI NE WE BS PR SR.",
    )
    parser.add_argument(
        "--types",
        nargs="+",
        default=DEFAULT_TYPES,
        help="Signal suffixes to keep, default: ERR CORR.",
    )
    parser.add_argument(
        "--save-all-channels",
        action="store_true",
        help="Also save the full raw channel list from FrChannels.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve the HPSS file path and exit without staging or listing channels.",
    )
    return parser.parse_args()


def run_command(
    cmd: list[str],
    *,
    check: bool = True,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, input=input_text, capture_output=True, text=True)
    if check and result.returncode != 0:
        stderr = result.stderr.strip()
        stdout = result.stdout.strip()
        detail = stderr or stdout or f"command failed with exit code {result.returncode}"
        raise RuntimeError(f"{' '.join(shlex.quote(x) for x in cmd)}\n{detail}")
    return result


def ssh_run(host: str, remote_command: str, *, check: bool = True) -> subprocess.CompletedProcess[str]:
    cmd = [
        "ssh",
        "-o",
        "BatchMode=yes",
        "-o",
        "ConnectTimeout=15",
        host,
        "bash",
        "-s",
        "--",
    ]
    return run_command(cmd, check=check, input_text=remote_command)


def frame_start_from_gps(gps: float) -> int:
    return int(math.floor(gps / 100.0) * 100)


def frame_name(frame_start: int) -> str:
    return f"V-raw-{frame_start}-100.gwf"


def hpss_path(run: str, frame_start: int, hpss_root: str) -> str:
    subdir = str(frame_start)[:4]
    return f"{hpss_root}/{run}/raw/{subdir}/{frame_name(frame_start)}"


def infer_frame_start_from_remote_path(remote_path: str) -> int:
    match = re.search(r"V-raw-(\d+)-100\.gwf$", remote_path)
    if not match:
        raise ValueError(f"cannot infer frame start from remote path: {remote_path}")
    return int(match.group(1))


def infer_run_from_remote_path(remote_path: str) -> str | None:
    match = re.search(r"/Run/(O4[bc])/raw/", remote_path)
    if not match:
        return None
    return match.group(1)


def remote_exists(host: str, remote_path: str) -> bool:
    quoted = shlex.quote(remote_path)
    probe = ssh_run(host, f"rfdir {quoted} >/dev/null 2>&1", check=False)
    return probe.returncode == 0


def resolve_remote_path(args: argparse.Namespace) -> tuple[str, str, int]:
    if args.remote_path:
        remote_path = args.remote_path
        run_name = infer_run_from_remote_path(remote_path) or args.run
        frame_start = infer_frame_start_from_remote_path(remote_path)
        if not remote_exists(args.host, remote_path):
            raise RuntimeError(f"remote file not found: {remote_path}")
        return run_name, remote_path, frame_start

    frame_start = args.raw_start if args.raw_start is not None else frame_start_from_gps(args.gps)
    runs = [args.run] if args.run != "auto" else ["O4c", "O4b"]
    tried: list[str] = []
    for run_name in runs:
        remote_path = hpss_path(run_name, frame_start, args.hpss_root)
        tried.append(remote_path)
        if remote_exists(args.host, remote_path):
            return run_name, remote_path, frame_start
    tried_msg = "\n  ".join(tried)
    raise RuntimeError(f"could not find raw frame on {args.host}:\n  {tried_msg}")


def list_remote_channels(host: str, remote_path: str) -> list[str]:
    quoted_path = shlex.quote(remote_path)
    local_name = frame_name(infer_frame_start_from_remote_path(remote_path))
    remote_script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_raw_asc_channels.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT

local_gwf="$stage_dir/{local_name}"
rfcp {quoted_path} "$local_gwf" >/dev/null

frchannels_bin="${{FRCHANNELS_BIN:-}}"
if [[ -z "$frchannels_bin" ]]; then
  if command -v FrChannels >/dev/null 2>&1; then
    frchannels_bin="$(command -v FrChannels)"
  else
    frchannels_bin="/cvmfs/software.igwn.org/conda/envs/igwn/bin/FrChannels"
  fi
fi

"$frchannels_bin" "$local_gwf"
"""
    result = ssh_run(host, remote_script)
    channels: list[str] = []
    seen: set[str] = set()
    for line in result.stdout.splitlines():
        match = re.search(r"(V1:[^\s]+)", line.strip())
        if not match:
            continue
        name = match.group(1)
        if name not in seen:
            channels.append(name)
            seen.add(name)
    if not channels:
        raise RuntimeError("FrChannels returned no V1:* channel names")
    return channels


def filter_channels(
    channels: list[str], towers: list[str], kinds: list[str]
) -> tuple[dict[str, dict[str, list[str]]], list[str], list[str]]:
    tower_set = {tower.upper() for tower in towers}
    kind_set = {kind.upper() for kind in kinds}
    grouped = {
        tower: {kind: [] for kind in kinds}
        for tower in towers
    }
    matched: list[str] = []
    other_asc: list[str] = []

    for channel in sorted(channels):
        match = ALL_ASC_KIND_RE.match(channel)
        if not match:
            continue
        scope = match.group("scope").upper()
        kind = match.group("kind").upper()
        if kind not in kind_set:
            continue
        if scope in tower_set:
            grouped[scope][kind].append(channel)
            matched.append(channel)
        else:
            other_asc.append(channel)

    return grouped, matched, other_asc


def render_report(
    *,
    host: str,
    run_name: str,
    remote_path: str,
    frame_start: int,
    total_channels: int,
    towers: list[str],
    kinds: list[str],
    grouped: dict[str, dict[str, list[str]]],
    matched: list[str],
    other_asc: list[str],
) -> str:
    lines = [
        "ASC raw ERR/CORR channel report",
        f"host: {host}",
        f"run: {run_name}",
        f"frame_start: {frame_start}",
        f"remote_path: {remote_path}",
        f"total_raw_channels: {total_channels}",
        f"matched_tower_channels: {len(matched)}",
        "",
    ]

    for tower in towers:
        lines.append(f"[{tower}]")
        for kind in kinds:
            chans = grouped[tower][kind]
            lines.append(f"{kind} ({len(chans)}):")
            if chans:
                lines.extend(chans)
            else:
                lines.append("(none)")
        lines.append("")

    if other_asc:
        lines.append("[OTHER_ASC_ERR_CORR]")
        lines.append(f"count: {len(other_asc)}")
        lines.extend(other_asc)
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()
    towers = [tower.upper() for tower in args.towers]
    kinds = [kind.upper() for kind in args.types]

    try:
        run_name, remote_path, frame_start = resolve_remote_path(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Resolved raw frame: {remote_path}")
    if args.dry_run:
        return 0

    print(f"Listing channels via ssh {args.host} ...", flush=True)
    try:
        channels = list_remote_channels(args.host, remote_path)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    grouped, matched, other_asc = filter_channels(channels, towers, kinds)

    OUTPUTS_DIR.mkdir(exist_ok=True)
    stem = f"raw_asc_channels_{run_name}_{frame_start}"
    report_path = OUTPUTS_DIR / f"{stem}.txt"
    json_path = OUTPUTS_DIR / f"{stem}.json"
    all_channels_path = OUTPUTS_DIR / f"{stem}_all_channels.txt"

    report = render_report(
        host=args.host,
        run_name=run_name,
        remote_path=remote_path,
        frame_start=frame_start,
        total_channels=len(channels),
        towers=towers,
        kinds=kinds,
        grouped=grouped,
        matched=matched,
        other_asc=other_asc,
    )
    report_path.write_text(report)

    payload = {
        "host": args.host,
        "run": run_name,
        "frame_start": frame_start,
        "remote_path": remote_path,
        "total_raw_channels": len(channels),
        "towers": towers,
        "types": kinds,
        "matched_tower_channels": matched,
        "grouped": grouped,
        "other_asc_err_corr": other_asc,
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    if args.save_all_channels:
        all_channels_path.write_text("\n".join(channels) + "\n")

    print(report, end="")
    print(f"Saved report: {report_path}")
    print(f"Saved JSON:   {json_path}")
    if args.save_all_channels:
        print(f"Saved full raw channel list: {all_channels_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
