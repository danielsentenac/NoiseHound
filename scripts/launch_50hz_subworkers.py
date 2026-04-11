#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path


def parse_range(text: str) -> tuple[int, int]:
    if ":" not in text:
        raise argparse.ArgumentTypeError(f"range must be start:end, got {text!r}")
    start_s, end_s = text.split(":", 1)
    start, end = int(start_s), int(end_s)
    if end <= start:
        raise argparse.ArgumentTypeError(f"invalid range {text!r}")
    return start, end


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--python", default="python")
    ap.add_argument("--extractor", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--chunk-size", type=int, default=3600)
    ap.add_argument("--checkpoint-every", type=int, default=1)
    ap.add_argument("--clear", action="store_true")
    ap.add_argument("--range", dest="ranges", action="append", type=parse_range, required=True)
    args = ap.parse_args()

    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    for idx, (start, end) in enumerate(args.ranges):
        csv_path = prefix.parent / f"{prefix.name}_part{idx}.csv"
        partial_path = prefix.parent / f"{prefix.name}_part{idx}.partial.csv"
        log_path = prefix.parent / f"{prefix.name}_part{idx}.log"
        if args.clear:
            for path in [csv_path, partial_path, log_path]:
                path.unlink(missing_ok=True)
        log_fh = open(log_path, "ab", buffering=0)
        cmd = [
            args.python,
            "-u",
            args.extractor,
            "--gps-start",
            str(start),
            "--gps-end",
            str(end),
            "--chunk-size",
            str(args.chunk_size),
            "--jobs",
            "1",
            "--checkpoint-every",
            str(args.checkpoint_every),
            "--output",
            str(csv_path),
        ]
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.DEVNULL,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        print(f"part{idx}:{proc.pid}:{start}:{end}:{csv_path}:{log_path}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
