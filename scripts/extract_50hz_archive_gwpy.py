#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import multiprocessing as mp
import socket
import time
from pathlib import Path

import pandas as pd

CHANNELS = [
    ("V1:DET_B1p_DC_50Hz", "det_b1p_dc"),
    ("V1:LSC_B1p_DC_50Hz", "b1p_dc"),
    ("V1:LSC_Etalon_NI_RH_SET", "ni_rh_set"),
    ("V1:LSC_Etalon_NI_RH_OUT", "ni_rh_out"),
    ("V1:LSC_Etalon_NI_RH_IN", "rh_ni_in"),
    ("V1:LSC_Etalon_NI_RH_ERR", "rh_ni_err"),
]


def gps_prefixes(start: int, end: int) -> list[str]:
    first = start // 10_000_000
    last = end // 10_000_000
    return [str(v) for v in range(first, last + 1)]


def list_files(base: str, gps_start: int, gps_end: int) -> list[str]:
    files: list[str] = []
    for prefix in gps_prefixes(gps_start, gps_end):
        d = Path(base) / prefix
        if not d.is_dir():
            continue
        for path in sorted(d.glob("V-RDS-*-2000.gwf")):
            try:
                frame_start = int(path.name.split("-")[2])
            except Exception:
                continue
            frame_end = frame_start + 2000
            if frame_end > gps_start and frame_start < gps_end:
                files.append(str(path))
    return files


def make_chunks(gps_start: int, gps_end: int, chunk_size: int) -> list[tuple[int, int]]:
    chunks: list[tuple[int, int]] = []
    cur = gps_start
    while cur < gps_end:
        nxt = min(cur + chunk_size, gps_end)
        chunks.append((cur, nxt))
        cur = nxt
    return chunks


def read_chunk(task: tuple[str, int, int, int]) -> pd.DataFrame:
    archive_root, start, end, bin_size = task
    from gwpy.timeseries import TimeSeriesDict

    files = list_files(archive_root, start, end)
    if not files:
        return pd.DataFrame(columns=["gps_bin", *(key for _, key in CHANNELS)])

    tsd = TimeSeriesDict.read(files, [ch for ch, _ in CHANNELS], start=start, end=end, nproc=1)
    rows: list[pd.DataFrame] = []
    for ch, key in CHANNELS:
        if ch not in tsd:
            continue
        ts = tsd[ch]
        df = pd.DataFrame(
            {
                "gps_bin": (ts.times.value.astype("int64") // bin_size) * bin_size,
                key: ts.value.astype(float),
            }
        )
        rows.append(df.groupby("gps_bin", as_index=False).mean(numeric_only=True))
    if not rows:
        return pd.DataFrame(columns=["gps_bin", *(key for _, key in CHANNELS)])
    return (
        pd.concat(rows, ignore_index=True)
        .groupby("gps_bin", as_index=False)
        .mean(numeric_only=True)
        .sort_values("gps_bin")
    )


def checkpoint_path(output: Path) -> Path:
    return output.with_name(f"{output.stem}.partial.csv")


def save_partial(parts: list[pd.DataFrame], output: Path) -> None:
    if not parts:
        return
    partial = (
        pd.concat(parts, ignore_index=True)
        .groupby("gps_bin", as_index=False)
        .mean(numeric_only=True)
        .sort_values("gps_bin")
    )
    partial.to_csv(checkpoint_path(output), index=False)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive-root", default="/data/archive_50Hz")
    ap.add_argument("--gps-start", type=int, required=True)
    ap.add_argument("--gps-end", type=int, required=True)
    ap.add_argument("--chunk-size", type=int, default=21600)
    ap.add_argument("--bin-size", type=int, default=3600)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--checkpoint-every", type=int, default=10)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    chunks = make_chunks(args.gps_start, args.gps_end, args.chunk_size)
    tasks = [(args.archive_root, start, end, args.bin_size) for start, end in chunks]
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    host = socket.gethostname()
    print(
        f"[{host}] start chunks={len(tasks)} gps=({args.gps_start},{args.gps_end}) "
        f"chunk_size={args.chunk_size} bin_size={args.bin_size} jobs={args.jobs}",
        flush=True,
    )
    t0 = time.monotonic()

    if args.jobs > 1:
        parts: list[pd.DataFrame] = []
        done = 0
        with mp.Pool(processes=args.jobs) as pool:
            for start, end, df in (
                (task[1], task[2], df) for task, df in zip(tasks, pool.imap(read_chunk, tasks))
            ):
                done += 1
                if not df.empty:
                    parts.append(df)
                elapsed = time.monotonic() - t0
                avg = elapsed / done
                remain = max(len(tasks) - done, 0)
                eta_s = avg * remain
                print(
                    f"[{host}] chunk {done}/{len(tasks)} gps=({start},{end}) "
                    f"rows={len(df)} elapsed={elapsed/60:.1f}m avg={avg:.1f}s "
                    f"eta={eta_s/3600:.1f}h",
                    flush=True,
                )
                if done % max(args.checkpoint_every, 1) == 0:
                    save_partial(parts, out_path)
    else:
        parts = []
        for done, task in enumerate(tasks, start=1):
            start, end = task[1], task[2]
            df = read_chunk(task)
            if not df.empty:
                parts.append(df)
            elapsed = time.monotonic() - t0
            avg = elapsed / done
            remain = max(len(tasks) - done, 0)
            eta_s = avg * remain
            print(
                f"[{host}] chunk {done}/{len(tasks)} gps=({start},{end}) "
                f"rows={len(df)} elapsed={elapsed/60:.1f}m avg={avg:.1f}s "
                f"eta={eta_s/3600:.1f}h",
                flush=True,
            )
            if done % max(args.checkpoint_every, 1) == 0:
                save_partial(parts, out_path)

    if not parts:
        raise SystemExit("No data extracted")

    out = (
        pd.concat(parts, ignore_index=True)
        .groupby("gps_bin", as_index=False)
        .mean(numeric_only=True)
        .sort_values("gps_bin")
    )
    out.to_csv(out_path, index=False)
    save_partial(parts, out_path)
    elapsed = time.monotonic() - t0
    print(f"[{host}] saved {len(out)} rows -> {out_path} in {elapsed/3600:.2f} h", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
