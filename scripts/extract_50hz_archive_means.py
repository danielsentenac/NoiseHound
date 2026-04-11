#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
BIN_SIZE = 3600

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


def list_files(base: Path, gps_start: int, gps_end: int) -> list[Path]:
    files: list[Path] = []
    for prefix in gps_prefixes(gps_start, gps_end):
        d = base / prefix
        if not d.is_dir():
            continue
        for path in sorted(d.glob("V-RDS-*-2000.gwf")):
            try:
                frame_start = int(path.name.split("-")[2])
            except Exception:
                continue
            frame_end = frame_start + 2000
            if frame_end > gps_start and frame_start < gps_end:
                files.append(path)
    return files


def run_copy(input_file: Path, output_file: Path, decimate: int) -> None:
    channel_list = " ".join(ch for ch, _ in CHANNELS)
    subprocess.run(
        [
            "FrCopy",
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "-t",
            channel_list,
            "-decimate",
            str(decimate),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def read_decimated(output_file: Path) -> pd.DataFrame:
    from gwpy.timeseries import TimeSeriesDict

    tsd = TimeSeriesDict.read(str(output_file), [ch for ch, _ in CHANNELS], nproc=1)
    rows: list[dict[str, float]] = []
    for channel, key in CHANNELS:
        if channel not in tsd:
            continue
        ts = tsd[channel]
        t_arr = [int(math.floor(t / BIN_SIZE)) * BIN_SIZE for t in ts.times.value]
        v_arr = [float(v) for v in ts.value]
        for gps_bin, val in zip(t_arr, v_arr):
            rows.append({"gps_bin": gps_bin, key: val})
    if not rows:
        return pd.DataFrame(columns=["gps_bin", *(key for _, key in CHANNELS)])
    return (
        pd.DataFrame(rows)
        .groupby("gps_bin", as_index=False)
        .mean(numeric_only=True)
        .sort_values("gps_bin")
    )


def process_file(task: tuple[int, int, str, int]) -> pd.DataFrame:
    idx, total, input_name, decimate = task
    input_file = Path(input_name)
    stage_dir = Path(tempfile.mkdtemp(prefix="nh_archive50_worker_"))
    out_file = stage_dir / f"{input_file.stem}_dec.gwf"
    try:
        run_copy(input_file, out_file, decimate)
        df = read_decimated(out_file)
        print(f"[{idx}/{total}] {input_file.name} OK", flush=True)
        return df
    except Exception as exc:
        print(f"[{idx}/{total}] {input_file.name} FAILED: {exc}", flush=True)
        return pd.DataFrame()
    finally:
        out_file.unlink(missing_ok=True)
        shutil.rmtree(stage_dir, ignore_errors=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive-root", default="/data/archive_50Hz")
    ap.add_argument("--gps-start", type=int, required=True)
    ap.add_argument("--gps-end", type=int, required=True)
    ap.add_argument("--decimate", type=int, default=100000)
    ap.add_argument("--bin-size", type=int, default=3600)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    global BIN_SIZE
    BIN_SIZE = int(args.bin_size)

    archive_root = Path(args.archive_root)
    files = list_files(archive_root, args.gps_start, args.gps_end)
    if not files:
        raise SystemExit("No archive 50Hz files found in the requested window")

    tasks = [(idx, len(files), str(path), args.decimate) for idx, path in enumerate(files, start=1)]
    if args.jobs > 1:
        with mp.Pool(processes=args.jobs) as pool:
            parts = [df for df in pool.imap_unordered(process_file, tasks) if not df.empty]
    else:
        parts = [df for df in map(process_file, tasks) if not df.empty]

    if not parts:
        raise SystemExit("No decimated data could be extracted")

    out = (
        pd.concat(parts, ignore_index=True)
        .groupby("gps_bin", as_index=False)
        .mean(numeric_only=True)
        .sort_values("gps_bin")
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"Saved {len(out)} rows -> {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
