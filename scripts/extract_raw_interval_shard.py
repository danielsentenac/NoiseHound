#!/usr/bin/env python3
"""
Extract one raw GWF file into a compact shard for interval rendering.

Each output shard stores per-channel time/value arrays for the requested
interval overlap, avoiding the huge sparse outer-joined CSV used by the
direct renderer.
"""
import argparse
import gzip
import os
import pickle
import re
import time
import warnings; warnings.filterwarnings("ignore")

from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-noisehound")

import numpy as np


RAW_CHANNELS = [
    ("V1:Hrec_hoft_16384Hz", "black"),
    ("V1:LSC_DARM_ERR", "dimgray"),
    ("V1:LSC_DARM_CORR", "tab:blue"),
    ("V1:LSC_MICH_CORR", "tab:purple"),
    ("V1:LSC_PRCL_CORR", "mediumorchid"),
    ("V1:LSC_SRCL_CORR", "plum"),
    ("V1:LSC_CARM_CORR", "indigo"),
    ("V1:LSC_BS_CORR", "tab:brown"),
    ("V1:LSC_PR_CORR", "sienna"),
    ("V1:LSC_SR_CORR", "peru"),
    ("V1:LSC_NI_CORR", "tab:olive"),
    ("V1:LSC_WI_CORR", "yellowgreen"),
    ("V1:LSC_NE_CORR", "tab:green"),
    ("V1:LSC_WE_CORR", "limegreen"),
]

LSC_SUPP_CHANNELS = [
    ("V1:LSC_DARM_ERR", "dimgray"),
    ("V1:LSC_DARM_CORR", "tab:blue"),
    ("V1:LSC_MICH_CORR", "tab:purple"),
    ("V1:LSC_PRCL_CORR", "mediumorchid"),
    ("V1:LSC_SRCL_CORR", "plum"),
    ("V1:LSC_CARM_CORR", "indigo"),
    ("V1:LSC_BS_CORR", "tab:brown"),
    ("V1:LSC_PR_CORR", "sienna"),
    ("V1:LSC_SR_CORR", "peru"),
    ("V1:LSC_NI_CORR", "tab:olive"),
    ("V1:LSC_WI_CORR", "yellowgreen"),
    ("V1:LSC_NE_CORR", "tab:green"),
    ("V1:LSC_WE_CORR", "limegreen"),
]


FP_AUX_CHANNELS = [
    ("V1:LSC_Etalon_INPOWER",   "tab:orange"),
    ("V1:LSC_Etalon_NI_RH_IN",  "sienna"),
    ("V1:LSC_Etalon_WI_RH_IN",  "peru"),
    ("V1:ENV_NEB_UPS_VOLT_R_2000Hz", "tab:red"),
    ("V1:ENV_CEB_UPS_VOLT_R_2000Hz", "crimson"),
    ("V1:ENV_CEB_UPS_CURR_R_2000Hz", "salmon"),
    ("V1:SBE_SNEB_GEO_GRWE_raw_500Hz", "tab:purple"),
    ("V1:TCS_NI_CO2_PWRLAS",    "deepskyblue"),
    ("V1:Sc_WI_FF50HZ_P_ERR",   "magenta"),
    ("V1:SQB1_Cam1_FitWaistY",  "teal"),
]


FP_AUX3_CHANNELS = [
    ("V1:LSC_DARM_PSTAB0_COUPLING_100Hz",       "tab:orange"),
    ("V1:ASC_SR_TY_DCP_mag_B1_I_10Hz",           "tab:red"),
    ("V1:ASC_SR_TX_DCP_mag_B1_I_10Hz",           "crimson"),
    ("V1:ASC_SR_TX",                              "salmon"),
    ("V1:ASC_SR_TY",                              "tab:purple"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_I_10Hz",        "deepskyblue"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_Q_10Hz",        "steelblue"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_phi0_10Hz",     "teal"),
    ("V1:LSC_DARM_BPC_TY_COUPLING_100Hz",        "magenta"),
]


FP_AUX4_CHANNELS = [
    ("V1:LSC_DARM_PSTAB0_COUPLING_100Hz",       "tab:orange"),
    ("V1:ASC_SR_TY_DCP_HF_B1p_I_10Hz",          "gold"),
    ("V1:ASC_SR_TY_DCP_HF_B1p_Q_10Hz",          "darkorange"),
    ("V1:ASC_SR_TX_DCP_HF_B1p_I_10Hz",          "orangered"),
    ("V1:ASC_SR_TX_DCP_HF_B1p_Q_10Hz",          "firebrick"),
    ("V1:ASC_SR_TY_DCP_mag_B1_I_10Hz",           "tab:red"),
    ("V1:ASC_SR_TX_DCP_mag_B1_I_10Hz",           "crimson"),
    ("V1:ASC_SR_TX",                              "salmon"),
    ("V1:ASC_SR_TY",                              "tab:purple"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_I_10Hz",        "deepskyblue"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_Q_10Hz",        "steelblue"),
    ("V1:LSC_DARM_BPC_TY_COUPLING_100Hz",        "magenta"),
]


FP_AUX2_CHANNELS = [
    ("V1:LSC_DARM_PSTAB0_COUPLING_100Hz",       "tab:orange"),
    ("V1:LSC_DARM_SIB1_LINE_mag_100Hz",          "sienna"),
    ("V1:Sc_WI_FF50HZ_PHASE",                    "magenta"),
    ("V1:ASC_SR_TY_DCP_mag_B1_I_10Hz",           "tab:red"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_I_10Hz",        "crimson"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_Q_10Hz",        "salmon"),
    ("V1:ASC_DIFFp_TY_DCP_mag_B1_phi0_10Hz",     "tab:purple"),
    ("V1:LSC_DARM_BPC_TY_COUPLING_100Hz",        "deepskyblue"),
    ("V1:Sc_WI_FF50HZ_G_ERR",                    "teal"),
]


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, seconds = divmod(seconds, 60.0)
    if minutes < 60:
        return f"{int(minutes)}m {seconds:.1f}s"
    hours, minutes = divmod(minutes, 60.0)
    return f"{int(hours)}h {int(minutes)}m {seconds:.0f}s"


def log_progress(message: str) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    print(f"[{stamp}] {message}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gwf", type=Path, required=True, help="Input raw GWF file.")
    parser.add_argument("--gps-start", type=float, required=True, help="Interval start GPS.")
    parser.add_argument("--gps-end", type=float, required=True, help="Interval end GPS.")
    parser.add_argument("--output", type=Path, required=True, help="Output shard path (.pkl.gz).")
    parser.add_argument("--channel-set", choices=("all", "lsc-supp", "fp-aux", "fp-aux2", "fp-aux3", "fp-aux4"), default="all",
                        help="Which channel set to extract: 'all' (default RAW_CHANNELS) "
                             "or 'lsc-supp' (LSC_SUPP_CHANNELS only).")
    return parser.parse_args()


def infer_file_bounds(gwf_path: Path) -> tuple[float, float]:
    match = re.search(r"V-raw-(\d+)-(\d+)\.gwf$", gwf_path.name)
    if not match:
        raise ValueError(f"Could not infer raw-file bounds from {gwf_path}")
    start = float(match.group(1))
    duration = float(match.group(2))
    return start, start + duration


def main():
    args = parse_args()
    if args.gps_end <= args.gps_start:
        raise SystemExit("--gps-end must be greater than --gps-start")
    if not args.gwf.exists():
        raise SystemExit(f"Missing GWF file: {args.gwf}")

    file_start, file_end = infer_file_bounds(args.gwf)
    seg_start = max(args.gps_start, file_start)
    seg_end = min(args.gps_end, file_end)

    log_progress(
        f"Shard extraction started for {args.gwf.name}: interval {args.gps_start:.0f}-{args.gps_end:.0f}, "
        f"file segment {file_start:.0f}-{file_end:.0f}"
    )

    from gwpy.io.gwf import iter_channel_names
    from gwpy.timeseries import TimeSeriesDict

    if seg_end <= seg_start:
        raise SystemExit(f"No overlap between interval and {args.gwf.name}")

    channel_names = set(iter_channel_names(str(args.gwf)))
    if args.channel_set == "lsc-supp":
        source_list = LSC_SUPP_CHANNELS
    elif args.channel_set == "fp-aux":
        source_list = FP_AUX_CHANNELS
    elif args.channel_set == "fp-aux2":
        source_list = FP_AUX2_CHANNELS
    elif args.channel_set == "fp-aux3":
        source_list = FP_AUX3_CHANNELS
    elif args.channel_set == "fp-aux4":
        source_list = FP_AUX4_CHANNELS
    else:
        source_list = RAW_CHANNELS
    selected = []
    for channel, _color in source_list:
        if channel in channel_names:
            selected.append(channel)
        else:
            log_progress(f"Skipping unavailable channel in {args.gwf.name}: {channel}")

    if not selected:
        raise RuntimeError(f"No requested channels available in {args.gwf.name}")

    payload = {
        "source_gwf": str(args.gwf),
        "segment_start": seg_start,
        "segment_end": seg_end,
        "channels": {},
    }

    read_started = time.perf_counter()
    log_progress(f"Reading {len(selected)} channels from {args.gwf.name} in one pass")
    tsd = TimeSeriesDict.read(str(args.gwf), selected, start=seg_start, end=seg_end)
    for channel in selected:
        ts = tsd[channel]
        payload["channels"][channel] = {
            "times": np.asarray(ts.times.value, dtype=np.float64),
            "values": np.asarray(ts.value, dtype=np.float32),
        }
    log_progress(
        f"Loaded {len(selected)} channels from {args.gwf.name} in "
        f"{format_elapsed(time.perf_counter() - read_started)}"
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_started = time.perf_counter()
    with gzip.open(args.output, "wb") as fp:
        pickle.dump(payload, fp, protocol=pickle.HIGHEST_PROTOCOL)
    log_progress(
        f"Wrote shard {args.output} in {format_elapsed(time.perf_counter() - write_started)} "
        f"({args.output.stat().st_size / 1e6:.1f} MB); total extraction time "
        f"{format_elapsed(time.perf_counter() - read_started)}"
    )


if __name__ == "__main__":
    main()
