#!/usr/bin/env python3
"""
Plot the onset of one Step-2 glitch for the top five ranked channels.

Important limitation:
Step 2 ranked 1 Hz trend products. For a 100 ms onset plot we need raw-frame
counterparts. Some of the top-ranked trend products do not exist verbatim in
raw frames; those panels are marked explicitly as unavailable instead of
guessing a proxy channel.
"""

from __future__ import annotations

import json
import subprocess
import warnings
from io import StringIO
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


WORKDIR = Path(__file__).resolve().parents[1]
CATALOG = WORKDIR / "data" / "full_25min_glitches_ER16-O4b.csv"
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"

STEP2_GPS_START = 1419724818
STEP2_GPS_END = 1419735618
WINDOW_BEFORE = 0.200
WINDOW_AFTER = 0.200
HPSS_ROOT = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
SSH_HOST = "cca.in2p3.fr"
REMOTE_PYTHON = "/cvmfs/software.igwn.org/conda/envs/igwn/bin/python3"

TOP5 = [
    {
        "rank": 1,
        "trend": "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean",
        "raw": "V1:LSC_DARM_IMC_LINE_mag_100Hz",
        "label": "LSC_DARM_IMC_LINE_mag_100Hz",
    },
    {
        "rank": 2,
        "trend": "V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean",
        "raw": "V1:LSC_DARM_PSTAB0_COUPLING_100Hz",
        "label": "LSC_DARM_PSTAB0_COUPLING_100Hz",
    },
    {
        "rank": 3,
        "trend": "V1:LSC_DARM_IMC_LINE_Q_100Hz_mean",
        "raw": "V1:LSC_DARM_IMC_LINE_Q_100Hz",
        "label": "LSC_DARM_IMC_LINE_Q_100Hz",
    },
    {
        "rank": 4,
        "trend": "V1:LSC_DARM_PSTAB0_I_FS_mean",
        "raw": "V1:LSC_DARM_PSTAB0_I_FS",
        "label": "LSC_DARM_PSTAB0_I_FS",
    },
    {
        "rank": 5,
        "trend": "V1:LSC_DARM_PSTAB0_I_100Hz_mean",
        "raw": "V1:LSC_DARM_PSTAB0_I_100Hz",
        "label": "LSC_DARM_PSTAB0_I_100Hz",
    },
]

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"


def choose_event() -> pd.Series:
    df = pd.read_csv(CATALOG)
    epoch = df[(df["time"] >= STEP2_GPS_START) & (df["time"] <= STEP2_GPS_END)].copy()
    if epoch.empty:
        raise RuntimeError("no glitches found in the Step-2 epoch")
    return epoch.sort_values(["snr", "time"], ascending=[False, True]).iloc[0]


def frame_start(gps: float) -> int:
    return int(gps // 100) * 100


def remote_path(gps: float) -> str:
    start = frame_start(gps)
    return f"{HPSS_ROOT}/{str(start)[:4]}/V-raw-{start}-100.gwf"


def run_remote_extract(event: pd.Series, csv_path: Path) -> dict:
    start = float(event["tstart"]) - WINDOW_BEFORE
    end = float(event["tstart"]) + WINDOW_AFTER
    remote_gwf = remote_path(float(event["time"]))
    channels = [REF_CHANNEL] + [item["raw"] for item in TOP5]
    channels_json = json.dumps(channels)

    remote_script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_step2_top5.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT

local_gwf="$stage_dir/V-raw-{frame_start(float(event["time"]))}-100.gwf"
rfcp {remote_gwf} "$local_gwf" >/dev/null

export GWF_PATH="$local_gwf"
export T0="{start:.9f}"
export T1="{end:.9f}"
export CHANNELS_JSON='{channels_json}'

"{REMOTE_PYTHON}" - <<'PY'
import json
import os
import sys
import pandas as pd
from gwpy.io.gwf import iter_channel_names
from gwpy.timeseries import TimeSeriesDict

gwf = os.environ["GWF_PATH"]
t0 = float(os.environ["T0"])
t1 = float(os.environ["T1"])
requested = json.loads(os.environ["CHANNELS_JSON"])

available = set(iter_channel_names(gwf))
found = [ch for ch in requested if ch in available]
missing = [ch for ch in requested if ch not in available]

rows = []
if found:
    tsd = TimeSeriesDict.read(gwf, found, start=t0, end=t1, nproc=1)
    for ch in found:
        ts = tsd[ch]
        frame = pd.DataFrame({{
            "channel": ch,
            "gps": ts.times.value,
            "value": ts.value.astype(float),
        }})
        rows.append(frame)

meta = {{
    "requested": requested,
    "found": found,
    "missing": missing,
}}
print("#META " + json.dumps(meta, sort_keys=True))
if rows:
    pd.concat(rows, ignore_index=True).to_csv(sys.stdout, index=False)
else:
    print("channel,gps,value")
PY
"""

    result = subprocess.run(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=15",
            SSH_HOST,
            "bash",
            "-s",
            "--",
        ],
        input=remote_script,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "remote extraction failed")

    stdout = result.stdout
    if not stdout.startswith("#META "):
        raise RuntimeError("remote extraction did not emit metadata header")
    meta_line, csv_body = stdout.split("\n", 1)
    meta = json.loads(meta_line[len("#META ") :])
    csv_path.write_text(csv_body)
    return meta


def get_channel_frame(df: pd.DataFrame, channel: str) -> pd.DataFrame:
    return df[df["channel"] == channel].copy().sort_values("gps")


def preprocess_reference(frame: pd.DataFrame, t_ref: float) -> tuple[np.ndarray, np.ndarray]:
    t_rel = frame["gps"].to_numpy() - t_ref
    values = frame["value"].to_numpy(dtype=float)
    if len(values) < 5:
        return t_rel, values
    fs = 1.0 / np.median(np.diff(frame["gps"].to_numpy()))
    if fs <= 2 * 55:
        return t_rel, values
    sos = butter(4, [30, 55], btype="band", fs=fs, output="sos")
    return t_rel, sosfiltfilt(sos, values)


def normalise_trace(t_rel: np.ndarray, values: np.ndarray) -> np.ndarray:
    baseline = (t_rel >= -WINDOW_BEFORE) & (t_rel <= -0.010)
    if baseline.sum() < 3:
        baseline = t_rel < 0
    if baseline.sum() < 2:
        baseline = np.ones_like(t_rel, dtype=bool)
    mu = float(np.mean(values[baseline]))
    sigma = float(np.std(values[baseline]))
    if not np.isfinite(sigma) or sigma == 0:
        sigma = max(float(np.max(np.abs(values - mu))), 1.0)
    return (values - mu) / sigma


def plot(meta: dict, df: pd.DataFrame, event: pd.Series, png_path: Path) -> None:
    t_start = float(event["tstart"])
    t_peak = float(event["time"])
    peak_offset_ms = (t_peak - t_start) * 1000.0

    fig, axes = plt.subplots(1 + len(TOP5), 1, figsize=(15, 11), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.07, left=0.13, right=0.97)

    ref = get_channel_frame(df, REF_CHANNEL)
    ax0 = axes[0]
    if not ref.empty:
        t_rel, values = preprocess_reference(ref, t_start)
        y = normalise_trace(t_rel, values)
        ax0.plot(t_rel * 1000.0, y, color="black", lw=0.9)
    else:
        ax0.text(
            0.5,
            0.5,
            f"{REF_CHANNEL} not available in raw frame",
            transform=ax0.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )
    ax0.axvline(0, color="crimson", lw=1.3, ls="--", alpha=0.8)
    ax0.axvline(peak_offset_ms, color="dimgray", lw=1.0, ls=":", alpha=0.8)
    ax0.set_ylabel("Hrec\nnorm.", fontsize=9)
    ax0.set_title("Reference: V1:Hrec_hoft_16384Hz  (30–55 Hz bandpass)", fontsize=9, loc="left")
    ax0.grid(axis="x", ls=":", alpha=0.3)

    for idx, item in enumerate(TOP5, start=1):
        ax = axes[idx]
        raw_name = item["raw"]
        panel = get_channel_frame(df, raw_name)
        if panel.empty:
            ax.text(
                0.5,
                0.5,
                "raw counterpart absent in frame\n(trend-only Step-2 product here)",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                color="firebrick",
            )
            ax.set_yticks([])
        else:
            t_rel = panel["gps"].to_numpy() - t_start
            values = panel["value"].to_numpy(dtype=float)
            y = normalise_trace(t_rel, values)
            ax.plot(t_rel * 1000.0, y, color="tab:blue", lw=1.0, marker="o", ms=3)
            ax.set_ylabel("norm.", fontsize=9)
        ax.axvline(0, color="crimson", lw=1.3, ls="--", alpha=0.8)
        ax.axvline(peak_offset_ms, color="dimgray", lw=1.0, ls=":", alpha=0.8)
        ax.grid(axis="x", ls=":", alpha=0.3)
        ax.set_title(
            f"Rank {item['rank']}: {item['trend']}  ->  {raw_name}",
            fontsize=8.5,
            loc="left",
        )

    x_min_ms = int(round(-WINDOW_BEFORE * 1000.0))
    x_max_ms = int(round(WINDOW_AFTER * 1000.0))
    xtick_step_ms = 50 if (WINDOW_BEFORE + WINDOW_AFTER) >= 0.4 else 10
    xticks = np.arange(x_min_ms, x_max_ms + xtick_step_ms, xtick_step_ms)
    for ax in axes:
        ax.set_xlim(x_min_ms, x_max_ms)
        ax.set_xticks(xticks)
    axes[-1].set_xlabel("Time relative to glitch start tstart (ms)", fontsize=9)

    fig.suptitle(
        "Step 2 top-ranked channels at glitch onset\n"
        f"Chosen event: peak={t_peak:.6f}, start={t_start:.6f}, "
        f"SNR={float(event['snr']):.1f}, f={float(event['frequency']):.2f} Hz  |  "
        f"window = {int(round((WINDOW_BEFORE + WINDOW_AFTER) * 1000.0))} ms",
        fontsize=10,
        fontweight="bold",
    )
    fig.text(
        0.13,
        0.015,
        "Crimson dashed line = Omicron tstart. Gray dotted line = catalog peak time.",
        fontsize=8,
    )
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    event = choose_event()
    event_id = int(float(event["time"]))
    csv_path = OUTPUTS / f"step2_top5_onset_{event_id}.csv"
    window_ms = int(round((WINDOW_BEFORE + WINDOW_AFTER) * 1000.0))
    png_path = OUT_DIR / f"step2_top5_onset_{event_id}_{window_ms}ms.png"
    meta_path = OUTPUTS / f"step2_top5_onset_{event_id}.json"

    meta = run_remote_extract(event, csv_path)
    meta.update(
        {
            "event_peak_gps": float(event["time"]),
            "event_start_gps": float(event["tstart"]),
            "event_end_gps": float(event["tend"]),
            "event_snr": float(event["snr"]),
            "event_frequency_hz": float(event["frequency"]),
            "window_before_s": WINDOW_BEFORE,
            "window_after_s": WINDOW_AFTER,
        }
    )
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")

    df = pd.read_csv(csv_path)
    plot(meta, df, event, png_path)

    print(f"Chosen Step-2 event peak:  {float(event['time']):.6f}")
    print(f"Chosen Step-2 event start: {float(event['tstart']):.6f}")
    print(f"Saved CSV:  {csv_path}")
    print(f"Saved meta: {meta_path}")
    print(f"Saved plot: {png_path}")
    if meta["missing"]:
        print("Missing raw counterparts:")
        for channel in meta["missing"]:
            print(f"  {channel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
