#!/usr/bin/env python3
"""
Plot a longer raw-window view for the Step-2 top-ranked channels.

Goal:
- use the same Step-2 event selection as the onset exercise
- show only literal raw counterparts that actually exist in the raw frame
- use a long enough window to make the reported +1.5 to +2 s delays visible
"""

from __future__ import annotations

import json
import subprocess
import warnings
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
WINDOW_BEFORE = 0.5
WINDOW_AFTER = 3.0
HPSS_ROOT = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4b/raw"
SSH_HOST = "cca.in2p3.fr"
REMOTE_PYTHON = "/cvmfs/software.igwn.org/conda/envs/igwn/bin/python3"

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"
TOP5 = [
    {
        "rank": 1,
        "trend": "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean",
        "raw": "V1:LSC_DARM_IMC_LINE_mag_100Hz",
        "label": "IMC line mag",
        "lag_s": 1.5,
        "color": "tab:blue",
    },
    {
        "rank": 2,
        "trend": "V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean",
        "raw": "V1:LSC_DARM_PSTAB0_COUPLING_100Hz",
        "label": "PSTAB0 coupling",
        "lag_s": 1.5,
        "color": "tab:orange",
    },
    {
        "rank": 3,
        "trend": "V1:LSC_DARM_IMC_LINE_Q_100Hz_mean",
        "raw": "V1:LSC_DARM_IMC_LINE_Q_100Hz",
        "label": "IMC line Q",
        "lag_s": 2.0,
        "color": "tab:green",
    },
    {
        "rank": 4,
        "trend": "V1:LSC_DARM_PSTAB0_I_FS_mean",
        "raw": "V1:LSC_DARM_PSTAB0_I_FS",
        "label": "PSTAB0 I FS",
        "lag_s": 2.0,
        "color": "tab:red",
    },
    {
        "rank": 5,
        "trend": "V1:LSC_DARM_PSTAB0_I_100Hz_mean",
        "raw": "V1:LSC_DARM_PSTAB0_I_100Hz",
        "label": "PSTAB0 I 100 Hz",
        "lag_s": 2.0,
        "color": "tab:purple",
    },
]


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
    t_ref = float(event["time"])
    start = t_ref - WINDOW_BEFORE
    end = t_ref + WINDOW_AFTER
    remote_gwf = remote_path(t_ref)
    channels = [REF_CHANNEL] + [item["raw"] for item in TOP5]
    channels_json = json.dumps(channels)

    remote_script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_step2_delay.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT

local_gwf="$stage_dir/V-raw-{frame_start(t_ref)}-100.gwf"
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
        rows.append(pd.DataFrame({{
            "channel": ch,
            "gps": ts.times.value,
            "value": ts.value.astype(float),
        }}))

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


def bandpass_hrec(frame: pd.DataFrame, t_ref: float) -> tuple[np.ndarray, np.ndarray]:
    t_rel = frame["gps"].to_numpy(dtype=float) - t_ref
    values = frame["value"].to_numpy(dtype=float)
    if len(values) < 5:
        return t_rel, values
    fs = 1.0 / np.median(np.diff(frame["gps"].to_numpy(dtype=float)))
    if fs <= 2 * 55:
        return t_rel, values
    sos = butter(4, [30, 55], btype="band", fs=fs, output="sos")
    return t_rel, sosfiltfilt(sos, values)


def z_norm(t_rel: np.ndarray, values: np.ndarray) -> np.ndarray:
    baseline = (t_rel >= -WINDOW_BEFORE) & (t_rel <= -0.1)
    if baseline.sum() < 3:
        baseline = t_rel < 0
    if baseline.sum() < 2:
        baseline = np.ones_like(t_rel, dtype=bool)
    mu = float(np.mean(values[baseline]))
    sigma = float(np.std(values[baseline]))
    if not np.isfinite(sigma) or sigma == 0:
        sigma = max(float(np.max(np.abs(values - mu))), 1.0)
    return (values - mu) / sigma


def plot(event: pd.Series, meta: dict, df: pd.DataFrame, png_path: Path) -> None:
    t_ref = float(event["time"])
    t_start = float(event["tstart"])
    start_offset = t_start - t_ref
    found_top = [item for item in TOP5 if item["raw"] in meta["found"]]

    n_rows = 1 + len(found_top)
    fig, axes = plt.subplots(n_rows, 1, figsize=(15, 2.7 * n_rows), sharex=True)
    if n_rows == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.08, top=0.90, bottom=0.11, left=0.11, right=0.97)

    ref_frame = get_channel_frame(df, REF_CHANNEL)
    ax0 = axes[0]
    if not ref_frame.empty:
        t_rel, values = bandpass_hrec(ref_frame, t_ref)
        ax0.plot(t_rel, z_norm(t_rel, values), color="black", lw=0.8)
    ax0.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax0.axvline(start_offset, color="dimgray", lw=1.0, ls=":", alpha=0.8)
    ax0.set_ylabel("Hrec\nz", fontsize=9)
    ax0.set_title("V1:Hrec_hoft_16384Hz  (30-55 Hz bandpass)", fontsize=9, loc="left")
    ax0.grid(axis="x", ls=":", alpha=0.3)

    for ax, item in zip(axes[1:], found_top):
        frame = get_channel_frame(df, item["raw"])
        t_rel = frame["gps"].to_numpy(dtype=float) - t_ref
        values = frame["value"].to_numpy(dtype=float)
        y = z_norm(t_rel, values)
        ax.plot(t_rel, y, color=item["color"], lw=1.1, marker="o", ms=3)
        ax.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
        ax.axvline(start_offset, color="dimgray", lw=1.0, ls=":", alpha=0.8)
        ax.axvline(item["lag_s"], color=item["color"], lw=1.0, ls=":", alpha=0.7)
        ax.set_ylabel("z", fontsize=9)
        ax.set_title(
            f"Rank {item['rank']}: {item['trend']} -> {item['raw']}  (table lag {item['lag_s']:+.1f} s)",
            fontsize=8.5,
            loc="left",
        )
        ax.grid(axis="x", ls=":", alpha=0.3)

    axes[-1].set_xlim(-WINDOW_BEFORE, WINDOW_AFTER)
    axes[-1].set_xticks(np.arange(-0.5, 3.01, 0.5))
    axes[-1].set_xlabel("Time relative to glitch peak GPS (s)", fontsize=9)

    missing_trend = [item["trend"] for item in TOP5 if item["raw"] in meta["missing"]]
    title = (
        "Step 2 top-ranked channels in raw data over a lag-resolving window\n"
        f"peak={t_ref:.6f}, start={t_start:.6f}, SNR={float(event['snr']):.1f}, "
        f"f={float(event['frequency']):.2f} Hz  |  window = -0.5 s to +3.0 s"
    )
    fig.suptitle(title, fontsize=10, fontweight="bold")
    if missing_trend:
        fig.text(
            0.11,
            0.04,
            "Trend-only in this comparison, omitted from raw plot: " + "; ".join(missing_trend),
            fontsize=8,
            color="firebrick",
        )
    fig.text(
        0.11,
        0.02,
        "Crimson dashed = glitch peak. Gray dotted = Omicron tstart. Colored dotted = Step 2 lag from table.",
        fontsize=8,
    )

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    event = choose_event()
    event_id = int(float(event["time"]))
    csv_path = OUTPUTS / f"step2_top_raw_delay_{event_id}.csv"
    meta_path = OUTPUTS / f"step2_top_raw_delay_{event_id}.json"
    png_path = OUT_DIR / f"step2_top_raw_delay_{event_id}_3p5s.png"

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
    plot(event, meta, df, png_path)

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
