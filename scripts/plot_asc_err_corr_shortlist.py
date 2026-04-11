#!/usr/bin/env python3
"""
Plot ASC ERR/CORR shortlist around a pre-Christmas 2025 glitch event.

Creates two figures from one raw extraction:
- 400 ms window  (-0.2 s .. +0.2 s)
- 3.5 s window   (-0.5 s .. +3.0 s)

Each figure has three panels:
1) Hrec reference
2) ERR channels overlay
3) CORR channels overlay
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"

SSH_HOST = "cca.in2p3.fr"
REMOTE_PYTHON = "/cvmfs/software.igwn.org/conda/envs/igwn/bin/python3"
HPSS_ROOT = "cchpss0:/hpss/in2p3.fr/group/virgo/Run/O4c/raw"

# pre-Christmas 2025 event used in the tower-onset work
WINDOW_BEFORE = 0.5
WINDOW_AFTER = 3.0

REF_CHANNEL = "V1:Hrec_hoft_16384Hz"
SHORTLIST = [
    "V1:ASC_NI_TX_CORR",
    "V1:ASC_NI_TY_CORR",
    "V1:ASC_WI_TX_CORR",
    "V1:ASC_WI_TY_CORR",
    "V1:ASC_NE_TX_CORR",
    "V1:ASC_NE_TY_CORR",
    "V1:ASC_WE_TX_CORR",
    "V1:ASC_WE_TY_CORR",
    "V1:ASC_BS_TX_ERR",
    "V1:ASC_BS_TY_ERR",
    "V1:ASC_BS_TX_CORR",
    "V1:ASC_BS_TY_CORR",
    "V1:ASC_PR_TX_ERR",
    "V1:ASC_PR_TY_ERR",
    "V1:ASC_PR_TX_CORR",
    "V1:ASC_PR_TY_CORR",
    "V1:ASC_SR_TX_ERR",
    "V1:ASC_SR_TY_ERR",
    "V1:ASC_SR_DOF_TX_CORR",
    "V1:ASC_SR_DOF_TY_CORR",
]

ERR_CH = [c for c in SHORTLIST if c.endswith("_ERR")]
CORR_CH = [c for c in SHORTLIST if c.endswith("_CORR")]

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot ASC ERR/CORR shortlist around one glitch event.")
    p.add_argument(
        "--event-gps",
        type=float,
        default=1446106739.098,
        help="Glitch peak GPS time (default: 1446106739.098).",
    )
    p.add_argument(
        "--run",
        default="O4c",
        help="Run folder name under /Run (default: O4c).",
    )
    return p.parse_args()


def frame_start(gps: float) -> int:
    return int(gps // 100) * 100


def remote_raw_path(gps: float, run_name: str) -> str:
    start = frame_start(gps)
    return f"cchpss0:/hpss/in2p3.fr/group/virgo/Run/{run_name}/raw/{str(start)[:4]}/V-raw-{start}-100.gwf"


def extract_raw(event_peak: float, run_name: str, csv_out: Path, meta_out: Path) -> dict:
    t0 = event_peak - WINDOW_BEFORE
    t1 = event_peak + WINDOW_AFTER
    remote = remote_raw_path(event_peak, run_name)
    channels = [REF_CHANNEL] + SHORTLIST
    channels_json = json.dumps(channels)

    script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_asc_shortlist.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT
local_gwf="$stage_dir/V-raw-{frame_start(event_peak)}-100.gwf"
rfcp {remote} "$local_gwf" >/dev/null

export GWF_PATH="$local_gwf"
export T0="{t0:.6f}"
export T1="{t1:.6f}"
export CHANNELS_JSON='{channels_json}'
"{REMOTE_PYTHON}" - <<'PY'
import json, os, sys
import pandas as pd
from gwpy.io.gwf import iter_channel_names
from gwpy.timeseries import TimeSeriesDict

gwf = os.environ["GWF_PATH"]
t0 = float(os.environ["T0"])
t1 = float(os.environ["T1"])
requested = json.loads(os.environ["CHANNELS_JSON"])
avail = set(iter_channel_names(gwf))
found = [c for c in requested if c in avail]
missing = [c for c in requested if c not in avail]
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
meta = {{"requested": requested, "found": found, "missing": missing}}
print("#META " + json.dumps(meta, sort_keys=True))
if rows:
    pd.concat(rows, ignore_index=True).to_csv(sys.stdout, index=False)
else:
    print("channel,gps,value")
PY
"""

    res = subprocess.run(
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
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip() or "raw extraction failed")
    if not res.stdout.startswith("#META "):
        raise RuntimeError("missing metadata header from remote extraction")
    meta_line, csv_body = res.stdout.split("\n", 1)
    meta = json.loads(meta_line[len("#META ") :])
    csv_out.write_text(csv_body)
    meta_out.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return meta


def series(df: pd.DataFrame, channel: str, event_peak: float) -> tuple[np.ndarray, np.ndarray]:
    d = df[df["channel"] == channel].sort_values("gps")
    t = d["gps"].to_numpy(dtype=float) - event_peak
    v = d["value"].to_numpy(dtype=float)
    return t, v


def z_norm(t: np.ndarray, v: np.ndarray) -> np.ndarray:
    base = (t >= -WINDOW_BEFORE) & (t <= -0.1)
    if not np.any(base):
        base = t < 0
    if not np.any(base):
        base = np.ones_like(t, dtype=bool)
    mu = np.mean(v[base])
    sig = np.std(v[base])
    if not np.isfinite(sig) or sig == 0:
        sig = 1.0
    return (v - mu) / sig


def hrec_bandpass(t: np.ndarray, v: np.ndarray, event_peak: float) -> np.ndarray:
    if len(v) < 5:
        return v
    fs = 1.0 / np.median(np.diff(t + event_peak))
    if fs <= 2 * 55:
        return v
    sos = butter(4, [30, 55], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def color_for(channel: str) -> str:
    if "_NI_" in channel:
        return "#1f77b4"
    if "_WI_" in channel:
        return "#ff7f0e"
    if "_NE_" in channel:
        return "#2ca02c"
    if "_WE_" in channel:
        return "#d62728"
    if "_BS_" in channel:
        return "#9467bd"
    if "_PR_" in channel:
        return "#8c564b"
    if "_SR_" in channel:
        return "#e377c2"
    return "black"


def short_label(channel: str) -> str:
    return channel.replace("V1:ASC_", "")


def plot_window(
    df: pd.DataFrame,
    meta: dict,
    win_lo: float,
    win_hi: float,
    out_png: Path,
    event_peak: float,
) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.92, bottom=0.08, left=0.09, right=0.98)

    # panel 1: Hrec
    ax0 = axes[0]
    if REF_CHANNEL in meta["found"]:
        t, v = series(df, REF_CHANNEL, event_peak)
        y = z_norm(t, hrec_bandpass(t, v, event_peak))
        m = (t >= win_lo) & (t <= win_hi)
        ax0.plot(t[m], y[m], color="black", lw=0.9)
    ax0.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax0.set_title("V1:Hrec_hoft_16384Hz (30-55 Hz bandpass)", loc="left", fontsize=10)
    ax0.set_ylabel("Hrec z")
    ax0.grid(axis="x", ls=":", alpha=0.35)

    # panel 2: ERR overlay
    ax1 = axes[1]
    err_found = [c for c in ERR_CH if c in meta["found"]]
    for ch in err_found:
        t, v = series(df, ch, event_peak)
        y = z_norm(t, v)
        m = (t >= win_lo) & (t <= win_hi)
        ax1.plot(t[m], y[m], lw=1.1, alpha=0.9, color=color_for(ch), label=short_label(ch))
    ax1.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax1.set_title("ASC ERR channels", loc="left", fontsize=10)
    ax1.set_ylabel("z")
    ax1.grid(axis="x", ls=":", alpha=0.35)
    if err_found:
        ax1.legend(fontsize=8, ncol=3, loc="upper right")

    # panel 3: CORR overlay
    ax2 = axes[2]
    corr_found = [c for c in CORR_CH if c in meta["found"]]
    for ch in corr_found:
        t, v = series(df, ch, event_peak)
        y = z_norm(t, v)
        m = (t >= win_lo) & (t <= win_hi)
        ax2.plot(t[m], y[m], lw=1.0, alpha=0.85, color=color_for(ch), label=short_label(ch))
    ax2.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax2.set_title("ASC CORR channels", loc="left", fontsize=10)
    ax2.set_ylabel("z")
    ax2.set_xlabel("Time relative to glitch peak GPS (s)")
    ax2.grid(axis="x", ls=":", alpha=0.35)
    if corr_found:
        ax2.legend(fontsize=7, ncol=4, loc="upper right")

    for ax in axes:
        ax.set_xlim(win_lo, win_hi)

    miss = [c for c in SHORTLIST if c not in meta["found"]]
    fig.suptitle(
        f"ASC ERR/CORR shortlist around glitch peak {event_peak:.3f} | window {win_lo:+.3f}..{win_hi:+.3f} s",
        fontsize=12,
        fontweight="bold",
    )
    if miss:
        fig.text(
            0.09,
            0.02,
            "Missing channels in this raw frame: " + "; ".join(miss),
            fontsize=8,
            color="firebrick",
        )
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    event_peak = float(args.event_gps)
    event_id = int(event_peak)
    csv_out = OUTPUTS / f"asc_err_corr_shortlist_{event_id}.csv"
    meta_out = OUTPUTS / f"asc_err_corr_shortlist_{event_id}.json"
    png_400ms = OUT_DIR / f"asc_err_corr_shortlist_{event_id}_400ms.png"
    png_3p5s = OUT_DIR / f"asc_err_corr_shortlist_{event_id}_3p5s.png"

    meta = extract_raw(event_peak, args.run, csv_out, meta_out)
    df = pd.read_csv(csv_out)
    plot_window(df, meta, -0.2, 0.2, png_400ms, event_peak)
    plot_window(df, meta, -0.5, 3.0, png_3p5s, event_peak)
    print(f"Saved CSV:  {csv_out}")
    print(f"Saved meta: {meta_out}")
    print(f"Saved plot: {png_400ms}")
    print(f"Saved plot: {png_3p5s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
