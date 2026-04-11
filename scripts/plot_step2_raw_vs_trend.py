#!/usr/bin/env python3
"""
Compare raw and trend representations of Step-2 top channels for one event.

Uses:
- local raw extraction already produced by plot_step2_top_raw_delay.py
- fresh trend extraction from the Jan 1, 2025 trend day file on cca.in2p3.fr
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"
OUTPUTS = WORKDIR / "outputs"

EVENT_PEAK = 1419726765.519531
EVENT_START = 1419726765.515625
WINDOW_BEFORE = 0.5
WINDOW_AFTER = 3.0

RAW_CSV = OUTPUTS / "step2_top_raw_delay_1419726765.csv"
TREND_CSV = OUTPUTS / "step2_top5_trend_window_1419726765.csv"
TREND_META = OUTPUTS / "step2_top5_trend_window_1419726765.json"
OUT_PNG = OUT_DIR / "step2_raw_vs_trend_1419726765_3p5s.png"

TREND_REMOTE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/2025/V-trend-1419724800-86400.gwf"
SSH_HOST = "cca.in2p3.fr"
REMOTE_PYTHON = "/cvmfs/software.igwn.org/conda/envs/igwn/bin/python3"

TOP5_TREND = [
    "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean",
    "V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean",
    "V1:LSC_DARM_IMC_LINE_Q_100Hz_mean",
    "V1:LSC_DARM_PSTAB0_I_FS_mean",
    "V1:LSC_DARM_PSTAB0_I_100Hz_mean",
]

RAW_EQUIV = {
    "V1:LSC_DARM_IMC_LINE_mag_100Hz_mean": "V1:LSC_DARM_IMC_LINE_mag_100Hz",
    "V1:LSC_DARM_PSTAB0_COUPLING_100Hz_mean": "V1:LSC_DARM_PSTAB0_COUPLING_100Hz",
    "V1:LSC_DARM_IMC_LINE_Q_100Hz_mean": "V1:LSC_DARM_IMC_LINE_Q_100Hz",
    "V1:LSC_DARM_PSTAB0_I_FS_mean": "V1:LSC_DARM_PSTAB0_I_FS",
    "V1:LSC_DARM_PSTAB0_I_100Hz_mean": "V1:LSC_DARM_PSTAB0_I_100Hz",
}


def znorm(t: np.ndarray, v: np.ndarray) -> np.ndarray:
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


def extract_trend_if_needed() -> None:
    if TREND_CSV.exists() and TREND_META.exists():
        return

    t0 = EVENT_PEAK - 5.0
    t1 = EVENT_PEAK + 5.0
    channels_json = json.dumps(TOP5_TREND)
    remote_script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_step2_trend_cmp.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT
local_gwf="$stage_dir/V-trend-1419724800-86400.gwf"
rfcp {TREND_REMOTE} "$local_gwf" >/dev/null
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
        input=remote_script,
        text=True,
        capture_output=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip() or "trend extraction failed")
    if not res.stdout.startswith("#META "):
        raise RuntimeError("trend extraction missing metadata header")
    meta_line, csv_body = res.stdout.split("\n", 1)
    meta = json.loads(meta_line[len("#META ") :])
    TREND_CSV.write_text(csv_body)
    TREND_META.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def main() -> int:
    if not RAW_CSV.exists():
        raise RuntimeError(f"raw CSV not found: {RAW_CSV}")
    extract_trend_if_needed()

    raw = pd.read_csv(RAW_CSV)
    trend = pd.read_csv(TREND_CSV)

    fig, axes = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.90, bottom=0.10, left=0.09, right=0.98)

    ax_raw, ax_trend = axes

    for ch_t in TOP5_TREND:
        ch_r = RAW_EQUIV[ch_t]
        d = raw[raw["channel"] == ch_r].sort_values("gps")
        if d.empty:
            continue
        t = d["gps"].to_numpy() - EVENT_PEAK
        z = znorm(t, d["value"].to_numpy(dtype=float))
        peak_t = t[int(np.argmax(z))]
        ax_raw.plot(t, z, lw=1.2, label=f"{ch_r}  (peak {peak_t:+.3f}s)")

    for ch_t in TOP5_TREND:
        d = trend[trend["channel"] == ch_t].sort_values("gps")
        if d.empty:
            continue
        t = d["gps"].to_numpy() - EVENT_PEAK
        z = znorm(t, d["value"].to_numpy(dtype=float))
        peak_t = t[int(np.argmax(z))]
        ax_trend.plot(t, z, marker="o", ms=4, lw=1.0, label=f"{ch_t}  (peak {peak_t:+.3f}s)")

    for ax in axes:
        ax.axvline(0.0, color="crimson", lw=1.2, ls="--", alpha=0.8)
        ax.axvline(EVENT_START - EVENT_PEAK, color="dimgray", lw=1.0, ls=":", alpha=0.8)
        ax.grid(axis="x", ls=":", alpha=0.35)
        ax.set_xlim(-WINDOW_BEFORE, WINDOW_AFTER)
        ax.set_xticks(np.arange(-0.5, 3.01, 0.5))
        ax.set_ylabel("z")

    ax_raw.set_title("Raw counterparts (only channels existing in raw frame)", loc="left", fontsize=10)
    ax_trend.set_title("Trend channels from Jan 1 day file (1 Hz products)", loc="left", fontsize=10)
    ax_raw.legend(fontsize=8, loc="upper right")
    ax_trend.legend(fontsize=8, loc="upper right")
    ax_trend.set_xlabel("Time relative to glitch peak GPS (s)")

    fig.suptitle(
        "Step 2 channel timing: raw vs trend on the same event\n"
        f"peak={EVENT_PEAK:.6f}, start={EVENT_START:.6f}, window=-0.5..+3.0 s",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.09,
        0.02,
        "Crimson dashed=peak. Gray dotted=tstart. Raw frame lacks verbatim IMC_LINE_Q/PSTAB0_I channels.",
        fontsize=8,
    )
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT_PNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
