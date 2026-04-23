#!/usr/bin/env python3
"""
Strict old-style glitch_response plot for O4c ASC CORR towers.

Reproduces the layout/logic of scripts/plot_glitch_response.py:
- Panel 1: V1:LSC_DARM_CORR bandpassed waveform (30-55 Hz)
- Panel 2: V1:LSC_DARM_CORR Hilbert envelope, baseline-normalized
- Panels 3-6: NI/WI/NE/WE ASC *_CORR Hilbert envelopes (TX/TY), baseline-normalized

Data are extracted directly from raw frames around one event.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import butter, hilbert, sosfiltfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt


WORKDIR = Path(__file__).resolve().parents[1]
OUTPUTS = WORKDIR / "outputs"
OUT_DIR = WORKDIR / "usecases" / "25-minute-glitch"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

SSH_HOST = "cca.in2p3.fr"
REMOTE_PYTHON = "/cvmfs/software.igwn.org/conda/envs/igwn/bin/python3"

DARM_COL = "V1:LSC_DARM_CORR"
ASC_TOWERS = [
    ("NI", [("V1:ASC_NI_TX_CORR", "tab:blue"), ("V1:ASC_NI_TY_CORR", "steelblue")]),
    ("WI", [("V1:ASC_WI_TX_CORR", "tab:orange"), ("V1:ASC_WI_TY_CORR", "goldenrod")]),
    ("NE", [("V1:ASC_NE_TX_CORR", "tab:green"), ("V1:ASC_NE_TY_CORR", "limegreen")]),
    ("WE", [("V1:ASC_WE_TX_CORR", "tab:red"), ("V1:ASC_WE_TY_CORR", "salmon")]),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Old-style glitch_response plot for O4c.")
    p.add_argument("--event-gps", type=float, required=True, help="Glitch peak GPS.")
    p.add_argument("--run", default="O4c", help="Run folder name under /Run (default O4c).")
    p.add_argument("--win", type=float, default=1.0, help="Half-window in seconds (default 1.0 => 2s total).")
    p.add_argument("--flo", type=float, default=30.0, help="Bandpass low cut (Hz).")
    p.add_argument("--fhi", type=float, default=55.0, help="Bandpass high cut (Hz).")
    p.add_argument("--env-smooth-ms", type=float, default=10.0, help="Hilbert envelope smoothing in ms.")
    p.add_argument("--baseline-start", type=float, default=-5.0, help="Baseline start relative to glitch (s).")
    p.add_argument("--baseline-end", type=float, default=-3.0, help="Baseline end relative to glitch (s).")
    p.add_argument("--refresh", action="store_true", help="Force raw re-extraction even if cache exists.")
    return p.parse_args()


def frame_start(gps: float) -> int:
    return int(gps // 100) * 100


def frame_starts_for_interval(t0: float, t1: float) -> list[int]:
    s0 = frame_start(t0)
    s1 = frame_start(t1 - 1e-9)
    return list(range(s0, s1 + 100, 100))


def extract_raw(
    event_gps: float,
    run_name: str,
    t0: float,
    t1: float,
    channels: list[str],
    csv_out: Path,
    meta_out: Path,
) -> dict:
    starts = frame_starts_for_interval(t0, t1)
    starts_shell = " ".join(str(s) for s in starts)
    channels_json = json.dumps(channels)
    script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_glitch_oldstyle.XXXXXX)"
cleanup() {{
  rm -rf "$stage_dir"
}}
trap cleanup EXIT
gwf_paths=()
for fs in {starts_shell}; do
  remote="cchpss0:/hpss/in2p3.fr/group/virgo/Run/{run_name}/raw/${{fs:0:4}}/V-raw-${{fs}}-100.gwf"
  local_gwf="$stage_dir/V-raw-${{fs}}-100.gwf"
  rfcp "$remote" "$local_gwf" >/dev/null
  gwf_paths+=("$local_gwf")
done
export GWF_PATHS="$(IFS=:; echo "${{gwf_paths[*]}}")"
export T0="{t0:.6f}"
export T1="{t1:.6f}"
export CHANNELS_JSON='{channels_json}'
"{REMOTE_PYTHON}" - <<'PY'
import json, os, sys
import pandas as pd
from gwpy.io.gwf import iter_channel_names
from gwpy.timeseries import TimeSeriesDict

gwf_paths = [p for p in os.environ["GWF_PATHS"].split(":") if p]
t0 = float(os.environ["T0"])
t1 = float(os.environ["T1"])
requested = json.loads(os.environ["CHANNELS_JSON"])
avail = set()
for gwf in gwf_paths:
    avail.update(iter_channel_names(gwf))
found = [c for c in requested if c in avail]
missing = [c for c in requested if c not in avail]
rows = []
if found:
    tsd = TimeSeriesDict.read(gwf_paths, found, start=t0, end=t1, nproc=1)
    for ch in found:
        ts = tsd[ch]
        rows.append(pd.DataFrame({{
            "channel": ch,
            "gps": ts.times.value,
            "value": ts.value.astype(float),
        }}))
meta = {{"requested": requested, "found": found, "missing": missing, "gwf_paths": gwf_paths}}
print("#META " + json.dumps(meta, sort_keys=True))
if rows:
    pd.concat(rows, ignore_index=True).to_csv(sys.stdout, index=False)
else:
    print("channel,gps,value")
PY
"""
    res = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", SSH_HOST, "bash", "-s", "--"],
        input=script,
        text=True,
        capture_output=True,
        check=False,
        timeout=1800,
    )
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip() or "raw extraction failed")
    if not res.stdout.startswith("#META "):
        raise RuntimeError("missing metadata header")
    meta_line, csv_body = res.stdout.split("\n", 1)
    meta = json.loads(meta_line[len("#META ") :])
    csv_out.write_text(csv_body)
    meta_out.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return meta


def bandpass(v: np.ndarray, fs: float, flo: float, fhi: float, order: int = 4) -> np.ndarray:
    if len(v) < 8 or fs <= 2.2 * fhi:
        return v
    sos = butter(order, [flo, fhi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def envelope(v: np.ndarray, fs: float, smooth_ms: float) -> np.ndarray:
    env = np.abs(hilbert(v))
    n = max(1, int(round(smooth_ms * 1e-3 * fs)))
    if n > 1:
        env = np.convolve(env, np.ones(n) / n, mode="same")
    return env


def get_fs_from_times(t: np.ndarray) -> float:
    if t.size < 3:
        return 1.0
    dt = np.median(np.diff(t[: min(200, t.size)]))
    if not np.isfinite(dt) or dt <= 0:
        return 1.0
    return float(1.0 / dt)


def load_series(df: pd.DataFrame, channel: str) -> tuple[np.ndarray, np.ndarray]:
    d = df[df["channel"] == channel].sort_values("gps")
    return d["gps"].to_numpy(dtype=float), d["value"].to_numpy(dtype=float)


def main() -> int:
    args = parse_args()
    event = float(args.event_gps)
    event_id = int(event)
    suffix = f"{int(args.win * 2)}s" if args.win >= 1 else f"{int(args.win * 2 * 1000)}ms"
    csv_out = OUTPUTS / f"glitch_response_oldstyle_{event_id}.csv"
    meta_out = OUTPUTS / f"glitch_response_oldstyle_{event_id}.json"
    out_png = OUT_DIR / f"glitch_response_{event_id}_{suffix}.png"

    t0 = event + min(args.baseline_start, -args.win) - 0.5
    t1 = event + max(args.baseline_end, args.win) + 0.5
    channels = [DARM_COL] + [c for _, group in ASC_TOWERS for c, _ in group]

    if args.refresh or (not csv_out.exists()) or (not meta_out.exists()):
        meta = extract_raw(event, args.run, t0, t1, channels, csv_out, meta_out)
    else:
        meta = json.loads(meta_out.read_text())

    df = pd.read_csv(csv_out).sort_values("gps")
    have = set(df["channel"])
    if DARM_COL not in have:
        raise RuntimeError(f"{DARM_COL} not found in extracted data")

    t_darm, v_darm = load_series(df, DARM_COL)
    t_rel_darm = t_darm - event
    fs_darm = get_fs_from_times(t_darm)
    darm_bp = bandpass(v_darm, fs_darm, args.flo, args.fhi)
    darm_env = envelope(darm_bp, fs_darm, args.env_smooth_ms)
    bl_mask = (t_rel_darm >= args.baseline_start) & (t_rel_darm <= args.baseline_end)
    if not np.any(bl_mask):
        bl_mask = t_rel_darm < 0
    if not np.any(bl_mask):
        bl_mask = np.ones_like(t_rel_darm, dtype=bool)
    darm_bl = np.mean(darm_env[bl_mask]) if np.any(bl_mask) else 1.0
    if not np.isfinite(darm_bl) or darm_bl <= 0:
        darm_bl = 1.0
    darm_env_norm = darm_env / darm_bl
    darm_rms = np.std(darm_bp[bl_mask]) if np.any(bl_mask) else np.std(darm_bp)
    if not np.isfinite(darm_rms) or darm_rms == 0:
        darm_rms = 1.0
    darm_bp_norm = darm_bp / darm_rms

    tower_traces: list[tuple[str, list[tuple[str, str, np.ndarray, np.ndarray]]]] = []
    for tower, chans in ASC_TOWERS:
        traces = []
        for ch, color in chans:
            if ch not in have:
                continue
            t, v = load_series(df, ch)
            t_rel = t - event
            fs = get_fs_from_times(t)
            bp = bandpass(v, fs, args.flo, args.fhi)
            env = envelope(bp, fs, args.env_smooth_ms)
            bl = (t_rel >= args.baseline_start) & (t_rel <= args.baseline_end)
            if not np.any(bl):
                bl = t_rel < 0
            if not np.any(bl):
                bl = np.ones_like(t_rel, dtype=bool)
            bl_mean = np.mean(env[bl]) if np.any(bl) else 1.0
            env_norm = env / bl_mean if bl_mean > 0 else env
            traces.append((ch, color, t_rel, env_norm))
        tower_traces.append((tower, traces))

    use_ms = args.win <= 0.5
    scale = 1e3 if use_ms else 1.0
    unit = "ms" if use_ms else "s"

    n_rows = 2 + len(tower_traces)
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, 3 * n_rows), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.06, left=0.13, right=0.97)

    m_d = (t_rel_darm >= -args.win) & (t_rel_darm <= args.win)
    t_plot = t_rel_darm[m_d] * scale

    ax1 = axes[0]
    ax1.plot(t_plot, darm_bp_norm[m_d], color="black", lw=0.8, alpha=0.9, label=f"{DARM_COL} ({args.flo:.0f}-{args.fhi:.0f} Hz)")
    ax1.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax1.set_ylabel("Amplitude\n(pre-glitch RMS)", fontsize=9)
    ax1.set_title(f"{DARM_COL} — bandpassed waveform", fontsize=9, loc="left")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(axis="x", ls=":", alpha=0.4)

    ax2 = axes[1]
    ax2.plot(t_plot, darm_env_norm[m_d], color="black", lw=1.5, alpha=0.9, label=f"{DARM_COL} envelope")
    ax2.axhline(1.0, color="black", lw=0.8, ls=":", alpha=0.5, label="baseline = 1")
    ax2.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
    ax2.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
    ax2.set_title(f"{DARM_COL} — Hilbert envelope", fontsize=9, loc="left")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(axis="x", ls=":", alpha=0.4)

    for i, (tower, traces) in enumerate(tower_traces):
        ax = axes[2 + i]
        for ch, color, t_rel, env_norm in traces:
            m = (t_rel >= -args.win) & (t_rel <= args.win)
            ax.plot(t_rel[m] * scale, env_norm[m], color=color, lw=1.2, alpha=0.9, label=ch)
        ax.axhline(1.0, color="gray", lw=0.8, ls=":", alpha=0.6)
        ax.axvline(0, color="crimson", lw=1.5, ls="--", alpha=0.8)
        ax.set_ylabel("Envelope\n(× pre-glitch)", fontsize=9)
        ax.set_title(f"{tower} — ASC CORR Hilbert envelope", fontsize=9, loc="left")
        if traces:
            ax.legend(fontsize=8, loc="upper right", handlelength=1.5)
        ax.grid(axis="x", ls=":", alpha=0.4)

    if use_ms:
        step = 50.0
    else:
        step = 0.2 if args.win <= 1.0 else 1.0
    xt = np.arange(-args.win * scale, args.win * scale + 1e-9, step)
    for ax in axes:
        ax.set_xticks(xt)
        if use_ms:
            ax.set_xticklabels([f"{v:+.0f}" for v in xt], fontsize=8)
        else:
            ax.set_xticklabels([f"{v:+.1f}" for v in xt], fontsize=8)
    axes[-1].set_xlabel(f"Time relative to catalog GPS {event:.6f} ({unit})", fontsize=9)

    utc = GPS_EPOCH + pd.to_timedelta(event, unit="s")
    fig.suptitle(
        "Glitch response style clone (old recipe)\n"
        f"GPS {event:.6f} ({utc.strftime('%Y-%m-%d %H:%M:%S UTC')})   "
        f"Window: ±{args.win*1e3:.0f} ms   Bandpass {args.flo:.0f}-{args.fhi:.0f} Hz   "
        f"Envelope smoothed {args.env_smooth_ms:.0f} ms",
        fontsize=10,
        fontweight="bold",
    )
    if meta.get("missing"):
        fig.text(
            0.13,
            0.012,
            "Missing channels: " + ", ".join(meta["missing"]),
            fontsize=8,
            color="dimgray",
        )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved CSV:  {csv_out}")
    print(f"Saved meta: {meta_out}")
    print(f"Saved -> {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
