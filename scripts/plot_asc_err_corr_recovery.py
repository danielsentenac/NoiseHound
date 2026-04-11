#!/usr/bin/env python3
"""
Plot ASC ERR/CORR recovery after a glitch over a long post-glitch window.

For one event:
- stage one raw frame
- extract ASC ERR/CORR shortlist + Hrec
- show long-window overlays
- estimate stabilization time from group spread returning to pre-glitch level
"""

from __future__ import annotations

import argparse
import json
import math
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
    p = argparse.ArgumentParser(description="Long-window ASC ERR/CORR recovery plot.")
    p.add_argument("--event-gps", type=float, required=True, help="Glitch peak GPS.")
    p.add_argument("--run", default="O4c", help="Run folder under /Run (default O4c).")
    p.add_argument("--pre", type=float, default=5.0, help="Seconds before event to extract.")
    p.add_argument("--post", type=float, default=30.0, help="Seconds after event to extract.")
    p.add_argument(
        "--hold",
        type=float,
        default=2.0,
        help="Seconds metric must remain near baseline to count as stabilized.",
    )
    p.add_argument(
        "--baseline-end",
        type=float,
        default=-1.0,
        help="End of pre-glitch baseline window (seconds, default -1.0).",
    )
    p.add_argument(
        "--fast-hp-hz",
        type=float,
        default=0.7,
        help="High-pass cutoff (Hz) for fast-component stabilization metric (default 0.7).",
    )
    p.add_argument(
        "--reuse-cached",
        action="store_true",
        help="Reuse existing CSV/meta outputs for this event if they exist.",
    )
    return p.parse_args()


def frame_start(gps: float) -> int:
    return int(gps // 100) * 100


def frame_starts_for_interval(t0: float, t1: float) -> list[int]:
    s0 = frame_start(t0)
    s1 = frame_start(t1 - 1e-9)
    return list(range(s0, s1 + 100, 100))


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


def robust_mad(x: np.ndarray) -> float:
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def stabilize_time(
    t: np.ndarray,
    metric: np.ndarray,
    baseline_mask: np.ndarray,
    hold_s: float,
    require_excursion: bool = False,
) -> tuple[float | None, float]:
    base = metric[baseline_mask]
    base_med = float(np.median(base))
    base_mad = float(robust_mad(base))
    threshold = base_med + max(3.0 * base_mad, 0.2)

    dt = float(np.median(np.diff(t)))
    n_hold = max(1, int(math.ceil(hold_s / dt)))
    valid = metric <= threshold
    start_idx = np.searchsorted(t, 0.0)
    search_start = start_idx
    if require_excursion:
        above_post = np.flatnonzero(metric[start_idx:] > threshold)
        if above_post.size == 0:
            return None, threshold
        search_start = start_idx + int(above_post[0]) + 1
    for i in range(search_start, len(t) - n_hold + 1):
        if np.all(valid[i : i + n_hold]):
            return float(t[i]), threshold
    return None, threshold


def extract_raw(event_peak: float, run_name: str, pre_s: float, post_s: float, csv_out: Path, meta_out: Path) -> dict:
    t0 = event_peak - pre_s
    t1 = event_peak + post_s
    starts = frame_starts_for_interval(t0, t1)
    starts_shell = " ".join(str(s) for s in starts)
    channels = [REF_CHANNEL] + SHORTLIST
    channels_json = json.dumps(channels)

    script = f"""
set -euo pipefail
stage_dir="$(mktemp -d /tmp/nh_asc_recovery.XXXXXX)"
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
        timeout=3600,
    )
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip() or "raw extraction failed")
    if not res.stdout.startswith("#META "):
        raise RuntimeError("missing metadata header from extraction")
    meta_line, csv_body = res.stdout.split("\n", 1)
    meta = json.loads(meta_line[len("#META ") :])
    csv_out.write_text(csv_body)
    meta_out.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return meta


def build_channel_matrix(df: pd.DataFrame, channels: list[str], t_common: np.ndarray) -> np.ndarray:
    mat = np.zeros((len(channels), len(t_common)), dtype=float)
    for i, ch in enumerate(channels):
        d = df[df["channel"] == ch].sort_values("t_rel")
        t = d["t_rel"].to_numpy()
        v = d["value"].to_numpy()
        mat[i, :] = np.interp(t_common, t, v)
    return mat


def z_normalize_matrix(mat: np.ndarray, baseline_mask: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat)
    for i in range(mat.shape[0]):
        v = mat[i]
        mu = float(np.mean(v[baseline_mask]))
        sig = float(np.std(v[baseline_mask]))
        if not np.isfinite(sig) or sig == 0:
            sig = 1.0
        out[i] = (v - mu) / sig
    return out


def bandpass_hrec(t: np.ndarray, v: np.ndarray) -> np.ndarray:
    if len(v) < 5:
        return v
    fs = 1.0 / np.median(np.diff(t))
    if fs <= 2 * 55:
        return v
    sos = butter(4, [30, 55], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, v)


def highpass_matrix(mat: np.ndarray, t: np.ndarray, cutoff_hz: float) -> np.ndarray:
    if mat.size == 0 or len(t) < 5 or cutoff_hz <= 0:
        return mat.copy()
    fs = 1.0 / np.median(np.diff(t))
    if not np.isfinite(fs) or fs <= 2.2 * cutoff_hz:
        return mat.copy()
    sos = butter(2, cutoff_hz, btype="highpass", fs=fs, output="sos")
    out = np.zeros_like(mat)
    for i in range(mat.shape[0]):
        out[i] = sosfiltfilt(sos, mat[i])
    return out


def plot_recovery(
    event_peak: float,
    pre_s: float,
    post_s: float,
    hold_s: float,
    baseline_end: float,
    fast_hp_hz: float,
    df_raw: pd.DataFrame,
    meta: dict,
    png_out: Path,
    summary_out: Path,
) -> None:
    # common grid from first ERR channel (all shortlist channels are 2000 Hz here)
    ref_grid_ch = next(ch for ch in ERR_CH if ch in meta["found"])
    grid = df_raw[df_raw["channel"] == ref_grid_ch].sort_values("t_rel")
    t = grid["t_rel"].to_numpy()

    baseline_mask = (t >= -pre_s) & (t <= baseline_end)
    if not np.any(baseline_mask):
        raise RuntimeError("empty baseline mask; adjust --pre/--baseline-end")

    err_found = [ch for ch in ERR_CH if ch in meta["found"]]
    corr_found = [ch for ch in CORR_CH if ch in meta["found"]]

    err_mat = build_channel_matrix(df_raw, err_found, t)
    corr_mat = build_channel_matrix(df_raw, corr_found, t)
    err_z = z_normalize_matrix(err_mat, baseline_mask)
    corr_z = z_normalize_matrix(corr_mat, baseline_mask)

    err_spread = np.std(err_z, axis=0)
    corr_spread = np.std(corr_z, axis=0)
    err_t_stab, err_thr = stabilize_time(t, err_spread, baseline_mask, hold_s)
    corr_t_stab, corr_thr = stabilize_time(t, corr_spread, baseline_mask, hold_s)

    err_fast = highpass_matrix(err_mat, t, fast_hp_hz)
    corr_fast = highpass_matrix(corr_mat, t, fast_hp_hz)
    err_fast_z = z_normalize_matrix(err_fast, baseline_mask)
    corr_fast_z = z_normalize_matrix(corr_fast, baseline_mask)
    err_fast_spread = np.std(err_fast_z, axis=0)
    corr_fast_spread = np.std(corr_fast_z, axis=0)
    err_fast_t_stab, err_fast_thr = stabilize_time(
        t, err_fast_spread, baseline_mask, hold_s, require_excursion=True
    )
    corr_fast_t_stab, corr_fast_thr = stabilize_time(
        t, corr_fast_spread, baseline_mask, hold_s, require_excursion=True
    )

    # Hrec panel on native high-rate data
    hrec = df_raw[df_raw["channel"] == REF_CHANNEL].sort_values("t_rel")
    th = hrec["t_rel"].to_numpy()
    vh = hrec["value"].to_numpy()
    vh_bp = bandpass_hrec(th, vh)
    h_base = (th >= -pre_s) & (th <= baseline_end)
    h_mu = float(np.mean(vh_bp[h_base]))
    h_sig = float(np.std(vh_bp[h_base]))
    if not np.isfinite(h_sig) or h_sig == 0:
        h_sig = 1.0
    h_z = (vh_bp - h_mu) / h_sig

    fig, axes = plt.subplots(4, 1, figsize=(17, 12), sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.93, bottom=0.08, left=0.08, right=0.98)

    ax0, ax1, ax2, ax3 = axes
    m_h = (th >= -pre_s) & (th <= post_s)
    ax0.plot(th[m_h], h_z[m_h], color="black", lw=0.8)
    ax0.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax0.set_ylabel("Hrec z")
    ax0.set_title("V1:Hrec_hoft_16384Hz (30-55 Hz bandpass)", loc="left", fontsize=10)
    ax0.grid(axis="x", ls=":", alpha=0.35)

    for i, ch in enumerate(err_found):
        ax1.plot(t, err_z[i], lw=1.0, alpha=0.85, color=color_for(ch), label=short_label(ch))
    ax1.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax1.set_ylabel("z")
    ax1.set_title("ASC ERR channels", loc="left", fontsize=10)
    ax1.grid(axis="x", ls=":", alpha=0.35)
    ax1.legend(fontsize=8, ncol=3, loc="upper right")

    for i, ch in enumerate(corr_found):
        ax2.plot(t, corr_z[i], lw=0.95, alpha=0.82, color=color_for(ch), label=short_label(ch))
    ax2.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax2.set_ylabel("z")
    ax2.set_title("ASC CORR channels", loc="left", fontsize=10)
    ax2.grid(axis="x", ls=":", alpha=0.35)
    ax2.legend(fontsize=7, ncol=4, loc="upper right")

    ax3.plot(t, err_spread, color="#1f77b4", lw=1.8, label=f"ERR spread (std), thr={err_thr:.2f}")
    ax3.plot(t, corr_spread, color="#2ca02c", lw=1.8, label=f"CORR spread (std), thr={corr_thr:.2f}")
    ax3.plot(
        t,
        err_fast_spread,
        color="#1f77b4",
        lw=1.4,
        ls="--",
        alpha=0.9,
        label=f"ERR fast spread (hp {fast_hp_hz:.2g}Hz), thr={err_fast_thr:.2f}",
    )
    ax3.plot(
        t,
        corr_fast_spread,
        color="#2ca02c",
        lw=1.4,
        ls="--",
        alpha=0.9,
        label=f"CORR fast spread (hp {fast_hp_hz:.2g}Hz), thr={corr_fast_thr:.2f}",
    )
    ax3.axvline(0, color="crimson", lw=1.4, ls="--", alpha=0.8)
    ax3.axhline(err_thr, color="#1f77b4", lw=1.0, ls=":", alpha=0.7)
    ax3.axhline(corr_thr, color="#2ca02c", lw=1.0, ls=":", alpha=0.7)
    ax3.axhline(err_fast_thr, color="#1f77b4", lw=0.9, ls="--", alpha=0.35)
    ax3.axhline(corr_fast_thr, color="#2ca02c", lw=0.9, ls="--", alpha=0.35)
    if err_t_stab is not None:
        ax3.axvline(err_t_stab, color="#1f77b4", lw=1.3, ls="--", alpha=0.9)
    if corr_t_stab is not None:
        ax3.axvline(corr_t_stab, color="#2ca02c", lw=1.3, ls="--", alpha=0.9)
    if err_fast_t_stab is not None:
        ax3.axvline(err_fast_t_stab, color="#1f77b4", lw=1.1, ls="-.", alpha=0.9)
    if corr_fast_t_stab is not None:
        ax3.axvline(corr_fast_t_stab, color="#2ca02c", lw=1.1, ls="-.", alpha=0.9)
    ax3.set_ylabel("spread")
    ax3.set_xlabel("Time relative to glitch peak GPS (s)")
    ax3.set_title(
        f"Group spread and stabilization (hold={hold_s:.1f}s, baseline=[{-pre_s:.1f},{baseline_end:.1f}]s)",
        loc="left",
        fontsize=10,
    )
    ax3.grid(axis="x", ls=":", alpha=0.35)
    ax3.legend(fontsize=9, loc="upper right")

    for ax in axes:
        ax.set_xlim(-pre_s, post_s)
    axes[-1].set_xticks(np.arange(-pre_s, post_s + 1e-9, 2.5))

    err_text = f"{err_t_stab:+.2f}s" if err_t_stab is not None else f"no return by +{post_s:.1f}s"
    corr_text = f"{corr_t_stab:+.2f}s" if corr_t_stab is not None else f"no return by +{post_s:.1f}s"
    err_fast_text = (
        f"{err_fast_t_stab:+.2f}s" if err_fast_t_stab is not None else f"no return by +{post_s:.1f}s"
    )
    corr_fast_text = (
        f"{corr_fast_t_stab:+.2f}s" if corr_fast_t_stab is not None else f"no return by +{post_s:.1f}s"
    )
    fig.suptitle(
        f"ASC ERR/CORR recovery around glitch peak {event_peak:.6f} | window {-pre_s:+.1f}..{post_s:+.1f}s\n"
        f"Stabilization (std spread): ERR={err_text}, CORR={corr_text} | "
        f"fast hp {fast_hp_hz:.2g}Hz: ERR={err_fast_text}, CORR={corr_fast_text}",
        fontsize=12,
        fontweight="bold",
    )

    png_out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_out, dpi=170, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "event_peak_gps": event_peak,
        "window_pre_s": pre_s,
        "window_post_s": post_s,
        "hold_s": hold_s,
        "baseline_end_s": baseline_end,
        "fast_hp_hz": fast_hp_hz,
        "err_stabilization_s": err_t_stab,
        "corr_stabilization_s": corr_t_stab,
        "err_threshold": err_thr,
        "corr_threshold": corr_thr,
        "err_fast_stabilization_s": err_fast_t_stab,
        "corr_fast_stabilization_s": corr_fast_t_stab,
        "err_fast_threshold": err_fast_thr,
        "corr_fast_threshold": corr_fast_thr,
        "err_found": err_found,
        "corr_found": corr_found,
        "missing": [c for c in SHORTLIST if c not in meta["found"]],
    }
    summary_out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main() -> int:
    args = parse_args()
    event_peak = float(args.event_gps)
    event_id = int(event_peak)

    csv_out = OUTPUTS / f"asc_err_corr_recovery_{event_id}.csv"
    meta_out = OUTPUTS / f"asc_err_corr_recovery_{event_id}.json"
    png_out = OUT_DIR / f"asc_err_corr_recovery_{event_id}_{int(args.pre + args.post)}s.png"
    summary_out = OUTPUTS / f"asc_err_corr_recovery_{event_id}_summary.json"

    if args.reuse_cached and csv_out.exists() and meta_out.exists():
        meta = json.loads(meta_out.read_text())
    else:
        meta = extract_raw(event_peak, args.run, args.pre, args.post, csv_out, meta_out)
    df = pd.read_csv(csv_out)
    df["t_rel"] = df["gps"] - event_peak
    if args.reuse_cached:
        t_min = float(df["t_rel"].min())
        t_max = float(df["t_rel"].max())
        if t_min > -args.pre + 0.2 or t_max < args.post - 0.2:
            raise RuntimeError(
                f"cached CSV does not cover requested window [-{args.pre}, +{args.post}]s "
                f"(found [{t_min:.3f}, {t_max:.3f}]s); rerun without --reuse-cached"
            )
    plot_recovery(
        event_peak=event_peak,
        pre_s=args.pre,
        post_s=args.post,
        hold_s=args.hold,
        baseline_end=args.baseline_end,
        fast_hp_hz=args.fast_hp_hz,
        df_raw=df,
        meta=meta,
        png_out=png_out,
        summary_out=summary_out,
    )

    print(f"Saved CSV:     {csv_out}")
    print(f"Saved meta:    {meta_out}")
    print(f"Saved summary: {summary_out}")
    print(f"Saved plot:    {png_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
