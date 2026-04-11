#!/usr/bin/env bash
#SBATCH -J nh-seg-scan
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-seg-scan-%j.out

# Scan ALL ~80k trend channels over UP/DOWN/UP segments to find channels
# that change consistently with the glitch-SNR dip.
#
# Two dips identified:
#   Dip 1 (Nov/Dec 2024):
#     sbatch --export=DIP=1,\
#       GPS_UP1_S=1411689600,GPS_UP1_E=1415318400,\
#       GPS_DN_S=1415318400,GPS_DN_E=1416528000,\
#       GPS_UP2_S=1416528000,GPS_UP2_E=1419552000 \
#       slurm/nh_segment_scan.slurm
#
#   Dip 2 (Apr/May 2025):
#     sbatch --export=DIP=2,\
#       GPS_UP1_S=1424995200,GPS_UP1_E=1429228800,\
#       GPS_DN_S=1429228800,GPS_DN_E=1431043200,\
#       GPS_UP2_S=1431043200,GPS_UP2_E=1435276800 \
#       slurm/nh_segment_scan.slurm
#
# Output: outputs/segment_scan_dip{DIP}.csv
#   columns: channel, mu_up1, mu_down, mu_up2, delta, sigma_up, z_score

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"

DIP="${DIP:-1}"
GPS_UP1_S="${GPS_UP1_S:-1411689600}"
GPS_UP1_E="${GPS_UP1_E:-1415318400}"
GPS_DN_S="${GPS_DN_S:-1415318400}"
GPS_DN_E="${GPS_DN_E:-1416528000}"
GPS_UP2_S="${GPS_UP2_S:-1416528000}"
GPS_UP2_E="${GPS_UP2_E:-1419552000}"

export STAGE_DIR="/tmp/nh_seg_scan_dip${DIP}_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== Segment scan DIP=${DIP} started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
echo "    UP1:  $GPS_UP1_S – $GPS_UP1_E"
echo "    DOWN: $GPS_DN_S  – $GPS_DN_E"
echo "    UP2:  $GPS_UP2_S – $GPS_UP2_E"

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << PYEOF
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from gwpy.timeseries import TimeSeriesDict
from gwpy.io.gwf import iter_channel_names

warnings.filterwarnings("ignore")

WORKDIR   = Path.home() / "NOISEHOUND"
STAGE_DIR = Path(os.environ["STAGE_DIR"])
STAGE_DIR.mkdir(exist_ok=True)

HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
DIP       = int(os.environ["DIP"])

GPS_UP1_S = int(os.environ["GPS_UP1_S"])
GPS_UP1_E = int(os.environ["GPS_UP1_E"])
GPS_DN_S  = int(os.environ["GPS_DN_S"])
GPS_DN_E  = int(os.environ["GPS_DN_E"])
GPS_UP2_S = int(os.environ["GPS_UP2_S"])
GPS_UP2_E = int(os.environ["GPS_UP2_E"])

BATCH_SIZE = 500
WINDOW_H   = 2     # hours of midday data to read per GWF
N_SAMPLE   = 7     # target number of days to sample per segment


def sample_days(gps_start, gps_end, n):
    """Sample n evenly-spaced GPS day starts within a segment."""
    days = list(range((gps_start // 86400) * 86400,
                      (gps_end   // 86400) * 86400,
                      86400))
    if len(days) <= n:
        return days
    idx = np.round(np.linspace(0, len(days) - 1, n)).astype(int)
    return [days[i] for i in idx]


def stage(gps_day):
    year   = (GPS_EPOCH + pd.to_timedelta(gps_day, unit="s")).year
    remote = f"{HPSS_BASE}/{year}/V-trend-{gps_day}-86400.gwf"
    local  = STAGE_DIR / f"V-trend-{gps_day}-86400.gwf"
    if not local.exists():
        r = subprocess.run(["rfcp", remote, str(local)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  rfcp failed {gps_day}: {r.stderr[:80]}", flush=True)
            return None
    return local


# ── Segment definitions ───────────────────────────────────────────────────────
segments = {
    "up1":  (GPS_UP1_S, GPS_UP1_E),
    "down": (GPS_DN_S,  GPS_DN_E),
    "up2":  (GPS_UP2_S, GPS_UP2_E),
}

# Denser sampling for the short DOWN segment
sample_counts = {"up1": N_SAMPLE, "down": max(4, N_SAMPLE // 2), "up2": N_SAMPLE}

sampled = {seg: sample_days(s, e, sample_counts[seg])
           for seg, (s, e) in segments.items()}

print("Sampled days per segment:", {k: len(v) for k, v in sampled.items()}, flush=True)

# ── Get channel list from first available file ────────────────────────────────
channels = None
for gps_day in sampled["up1"] + sampled["down"] + sampled["up2"]:
    gwf = stage(gps_day)
    if gwf is None:
        continue
    channels = list(iter_channel_names(str(gwf)))
    print(f"Channel list: {len(channels)} channels from day {gps_day}", flush=True)
    gwf.unlink(missing_ok=True)
    break

if channels is None:
    print("Could not obtain channel list!", flush=True)
    sys.exit(1)

n_ch = len(channels)
ch_index = {c: i for i, c in enumerate(channels)}

# ── Read data: data_mat[seg][day_i, ch_j] = daily mean ───────────────────────
# Store as dict of lists of arrays
seg_daily = {seg: [] for seg in segments}

all_days = [(seg, d) for seg, days in sampled.items() for d in days]

for seg, gps_day in all_days:
    gwf = stage(gps_day)
    if gwf is None:
        continue
    t0 = gps_day + 11 * 3600   # 11h UTC midday start
    t1 = gps_day + (11 + WINDOW_H) * 3600
    row = np.full(n_ch, np.nan, dtype=np.float32)
    n_ok = 0
    for b0 in range(0, n_ch, BATCH_SIZE):
        batch = channels[b0:b0 + BATCH_SIZE]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tsd = TimeSeriesDict.read(str(gwf), batch, start=t0, end=t1,
                                          nproc=1)
            for ch, ts in tsd.items():
                v = np.nanmean(np.asarray(ts.value, float))
                if np.isfinite(v):
                    row[ch_index[ch]] = v
                    n_ok += 1
        except Exception:
            pass
    seg_daily[seg].append(row)
    gwf.unlink(missing_ok=True)
    pct = 100 * n_ok / n_ch
    print(f"  [{seg}] day {gps_day}: {n_ok}/{n_ch} channels ({pct:.0f}%)", flush=True)

# ── Compute per-segment medians ───────────────────────────────────────────────
def seg_median(daily_list):
    if not daily_list:
        return np.full(n_ch, np.nan)
    return np.nanmedian(np.stack(daily_list, axis=0), axis=0)

mu_up1 = seg_median(seg_daily["up1"])
mu_dn  = seg_median(seg_daily["down"])
mu_up2 = seg_median(seg_daily["up2"])

# ── Scoring ───────────────────────────────────────────────────────────────────
# Baseline: average of both UP segments
mu_up_avg = (mu_up1 + mu_up2) / 2.0

# Delta: how much the channel changes during the DOWN segment
delta = mu_dn - mu_up_avg

# Noise: std of the UP daily samples (natural variability when SNR is "normal")
all_up = seg_daily["up1"] + seg_daily["up2"]
if all_up:
    up_mat = np.stack(all_up, axis=0)   # [n_up_days, n_ch]
    sigma_up = np.nanstd(up_mat, axis=0, ddof=1)
else:
    sigma_up = np.ones(n_ch)

# Z-score: normalised change. |z| > 2 → significant change during dip.
# Positive z = channel goes UP when SNR goes DOWN.
# Negative z = channel goes DOWN when SNR goes DOWN (anti-correlated with dip).
epsilon = np.nanmedian(np.abs(mu_up_avg[np.isfinite(mu_up_avg)])) * 1e-4
epsilon = max(epsilon, 1e-10)
z_score = delta / (sigma_up + epsilon)

# ── Build output ──────────────────────────────────────────────────────────────
valid = (np.isfinite(mu_up1) & np.isfinite(mu_dn) &
         np.isfinite(mu_up2) & np.isfinite(z_score))

results = []
for j in np.where(valid)[0]:
    if abs(z_score[j]) > 0.5:   # loose threshold — keep borderline cases
        results.append({
            "channel":  channels[j],
            "mu_up1":   round(float(mu_up1[j]),  6),
            "mu_down":  round(float(mu_dn[j]),   6),
            "mu_up2":   round(float(mu_up2[j]),  6),
            "delta":    round(float(delta[j]),    6),
            "sigma_up": round(float(sigma_up[j]),6),
            "z_score":  round(float(z_score[j]), 4),
        })

results.sort(key=lambda r: abs(r["z_score"]), reverse=True)
print(f"\nChannels with |z|>0.5: {len(results)}", flush=True)
print("Top 30:", flush=True)
for r in results[:30]:
    print(f"  z={r['z_score']:+7.2f}  delta={r['delta']:+.4g}  {r['channel']}", flush=True)

out = WORKDIR / "outputs" / f"segment_scan_dip{DIP}.csv"
out.parent.mkdir(exist_ok=True)
pd.DataFrame(results).to_csv(out, index=False)
print(f"\nSaved {len(results)} entries -> {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
