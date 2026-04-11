#!/usr/bin/env bash
#SBATCH -J nh-snr-corr
#SBATCH -p htc
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH -o /pbs/home/s/sentenac/NOISEHOUND/slurm-nh-snr-corr-%j.out

# Correlate ALL physical _mean channels from trend GWF files against the
# smoothed glitch SNR trend over full O4.
#
# Strategy:
#   - Sample O4 at one GWF per 2 weeks (~65 files)
#   - Read all physical _mean channels from each file, compute daily mean
#   - Build channel matrix [N_days x N_channels]
#   - Compute Spearman rank correlation vs SNR trend with lag search ±28 days
#   - Output ranked channel list (top 1000 by |rho|)
#
# Output: outputs/snr_channel_correlations.csv

set -euo pipefail

WORKDIR="$HOME/NOISEHOUND"
IGWN_PYTHON="${IGWN_PYTHON:-/cvmfs/software.igwn.org/conda/envs/igwn/bin/python}"
export STAGE_DIR="/tmp/nh_snr_corr_${SLURM_JOB_ID:-local}"
mkdir -p "$STAGE_DIR"

echo "=== SNR channel correlation started $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="

export PYTHONPATH="$WORKDIR${PYTHONPATH:+:$PYTHONPATH}"

"$IGWN_PYTHON" - << 'PYEOF'
import sys, warnings, subprocess, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from gwpy.timeseries import TimeSeriesDict
from gwpy.io.gwf import iter_channel_names

warnings.filterwarnings("ignore")

WORKDIR   = Path.home() / "NOISEHOUND"
STAGE_DIR = Path(os.environ["STAGE_DIR"])
HPSS_BASE = "cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend"
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")

# O4 GPS range, sample every 14 days
GPS_O4_START = 1364774400
GPS_O4_END   = 1459468800
STEP_S       = 14 * 86400
BATCH_SIZE   = 500   # channels per gwpy read call
WINDOW_H     = 2     # hours of midday data to read per GWF (instead of full day)

# ── Load SNR trend ───────────────────────────────────────────────────────────
snr_df = pd.read_csv(WORKDIR / "outputs" / "snr_trend_daily.csv")
snr_df = snr_df.dropna().sort_values("gps_bin").reset_index(drop=True)
snr_gps  = snr_df["gps_bin"].values.astype(int)
snr_vals = snr_df["snr_trend"].values
print(f"SNR trend: {len(snr_gps)} daily bins, GPS {snr_gps[0]}–{snr_gps[-1]}", flush=True)

# ── Get channel list from a staged GWF ──────────────────────────────────────
SKIP_PREFIXES = {"DQ_", "Qc_", "Sc_", "Nds", "nds"}
# Suffixes that indicate non-analog channels (flags, counters, etc.)
SKIP_SUFFIXES = {"_min", "_max", "_rms",
                 "_ST", "_FLAG", "_ALRM", "_ALARM", "_LOCK",
                 "_ENBL", "_ON", "_OFF", "_RBT", "_COM_ST"}
# Name fragments that indicate non-physical integer/counter channels
SKIP_FRAGMENTS = {"cntLoad", "nApply", "nLoadApply", "nLoadOnly", "nWkLoad",
                  "restart", "runTime", "setTime", "lTmo", "NUMITERS",
                  "CONTIG_RT", "CTRL_MAN_SIM", "CTRL_SOFT_SIM",
                  "CfgChange", "elapsed_time"}

def get_channel_list(gwf_path: str) -> list[str]:
    all_ch = list(iter_channel_names(gwf_path))
    result = []
    for c in all_ch:
        if any(c.startswith(s) for s in SKIP_PREFIXES):
            continue
        if any(c.endswith(s) for s in SKIP_SUFFIXES):
            continue
        if any(f in c for f in SKIP_FRAGMENTS):
            continue
        result.append(c)
    return result

def stage(gps_day: int) -> Path | None:
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

# ── Sample GPS days across O4 ────────────────────────────────────────────────
sample_days = list(range(GPS_O4_START, GPS_O4_END, STEP_S))
print(f"Sampling {len(sample_days)} GWF files across O4", flush=True)

# First file: get channel list
channels = None
for gps_day in sample_days:
    gwf = stage(gps_day)
    if gwf is None:
        continue
    channels = get_channel_list(str(gwf))
    print(f"Channel list: {len(channels)} physical _mean channels", flush=True)
    gwf.unlink(missing_ok=True)
    break

if channels is None:
    print("Could not get channel list!", flush=True)
    sys.exit(1)

# ── Build daily-mean matrix ──────────────────────────────────────────────────
# data_matrix[i, j] = daily mean of channel j on day i
n_days = len(sample_days)
n_ch   = len(channels)
data_matrix = np.full((n_days, n_ch), np.nan, dtype=np.float32)
day_gps     = np.array(sample_days, dtype=int)

ch_index = {c: i for i, c in enumerate(channels)}

for day_i, gps_day in enumerate(sample_days):
    gwf = stage(gps_day)
    if gwf is None:
        continue
    # Read a 2h midday window only (much faster than full day)
    t0 = gps_day + 11 * 3600   # 11h UTC
    t1 = gps_day + 13 * 3600   # 13h UTC
    n_ok = 0
    # Read in batches
    for b_start in range(0, n_ch, BATCH_SIZE):
        batch = channels[b_start:b_start + BATCH_SIZE]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tsd = TimeSeriesDict.read(str(gwf), batch, start=t0, end=t1,
                                          nproc=1)
            for ch, ts in tsd.items():
                v = np.nanmean(np.asarray(ts.value, float))
                if np.isfinite(v):
                    data_matrix[day_i, ch_index[ch]] = v
                    n_ok += 1
        except Exception as e:
            pass  # silent: individual channels may be absent
    gwf.unlink(missing_ok=True)
    pct = 100 * n_ok / n_ch
    print(f"  day {gps_day} ({day_i+1}/{n_days}): {n_ok}/{n_ch} channels ({pct:.0f}%)", flush=True)

print(f"Matrix shape: {data_matrix.shape}", flush=True)

# ── Align SNR trend to sampled days ─────────────────────────────────────────
# Interpolate SNR trend to sample_days (many sample days won't have glitches =
# NaN in snr_trend; keep only days where both exist)
snr_interp = np.interp(day_gps, snr_gps, snr_vals,
                        left=np.nan, right=np.nan)

# ── Spearman correlation with lag search ─────────────────────────────────────
MAX_LAG_DAYS = 4   # in units of STEP_S (14d steps) → ±56 days
lags = list(range(-MAX_LAG_DAYS, MAX_LAG_DAYS + 1))  # in sample steps

results = []
valid_mask_snr = np.isfinite(snr_interp)

print(f"Computing Spearman correlations for {n_ch} channels, {len(lags)} lags...", flush=True)

for j in range(n_ch):
    col = data_matrix[:, j]
    best_rho = 0.0
    best_lag = 0
    best_p   = 1.0
    for lag in lags:
        if lag >= 0:
            x = col[:n_days - lag] if lag > 0 else col
            y = snr_interp[lag:] if lag > 0 else snr_interp
        else:
            x = col[-lag:]
            y = snr_interp[:n_days + lag]
        mask = valid_mask_snr[lag:] if lag >= 0 else valid_mask_snr[:n_days + lag]
        mask = mask & np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 10:
            continue
        try:
            rho, p = spearmanr(x[mask], y[mask])
        except Exception:
            continue
        if abs(rho) > abs(best_rho):
            best_rho = rho
            best_lag = lag * STEP_S // 86400  # convert to days
            best_p   = p
    if abs(best_rho) > 0.1:
        results.append({
            "channel":  channels[j],
            "rho":      round(best_rho, 4),
            "lag_days": best_lag,
            "p_value":  round(best_p, 6),
        })

results.sort(key=lambda r: abs(r["rho"]), reverse=True)
print(f"Channels with |rho| > 0.1: {len(results)}", flush=True)
print("Top 20:", flush=True)
for r in results[:20]:
    print(f"  {r['rho']:+.3f}  lag={r['lag_days']:+3d}d  {r['channel']}", flush=True)

out = WORKDIR / "outputs" / "snr_channel_correlations.csv"
pd.DataFrame(results).to_csv(out, index=False)
print(f"\nSaved {len(results)} entries -> {out}", flush=True)
PYEOF

echo "=== Done $(date -u '+%Y-%m-%d %H:%M:%SZ') ==="
rm -rf "$STAGE_DIR"
