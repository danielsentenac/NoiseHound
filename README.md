# NOISEHOUND

A toolkit for investigating long-lived periodic broadband glitches typically seen in Virgo. It replaces one-by-one coherence browsing with a repeatable multi-step workflow based on 1 Hz trend GWF files staged from the CC-IN2P3 HPSS tape archive.

## Workflow overview

```
Step 1  detect      Identify glitch times in Hrec_hoft (short epoch, raw or HoftOnline frames)
Step 2  rank        Rank auxiliary channels by event-synchronous response (1 Hz trend GWFs)
Step 3  causality   Granger / Transfer Entropy on top-ranked channels (same epoch)
Step 4  correlate   Pearson/Spearman rate correlation over full run baseline (1h bins, trend GWFs)
Step 5  refine      Lag refinement with finer bins (15 min) on priority channels
Step 6  ccm         Convergent Cross Mapping if causality remains inconclusive
Step 7  control     Repeat Steps 4–5 on a glitch-free epoch for negative control
```

Steps 1–3 work on a short epoch (hours). Steps 4–7 use the full run baseline (months to years) via SLURM array jobs on CC-IN2P3.

## Data source

All auxiliary analysis uses **1 Hz daily trend GWF files** from the CC-IN2P3 HPSS archive:

```
cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend/{YEAR}/V-trend-{GPS}-86400.gwf
```

Files are staged with `rfcp` inside SLURM jobs. Recently staged files stay in the HPSS disk buffer; successive jobs on overlapping time ranges benefit significantly from this cache.

Detection (Step 1) uses `Hrec_hoft` frames (HoftOnline or raw), typically a short epoch submitted interactively or as a single SLURM job.

## Install

```bash
python3 -m pip install -e .
```

Or run without installation:

```bash
python3 -m noisehound.cli --help
```

The recommended Python on CC-IN2P3 is the IGWN CVMFS environment:

```bash
/cvmfs/software.igwn.org/conda/envs/igwn/bin/python
```

## CC-IN2P3 / SLURM workflow

### Sync repo

```bash
# On CC-IN2P3 login node
cd ~
git clone git@github.com:danielsentenac/NoiseHound.git NOISEHOUND
# or, after initial clone:
cd ~/NOISEHOUND && git pull origin main
```

### Step 1 — Detect glitch times

Run `noisehound detect` on a short Hrec_hoft epoch:

```bash
noisehound detect \
  --channel V1:Hrec_hoft_16384Hz \
  --gps-start 1419724818 \
  --gps-end   1419735618 \
  --output    outputs/triggers.csv \
  /path/to/V-HoftOnline-*.gwf
```

Or via SLURM:

```bash
sbatch slurm/noisehound_detect.slurm
```

Output: `outputs/detect_{GPS_START}_{GPS_END}/triggers.csv`

### Step 2 — Rank auxiliary channels

Ranks all 1 Hz channels by z-score, hit fraction, and cross-correlation lag over ±5 s windows around each trigger.

Local (trend GWF already staged):

```bash
python -m noisehound.cli rank \
    --triggers  outputs/detect_.../triggers.csv \
    --include   "V1:(ENV|TCS|INF|SBE).*" \
    --top       50 \
    --output    outputs/ranking.csv \
    /path/to/V-trend-{GPS}-86400.gwf
```

Via SLURM (stages the trend GWF from HPSS then ranks):

```bash
HPSS_TREND=/hpss/in2p3.fr/group/virgo/DATA/trend/{YEAR}/V-trend-{GPS}-86400.gwf \
TRIGGERS=$HOME/NOISEHOUND/outputs/detect_.../triggers.csv \
sbatch slurm/noisehound_rank_trend.slurm
```

Output: `outputs/rank_trend_{JOBID}/ranking.csv`

### Step 3 — Causality analysis

Granger causality (VAR F-test) and Transfer Entropy (Schreiber 2000) applied to the top-ranked channels from Step 2, over ±5 s windows around each trigger.

Local:

```bash
python scripts/causality_analysis.py \
    --trend     /path/to/V-trend-{GPS}-86400.gwf \
    --hrec-dir  /path/to/hoftonline_frames/ \
    --triggers  outputs/detect_.../triggers.csv \
    --ranking   outputs/rank_trend_.../ranking.csv \
    --gps-start 1419724818 \
    --gps-end   1419735618 \
    --output    outputs/causality/
```

No dedicated SLURM script — runs in a few minutes on the login node for a short epoch.

Output: `outputs/causality/causality_summary.csv` and `causality_report.txt`.

### Step 4 — Rate correlation (1-hour bins, full baseline)

Pearson and Spearman correlation between the hourly glitch rate and each slow channel, plus a cross-correlation lag scan over ±7 days.

Local (single time block, trend GWF already staged):

```bash
python scripts/rate_correlation_direct.py \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --gps-start  1364774400 \
    --gps-end    1372636800 \
    --output     outputs/rate_correlation/partial_0.csv \
    --hpss-base  cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend \
    --bin-size-s 3600
```

Via SLURM (12 tasks × 3-month blocks, stages GWFs from HPSS):

```bash
sbatch --array=0-11 slurm/noisehound_rate_correlation.slurm
```

After all tasks finish, pull results locally and merge:

```bash
rsync -av sentenac@cca.in2p3.fr:~/NOISEHOUND/outputs/rate_correlation/ \
    outputs/rate_correlation/

python scripts/rate_correlation_direct.py \
    --merge-only \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --output-dir outputs/rate_correlation
```

Output: `outputs/rate_correlation/binned_summary.csv`, `rate_correlation_table.csv`, plots.

### Step 5 — Lag refinement (15-minute bins)

Same pipeline as Step 4, 15-minute bins on the priority channels identified in Step 4.

Local (single time block):

```bash
python scripts/rate_correlation_direct.py \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --gps-start  1364774400 \
    --gps-end    1372636800 \
    --output     outputs/rate_correlation_step4/partial_0.csv \
    --hpss-base  cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend \
    --bin-size-s 900 \
    --channels   ni_bottom_te1 wi_bottom_te1 ni_co2_env_te wi_co2_env_te
```

Via SLURM:

```bash
sbatch --array=0-11 slurm/noisehound_rate_corr_step4.slurm
```

After all tasks finish, pull and merge:

```bash
rsync -av sentenac@cca.in2p3.fr:~/NOISEHOUND/outputs/rate_correlation_step4/ \
    outputs/rate_correlation_step4/

python scripts/rate_correlation_direct.py \
    --merge-only \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --output-dir outputs/rate_correlation_step4 \
    --channels   ni_bottom_te1 wi_bottom_te1 ni_co2_env_te wi_co2_env_te
```

Output: `outputs/rate_correlation_step4/binned_summary.csv`, lag scan at 15-min resolution.

## Key scripts

| Script | Purpose |
|--------|---------|
| `scripts/rate_correlation_direct.py` | Stage trend GWFs and compute rate correlation + lag scan |
| `scripts/verify_channel_names.py` | Verify channel names against a GWF file |
| `scripts/probe_gwf_channels.py` | Stage one GWF and list channels matching a pattern |
| `scripts/fetch_logbook_entries.py` | Download Virgo logbook search results (requires SSO cookie) |

## `rate_correlation_direct.py` reference

```
# Extraction (one array task):
python scripts/rate_correlation_direct.py \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --gps-start  1364774400 \
    --gps-end    1372636800 \
    --output     outputs/rate_correlation/partial_0.csv \
    --hpss-base  cchpss0:/hpss/in2p3.fr/group/virgo/DATA/trend \
    --bin-size-s 3600

# Merge + analyse (after all tasks done):
python scripts/rate_correlation_direct.py \
    --merge-only \
    --triggers   data/full_25min_glitches_ER16-O4b.csv \
    --output-dir outputs/rate_correlation

# Restrict to specific channels (short label names from CHANNELS list):
    --channels ni_bottom_te1 wi_bottom_te1 ni_co2_env_te
```

## Use cases

| Directory | Glitch | Status |
|-----------|--------|--------|
| `usecases/25-minute-glitch/` | Periodic broadband ~25 min, SNR ~400, 40–50 Hz | Steps 1–5 ongoing |

## Analysis logic

### Detection

`noisehound detect` processes frames sequentially: computes a spectrogram per frame, normalises each frequency bin by its median over time, collapses the 20–210 Hz band into a broadband excess score, and finds robust outliers. Designed for thick full-band vertical glitches without requiring a full day of `Hrec_hoft` in memory.

### Ranking

`noisehound rank` reads one long time series per channel covering all trigger times, scores the maximum robust excursion in a short window around each glitch, compares it to matched control times offset away from the glitch, and ranks by response strength, enrichment over controls, and lag stability. This is appropriate for periodic, strong, possibly control-related disturbances that may not appear as stationary coherence excess.

### Rate correlation

`rate_correlation_direct.py` bins the trigger rate into equal-width time bins, reads the median value of each slow channel per bin from 1 Hz trend GWFs, and computes Pearson and Spearman r between rate and each channel. A cross-correlation lag scan over ±7 days identifies thermal lead times. Runs as a 12-task SLURM array (one task per 3-month block); tasks are resumable if interrupted.
