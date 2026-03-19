# NOISEHOUND

`NOISEHOUND` is a small local toolkit for investigating the long-lived periodic broadband Virgo glitch you described. It is built to replace one-by-one coherence browsing with a repeatable workflow:

1. discover frame coverage and URLs with `gwdatafind`
2. detect the broadband glitch times directly in `Hrec_hoft`
3. stack many glitches together
4. rank auxiliary channels by event-synchronous response, not by stationary coherence alone

## Current data-access findings

These are already verified from this workspace.

- Discovery works through the Python `gwdatafind` API with a token that has `gwdatafind.read`.
- If your computing centre exports `GWDATAFIND_SERVER`, `NOISEHOUND` can use that automatically.
- Virgo private frame types visible there include `raw`, `HoftOnline`, `HoftAR1`, `HoftAR1U02`, and `envDataTransfer`.
- For Virgo, the useful URL types are `osdf` and `https`. The default `file` URL query returns nothing.
- Example URL resolution for the target day works for non-raw streams:
  - `HoftOnline` around GPS `1424300399` resolves to `osdf:///igwn/virgo/frames/O4/HoftOnline/...`
  - `envDataTransfer` around GPS `1424300399` resolves to `osdf:///igwn/virgo/frames/envmon/envDataTransfer/...`
- The current token only has `gwdatafind.read`, so discovery works but OSDF download fails with `403 Forbidden`.
- The current working assumption is split-source:
  - use `HoftOnline` to detect the glitch in `Hrec_hoft`
  - use `raw` for the auxiliary/control-channel investigation
- `raw` is still the preferred auxiliary source because it should contain the broad channel space needed for NI/NE/WI/WE, control, thermal, suspension, and environmental follow-up.
- On CC, `raw` should be treated as an HPSS archive workflow, not an IGWN datafind workflow:
  - list archive content with `rfdir cchpss0:/hpss/in2p3.fr/group/virgo/DATA/archive_raw`
  - on `cca020`, copy a file with `rfcp cchpss0:/hpss/in2p3.fr/group/virgo/<path-to-file> /scratch/$USER`

## Install

From this directory:

```bash
python3 -m pip install -e .
```

Or run directly without installation:

```bash
python3 -m noisehound.cli --help
```

## CC / SLURM workflow

The intended production workflow is on the CC login node and Slurm cluster, not on a laptop.

1. Sync this repo to `cca020`:

```bash
./scripts/cc_sync.sh sentenac@cca020.in2p3.fr /pbs/home/s/sentenac/NOISEHOUND
```

2. Log in to CC and create the virtual environment:

```bash
cd ~/NOISEHOUND
./scripts/cc_bootstrap_env.sh
```

If your centre config exports a local GWDataFind server, check it:

```bash
echo "$GWDATAFIND_SERVER"
```

3. Prepare the detection and auxiliary frame-list files:

```bash
GPS_START=1424300000 GPS_END=1424388000 ./scripts/cc_make_frame_lists.sh
```

This writes two frame lists, by default:

- `frame_lists/hrec_HoftOnline_1424300000_1424388000_osdf.txt`
- `frame_lists/aux_raw_1424300000_1424388000_osdf.txt`

The intended workflow is:

- detect glitch times from `HoftOnline`
- rank candidate auxiliary channels from `raw`

The current Slurm scripts still expect local `.gwf` paths for the actual analysis step. For `HoftOnline`, that means local staging if your list contains `osdf://...`. For `raw`, the expected path is HPSS staging to a compute-node-local filesystem.

For `.gwf` frame I/O, the Slurm wrappers now prefer the IGWN CVMFS Python at `/cvmfs/software.igwn.org/conda/envs/igwn/bin/python`, because it already includes working `LDAStools.frameCPP` / `framel` / `lalframe` backends.

4. Stage HPSS raw files to local scratch once you have identified the raw frame paths.

For long HPSS recalls, prefer a Slurm job instead of a login-node background process:

```bash
SOURCE_LIST=$HOME/NOISEHOUND/frame_lists/aux_raw_hpss_sources.txt \
OUTPUT_LIST=$HOME/NOISEHOUND/frame_lists/aux_raw_1424300000_1424388000_local.txt \
sbatch $HOME/NOISEHOUND/slurm/noisehound_stage_hpss.slurm
```

By default the job resolves a writable staging directory at runtime on the compute node. Only set `DEST_DIR` explicitly if you already know a valid compute-node-local path.

For a short interactive test, you can still run the staging helper directly:

```bash
SOURCE_LIST=$HOME/NOISEHOUND/frame_lists/aux_raw_hpss_sources.txt \
OUTPUT_LIST=$HOME/NOISEHOUND/frame_lists/aux_raw_1424300000_1424388000_local.txt \
./scripts/cc_stage_hpss.sh
```

The source list can contain any of these forms, one per line:

- `chpss0:/hpss/in2p3.fr/group/virgo/DATA/archive_raw/.../V-raw-....gwf`
- `cchpss0:/hpss/in2p3.fr/group/virgo/DATA/archive_raw/.../V-raw-....gwf`
- `/hpss/in2p3.fr/group/virgo/DATA/archive_raw/.../V-raw-....gwf`
- `DATA/archive_raw/.../V-raw-....gwf`

If you are staging from `cca020`, the helper defaults to `HPSS_HOST=cchpss0` and picks the first writable local destination among `SLURM_TMPDIR`, `TMPDIR`, `/tmp/$USER`, and a last-resort fallback under `$HOME/NOISEHOUND/staging`.

The Slurm job logs to `slurm-nh-stage-<jobid>.out` and writes detailed copy progress through `scripts/cc_stage_hpss.sh`.

If the runtime staging directory resolves to node-local `/tmp`, the copied raw files are only useful inside that same job. In that case, prefer a combined stage-and-rank job:

```bash
AUX_SOURCE_LIST=$HOME/NOISEHOUND/frame_lists/aux_raw_hpss_sources.txt \
TRIGGERS=$HOME/NOISEHOUND/outputs/detect_1420214800_1420218400/triggers.csv \
CHANNEL_FILE=$HOME/NOISEHOUND/candidate_channels.txt \
MAX_SAMPLE_RATE=64 \
sbatch $HOME/NOISEHOUND/slurm/noisehound_rank_hpss.slurm
```

This stages the HPSS raw frames and runs `noisehound rank` in the same Slurm allocation, so node-local `/tmp` is no longer a problem.

For interactive frame inspection on CC, use the same interpreter directly:

```bash
PYTHONPATH=$HOME/NOISEHOUND \
/cvmfs/software.igwn.org/conda/envs/igwn/bin/python \
  -m noisehound.cli channels /path/to/frame.gwf --pattern 'V1:Hrec_hoft.*'
```

5. Submit detection on the `htc` partition:

```bash
FRAME_LIST_FILE=$HOME/NOISEHOUND/frame_lists/hrec_HoftOnline_1424300000_1424388000_local.txt \
GPS_START=1424300000 \
GPS_END=1424388000 \
sbatch $HOME/NOISEHOUND/slurm/noisehound_detect.slurm
```

6. Submit auxiliary ranking against raw frames:

```bash
AUX_FRAME_LIST_FILE=$HOME/NOISEHOUND/frame_lists/aux_raw_1424300000_1424388000_local.txt \
TRIGGERS=$HOME/NOISEHOUND/outputs/detect_1424300000_1424388000/triggers.csv \
CHANNEL_FILE=$HOME/NOISEHOUND/candidate_channels.txt \
MAX_SAMPLE_RATE=64 \
sbatch $HOME/NOISEHOUND/slurm/noisehound_rank.slurm
```

If your raw frames are still in HPSS rather than already-local shared storage, use `slurm/noisehound_rank_hpss.slurm` instead of `slurm/noisehound_rank.slurm`.

On `cca020`, `sbatch`, `squeue`, and `sinfo` are available, and the visible partitions include `htc`, `htc_highmem`, and GPU queues. The current scripts default to `htc`.

## Commands

Show Virgo frame types:

```bash
noisehound types --server datafind.igwn.org --observatory V --token-file /path/to/token
```

If `GWDATAFIND_SERVER` is set, you can omit `--server`.

Show availability segments for a frame type:

```bash
noisehound segments \
  --server datafind.igwn.org \
  --observatory V \
  --frame-type HoftOnline \
  --gps-start 1424300000 \
  --gps-end 1424388000 \
  --token-file /path/to/token
```

Resolve `HoftOnline` URLs for `Hrec_hoft`:

```bash
noisehound urls \
  --server datafind.igwn.org \
  --observatory V \
  --frame-type HoftOnline \
  --gps-start 1424300399 \
  --gps-end 1424300499 \
  --url-type osdf \
  --token-file /path/to/token
```

Resolve `raw` URLs for the auxiliary/control investigation if your centre exposes them via datafind:

```bash
noisehound urls \
  --server datafind.igwn.org \
  --observatory V \
  --frame-type raw \
  --gps-start 1424300399 \
  --gps-end 1424300499 \
  --url-type osdf \
  --token-file /path/to/token
```

List channels from a local raw frame:

```bash
noisehound channels /path/to/V-raw-*.gwf --pattern 'V1:(ENV|TCS|SUSP|SBE|.*NI.*|.*NE.*|.*WI.*|.*WE.*|Hrec_hoft.*)'
```

Detect the periodic broadband glitch directly from local `HoftOnline` frames:

```bash
noisehound detect \
  --channel V1:Hrec_hoft_16384Hz \
  --gps-start 1424300000 \
  --gps-end 1424388000 \
  --output triggers.csv \
  --plot triggers.png \
  /path/to/V-HoftOnline-*.gwf
```

Rank candidate auxiliary channels against those triggers:

```bash
noisehound rank \
  --triggers triggers.csv \
  --channel-file candidate_channels.txt \
  --output ranking.csv \
  /path/to/V-raw-*.gwf
```

## Analysis logic

### Detection

`noisehound detect` does not look for a narrow line. It processes frames sequentially, computes a spectrogram for each frame, normalizes each frequency bin by its own median over time, collapses the `20-210 Hz` band into a broadband excess score, then finds robust outliers in that score. This is meant for the thick full-band vertical glitches visible in your spectrogram without requiring one full day of `Hrec_hoft` to sit in memory at once.

### Ranking

`noisehound rank` is intentionally not a coherence replacement. For each candidate channel, it:

- reads one long time series covering all trigger times
- scores the maximum robust excursion in a short window around each glitch
- compares that to matched control times offset away from the glitch
- ranks channels by repeated response strength, enrichment over controls, and lag stability

This is the right direction for a periodic, strong, possibly control-related disturbance that may not appear as a simple stationary coherence excess.

## Suggested first-pass candidate groups

Once you can read raw locally on CC, start by building channel lists around:

- NI / NE / WI / WE cavity controls
- longitudinal control channels first
- thermal and temperature regulation channels
- TCS channels
- suspension and bench electronics channels
- environmental electric and magnetic monitors near the towers

## Important limitation

The toolkit works on local frame files. At the moment the remaining blockers are frame access on CC:

- discovery token: already enough for `gw_data_find`
- `HoftOnline` analysis still needs either staging or a local readable path
- `raw` auxiliary analysis should use the HPSS raw archive on CC, staged locally with `rfcp`

Once that is in place, this directory is ready for actual data analysis rather than only URL discovery.
