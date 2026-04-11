#!/usr/bin/env python3
"""Merge 4 quarterly logbook_channels CSVs into full O4 and run the plot."""
import glob, subprocess, sys
import pandas as pd
from pathlib import Path

WORKDIR = Path(__file__).parent.parent
files = sorted(glob.glob(str(WORKDIR / "outputs" / "logbook_channels_*_*.csv")))

if len(files) < 4:
    print(f"Only {len(files)}/4 quarterly files found — aborting.")
    sys.exit(1)

dfs = [pd.read_csv(f) for f in files]
print("Shapes:", [d.shape for d in dfs])
df = pd.concat(dfs).groupby("gps_bin").mean().reset_index().sort_values("gps_bin")
out = WORKDIR / "outputs" / "logbook_channels_full_o4.csv"
df.to_csv(out, index=False)
print(f"Merged {len(df)} bins → {out}")
print("Columns:", list(df.columns))

subprocess.run([sys.executable,
                str(WORKDIR / "scripts" / "plot_logbook_channels_timeseries.py")],
               check=True)
