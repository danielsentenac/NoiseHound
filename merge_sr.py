import pandas as pd
from pathlib import Path

out = Path("outputs/rate_correlation")
base = pd.read_csv(out / "binned_summary.csv")
sr   = pd.read_csv(out / "partial_sr.csv")

# Only join the new SR columns (drop n_triggers from sr)
sr_cols = [c for c in sr.columns if c not in base.columns]
merged = base.merge(sr[["gps_bin"] + sr_cols], on="gps_bin", how="left")
merged.to_csv(out / "binned_summary.csv", index=False, float_format="%.6g")
print(f"binned_summary: {len(merged)} bins, {len(merged.columns)} columns")
print("New SR columns:", sr_cols)
