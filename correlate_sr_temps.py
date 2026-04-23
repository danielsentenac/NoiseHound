import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
EPOCH_START = 1443657600  # 2025-10-01
LAST_TRIGGER = 1449582668  # 2025-12-12 — only use active glitch period

df = pd.read_csv("outputs/rate_correlation/binned_summary.csv")
df = df[(df.gps_bin >= EPOCH_START) & (df.gps_bin <= LAST_TRIGGER)]

rate = df["n_triggers"].values.astype(float)

sr_cols = [
    "V1:INF_SR_MIR_COIL_UL_TE",
    "V1:INF_TCS_SR_RH_TE",
]
ni_refs = [
    "V1:INF_NI_BOTTOM_TE1",
    "V1:ENV_TCS_CO2_NI_TE",
    "V1:ENV_TCS_CO2_WI_TE",
]

print(f"Active glitch period bins: {len(df)}  (bins with triggers: {(rate>0).sum()})\n")
print(f"{'Channel':<45} {'Pearson r':>10} {'Spearman r':>11} {'n':>6}")
print("-" * 76)
for col in sr_cols + ni_refs:
    if col not in df.columns:
        print(f"  {col}: NOT IN DATA")
        continue
    v = df[col].values.astype(float)
    mask = np.isfinite(v) & np.isfinite(rate)
    if mask.sum() < 10:
        print(f"  {col}: too few valid bins ({mask.sum()})")
        continue
    pr, pp = pearsonr(rate[mask], v[mask])
    sr, sp = spearmanr(rate[mask], v[mask])
    print(f"  {col:<43} {pr:>+10.3f} {sr:>+11.3f} {mask.sum():>6}")
