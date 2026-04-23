import pandas as pd
from pathlib import Path

out = Path("outputs/rate_correlation")
df = pd.read_csv(out / "binned_summary.csv")

print("Before:", df.shape)

# Merge ni_co2_tc into V1:TCS_NI_TE_CO2Laser (fill gaps)
for alias, v1 in [("ni_co2_tc", "V1:TCS_NI_TE_CO2Laser"),
                  ("wi_co2_tc", "V1:TCS_WI_TE_CO2Laser")]:
    if alias in df.columns and v1 in df.columns:
        filled = df[alias].notna() & df[v1].isna()
        print(f"  Filling {filled.sum()} rows from {alias!r} into {v1!r}")
        df[v1] = df[v1].fillna(df[alias])
        df.drop(columns=[alias], inplace=True)
    elif alias in df.columns:
        print(f"  Renaming {alias!r} → {v1!r}")
        df.rename(columns={alias: v1}, inplace=True)

print("After:", df.shape)
print("GPS range:", int(df.gps_bin.min()), "-", int(df.gps_bin.max()))
v1 = "V1:TCS_NI_TE_CO2Laser"
nn = df[v1].notna()
print(f"{v1} non-null: {nn.sum()}, GPS {int(df.loc[nn,'gps_bin'].min())} - {int(df.loc[nn,'gps_bin'].max())}")
df.to_csv(out / "binned_summary.csv", index=False, float_format="%.6g")
print("Written.")
