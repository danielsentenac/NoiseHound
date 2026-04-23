import pandas as pd
import numpy as np

df = pd.read_csv("outputs/rate_correlation/partial_sr.csv")

GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
def gps_to_dt(g):
    return GPS_EPOCH + pd.to_timedelta(g, unit="s")

# Step7 epoch only
df = df[(df.gps_bin >= 1443657600) & (df.gps_bin < 1459468800)]

for col in df.columns[2:]:
    v = df[col].dropna()
    dt_min = gps_to_dt(df.loc[df[col].notna(), "gps_bin"].min())
    dt_max = gps_to_dt(df.loc[df[col].notna(), "gps_bin"].max())
    print(f"{col}")
    print(f"  min={v.min():.3f}  max={v.max():.3f}  mean={v.mean():.3f}  std={v.std():.3f}")
    print(f"  range: {dt_min.date()} – {dt_max.date()}")
