import pandas as pd
df = pd.read_csv("outputs/rate_correlation/partial_sr.csv")
print("Columns:", list(df.columns))
print("Bins:", len(df))
for col in df.columns[2:]:
    nn = df[col].notna().sum()
    print(f"  {col}: {nn} non-null ({nn/len(df)*100:.0f}%)")
