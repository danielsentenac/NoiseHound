import pandas as pd
df = pd.read_csv("outputs/rate_correlation/binned_summary.csv")
laser = "V1:TCS_NI_TE_CO2Laser"
old = "ni_co2_tc"
nn_l = df[laser].notna()
nn_o = df[old].notna()
print("V1:TCS_NI_TE_CO2Laser non-null:", nn_l.sum(),
      "GPS", int(df.loc[nn_l,"gps_bin"].min()), "-", int(df.loc[nn_l,"gps_bin"].max()))
print("ni_co2_tc non-null:", nn_o.sum(),
      "GPS", int(df.loc[nn_o,"gps_bin"].min()), "-", int(df.loc[nn_o,"gps_bin"].max()))
print("both non-null:", (nn_l & nn_o).sum())
print("only ni_co2_tc:", (nn_o & ~nn_l).sum())
print("only laser:", (nn_l & ~nn_o).sum())
