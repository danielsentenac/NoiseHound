#!/usr/bin/env python3
"""
Merge segment_scan_dip1.csv + segment_scan_dip2.csv and rank channels
by how strongly and proportionally they vary with the glitch-SNR dip.

Three checks, all required for a strong candidate:

  1. V-shape:        both UP segments sit on the same side of DOWN (within each dip)
  2. Consistency:    same direction of change in BOTH dips
  3. Proportionality: channel delta scales with SNR delta depth across dips
                     sensitivity = channel_delta / snr_delta
                     prop_score  = min(|s1|,|s2|) / max(|s1|,|s2|)  ∈ [0,1]
                     (= 1 when perfectly proportional, 0 when one dip shows no change)

Final score:
  combined_z   = sign * geom_mean(|z1|, |z2|)    [strength]
  final_score  = |combined_z| * prop_score         [strength × proportionality]

Output: outputs/segment_scan_combined.csv
"""
import numpy as np
import pandas as pd
from pathlib import Path

WORKDIR     = Path(__file__).parent.parent
GLITCH_FILE = WORKDIR / "data" / "full_25min_glitches_ER16-O4b.csv"

# ── SNR segment medians (from glitch catalog) ─────────────────────────────────
GPS_EPOCH = pd.Timestamp("1980-01-06", tz="UTC")
glitch = pd.read_csv(GLITCH_FILE)
glitch["time"] = GPS_EPOCH + pd.to_timedelta(glitch["time"], unit="s")
glitch["gps"]  = glitch["time"].apply(lambda t: (t - GPS_EPOCH).total_seconds())

def snr_median(gps_s, gps_e):
    m = glitch[(glitch["gps"] >= gps_s) & (glitch["gps"] < gps_e)]["snr"]
    return float(m.median()) if len(m) > 0 else np.nan

# Dip 1 segment GPS boundaries
D1_UP1_S, D1_UP1_E = 1411689600, 1415318400
D1_DN_S,  D1_DN_E  = 1415318400, 1416528000
D1_UP2_S, D1_UP2_E = 1416528000, 1419552000

# Dip 2 segment GPS boundaries
D2_UP1_S, D2_UP1_E = 1424995200, 1429228800
D2_DN_S,  D2_DN_E  = 1429228800, 1431043200
D2_UP2_S, D2_UP2_E = 1431043200, 1435276800

snr_d1_up  = (snr_median(D1_UP1_S, D1_UP1_E) + snr_median(D1_UP2_S, D1_UP2_E)) / 2
snr_d1_dn  =  snr_median(D1_DN_S,  D1_DN_E)
snr_d2_up  = (snr_median(D2_UP1_S, D2_UP1_E) + snr_median(D2_UP2_S, D2_UP2_E)) / 2
snr_d2_dn  =  snr_median(D2_DN_S,  D2_DN_E)

snr_delta1 = snr_d1_dn - snr_d1_up   # negative (SNR drops during dip)
snr_delta2 = snr_d2_dn - snr_d2_up

print(f"SNR dip 1: UP={snr_d1_up:.1f}  DOWN={snr_d1_dn:.1f}  delta={snr_delta1:.1f}")
print(f"SNR dip 2: UP={snr_d2_up:.1f}  DOWN={snr_d2_dn:.1f}  delta={snr_delta2:.1f}")
print(f"Expected depth ratio dip1/dip2 = {snr_delta1/snr_delta2:.3f}")

# ── Load channel scan results ─────────────────────────────────────────────────
d1 = pd.read_csv(WORKDIR / "outputs" / "segment_scan_dip1.csv").set_index("channel")
d2 = pd.read_csv(WORKDIR / "outputs" / "segment_scan_dip2.csv").set_index("channel")

common = d1.index.intersection(d2.index)
print(f"\nChannels in dip1: {len(d1)}  dip2: {len(d2)}  common: {len(common)}")

d1c = d1.loc[common]
d2c = d2.loc[common]

z1     = d1c["z_score"].values
z2     = d2c["z_score"].values
delta1 = d1c["delta"].values    # channel_delta for dip 1
delta2 = d2c["delta"].values    # channel_delta for dip 2

# ── Check 1: V-shape within each dip ─────────────────────────────────────────
v1 = (np.sign(d1c["mu_down"].values - d1c["mu_up1"].values) ==
      np.sign(d1c["mu_down"].values - d1c["mu_up2"].values))
v2 = (np.sign(d2c["mu_down"].values - d2c["mu_up1"].values) ==
      np.sign(d2c["mu_down"].values - d2c["mu_up2"].values))
v_shape_both = v1 & v2

# ── Check 2: Direction consistency across dips ────────────────────────────────
consistent = np.sign(z1) == np.sign(z2)

# Combined z (strength, zero if inconsistent)
combined_z = np.where(consistent,
                      np.sign(z1) * np.sqrt(np.abs(z1) * np.abs(z2)),
                      0.0)

# ── Check 3: Proportionality ──────────────────────────────────────────────────
# sensitivity = channel_delta / snr_delta  (units: channel_unit / SNR_unit)
# For a truly proportional channel: sensitivity_dip1 ≈ sensitivity_dip2
# We use the absolute ratio of sensitivities — should be ~1 for proportional channels.

eps = 1e-30
s1 = delta1 / (snr_delta1 + eps)   # sensitivity dip 1
s2 = delta2 / (snr_delta2 + eps)   # sensitivity dip 2

# prop_score in [0, 1]:
#   1.0 = perfectly proportional (s1 == s2)
#   0.5 = one sensitivity is twice the other
#   0.0 = one sensitivity is zero while the other is not
abs_s1 = np.abs(s1)
abs_s2 = np.abs(s2)
prop_score = np.where(
    (abs_s1 > 0) & (abs_s2 > 0),
    np.minimum(abs_s1, abs_s2) / np.maximum(abs_s1, abs_s2),
    0.0
)
# Only meaningful when consistent
prop_score = np.where(consistent, prop_score, 0.0)

# ── Final score: strength × proportionality ───────────────────────────────────
final_score = np.abs(combined_z) * prop_score

direction = np.where(combined_z > 0, "UP_during_dip",
            np.where(combined_z < 0, "DOWN_during_dip", "inconsistent"))

# ── Output ────────────────────────────────────────────────────────────────────
out = pd.DataFrame({
    "channel":       list(common),
    "final_score":   final_score.round(3),
    "combined_z":    combined_z.round(3),
    "prop_score":    prop_score.round(3),
    "direction":     direction,
    "v_shape_both":  v_shape_both,
    "consistent":    consistent,
    "z_dip1":        z1.round(3),
    "z_dip2":        z2.round(3),
    "sensitivity1":  s1.round(6),
    "sensitivity2":  s2.round(6),
    "mu_up1_d1":     d1c["mu_up1"].values,
    "mu_down_d1":    d1c["mu_down"].values,
    "mu_up2_d1":     d1c["mu_up2"].values,
    "mu_up1_d2":     d2c["mu_up1"].values,
    "mu_down_d2":    d2c["mu_down"].values,
    "mu_up2_d2":     d2c["mu_up2"].values,
    "delta_d1":      d1c["delta"].values,
    "delta_d2":      d2c["delta"].values,
    "sigma_up_d1":   d1c["sigma_up"].values,
    "sigma_up_d2":   d2c["sigma_up"].values,
})
out = out.sort_values("final_score", ascending=False)

out_file = WORKDIR / "outputs" / "segment_scan_combined.csv"
out.to_csv(out_file, index=False)
print(f"Saved {len(out)} combined entries → {out_file}")

# ── Summary ───────────────────────────────────────────────────────────────────
mask_best = out["consistent"] & out["v_shape_both"]
print(f"\nConsistent + v-shape:          {mask_best.sum():5d} / {len(out)}")
print(f"  of which prop_score > 0.8:   {(mask_best & (out['prop_score'] > 0.8)).sum():5d}")
print(f"  of which prop_score > 0.5:   {(mask_best & (out['prop_score'] > 0.5)).sum():5d}")
print(f"final_score > 3:               {(out['final_score'] > 3).sum():5d}")
print(f"final_score > 5:               {(out['final_score'] > 5).sum():5d}")

# ── Top ranked ───────────────────────────────────────────────────────────────
print("\n=== Top 50: consistent + v-shape, ranked by final_score = |combined_z| × prop_score ===")
top = out[mask_best].head(50)
for _, r in top.iterrows():
    print(f"  final={r.final_score:6.2f}  z={r.combined_z:+6.2f}  prop={r.prop_score:.2f}"
          f"  [{r.direction:16s}]  z1={r.z_dip1:+5.2f} z2={r.z_dip2:+5.2f}"
          f"   {r.channel}")

print("\n=== Top 20 UP_during_dip ===")
for _, r in out[out["direction"] == "UP_during_dip"].head(20).iterrows():
    print(f"  final={r.final_score:6.2f}  z={r.combined_z:+6.2f}  prop={r.prop_score:.2f}"
          f"   s1={r.sensitivity1:+.4g}  s2={r.sensitivity2:+.4g}   {r.channel}")

print("\n=== Top 20 DOWN_during_dip ===")
for _, r in out[out["direction"] == "DOWN_during_dip"].head(20).iterrows():
    print(f"  final={r.final_score:6.2f}  z={r.combined_z:+6.2f}  prop={r.prop_score:.2f}"
          f"   s1={r.sensitivity1:+.4g}  s2={r.sensitivity2:+.4g}   {r.channel}")
