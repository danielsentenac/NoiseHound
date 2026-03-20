"""
migrate_csv_columns.py — rename alias column names to V1: channel names in
existing partial and summary CSV files produced by rate_correlation_direct.py.

Run this once after updating the scripts to use V1: names as column names.

Usage:
    python scripts/migrate_csv_columns.py [--dry-run]
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import pandas as pd

# alias → V1: channel name mapping
ALIAS_TO_V1 = {
    "ni_bottom_te1":  "V1:INF_NI_BOTTOM_TE1",
    "wi_bottom_te1":  "V1:INF_WI_BOTTOM_TE1",
    "ni_mir_coil_te": "V1:INF_NI_MIR_COIL_UL_TE",
    "wi_mir_coil_te": "V1:INF_WI_MIR_COIL_DR_TE",
    "ni_hws_te1":     "V1:TCS_HWS_NI_TE1_mean",
    "ni_hws_te2":     "V1:TCS_HWS_NI_TE2_mean",
    "ni_co2_env_te":  "V1:ENV_TCS_CO2_NI_TE",
    "ni_co2_laser_te":"V1:TCS_NI_TE_CO2Laser",
    "wi_hws_te1":     "V1:TCS_HWS_WI_TE1_mean",
    "wi_hws_te2":     "V1:TCS_HWS_WI_TE2_mean",
    "wi_co2_env_te":  "V1:ENV_TCS_CO2_WI_TE",
    "wi_co2_laser_te":"V1:TCS_WI_TE_CO2Laser",
    "ne_hws_te1":     "V1:TCS_HWS_NE_TE1_mean",
    "ni_co2_pwr":     "V1:TCS_NI_CO2_PWRLAS_mean",
    "ni_rh_set":      "V1:LSC_Etalon_NI_RH_SET_mean",
    "ni_rh_out":      "V1:LSC_Etalon_NI_RH_OUT_mean",
    "ni_rh_in":       "V1:LSC_Etalon_NI_RH_IN_mean",
    "ni_rh_err":      "V1:LSC_Etalon_NI_RH_ERR_mean",
    "wi_rh_set":      "V1:LSC_Etalon_WI_RH_SET_mean",
    "wi_rh_out":      "V1:LSC_Etalon_WI_RH_OUT_mean",
    "wi_rh_err":      "V1:LSC_Etalon_WI_RH_ERR_mean",
    "neb_ups_volt_r": "V1:ENV_NEB_UPS_VOLT_R_mean",
    "ceb_ups_volt_r": "V1:ENV_CEB_UPS_VOLT_R_mean",
    "web_ups_volt_r": "V1:ENV_WEB_UPS_VOLT_R_mean",
    "mcb_ips_curr_t": "V1:ENV_MCB_IPS_CURR_T_mean",
    "ceb_ups_curr_r": "V1:ENV_CEB_UPS_CURR_R_mean",
    "sneb_geo_we":    "V1:SBE_SNEB_GEO_GRWE_raw_mean",
    "sneb_geo_ns":    "V1:SBE_SNEB_GEO_GRNS_raw_mean",
    "sweb_geo_we":    "V1:SBE_SWEB_GEO_GRWE_raw_mean",
    "sweb_geo_ns":    "V1:SBE_SWEB_GEO_GRNS_raw_mean",
    "ni_rh_te":       "V1:INF_TCS_NI_RH_TE",
    "wi_rh_te":       "V1:INF_TCS_WI_RH_TE",
    "ceb_n_te":       "V1:ENV_CEB_N_TE",
}

# Repo root (two levels up from this script)
REPO_ROOT = Path(__file__).resolve().parent.parent

# CSV files to migrate
CSV_PATTERNS = [
    "outputs/rate_correlation/binned_summary.csv",
    "outputs/rate_correlation/partial_gap.csv",
    "outputs/rate_correlation/partial_janapr.csv",
    "outputs/rate_correlation_step4/partial_*.csv",
]


def migrate_file(path: Path, dry_run: bool) -> None:
    df = pd.read_csv(path)
    rename_map = {alias: v1 for alias, v1 in ALIAS_TO_V1.items() if alias in df.columns}
    if not rename_map:
        print(f"  SKIP (no alias columns found): {path}")
        return
    print(f"  {'DRY-RUN ' if dry_run else ''}RENAME {len(rename_map)} columns in {path}:")
    for alias, v1 in rename_map.items():
        print(f"    {alias!r} → {v1!r}")
    if not dry_run:
        df = df.rename(columns=rename_map)
        df.to_csv(path, index=False, float_format="%.6g")
        print(f"  Written: {path}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would be done without writing any files")
    args = p.parse_args()

    csv_files: list[Path] = []
    for pattern in CSV_PATTERNS:
        matches = sorted(REPO_ROOT.glob(pattern))
        csv_files.extend(matches)

    if not csv_files:
        print("No CSV files found matching the expected patterns.")
        return

    print(f"Found {len(csv_files)} CSV file(s) to check.")
    for path in csv_files:
        migrate_file(path, dry_run=args.dry_run)

    if args.dry_run:
        print("\nDry-run complete — no files were modified.")
    else:
        print("\nMigration complete.")


if __name__ == "__main__":
    main()
