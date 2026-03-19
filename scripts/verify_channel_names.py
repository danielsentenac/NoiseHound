"""
Verify that every channel name in CHANNELS exists verbatim in a GWF file.
For any that are missing, print close matches found in the file.

Usage (as a SLURM --wrap or standalone):
    python verify_channel_names.py --gwf /path/to/day.gwf
"""
import argparse
import sys

CHANNELS = [
    ("V1:INF_NI_BOTTOM_TE1",           "ni_bottom_te1"),
    ("V1:INF_WI_BOTTOM_TE1",           "wi_bottom_te1"),
    ("V1:INF_NI_MIR_COIL_UL_TE",      "ni_mir_coil_te"),
    ("V1:INF_WI_MIR_COIL_DR_TE",      "wi_mir_coil_te"),
    ("V1:TCS_HWS_NI_TE1_mean",         "ni_hws_te1"),
    ("V1:TCS_HWS_NI_TE2_mean",         "ni_hws_te2"),
    ("V1:ENV_TCS_CO2_NI_TE",            "ni_co2_env_te"),
    ("V1:TCS_NI_TE_CO2Laser",          "ni_co2_laser_te"),
    ("V1:TCS_HWS_WI_TE1_mean",         "wi_hws_te1"),
    ("V1:TCS_HWS_WI_TE2_mean",         "wi_hws_te2"),
    ("V1:ENV_TCS_CO2_WI_TE",           "wi_co2_env_te"),
    ("V1:TCS_WI_TE_CO2Laser",          "wi_co2_laser_te"),
    ("V1:TCS_HWS_NE_TE1_mean",         "ne_hws_te1"),
    ("V1:TCS_NI_CO2_PWRLAS_mean",      "ni_co2_pwr"),
    ("V1:LSC_Etalon_NI_RH_SET_mean",   "ni_rh_set"),
    ("V1:LSC_Etalon_NI_RH_OUT_mean",   "ni_rh_out"),
    ("V1:LSC_Etalon_NI_RH_IN_mean",    "ni_rh_in"),
    ("V1:LSC_Etalon_NI_RH_ERR_mean",   "ni_rh_err"),
    ("V1:LSC_Etalon_WI_RH_SET_mean",   "wi_rh_set"),
    ("V1:LSC_Etalon_WI_RH_OUT_mean",   "wi_rh_out"),
    ("V1:LSC_Etalon_WI_RH_ERR_mean",   "wi_rh_err"),
    ("V1:ENV_NEB_UPS_VOLT_R_mean",     "neb_ups_volt_r"),
    ("V1:ENV_CEB_UPS_VOLT_R_mean",     "ceb_ups_volt_r"),
    ("V1:ENV_WEB_UPS_VOLT_R_mean",     "web_ups_volt_r"),
    ("V1:ENV_MCB_IPS_CURR_T_mean",     "mcb_ips_curr_t"),
    ("V1:ENV_CEB_UPS_CURR_R_mean",     "ceb_ups_curr_r"),
    ("V1:SBE_SNEB_GEO_GRWE_raw_mean", "sneb_geo_we"),
    ("V1:SBE_SNEB_GEO_GRNS_raw_mean", "sneb_geo_ns"),
    ("V1:SBE_SWEB_GEO_GRWE_raw_mean", "sweb_geo_we"),
    ("V1:SBE_SWEB_GEO_GRNS_raw_mean", "sweb_geo_ns"),
]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gwf", required=True)
    args = p.parse_args()

    from gwpy.io.gwf import iter_channel_names
    all_channels = set(iter_channel_names(args.gwf))
    print(f"Total channels in file: {len(all_channels)}\n")

    ok, missing = [], []
    for ch_name, label in CHANNELS:
        if ch_name in all_channels:
            ok.append((ch_name, label))
        else:
            missing.append((ch_name, label))

    print(f"=== FOUND ({len(ok)}/{len(CHANNELS)}) ===")
    for ch_name, label in ok:
        print(f"  OK  {ch_name}")

    if missing:
        print(f"\n=== MISSING ({len(missing)}/{len(CHANNELS)}) ===")
        for ch_name, label in missing:
            # strip V1: prefix and any trailing _mean/_min/_max/_rms/_raw
            base = ch_name.split(":")[1]
            for suffix in ("_mean", "_min", "_max", "_rms", "_raw"):
                if base.endswith(suffix):
                    base = base[: -len(suffix)]
                    break
            close = sorted(c for c in all_channels if base.upper() in c.upper())
            print(f"\n  MISS {ch_name}  [{label}]")
            if close:
                print(f"       close matches in GWF:")
                for c in close[:10]:
                    print(f"         {c}")
            else:
                print(f"       no match for base '{base}' in GWF")


if __name__ == "__main__":
    main()
