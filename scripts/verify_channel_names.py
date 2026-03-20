"""
Verify that every channel name in CHANNELS exists verbatim in a GWF file.
For any that are missing, print close matches found in the file.

Usage (as a SLURM --wrap or standalone):
    python verify_channel_names.py --gwf /path/to/day.gwf
"""
import argparse
import sys

CHANNELS = [
    ("V1:INF_NI_BOTTOM_TE1",           "NI tower bottom TE1 [°C]  ★"),
    ("V1:INF_WI_BOTTOM_TE1",           "WI tower bottom TE1 [°C]"),
    ("V1:INF_NI_MIR_COIL_UL_TE",      "NI mirror coil UL TE [°C]"),
    ("V1:INF_WI_MIR_COIL_DR_TE",      "WI mirror coil DR TE [°C]"),
    ("V1:TCS_HWS_NI_TE1_mean",         "NI HWS TE1 [°C]"),
    ("V1:TCS_HWS_NI_TE2_mean",         "NI HWS TE2 [°C]"),
    ("V1:ENV_TCS_CO2_NI_TE",           "NI CO2 bench ambient TE [°C]"),
    ("V1:TCS_NI_TE_CO2Laser",          "NI CO2 laser body TE [°C]"),
    ("V1:TCS_HWS_WI_TE1_mean",         "WI HWS TE1 [°C]"),
    ("V1:TCS_HWS_WI_TE2_mean",         "WI HWS TE2 [°C]"),
    ("V1:ENV_TCS_CO2_WI_TE",           "WI CO2 bench ambient TE [°C]"),
    ("V1:TCS_WI_TE_CO2Laser",          "WI CO2 laser body TE [°C]"),
    ("V1:TCS_HWS_NE_TE1_mean",         "NE HWS TE1 [°C]"),
    ("V1:TCS_NI_CO2_PWRLAS_mean",      "NI CO2 laser power [W]"),
    ("V1:LSC_Etalon_NI_RH_SET_mean",   "NI ring heater setpoint [W]"),
    ("V1:LSC_Etalon_NI_RH_OUT_mean",   "NI ring heater output [W]"),
    ("V1:LSC_Etalon_NI_RH_IN_mean",    "NI ring heater input [W]"),
    ("V1:LSC_Etalon_NI_RH_ERR_mean",   "NI ring heater error [W]"),
    ("V1:LSC_Etalon_WI_RH_SET_mean",   "WI ring heater setpoint [W]"),
    ("V1:LSC_Etalon_WI_RH_OUT_mean",   "WI ring heater output [W]"),
    ("V1:LSC_Etalon_WI_RH_ERR_mean",   "WI ring heater error [W]"),
    ("V1:ENV_NEB_UPS_VOLT_R_mean",     "NEB UPS voltage R [V]"),
    ("V1:ENV_CEB_UPS_VOLT_R_mean",     "CEB UPS voltage R [V]"),
    ("V1:ENV_WEB_UPS_VOLT_R_mean",     "WEB UPS voltage R [V]"),
    ("V1:ENV_MCB_IPS_CURR_T_mean",     "MCB IPS current T [A]"),
    ("V1:ENV_CEB_UPS_CURR_R_mean",     "CEB UPS current R [A]"),
    ("V1:SBE_SNEB_GEO_GRWE_raw_mean", "SNEB geophone W-E [raw]"),
    ("V1:SBE_SNEB_GEO_GRNS_raw_mean", "SNEB geophone N-S [raw]"),
    ("V1:SBE_SWEB_GEO_GRWE_raw_mean", "SWEB geophone W-E [raw]"),
    ("V1:SBE_SWEB_GEO_GRNS_raw_mean", "SWEB geophone N-S [raw]"),
    # Logbook-backed channels confirmed in GWF (Apr 2023 trend file)
    ("V1:INF_TCS_NI_RH_TE",          "NI ring heater thermistor [°C]"),
    ("V1:INF_TCS_WI_RH_TE",          "WI ring heater thermistor [°C]"),
    ("V1:ENV_CEB_N_TE",              "CEB north ambient TE [°C]"),
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
