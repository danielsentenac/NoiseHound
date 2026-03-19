from __future__ import annotations

import argparse
import re
from pathlib import Path


def load_channels(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def filter_channels(channels: list[str], include: str, exclude: str | None = None) -> list[str]:
    include_re = re.compile(include)
    exclude_re = re.compile(exclude) if exclude else None
    selected = [channel for channel in channels if include_re.search(channel)]
    if exclude_re:
        selected = [channel for channel in selected if not exclude_re.search(channel)]
    return sorted(dict.fromkeys(selected))


def write_list(path: Path, channels: list[str]) -> None:
    path.write_text("\n".join(channels) + "\n")


def build_lists(channels: list[str]) -> dict[str, list[str]]:
    return {
        "ln2ni_ne_core": filter_channels(
            channels,
            include=(
                r"^V1:VAC_LN2(NI|NE)_("
                r"CLLEVEL_(AUTO|DERIV_GAIN|INT_GAIN|MAN_OUTPUT|MAX_OUTPUT|MAX_RATE|MIN_OUTPUT|PROP_GAIN|RESET|SETPOINT)"
                r"|HEATER_(AUTO|DERIV_GAIN|INT_GAIN|MAN_OUTPUT|MAX_OUTPUT|MAX_RATE|MIN_OUTPUT|PROP_GAIN|RESET|SETPOINT)"
                r"|PSLEVEL_(AUTO|DERIV_GAIN|INT_GAIN|MAN_OUTPUT|MAX_OUTPUT|MAX_RATE|MIN_OUTPUT|PROP_GAIN|RESET|SETPOINT)"
                r"|CL_(CTRLSTATE|LEVEL_ALRM|NONLEVEL_ALRM|PIDINPUT|PIDOUTPUT)"
                r"|HEATER_PID(INPUT|OUTPUT)"
                r"|PS_PID(INPUT|OUTPUT)"
                r"|DT0[12]_DISPL0[12]"
                r"|CV0[12]"
                r"|GR0[123]_.+"
                r")$"
            ),
        ),
        "hvac_inf_conditioning_core": filter_channels(
            channels,
            include=r"^V1:(HVAC|INF)_.*(TE|TEMP|HUM|PRES|FLOW|FLUX|CHILL|COOL|AIRQ|RH)",
            exclude=r"(ALARM|_ERR$|_CORR$|_SET$|State$|_ST$|_CMD$|_ON$|_OFF$)",
        ),
        "lsc_etalon_key": filter_channels(
            channels,
            include=(
                r"^V1:(LSC_Etalon_(NI|WI)_(HB_cmd|DOWN_OFFSET|RH_(ERR|IN|INPUT|OUT_(BOOST|CLIP|NORMAL)|SET))"
                r"|LSC_(NI|NE|WI|WE)_CORR_(CARM_CORR|DARM_CORR|NArm_CORR|WArm_CORR|MICH_CORR|PRCL_CORR|SRCL_CORR)"
                r"|LSC_(NI|NE|WI|WE)_LOCK_.*LOCK_ON)$"
            ),
            exclude=r"(_cnd$|_trig$|_NOISE$|_Locked$|_ENBL$|_SET_NOISE$)",
        ),
        "env_sensors_core": filter_channels(
            channels,
            include=r"^V1:ENV_(NI|NE|WI|WE|BS|CEB|DET)_.*(ACC|MIC|MAG|SEIS|IPS_(CURR|VOLT))",
            exclude=r"(_count$|_Tcomm$|_Tconv$)",
        ),
        "inf_power_core": filter_channels(
            channels,
            include=r"^V1:INF_.*(UPS|BATTERY|POWER|VOLT|CURR|CURRENT|FREQ|LOAD)",
            exclude=r"(_ST$|_ERR$|_COM_ST$|_LOCKED_ST$|_FAILURE_ST$|_INPUT_MAIN_ST$|_NORMAL_FUNC_ST$|_BYPASS_ON_ST$|_BATTERY_LOW_ST$|_BATTERY_WORK_ST$)",
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build smaller physics-driven raw candidate lists from a channel inventory.")
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()

    inventory_path = Path(args.inventory)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = f"_{args.tag}" if args.tag else ""

    channels = load_channels(inventory_path)
    groups = build_lists(channels)
    for name, items in groups.items():
        path = output_dir / f"candidate_channels_{name}{tag}.txt"
        write_list(path, items)
        print(f"{name}: {len(items)} -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
