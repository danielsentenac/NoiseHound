#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image


WORKDIR = Path("/home/sentenac/NOISEHOUND/usecases/25-minute-glitch")
TEMPLATE = WORKDIR / "legend_crops" / "legend_template_1hz_bigmarkers.png"

# Legend rectangles in the template coordinate system.
TEMPLATE_LEGEND_BOXES = [
    (45, 176, 500, 256),
    (48, 579, 290, 650),
    (48, 854, 635, 944),
    (48, 1126, 760, 1216),
    (48, 1463, 360, 1537),
    (48, 1837, 650, 1918),
]


def big_point(color: str, label: str) -> Line2D:
    return Line2D([], [], ls="None", marker="o", ms=9.0, color=color, label=label)


def line_handle(color: str, label: str) -> Line2D:
    return Line2D([], [], ls="--", lw=1.2, color=color, label=label)


def legend(ax: plt.Axes, handles, **kwargs):
    return ax.legend(handles=handles, labels=[h.get_label() for h in handles], **kwargs)


def build_template(path: Path) -> tuple[int, int]:
    fig, axes = plt.subplots(
        6,
        1,
        figsize=(16, 16),
        sharex=True,
        gridspec_kw={"height_ratios": [1.5, 1.0, 1.2, 1.0, 0.9, 1.1], "hspace": 0.08},
    )
    start = pd.Timestamp("2025-10-01", tz="UTC")
    end = pd.Timestamp("2026-02-28 23:59:59", tz="UTC")
    first_diaphragm = pd.Timestamp("2025-12-25", tz="UTC")
    second_diaphragm = pd.Timestamp("2026-02-07", tz="UTC")

    legends = [
        [
            big_point("#1f77b4", "V1:DET_B1p_DC_mean"),
            line_handle("#666666", "70% line (145.0 mW)"),
            line_handle("#999999", "50% line (75.0 mW)"),
        ],
        [big_point("#2c7fb8", "V1:Hrec_Range_BNS")],
        [
            big_point("tab:red", "V1:INF_NI_BOTTOM_TE1"),
            big_point("tab:orange", "V1:INF_NI_MIR_COIL_UL_TE"),
            big_point("tab:purple", "V1:INF_TCS_NI_RH_TE"),
            big_point("darkgreen", "V1:ENV_TCS_CO2_NI_TE"),
        ],
        [
            big_point("tab:brown", "V1:ENV_NI_F0_TE1"),
            big_point("peru", "V1:ENV_NI_F0_TE2"),
            big_point("tab:green", "V1:ENV_NI_F4_TE1"),
            big_point("limegreen", "V1:ENV_NI_F4_TE2"),
            big_point("tab:cyan", "V1:ENV_NI_F7_TE1"),
            big_point("deepskyblue", "V1:ENV_NI_F7_TE2"),
        ],
        [
            big_point("tab:blue", "V1:TCS_NI_TE_CO2Laser"),
            big_point("tab:olive", "V1:TCS_NI_TE_AUXLaser"),
        ],
        [
            big_point("tab:blue", "V1:LSC_Etalon_NI_RH_OUT_mean"),
            big_point("tab:green", "V1:LSC_Etalon_NI_RH_SET_mean"),
            big_point("goldenrod", "V1:LSC_Etalon_NI_RH_IN_mean"),
            big_point("crimson", "V1:LSC_Etalon_NI_RH_ERR_mean"),
        ],
    ]
    locs = [
        dict(loc="upper left", fontsize=8, ncol=2, framealpha=0.9),
        dict(loc="lower left", fontsize=8, framealpha=0.9),
        dict(loc="lower left", fontsize=8, ncol=2, framealpha=0.9),
        dict(loc="lower left", fontsize=8, ncol=3, framealpha=0.9),
        dict(loc="center left", fontsize=8, framealpha=0.9),
        dict(loc="upper left", fontsize=8, ncol=2, framealpha=0.9),
    ]

    for idx, ax in enumerate(axes):
        ax.set_xlim(start, end)
        ax.grid(alpha=0.25)
        ax.axvline(first_diaphragm, color="black", lw=1.0, ls="--", alpha=0.65)
        ax.axvline(second_diaphragm, color="black", lw=1.0, ls="--", alpha=0.65)
        legend(ax, legends[idx], **locs[idx])

    axes[0].set_ylabel("B1p carrier\n[mW]")
    axes[0].set_ylim(0, 500)
    axes[1].set_ylabel("Hrec range BNS\n[Mpc]")
    axes[2].set_ylabel("Temperature\n[°C]")
    axes[3].set_ylabel("Additional NI thermal\n[degC]")
    axes[4].set_ylabel("NI laser TE\n[degC]")
    axes[5].set_ylabel("NI etalon\nraw value")
    axes[5].set_ylim(0, 3000)
    axes[-1].xaxis.set_major_locator(mdates.MonthLocator())
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    axes[-1].tick_params(axis="x", rotation=25, labelsize=9)
    for axis in axes[:-1]:
        axis.tick_params(labelsize=9)
    fig.suptitle(
        "Slide-6 cryogenic-trap-leak focus from 1 Hz trend points\n"
        "pre-Christmas lock_state >= 134; post-Christmas proxy mask: lock_state >= 100",
        fontsize=12,
        y=0.995,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return Image.open(path).size


def patch_file(target: Path, template_img: Image.Image, template_size: tuple[int, int]) -> None:
    im = Image.open(target).convert("RGBA")
    if not target.with_name(target.stem + "_before_legendboost.png").exists():
        target.with_name(target.stem + "_before_legendboost.png").write_bytes(target.read_bytes())
    scaled_template = template_img.resize(im.size, Image.Resampling.BICUBIC)
    sx = im.size[0] / template_size[0]
    sy = im.size[1] / template_size[1]
    for x0, y0, x1, y1 in TEMPLATE_LEGEND_BOXES:
        box = (
            int(round(x0 * sx)),
            int(round(y0 * sy)),
            int(round(x1 * sx)),
            int(round(y1 * sy)),
        )
        im.paste(scaled_template.crop(box), box)
    im.save(target)
    print(f"Updated {target}")


def main() -> None:
    template_size = build_template(TEMPLATE)
    template_img = Image.open(TEMPLATE).convert("RGBA")
    for name in [
        "slide6_cryotrap_focus_proxy_1hz.png",
        "slide6_cryotrap_focus_strict_1hz.png",
    ]:
        patch_file(WORKDIR / name, template_img, template_size)


if __name__ == "__main__":
    main()
