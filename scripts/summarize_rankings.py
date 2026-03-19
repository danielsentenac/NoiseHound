from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def classify_family(channel: str) -> str:
    if channel.startswith("V1:VAC_LN2"):
        return "LN2"
    if channel.startswith("V1:HVAC_"):
        return "HVAC"
    if channel.startswith("V1:INF_"):
        return "INF"
    if channel.startswith("V1:ENV_"):
        return "ENV"
    if channel.startswith("V1:LSC_Etalon_"):
        return "LSC_Etalon"
    if channel.startswith("V1:LSC_"):
        return "LSC"
    if channel.startswith("V1:ASC_"):
        return "ASC"
    if channel.startswith("V1:SAT_"):
        return "SAT"
    if channel.startswith("V1:Sa_"):
        return "Sa"
    if "_LC_X" in channel:
        return "LocalControlX"
    if "TCS" in channel:
        return "TCS"
    if "SBE" in channel:
        return "SBE"
    if "SUSP" in channel:
        return "SUSP"
    if "TEMP" in channel or "THERM" in channel:
        return "Thermal"
    return "Other"


def load_rankings(result_dirs: list[Path]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for result_dir in result_dirs:
        ranking_path = result_dir / "ranking.csv"
        if not ranking_path.exists():
            continue
        frame = pd.read_csv(ranking_path)
        if frame.empty:
            continue
        frame["result_dir"] = str(result_dir)
        frame["family"] = frame["channel"].astype(str).map(classify_family)
        rows.append(frame)
    if not rows:
        return pd.DataFrame()
    merged = pd.concat(rows, ignore_index=True)
    return merged.sort_values("rank_score", ascending=False, na_position="last").reset_index(drop=True)


def write_markdown_summary(frame: pd.DataFrame, summary_path: Path) -> None:
    lines: list[str] = ["# January 1, 2025 raw ranking summary", ""]
    if frame.empty:
        lines.append("No ranking outputs were available.")
        summary_path.write_text("\n".join(lines) + "\n")
        return

    available_dirs = frame["result_dir"].drop_duplicates().tolist()
    lines.append("## Available ranking outputs")
    for result_dir in available_dirs:
        lines.append(f"- {result_dir}")
    lines.append("")

    top_overall = frame[["channel", "family", "rank_score", "median_peak_z", "hit_fraction", "lag_median_s", "result_dir"]]
    lines.append("## Top overall")
    lines.append("```text")
    lines.append(top_overall.head(25).to_string(index=False))
    lines.append("```")
    lines.append("")

    lines.append("## Top by family")
    family_top = (
        frame.sort_values("rank_score", ascending=False, na_position="last")
        .groupby("family", as_index=False)
        .head(5)[["family", "channel", "rank_score", "median_peak_z", "hit_fraction", "lag_median_s", "result_dir"]]
    )
    lines.append("```text")
    lines.append(family_top.to_string(index=False))
    lines.append("```")
    lines.append("")

    for result_dir, group in frame.groupby("result_dir", sort=False):
        lines.append(f"## {result_dir}")
        view = group[["channel", "family", "rank_score", "median_peak_z", "hit_fraction", "lag_median_s"]]
        lines.append("```text")
        lines.append(view.head(20).to_string(index=False))
        lines.append("```")
        lines.append("")

    summary_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize one or more ranking outputs.")
    parser.add_argument("--summary-path", required=True)
    parser.add_argument("--merged-csv")
    parser.add_argument("result_dirs", nargs="+")
    args = parser.parse_args()

    result_dirs = [Path(item) for item in args.result_dirs]
    summary_path = Path(args.summary_path)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    merged = load_rankings(result_dirs)
    write_markdown_summary(merged, summary_path)

    if args.merged_csv:
        merged_path = Path(args.merged_csv)
        merged_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(merged_path, index=False)

    print(f"wrote summary to {summary_path}")
    if args.merged_csv:
        print(f"wrote merged csv to {args.merged_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
