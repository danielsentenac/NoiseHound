from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

from .datafind import find_urls, list_segments, list_types
from .detect import detect_broadband_glitches, save_detection_plot
from .frames import list_channels
from .rank import load_candidate_channels, load_trigger_times, rank_channels_against_triggers


def add_datafind_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--server",
        default=None,
        help="GWDataFind server. If omitted, use GWDATAFIND_SERVER/LIGO_DATAFIND_SERVER from the environment.",
    )
    parser.add_argument("--observatory", default="V")
    parser.add_argument("--token-file")


def command_types(arguments: argparse.Namespace) -> int:
    for frame_type in list_types(
        server=arguments.server,
        observatory=arguments.observatory,
        token_file=arguments.token_file,
    ):
        print(frame_type)
    return 0


def command_segments(arguments: argparse.Namespace) -> int:
    rows = list_segments(
        server=arguments.server,
        observatory=arguments.observatory,
        frame_type=arguments.frame_type,
        gps_start=arguments.gps_start,
        gps_end=arguments.gps_end,
        token_file=arguments.token_file,
    )
    print("start,stop,duration")
    for row in rows:
        print(f"{row.start},{row.stop},{row.duration}")
    return 0


def command_urls(arguments: argparse.Namespace) -> int:
    for url in find_urls(
        server=arguments.server,
        observatory=arguments.observatory,
        frame_type=arguments.frame_type,
        gps_start=arguments.gps_start,
        gps_end=arguments.gps_end,
        url_type=arguments.url_type,
        token_file=arguments.token_file,
        match=arguments.match,
    ):
        print(url)
    return 0


def command_channels(arguments: argparse.Namespace) -> int:
    channels = list_channels(arguments.frames, pattern=arguments.pattern)
    if arguments.output:
        Path(arguments.output).write_text("\n".join(channels) + "\n")
    for channel in channels:
        print(channel)
    return 0


def command_detect(arguments: argparse.Namespace) -> int:
    result = detect_broadband_glitches(
        arguments.frames,
        channel=arguments.channel,
        gps_start=arguments.gps_start,
        gps_end=arguments.gps_end,
        fmin=arguments.fmin,
        fmax=arguments.fmax,
        stride=arguments.stride,
        fftlength=arguments.fftlength,
        overlap=arguments.overlap,
        summary_stat=arguments.summary_stat,
        threshold=arguments.threshold,
        background_window_s=arguments.background_window_s,
        min_separation_s=arguments.min_separation_s,
        keep_spectrogram=bool(arguments.plot),
    )
    result.events.to_csv(arguments.output, index=False)
    if arguments.plot:
        save_detection_plot(result, arguments.plot, title=arguments.channel)

    print(f"wrote {len(result.events)} triggers to {arguments.output}")
    if len(result.events) > 1 and "delta_prev_s" in result.events.columns:
        median_spacing = result.events["delta_prev_s"].dropna().median()
        print(f"median spacing: {median_spacing:.1f} s")
    return 0


def command_rank(arguments: argparse.Namespace) -> int:
    trigger_times = load_trigger_times(arguments.triggers)
    channels = load_candidate_channels(
        arguments.frames,
        channel_file=arguments.channel_file,
        include_pattern=arguments.include,
        exclude_pattern=arguments.exclude,
        limit=arguments.limit,
    )
    ranking = rank_channels_against_triggers(
        arguments.frames,
        trigger_times=trigger_times,
        candidate_channels=channels,
        gps_start=arguments.gps_start,
        gps_end=arguments.gps_end,
        window_before_s=arguments.window_before_s,
        search_before_s=arguments.search_before_s,
        search_after_s=arguments.search_after_s,
        control_offset_s=arguments.control_offset_s,
        hit_threshold=arguments.hit_threshold,
        max_sample_rate=arguments.max_sample_rate,
        checkpoint_path=arguments.output,
    )
    ranking.to_csv(arguments.output, index=False)
    print(f"wrote ranking for {len(ranking)} channels to {arguments.output}")
    if len(ranking):
        preview = ranking.head(arguments.top)
        print(preview.to_string(index=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    default_server = os.environ.get("GWDATAFIND_SERVER") or os.environ.get("LIGO_DATAFIND_SERVER")
    parser = argparse.ArgumentParser(
        prog="noisehound",
        description="Virgo glitch triage toolkit built around broadband triggers and event-synchronous channel ranking.",
    )
    if default_server:
        parser.epilog = f"Default datafind server from environment: {default_server}"
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_types = subparsers.add_parser("types", help="list frame types available through gw_data_find")
    add_datafind_arguments(parser_types)
    parser_types.set_defaults(func=command_types)

    parser_segments = subparsers.add_parser("segments", help="show availability segments for a frame type")
    add_datafind_arguments(parser_segments)
    parser_segments.add_argument("--frame-type", required=True)
    parser_segments.add_argument("--gps-start", required=True, type=int)
    parser_segments.add_argument("--gps-end", required=True, type=int)
    parser_segments.set_defaults(func=command_segments)

    parser_urls = subparsers.add_parser("urls", help="resolve data URLs for a frame type")
    add_datafind_arguments(parser_urls)
    parser_urls.add_argument("--frame-type", required=True)
    parser_urls.add_argument("--gps-start", required=True, type=int)
    parser_urls.add_argument("--gps-end", required=True, type=int)
    parser_urls.add_argument("--url-type", default="osdf")
    parser_urls.add_argument("--match")
    parser_urls.set_defaults(func=command_urls)

    parser_channels = subparsers.add_parser("channels", help="list channels available in local frame files")
    parser_channels.add_argument("frames", nargs="+")
    parser_channels.add_argument("--pattern")
    parser_channels.add_argument("--output")
    parser_channels.set_defaults(func=command_channels)

    parser_detect = subparsers.add_parser("detect", help="find broadband glitches in a target strain channel")
    parser_detect.add_argument("frames", nargs="+")
    parser_detect.add_argument("--channel", default="V1:Hrec_hoft_20000Hz")
    parser_detect.add_argument("--gps-start", type=int)
    parser_detect.add_argument("--gps-end", type=int)
    parser_detect.add_argument("--fmin", type=float, default=20.0)
    parser_detect.add_argument("--fmax", type=float, default=210.0)
    parser_detect.add_argument("--stride", type=float, default=4.0)
    parser_detect.add_argument("--fftlength", type=float, default=4.0)
    parser_detect.add_argument("--overlap", type=float, default=2.0)
    parser_detect.add_argument(
        "--summary-stat",
        default="p95",
        choices=["median", "p95", "p98", "top5pct_mean", "top2pct_mean"],
        help="How to collapse the spectrogram across frequency bins before peak finding.",
    )
    parser_detect.add_argument("--threshold", type=float, default=6.0)
    parser_detect.add_argument("--background-window-s", type=float, default=1800.0)
    parser_detect.add_argument("--min-separation-s", type=float, default=600.0)
    parser_detect.add_argument("--output", default="triggers.csv")
    parser_detect.add_argument("--plot")
    parser_detect.set_defaults(func=command_detect)

    parser_rank = subparsers.add_parser("rank", help="rank auxiliary channels against glitch trigger times")
    parser_rank.add_argument("frames", nargs="+")
    parser_rank.add_argument("--triggers", required=True)
    parser_rank.add_argument("--channel-file")
    parser_rank.add_argument("--include")
    parser_rank.add_argument("--exclude")
    parser_rank.add_argument("--limit", type=int)
    parser_rank.add_argument("--gps-start", type=float)
    parser_rank.add_argument("--gps-end", type=float)
    parser_rank.add_argument("--window-before-s", type=float, default=60.0)
    parser_rank.add_argument("--search-before-s", type=float, default=5.0)
    parser_rank.add_argument("--search-after-s", type=float, default=5.0)
    parser_rank.add_argument("--control-offset-s", type=float, default=600.0)
    parser_rank.add_argument("--hit-threshold", type=float, default=6.0)
    parser_rank.add_argument("--max-sample-rate", type=float)
    parser_rank.add_argument("--top", type=int, default=20)
    parser_rank.add_argument("--output", default="ranking.csv")
    parser_rank.set_defaults(func=command_rank)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        return int(arguments.func(arguments))
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
