from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Optional

from gwpy.io.gwf import get_channel_names
from gwpy.timeseries import TimeSeries, TimeSeriesDict


def normalize_frame_paths(frame_paths: Iterable[str]) -> list[str]:
    normalized = [path for path in frame_paths if path]
    if not normalized:
        raise ValueError("at least one frame path is required")
    return normalized


def list_channels(frame_paths: Iterable[str], pattern: Optional[str] = None) -> list[str]:
    compiled = re.compile(pattern) if pattern else None
    channels: set[str] = set()
    for frame_path in normalize_frame_paths(frame_paths):
        for channel in get_channel_names(frame_path):
            if compiled and not compiled.search(channel):
                continue
            channels.add(channel)
    return sorted(channels)


def read_timeseries(
    frame_paths: Iterable[str],
    *,
    channel: str,
    gps_start: Optional[float] = None,
    gps_end: Optional[float] = None,
) -> TimeSeries:
    return TimeSeries.read(
        normalize_frame_paths(frame_paths),
        channel,
        start=gps_start,
        end=gps_end,
    )


def read_timeseries_dict(
    frame_paths: Iterable[str],
    *,
    channels: Iterable[str],
    gps_start: Optional[float] = None,
    gps_end: Optional[float] = None,
) -> TimeSeriesDict:
    return TimeSeriesDict.read(
        normalize_frame_paths(frame_paths),
        list(channels),
        start=gps_start,
        end=gps_end,
    )
