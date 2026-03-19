from __future__ import annotations

import os
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional

from gwdatafind.ui import find_times as ui_find_times
from gwdatafind.ui import find_types as ui_find_types
from gwdatafind.ui import find_urls as ui_find_urls


@dataclass(frozen=True)
class Segment:
    start: float
    stop: float
    duration: float


@contextmanager
def with_token(token_file: Optional[str]):
    previous_file = os.environ.get("BEARER_TOKEN_FILE")
    previous_token = os.environ.get("BEARER_TOKEN")
    if token_file:
        os.environ["BEARER_TOKEN_FILE"] = token_file
        # Ensure an unrelated inline bearer token doesn't override the explicit file.
        os.environ.pop("BEARER_TOKEN", None)
    try:
        yield
    finally:
        if token_file:
            if previous_file is None:
                os.environ.pop("BEARER_TOKEN_FILE", None)
            else:
                os.environ["BEARER_TOKEN_FILE"] = previous_file
            if previous_token is None:
                os.environ.pop("BEARER_TOKEN", None)
            else:
                os.environ["BEARER_TOKEN"] = previous_token


def list_types(
    *,
    server: Optional[str],
    observatory: Optional[str],
    token_file: Optional[str],
) -> list[str]:
    with with_token(token_file):
        frame_types = ui_find_types(site=observatory, host=server)
    return sorted(frame_types)


def list_segments(
    *,
    server: Optional[str],
    observatory: str,
    frame_type: str,
    gps_start: int,
    gps_end: int,
    token_file: Optional[str],
) -> list[Segment]:
    with with_token(token_file):
        segments = ui_find_times(
            observatory,
            frame_type,
            gpsstart=gps_start,
            gpsend=gps_end,
            host=server,
        )
    return [Segment(start=float(seg[0]), stop=float(seg[1]), duration=float(abs(seg))) for seg in segments]


def find_urls(
    *,
    server: Optional[str],
    observatory: str,
    frame_type: str,
    gps_start: int,
    gps_end: int,
    url_type: str,
    token_file: Optional[str],
    match: Optional[str] = None,
) -> list[str]:
    with with_token(token_file):
        urls = ui_find_urls(
            observatory,
            frame_type,
            gpsstart=gps_start,
            gpsend=gps_end,
            match=match,
            urltype=url_type,
            host=server,
        )
    return list(urls)
