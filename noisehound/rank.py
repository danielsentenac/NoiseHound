from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Optional, Union

import numpy as np
import pandas as pd

from .frames import list_channels, read_timeseries_dict


TRIGGER_COLUMNS = ("gps_peak", "gps", "time", "trigger_gps")
FRAME_NAME_PATTERN = re.compile(r"-(\d+)-(\d+)\.gwf$")


@dataclass(frozen=True)
class EventScore:
    peak_abs_z: float
    signed_peak_z: float
    lag_s: float


@dataclass(frozen=True)
class FrameSpan:
    path: str
    start: float
    stop: float


def load_trigger_times(trigger_path: str) -> np.ndarray:
    frame = pd.read_csv(trigger_path)
    for column in TRIGGER_COLUMNS:
        if column in frame.columns:
            values = frame[column].dropna().to_numpy(dtype=float)
            if len(values):
                return np.sort(values)
    raise ValueError(f"could not find a trigger time column in {trigger_path}")


def load_candidate_channels(
    frame_paths: Iterable[str],
    *,
    channel_file: Optional[str],
    include_pattern: Optional[str],
    exclude_pattern: Optional[str],
    limit: Optional[int],
) -> list[str]:
    if channel_file:
        channels = []
        for line in Path(channel_file).read_text().splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            channels.append(stripped)
    else:
        channels = list_channels(frame_paths, pattern=include_pattern)

    if exclude_pattern:
        excluded = pd.Series(channels, dtype="string").str.contains(exclude_pattern, regex=True, na=False)
        channels = [channel for channel, reject in zip(channels, excluded) if not reject]
    if limit:
        channels = channels[:limit]
    return channels


def build_frame_index(frame_paths: Iterable[str]) -> list[FrameSpan]:
    spans: list[FrameSpan] = []
    for frame_path in frame_paths:
        match = FRAME_NAME_PATTERN.search(Path(frame_path).name)
        if not match:
            raise ValueError(f"could not parse frame start/duration from {frame_path}")
        start = float(match.group(1))
        duration = float(match.group(2))
        spans.append(FrameSpan(path=frame_path, start=start, stop=start + duration))
    return sorted(spans, key=lambda item: item.start)


def select_frame_paths(frame_index: Iterable[FrameSpan], start: float, stop: float) -> list[str]:
    return [item.path for item in frame_index if item.stop > start and item.start < stop]


def score_single_window(
    *,
    times: np.ndarray,
    values: np.ndarray,
    trigger: float,
    window_before_s: float,
    search_before_s: float,
    search_after_s: float,
    min_baseline_samples: int = 16,
    floor: float = 1e-12,
) -> Optional[EventScore]:
    baseline_mask = (times >= trigger - window_before_s) & (times < trigger - search_before_s)
    search_mask = (times >= trigger - search_before_s) & (times <= trigger + search_after_s)
    baseline = values[baseline_mask]
    search = values[search_mask]
    search_times = times[search_mask]
    if len(baseline) < min_baseline_samples or len(search) == 0:
        return None

    center = float(np.nanmedian(baseline))
    scale = float(1.4826 * np.nanmedian(np.abs(baseline - center)))
    if not np.isfinite(scale) or scale <= floor:
        scale = float(np.nanstd(baseline))
    scale = max(scale, floor)

    z = (search - center) / scale
    peak_index = int(np.nanargmax(np.abs(z)))
    return EventScore(
        peak_abs_z=float(abs(z[peak_index])),
        signed_peak_z=float(z[peak_index]),
        lag_s=float(search_times[peak_index] - trigger),
    )


def summarize_channel(
    channel: str,
    true_scores: list[EventScore],
    control_scores: list[EventScore],
    *,
    hit_threshold: float,
) -> dict[str, Union[float, str]]:
    true_peaks = np.asarray([s.peak_abs_z for s in true_scores], dtype=float)
    true_lags = np.asarray([s.lag_s for s in true_scores], dtype=float)
    control_peaks = np.asarray([s.peak_abs_z for s in control_scores], dtype=float)
    control_median = float(np.nanmedian(control_peaks)) if len(control_peaks) else 0.0
    lag_iqr = float(np.subtract(*np.quantile(true_lags, [0.75, 0.25]))) if len(true_lags) else np.nan
    median_true = float(np.nanmedian(true_peaks)) if len(true_peaks) else 0.0
    hit_fraction = float(np.mean(true_peaks >= hit_threshold)) if len(true_peaks) else 0.0
    enrichment = median_true / max(control_median, 1.0)
    rank_score = hit_fraction * np.log1p(median_true) * enrichment / (1.0 + (lag_iqr if np.isfinite(lag_iqr) else 10.0))
    return {
        "channel": channel,
        "events_used": int(len(true_peaks)),
        "controls_used": int(len(control_peaks)),
        "median_peak_z": median_true,
        "p95_peak_z": float(np.quantile(true_peaks, 0.95)) if len(true_peaks) else 0.0,
        "median_control_peak_z": control_median,
        "hit_fraction": hit_fraction,
        "lag_median_s": float(np.nanmedian(true_lags)) if len(true_lags) else np.nan,
        "lag_iqr_s": lag_iqr,
        "rank_score": rank_score,
    }


def rank_channels_against_triggers(
    frame_paths: list[str],
    *,
    trigger_times: np.ndarray,
    candidate_channels: list[str],
    gps_start: Optional[float],
    gps_end: Optional[float],
    window_before_s: float,
    search_before_s: float,
    search_after_s: float,
    control_offset_s: float,
    hit_threshold: float,
    max_sample_rate: Optional[float],
    checkpoint_path: Optional[str] = None,
) -> pd.DataFrame:
    """Read each frame file once, score all channels against all triggers."""
    if len(trigger_times) == 0:
        raise ValueError("no trigger times were provided")

    control_times = trigger_times + control_offset_s
    if gps_start is not None or gps_end is not None:
        lower = float(gps_start) if gps_start is not None else -np.inf
        upper = float(gps_end) if gps_end is not None else np.inf
        trigger_times = trigger_times[(trigger_times >= lower) & (trigger_times <= upper)]
        control_times = control_times[(control_times >= lower) & (control_times <= upper)]

    frame_index = build_frame_index(frame_paths)
    all_times = np.concatenate([trigger_times, control_times])
    read_start = float(np.min(all_times) - window_before_s)
    read_end = float(np.max(all_times) + search_after_s)
    interval_paths = select_frame_paths(frame_index, read_start, read_end)

    if not interval_paths:
        raise RuntimeError("no frame files cover the required time range")

    channel_set = set(candidate_channels)
    true_scores: dict[str, list[EventScore]] = {ch: [] for ch in candidate_channels}
    control_scores: dict[str, list[EventScore]] = {ch: [] for ch in candidate_channels}
    errors: dict[str, str] = {}

    n_frames = len(interval_paths)
    interval_path_set = set(interval_paths)
    print(f"[rank] {len(candidate_channels)} channels, {len(trigger_times)} triggers, {n_frames} frames to read")

    # Read each frame once, score all channels against all overlapping triggers.
    for frame_idx, frame_span in enumerate(frame_index):
        if frame_span.path not in interval_path_set:
            continue

        relevant_true = [
            t for t in trigger_times
            if frame_span.stop > t - window_before_s and frame_span.start < t + search_after_s
        ]
        relevant_control = [
            t for t in control_times
            if frame_span.stop > t - window_before_s and frame_span.start < t + search_after_s
        ]
        if not relevant_true and not relevant_control:
            continue

        active_channels = [ch for ch in candidate_channels if ch not in errors]
        print(f"[rank] frame {frame_idx + 1}/{n_frames}: {Path(frame_span.path).name} "
              f"({len(relevant_true)} triggers, {len(relevant_control)} controls, {len(active_channels)} channels)", flush=True)

        try:
            series_dict = read_timeseries_dict(
                [frame_span.path],
                channels=active_channels,
                gps_start=None,
                gps_end=None,
            )
        except Exception as e:
            print(f"[rank] batch read failed ({e}), falling back to batches of 64", flush=True)
            series_dict = {}
            batch_size = 64
            for i in range(0, len(active_channels), batch_size):
                batch = active_channels[i:i + batch_size]
                try:
                    series_dict.update(read_timeseries_dict(
                        [frame_span.path], channels=batch, gps_start=None, gps_end=None,
                    ))
                except Exception:
                    from .frames import read_timeseries
                    for ch in batch:
                        try:
                            series_dict[ch] = read_timeseries([frame_span.path], channel=ch, gps_start=None, gps_end=None)
                        except Exception as ch_e:
                            errors[ch] = str(ch_e)

        print(f"[rank]   read {len(series_dict)} channels from frame", flush=True)

        for ch_idx, (channel, series) in enumerate(series_dict.items()):
            if channel not in channel_set:
                continue
            if max_sample_rate:
                current = float(series.sample_rate.value)
                if current > max_sample_rate:
                    series = series.resample(max_sample_rate)
            times = np.asarray(series.times.value, dtype=float)
            values = np.asarray(series.value, dtype=float)

            for trigger in relevant_true:
                score = score_single_window(
                    times=times, values=values, trigger=trigger,
                    window_before_s=window_before_s,
                    search_before_s=search_before_s,
                    search_after_s=search_after_s,
                )
                if score is not None:
                    true_scores[channel].append(score)

            for trigger in relevant_control:
                score = score_single_window(
                    times=times, values=values, trigger=trigger,
                    window_before_s=window_before_s,
                    search_before_s=search_before_s,
                    search_after_s=search_after_s,
                )
                if score is not None:
                    control_scores[channel].append(score)

            if (ch_idx + 1) % 500 == 0:
                print(f"[rank]   scored {ch_idx + 1}/{len(series_dict)} channels", flush=True)

        print(f"[rank] frame {frame_idx + 1}/{n_frames} done", flush=True)

    rows = []
    for channel in candidate_channels:
        if channel in errors:
            rows.append({
                "channel": channel, "events_used": 0, "controls_used": 0,
                "median_peak_z": np.nan, "p95_peak_z": np.nan,
                "median_control_peak_z": np.nan, "hit_fraction": np.nan,
                "lag_median_s": np.nan, "lag_iqr_s": np.nan,
                "rank_score": np.nan, "error": errors[channel],
            })
        else:
            rows.append(summarize_channel(
                channel, true_scores[channel], control_scores[channel],
                hit_threshold=hit_threshold,
            ))

    if checkpoint_path:
        pd.DataFrame(rows).sort_values("rank_score", ascending=False, na_position="last").to_csv(checkpoint_path, index=False)

    result = pd.DataFrame(rows)
    return result.sort_values("rank_score", ascending=False, na_position="last").reset_index(drop=True)
