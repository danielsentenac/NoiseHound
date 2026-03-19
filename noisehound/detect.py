from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import find_peaks, peak_widths

from .frames import read_timeseries


@dataclass
class DetectionResult:
    events: pd.DataFrame
    times: np.ndarray
    frequencies: Optional[np.ndarray]
    relative_excess: Optional[np.ndarray]
    broadband_score: np.ndarray
    robust_zscore: np.ndarray
    threshold: float


def summarize_relative_excess(relative_excess: np.ndarray, statistic: str) -> np.ndarray:
    if statistic == "median":
        return np.nanmedian(relative_excess, axis=1)
    if statistic == "p95":
        return np.nanpercentile(relative_excess, 95, axis=1)
    if statistic == "p98":
        return np.nanpercentile(relative_excess, 98, axis=1)
    if statistic == "top5pct_mean":
        cutoff = max(0, int(relative_excess.shape[1] * 0.95))
        return np.nanmean(np.sort(relative_excess, axis=1)[:, cutoff:], axis=1)
    if statistic == "top2pct_mean":
        cutoff = max(0, int(relative_excess.shape[1] * 0.98))
        return np.nanmean(np.sort(relative_excess, axis=1)[:, cutoff:], axis=1)
    raise ValueError(f"unsupported summary statistic: {statistic}")


def rolling_median(values: np.ndarray, window_bins: int) -> np.ndarray:
    series = pd.Series(values)
    return (
        series.rolling(window=window_bins, center=True, min_periods=max(5, window_bins // 4))
        .median()
        .bfill()
        .ffill()
        .to_numpy()
    )


def rolling_robust_zscore(values: np.ndarray, window_bins: int, floor: float = 1e-12) -> np.ndarray:
    baseline = rolling_median(values, window_bins)
    residual = np.abs(values - baseline)
    scale = (
        pd.Series(residual)
        .rolling(window=window_bins, center=True, min_periods=max(5, window_bins // 4))
        .median()
        .bfill()
        .ffill()
        .to_numpy()
    )
    return (values - baseline) / np.maximum(1.4826 * scale, floor)


def broadband_excess_from_matrix(
    matrix: np.ndarray,
    *,
    statistic: str = "p95",
    floor: Optional[float] = None,
) -> tuple[np.ndarray, np.ndarray]:
    finite_positive = matrix[np.isfinite(matrix) & (matrix > 0)]
    if floor is None:
        if finite_positive.size:
            floor = max(float(np.nanmin(finite_positive)) * 0.1, np.finfo(float).tiny)
        else:
            floor = np.finfo(float).tiny
    reference = np.nanmedian(matrix, axis=0)
    relative_excess = np.log10(np.maximum(matrix, floor) / np.maximum(reference, floor))
    broadband_score = summarize_relative_excess(relative_excess, statistic)
    return relative_excess, broadband_score


def detect_broadband_glitches(
    frame_paths: list[str],
    *,
    channel: str,
    gps_start: Optional[int],
    gps_end: Optional[int],
    fmin: float,
    fmax: float,
    stride: float,
    fftlength: float,
    overlap: float,
    summary_stat: str,
    threshold: float,
    background_window_s: float,
    min_separation_s: float,
    keep_spectrogram: bool = False,
) -> DetectionResult:
    time_blocks: list[np.ndarray] = []
    score_blocks: list[np.ndarray] = []
    matrix_blocks: list[np.ndarray] = []
    frequencies: Optional[np.ndarray] = None
    last_error: Optional[Exception] = None

    for frame_path in frame_paths:
        try:
            target = read_timeseries([frame_path], channel=channel, gps_start=None, gps_end=None)
            block_start = float(target.span[0])
            block_end = float(target.span[1])
            if gps_start is not None:
                block_start = max(block_start, float(gps_start))
            if gps_end is not None:
                block_end = min(block_end, float(gps_end))
            if block_end <= block_start:
                continue
            target = target.crop(block_start, block_end)
            spectrogram = target.spectrogram(
                stride=stride,
                fftlength=fftlength,
                overlap=overlap,
                method="median",
            ).crop_frequencies(fmin, fmax)
        except Exception as error:
            last_error = error
            continue

        matrix = np.asarray(spectrogram.value, dtype=float)
        if matrix.size == 0:
            continue

        block_times = np.asarray(spectrogram.times.value, dtype=float)
        block_frequencies = np.asarray(spectrogram.frequencies.value, dtype=float)
        if frequencies is None:
            frequencies = block_frequencies
        elif not np.allclose(frequencies, block_frequencies):
            raise ValueError("frequency grid changed between frame files")

        relative_excess, broadband_score = broadband_excess_from_matrix(matrix, statistic=summary_stat)
        time_blocks.append(block_times)
        score_blocks.append(broadband_score.astype(np.float32))
        if keep_spectrogram:
            matrix_blocks.append(relative_excess.astype(np.float32))

    if not time_blocks:
        raise RuntimeError(str(last_error) if last_error else "no usable spectrogram blocks were produced")

    times = np.concatenate(time_blocks)
    broadband_score = np.concatenate(score_blocks)
    order = np.argsort(times)
    times = times[order]
    broadband_score = broadband_score[order]
    unique_mask = np.concatenate(([True], np.diff(times) > 0))
    times = times[unique_mask]
    broadband_score = broadband_score[unique_mask]

    relative_excess = None
    if keep_spectrogram:
        relative_excess = np.vstack(matrix_blocks)[order][unique_mask]

    background_window_bins = max(5, int(round(background_window_s / stride)))
    zscore = rolling_robust_zscore(broadband_score, background_window_bins)
    distance_bins = max(1, int(round(min_separation_s / stride)))
    peaks, properties = find_peaks(zscore, height=threshold, distance=distance_bins)

    if len(peaks):
        widths, _, left_ips, right_ips = peak_widths(zscore, peaks, rel_height=0.5)
        rows = []
        for index, peak in enumerate(peaks):
            rows.append(
                {
                    "gps_peak": times[peak],
                    "score_z": float(properties["peak_heights"][index]),
                    "broadband_score": float(broadband_score[peak]),
                    "width_s": float(widths[index] * stride),
                    "start_gps": float(np.interp(left_ips[index], np.arange(len(times)), times)),
                    "stop_gps": float(np.interp(right_ips[index], np.arange(len(times)), times)),
                }
            )
        events = pd.DataFrame(rows).sort_values("gps_peak").reset_index(drop=True)
        if len(events) > 1:
            events["delta_prev_s"] = events["gps_peak"].diff()
    else:
        events = pd.DataFrame(
            columns=["gps_peak", "score_z", "broadband_score", "width_s", "start_gps", "stop_gps"]
        )

    return DetectionResult(
        events=events,
        times=times,
        frequencies=frequencies,
        relative_excess=relative_excess,
        broadband_score=broadband_score,
        robust_zscore=zscore,
        threshold=threshold,
    )


def save_detection_plot(result: DetectionResult, output_path: str, *, title: str) -> None:
    if result.relative_excess is None or result.frequencies is None:
        raise ValueError("spectrogram data were not retained; rerun detection with plotting enabled")

    figure = plt.figure(figsize=(14, 8))
    grid = figure.add_gridspec(
        2,
        2,
        height_ratios=[3, 1],
        width_ratios=[40, 2],
        hspace=0.10,
        wspace=0.10,
    )
    ax_spec = figure.add_subplot(grid[0, 0])
    ax_score = figure.add_subplot(grid[1, 0], sharex=ax_spec)
    ax_color = figure.add_subplot(grid[0, 1])

    image = ax_spec.imshow(
        result.relative_excess.T,
        aspect="auto",
        origin="lower",
        extent=[
            result.times[0],
            result.times[-1],
            result.frequencies[0],
            result.frequencies[-1],
        ],
        cmap="jet",
        vmin=-1.5,
        vmax=2.0,
    )
    ax_spec.set_ylabel("Frequency [Hz]")
    ax_spec.set_title(title)
    ax_spec.grid(False)
    figure.colorbar(image, cax=ax_color, label="log10(power / median)")
    ax_spec.tick_params(axis="x", labelbottom=False)

    ax_score.plot(result.times, result.robust_zscore, color="tab:blue", linewidth=1.2)
    ax_score.axhline(result.threshold, color="tab:red", linestyle="--", linewidth=1.0)
    if not result.events.empty:
        for event in result.events.itertuples(index=False):
            start_gps = getattr(event, "start_gps", None)
            stop_gps = getattr(event, "stop_gps", None)
            gps_peak = getattr(event, "gps_peak", None)
            if start_gps is not None and stop_gps is not None:
                ax_score.axvspan(start_gps, stop_gps, color="tab:red", alpha=0.08, linewidth=0)
            if gps_peak is not None:
                ax_score.axvline(gps_peak, color="tab:red", alpha=0.25, linewidth=0.8)
    ax_score.set_xlabel("GPS time [s]")
    ax_score.set_ylabel("MAD-normalized excess")
    ax_score.grid(alpha=0.2)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
