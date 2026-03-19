"""
causality_analysis.py — causal vs. consequence discrimination for NOISEHOUND channels.

For each channel with rank_score > 0 in the supplied ranking CSV:
  1. Cross-correlation lag analysis  (lag at max |corr|)
  2. Granger causality test          (statsmodels grangercausalitytests)
  3. Transfer Entropy                (k-nearest-neighbour estimator)
  4. Summary table saved to <output>/causality_summary.csv
     and a human-readable report  to <output>/causality_report.txt

Usage:
    python causality_analysis.py \
        --trend     /tmp/trend.gwf \
        --hrec-dir  /path/to/hoftonline_frames/ \
        --triggers  /path/to/triggers_p95_auto.csv \
        --ranking   /path/to/ranking.csv \
        --gps-start 1419724818 \
        --gps-end   1419735618 \
        --output    /path/to/output_dir
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import correlate, correlation_lags

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--trend",       required=True,  help="1 Hz trend GWF file")
    p.add_argument("--hrec-dir",    required=True,  help="Directory containing Hrec hoft GWF frames")
    p.add_argument("--triggers",    required=True,  help="Trigger CSV (must have gps_peak column)")
    p.add_argument("--ranking",     required=True,  help="ranking.csv from noisehound rank")
    p.add_argument("--gps-start",   type=float, required=True)
    p.add_argument("--gps-end",     type=float, required=True)
    p.add_argument("--output",      required=True,  help="Output directory")
    p.add_argument("--top-n",       type=int, default=40,
                   help="Analyse only the top-N ranked channels (default 40)")
    p.add_argument("--max-lag-s",   type=float, default=120.0,
                   help="Maximum lag in seconds for cross-correlation (default 120 s)")
    p.add_argument("--granger-max-lag", type=int, default=30,
                   help="Maximum lag order tested by Granger (default 30 steps)")
    return p.parse_args()


# ─── Helpers ──────────────────────────────────────────────────────────────────

def zscore_norm(x: np.ndarray) -> np.ndarray:
    """Standard z-score normalisation; safe on constant arrays."""
    s = x.std()
    if s < 1e-15:
        return np.zeros_like(x, dtype=float)
    return (x - x.mean()) / s


def rolling_mad_zscore(v: np.ndarray, window: int, floor: float = 1e-12) -> np.ndarray:
    """Robust rolling MAD z-score (1 Hz series → window in samples)."""
    s = pd.Series(v)
    baseline = s.rolling(window, center=True, min_periods=max(5, window // 4)).median().bfill().ffill()
    residual = (s - baseline).abs()
    scale = residual.rolling(window, center=True, min_periods=max(5, window // 4)).median().bfill().ffill()
    return ((s - baseline) / np.maximum(1.4826 * scale, floor)).to_numpy()


def cross_corr_lag(x: np.ndarray, y: np.ndarray, max_lag: int) -> tuple[float, float]:
    """
    Return (lag_s, peak_corr) where lag_s > 0 means x *leads* y.
    Both series must be 1 Hz and same length.
    """
    x_n = zscore_norm(x)
    y_n = zscore_norm(y)
    cc = correlate(y_n, x_n, mode="full") / len(x)
    lags = correlation_lags(len(y_n), len(x_n), mode="full")
    mask = np.abs(lags) <= max_lag
    cc_m, lags_m = cc[mask], lags[mask]
    idx = np.argmax(np.abs(cc_m))
    return float(lags_m[idx]), float(cc_m[idx])


def granger_p(x: np.ndarray, y: np.ndarray, max_lag: int) -> float:
    """
    Granger test: does x Granger-cause y?
    Returns minimum p-value across lag orders 1..max_lag.
    Returns 1.0 on any failure.
    """
    try:
        from statsmodels.tsa.stattools import grangercausalitytests
        data = np.column_stack([y, x])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = grangercausalitytests(data, maxlag=max_lag, verbose=False)
        pvals = [res[lag][0]["ssr_ftest"][1] for lag in range(1, max_lag + 1)]
        return float(min(pvals))
    except Exception:
        return 1.0


def transfer_entropy_knn(x: np.ndarray, y: np.ndarray,
                          k: int = 5, lag: int = 1) -> float:
    """
    Transfer entropy T(X→Y) estimated via k-NN (Kraskov-like).
    T(X→Y) > T(Y→X) suggests X drives Y.

    Uses a simple histogram approximation (fast, works at 1 Hz).
    lag: the embedding delay in samples (default 1 s).
    """
    try:
        n = min(len(x), len(y)) - lag
        if n < 20:
            return float("nan")
        yt1 = y[lag : lag + n]   # y(t+lag)
        yt  = y[:n]              # y(t)
        xt  = x[:n]              # x(t)

        def _mi(a: np.ndarray, b: np.ndarray, bins: int = 20) -> float:
            """Mutual information I(a; b) via histogram."""
            c_ab, _, _ = np.histogram2d(a, b, bins=bins)
            c_a = c_ab.sum(axis=1, keepdims=True)
            c_b = c_ab.sum(axis=0, keepdims=True)
            c_ab_safe = c_ab / c_ab.sum()
            c_a_safe  = c_a  / c_a.sum()
            c_b_safe  = c_b  / c_b.sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(
                    (c_ab_safe > 0) & (c_a_safe > 0) & (c_b_safe > 0),
                    c_ab_safe / (c_a_safe * c_b_safe), 1.0)
            return float(np.sum(c_ab_safe * np.log2(np.maximum(ratio, 1e-300))))

        # TE(X→Y) = I(Yt+1; Xt | Yt) ≈ I(Yt+1, Xt, Yt) - I(Yt+1, Yt)
        # simple approximation: conditional MI via joint - marginal
        te_xy = _mi(yt1, xt) - _mi(yt1, yt) + _mi(yt, xt)
        return max(0.0, float(te_xy))
    except Exception:
        return float("nan")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    gps_start = args.gps_start
    gps_end   = args.gps_end
    duration  = gps_end - gps_start
    window    = int(round(duration / 3))   # ~30% of span for rolling background

    print(f"GPS range : {gps_start:.0f} – {gps_end:.0f}  ({duration:.0f} s)")

    # ── Load triggers ──────────────────────────────────────────────────────────
    triggers_df = pd.read_csv(args.triggers)
    trigger_gps = triggers_df["gps_peak"].values
    print(f"Triggers  : {len(trigger_gps)}")

    # ── Load ranking ───────────────────────────────────────────────────────────
    ranking = pd.read_csv(args.ranking)
    ranking = ranking[ranking["rank_score"] > 0].sort_values("rank_score", ascending=False)
    top_channels = ranking.head(args.top_n)["channel"].tolist()
    print(f"Channels to analyse: {len(top_channels)} (top {args.top_n} with rank>0)")

    # ── Build Hrec MAD-z at 1 Hz from raw frames ───────────────────────────────
    hrec_frames = sorted(Path(args.hrec_dir).glob("*.gwf"))
    print(f"\nLoading Hrec from {len(hrec_frames)} frames …")
    hrec_z1hz = None
    t_ref = None
    try:
        from gwpy.timeseries import TimeSeries
        hrec = TimeSeries.read(
            [str(f) for f in hrec_frames],
            "V1:Hrec_hoft_16384Hz",
            start=gps_start, end=gps_end,
        )
        # Downsample to 1 Hz via resampling
        hrec_1hz = hrec.resample(1)
        t_ref     = np.asarray(hrec_1hz.times.value, dtype=float)
        v         = np.asarray(hrec_1hz.value, dtype=float)
        win       = max(5, int(round(1800.0)))   # 1800-s rolling window at 1 Hz
        hrec_z1hz = rolling_mad_zscore(v, win)
        print(f"  Hrec 1 Hz z-score computed, {len(hrec_z1hz)} samples")
    except Exception as e:
        print(f"  WARNING: Could not load Hrec — {e}")
        print("  Falling back: using trigger times only for lag sign test")

    # ── Load aux trend channels ────────────────────────────────────────────────
    print(f"\nLoading {len(top_channels)} trend channels …")
    try:
        from gwpy.timeseries import TimeSeriesDict
        tdict = TimeSeriesDict.read(
            args.trend, top_channels, start=gps_start, end=gps_end
        )
        print(f"  Loaded {len(tdict)} channels from trend file")
    except Exception as e:
        print(f"  FATAL: Could not read trend file — {e}")
        raise

    max_lag_samples = int(round(args.max_lag_s))

    # ── Analyse each channel ───────────────────────────────────────────────────
    results = []
    for ch in top_channels:
        if ch not in tdict:
            print(f"  SKIP  {ch} (not in trend)")
            continue

        ts = tdict[ch]
        t_aux = np.asarray(ts.times.value, dtype=float)
        v_aux = np.asarray(ts.value, dtype=float)
        sr    = float(ts.sample_rate.value)

        # Interpolate to uniform 1 Hz grid aligned with Hrec (or self-defined)
        n_expected = int(round(duration * sr))
        t_uniform  = np.linspace(gps_start, gps_end - 1.0 / sr, n_expected)
        v_interp   = np.interp(t_uniform, t_aux, v_aux)

        win_ch = max(5, int(round(1800.0 * sr)))
        z_aux  = rolling_mad_zscore(v_interp, win_ch)

        # ── 1. Cross-correlation lag (aux vs Hrec) ─────────────────────────────
        lag_cc = float("nan")
        peak_cc = float("nan")
        if hrec_z1hz is not None:
            # Resample Hrec to same sr if needed
            if sr != 1.0:
                t_hrec_r = np.linspace(gps_start, gps_end - 1.0 / sr, n_expected)
                hrec_r   = np.interp(t_hrec_r, t_ref, hrec_z1hz)
            else:
                hrec_r = hrec_z1hz[: len(z_aux)]
            n_common = min(len(z_aux), len(hrec_r))
            lag_cc, peak_cc = cross_corr_lag(
                z_aux[:n_common], hrec_r[:n_common],
                max_lag=int(round(args.max_lag_s * sr))
            )
            lag_cc = lag_cc / sr   # convert to seconds

        # ── 2. Granger causality (aux → Hrec at 1 Hz) ──────────────────────────
        p_aux_causes_hrec = 1.0
        p_hrec_causes_aux = 1.0
        if hrec_z1hz is not None and sr <= 1.0:
            hrec_1 = hrec_z1hz[: len(z_aux)]
            n1     = min(len(z_aux), len(hrec_1))
            gl     = min(args.granger_max_lag, n1 // 10)
            if gl >= 1:
                p_aux_causes_hrec = granger_p(z_aux[:n1], hrec_1[:n1], gl)
                p_hrec_causes_aux = granger_p(hrec_1[:n1], z_aux[:n1], gl)

        # ── 3. Transfer entropy (mutual information approach) ──────────────────
        te_aux_to_hrec = float("nan")
        te_hrec_to_aux = float("nan")
        if hrec_z1hz is not None and sr <= 1.0:
            hrec_1 = hrec_z1hz[: len(z_aux)]
            n1     = min(len(z_aux), len(hrec_1))
            te_aux_to_hrec = transfer_entropy_knn(z_aux[:n1], hrec_1[:n1])
            te_hrec_to_aux = transfer_entropy_knn(hrec_1[:n1], z_aux[:n1])

        # ── 4. Retrieve noisehound lag from ranking ────────────────────────────
        row = ranking[ranking["channel"] == ch]
        rank_score  = float(row["rank_score"].values[0]) if len(row) else float("nan")
        lag_nh      = float(row["lag_median_s"].values[0]) if "lag_median_s" in row.columns and len(row) else float("nan")
        hit_frac    = float(row["hit_fraction"].values[0]) if "hit_fraction" in row.columns and len(row) else float("nan")
        med_z       = float(row["median_peak_z"].values[0]) if "median_peak_z" in row.columns and len(row) else float("nan")

        # ── Verdict heuristic ──────────────────────────────────────────────────
        # Causal candidates: lag_nh < 0 (aux fires BEFORE glitch) AND
        #   cross-correlation lag < -2 s AND Granger p(aux→hrec) < 0.05
        # Consequence candidates: lag_nh > 0 OR (Granger p(hrec→aux) < p(aux→hrec))
        verdict = "ambiguous"
        if (not np.isnan(lag_nh)) and lag_nh < -2.0:
            if p_aux_causes_hrec < 0.05 and p_aux_causes_hrec < p_hrec_causes_aux:
                verdict = "likely_causal"
            elif p_hrec_causes_aux < p_aux_causes_hrec:
                verdict = "likely_consequence"
            else:
                verdict = "precedes_glitch"
        elif (not np.isnan(lag_nh)) and lag_nh > 2.0:
            verdict = "likely_consequence"
        elif p_aux_causes_hrec < 0.05 and p_aux_causes_hrec < p_hrec_causes_aux:
            verdict = "possible_cause"

        te_direction = "→hrec" if (not np.isnan(te_aux_to_hrec)) and \
                                   (not np.isnan(te_hrec_to_aux)) and \
                                   te_aux_to_hrec > te_hrec_to_aux else "←hrec"

        results.append(dict(
            channel             = ch,
            rank_score          = rank_score,
            hit_fraction        = hit_frac,
            median_peak_z       = med_z,
            lag_nh_s            = lag_nh,
            lag_xcorr_s         = lag_cc,
            peak_xcorr          = peak_cc,
            granger_p_aux2hrec  = p_aux_causes_hrec,
            granger_p_hrec2aux  = p_hrec_causes_aux,
            te_aux2hrec         = te_aux_to_hrec,
            te_hrec2aux         = te_hrec_to_aux,
            te_direction        = te_direction,
            verdict             = verdict,
        ))

        print(
            f"  {ch[:55]:<55}  "
            f"lag_nh={lag_nh:+6.1f}s  "
            f"lag_cc={lag_cc:+6.1f}s  "
            f"Gp={p_aux_causes_hrec:.3f}  "
            f"TE={te_direction}  "
            f"→ {verdict}"
        )

    # ── Save outputs ───────────────────────────────────────────────────────────
    df = pd.DataFrame(results).sort_values("rank_score", ascending=False)
    csv_path = out_dir / "causality_summary.csv"
    df.to_csv(csv_path, index=False, float_format="%.4f")
    print(f"\nSaved {csv_path}")

    # Text report
    report_path = out_dir / "causality_report.txt"
    with open(report_path, "w") as f:
        f.write("NOISEHOUND — Causality Analysis Report\n")
        f.write(f"GPS range: {gps_start:.0f} – {gps_end:.0f}\n")
        f.write(f"Channels analysed: {len(df)}\n\n")

        for verdict_label, title in [
            ("likely_causal",     "=== LIKELY CAUSAL (aux precedes Hrec, Granger p<0.05) ==="),
            ("precedes_glitch",   "=== PRECEDES GLITCH (lag_nh<-2s, Granger ambiguous) ==="),
            ("possible_cause",    "=== POSSIBLE CAUSE (Granger only) ==="),
            ("ambiguous",         "=== AMBIGUOUS ==="),
            ("likely_consequence","=== LIKELY CONSEQUENCE (follows glitch) ==="),
        ]:
            sub = df[df["verdict"] == verdict_label]
            if sub.empty:
                continue
            f.write(f"\n{title}\n")
            f.write(f"{'Channel':<60} {'rank':>8} {'lag_nh':>8} {'lag_cc':>8} "
                    f"{'Gp_fwd':>8} {'Gp_rev':>8} {'TE_dir':>8}\n")
            f.write("-" * 110 + "\n")
            for _, row in sub.iterrows():
                f.write(
                    f"{row['channel']:<60} "
                    f"{row['rank_score']:>8.2f} "
                    f"{row['lag_nh_s']:>+8.1f} "
                    f"{row['lag_xcorr_s']:>+8.1f} "
                    f"{row['granger_p_aux2hrec']:>8.3f} "
                    f"{row['granger_p_hrec2aux']:>8.3f} "
                    f"{row['te_direction']:>8}\n"
                )

    print(f"Saved {report_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
