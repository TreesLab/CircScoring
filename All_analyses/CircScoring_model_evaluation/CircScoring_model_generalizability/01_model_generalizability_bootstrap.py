#!/usr/bin/env python3
"""Reproduce the CircScoring model-generalizability bootstrap analysis.

The four finalized CircScoring score vectors are evaluated without refitting in
three independent benchmark pairs: P1N1, P2N2 and P2'N2'. Confidence intervals
are class-stratified nonparametric bootstrap percentile intervals.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260803
BOOTSTRAP_ITERATIONS = 1000

BENCHMARKS = [
    ("P1N1", "P1", "N1", 797, 504),
    ("P2N2", "P2", "N2", 805, 484),
    ("P2'N2'", "P2'", "N2'", 203, 162),
]

# This order is part of the historical random-number stream and must not change.
SCORE_MODELS = [
    ("CS-R", "XGBoost", "XGBoost_CS-R_score"),
    ("CS-R", "Elastic Net", "ElasticNet_CS-R_score"),
    ("CS-C", "XGBoost", "XGBoost_CS-C_score"),
    ("CS-C", "Elastic Net", "ElasticNet_CS-C_score"),
]


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def auroc(y: np.ndarray, score: np.ndarray) -> float:
    """Mann-Whitney AUROC with average ranks for tied scores."""
    ranks = pd.Series(score).rank(method="average").to_numpy(dtype=float)
    n_positive = int(y.sum())
    n_negative = int(len(y) - n_positive)
    if n_positive == 0 or n_negative == 0:
        raise ValueError("AUROC requires both classes.")
    numerator = ranks[y == 1].sum() - n_positive * (n_positive + 1) / 2
    return float(numerator / (n_positive * n_negative))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    """Non-interpolated precision-recall area (average precision).

    AP = sum_n (R_n - R_{n-1}) * P_n, evaluated at the end of each distinct
    score threshold. This is the computation labeled AUPRC in the retained
    results and is not trapezoidal integration of the PR curve.
    """
    order = np.argsort(-score, kind="stable")
    sorted_y = y[order]
    sorted_score = score[order]
    true_positive = np.cumsum(sorted_y)
    false_positive = np.cumsum(1 - sorted_y)
    threshold_ends = np.r_[
        np.flatnonzero(np.diff(sorted_score) != 0),
        len(sorted_score) - 1,
    ]
    recall = true_positive[threshold_ends] / true_positive[-1]
    precision = true_positive[threshold_ends] / (
        true_positive[threshold_ends] + false_positive[threshold_ends]
    )
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def prepare_benchmark(
    data: pd.DataFrame,
    positive_column: str,
    negative_column: str,
    expected_positive: int,
    expected_negative: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    positive = pd.to_numeric(data[positive_column], errors="raise").to_numpy()
    negative = pd.to_numeric(data[negative_column], errors="raise").to_numpy()
    keep = (positive == 1) | (negative == 1)
    subset = data.loc[keep].reset_index(drop=True).copy()
    if ((subset[positive_column] == 1) & (subset[negative_column] == 1)).any():
        raise ValueError(f"{positive_column} and {negative_column} must be mutually exclusive.")
    y = (pd.to_numeric(subset[positive_column], errors="raise") == 1).astype(int).to_numpy()
    observed_positive = int(y.sum())
    observed_negative = int(len(y) - observed_positive)
    if (observed_positive, observed_negative) != (expected_positive, expected_negative):
        raise ValueError(
            f"Unexpected {positive_column}/{negative_column} composition: "
            f"observed {observed_positive}/{observed_negative}, expected "
            f"{expected_positive}/{expected_negative}."
        )
    return subset, y


def bootstrap_percentile_interval(
    y: np.ndarray,
    score: np.ndarray,
    statistic,
    rng: np.random.Generator,
) -> tuple[float, float]:
    positive_indices = np.flatnonzero(y == 1)
    negative_indices = np.flatnonzero(y == 0)
    values = np.empty(BOOTSTRAP_ITERATIONS, dtype=float)
    for iteration in range(BOOTSTRAP_ITERATIONS):
        sampled_positive = rng.choice(
            positive_indices, size=len(positive_indices), replace=True
        )
        sampled_negative = rng.choice(
            negative_indices, size=len(negative_indices), replace=True
        )
        sampled_indices = np.r_[sampled_positive, sampled_negative]
        values[iteration] = statistic(y[sampled_indices], score[sampled_indices])
    lower, upper = np.quantile(values, [0.025, 0.975])
    return float(lower), float(upper)


def run_analysis(data: pd.DataFrame) -> pd.DataFrame:
    required = {
        *(item for _, p, n, _, _ in BENCHMARKS for item in (p, n)),
        *(column for _, _, column in SCORE_MODELS),
    }
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # A single generator and a fixed loop order reproduce the retained stream.
    rng = np.random.default_rng(SEED)
    rows: list[dict[str, object]] = []
    for benchmark, positive_col, negative_col, n_positive, n_negative in BENCHMARKS:
        subset, y = prepare_benchmark(
            data, positive_col, negative_col, n_positive, n_negative
        )
        for score_name, model, score_column in SCORE_MODELS:
            score = pd.to_numeric(subset[score_column], errors="raise").to_numpy(dtype=float)
            if not np.isfinite(score).all():
                raise ValueError(f"Non-finite values in {score_column} for {benchmark}.")
            for metric_name, statistic in (
                ("AUROC", auroc),
                ("AUPRC", average_precision),
            ):
                estimate = statistic(y, score)
                ci_lower, ci_upper = bootstrap_percentile_interval(
                    y, score, statistic, rng
                )
                rows.append(
                    {
                        "Benchmark": benchmark,
                        "Score": score_name,
                        "Model": model,
                        "Metric": metric_name,
                        "Estimate": estimate,
                        "CI_lower": ci_lower,
                        "CI_upper": ci_upper,
                        "n_positive": n_positive,
                        "n_negative": n_negative,
                        "bootstrap_iterations": BOOTSTRAP_ITERATIONS,
                        "seed": SEED,
                    }
                )
    return pd.DataFrame(rows)


def verify_against_expected(observed: pd.DataFrame, expected_path: Path) -> None:
    expected = pd.read_csv(expected_path)
    if list(observed.columns) != list(expected.columns):
        raise AssertionError("Observed and expected column names/order differ.")
    if observed.shape != expected.shape:
        raise AssertionError(f"Shape differs: {observed.shape} versus {expected.shape}.")

    categorical = ["Benchmark", "Score", "Model", "Metric"]
    integer = ["n_positive", "n_negative", "bootstrap_iterations", "seed"]
    numeric = ["Estimate", "CI_lower", "CI_upper"]
    for column in categorical + integer:
        if not np.array_equal(observed[column].to_numpy(), expected[column].to_numpy()):
            raise AssertionError(f"Exact verification failed for {column}.")
    maximum_difference = float(
        np.max(np.abs(observed[numeric].to_numpy() - expected[numeric].to_numpy()))
    )
    if maximum_difference > 5e-14:
        raise AssertionError(
            "Numerical verification failed; maximum absolute difference="
            f"{maximum_difference:.3g}."
        )
    print(
        "Verification passed: maximum absolute estimate/CI difference="
        f"{maximum_difference:.3g}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = read_table(args.input)
    result = run_analysis(data)
    if args.expected:
        verify_against_expected(result, args.expected)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result)} rows to {args.output.resolve()}")


if __name__ == "__main__":
    main()
