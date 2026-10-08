#!/usr/bin/env python3
"""TCCIA circRNA-prevalence validation for CircScoring.

Numerical analysis only: BSJ-level Spearman correlations, paired nonparametric
bootstrap confidence intervals, and cohort-group median/IQR summaries.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260917
DEFAULT_BOOTSTRAPS = 10_000

BENCHMARK_COLUMNS = [
    "positive1_P1N1",
    "positive2_P4N4_L",
    "positive3_P4N4_H",
    "negative1_P1N1",
    "negative2_P4N4_L",
    "negative3_P4N4_H",
]

RANK_COLUMNS = {
    "XGBoost CS-R": "XGBoost_CS-R_percentile_rank",
    "XGBoost CS-C": "XGBoost_CS-C_percentile_rank",
    "Elastic Net CS-R": "ElasticNet_CS-R_percentile_rank",
    "Elastic Net CS-C": "ElasticNet_CS-C_percentile_rank",
}

EXPECTED_RHO = {
    "XGBoost CS-R": 0.581146058659086,
    "XGBoost CS-C": 0.302669663146107,
    "Elastic Net CS-R": 0.570647894746137,
    "Elastic Net CS-C": 0.299908499595978,
}

EXPECTED_GROUP_COUNTS = {
    "1": 247_586,
    "2": 109_862,
    "3": 35_056,
    "4": 20_642,
    "5": 14_840,
    "6": 11_066,
    "7": 9_105,
    "8": 7_903,
    "9": 7_105,
    "10": 6_545,
    "11-15": 24_597,
    ">=16": 4_415,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "BSJ-level TCCIA cohort-count versus CircScoring percentile-rank "
            "Spearman and bootstrap analysis."
        )
    )
    parser.add_argument("input", type=Path, help="TCCIA analysis TSV[.gz].")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--relaxed-qc",
        action="store_true",
        help="Disable fixed-population controls; intended only for testing.",
    )
    parser.add_argument(
        "--no-bootstrap-output",
        action="store_true",
        help="Do not save the iteration-level bootstrap distribution.",
    )
    args = parser.parse_args()
    if args.bootstraps < 1:
        parser.error("--bootstraps must be a positive integer")
    return args


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("rt", encoding="utf-8-sig", newline="")


def detect_header(path: Path) -> tuple[int, str]:
    with open_text(path) as handle:
        first = handle.readline().rstrip("\r\n")
        second = handle.readline().rstrip("\r\n")
    delimiter = "\t" if "\t" in first else ","
    h1 = first.split(delimiter)
    h2 = second.split(delimiter)
    required = {"Cohort_count", *RANK_COLUMNS.values(), *BENCHMARK_COLUMNS}
    if required.issubset(h1):
        return 0, delimiter
    if required.issubset(h2):
        return 1, delimiter
    missing = sorted(required.difference(set(h1) | set(h2)))
    raise ValueError(f"Required columns were not found: {', '.join(missing)}")


def numeric_blank_zero(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float)


def load_data(path: Path) -> tuple[pd.DataFrame, int, int, int]:
    skiprows, delimiter = detect_header(path)
    usecols = ["Cohort_count", *RANK_COLUMNS.values(), *BENCHMARK_COLUMNS]
    data = pd.read_csv(
        path,
        sep=delimiter,
        skiprows=skiprows,
        usecols=usecols,
        low_memory=False,
    )
    source_rows = len(data)

    benchmark_overlap = np.zeros(source_rows, dtype=bool)
    for column in BENCHMARK_COLUMNS:
        benchmark_overlap |= numeric_blank_zero(data[column]) != 0
    benchmark_excluded = int(benchmark_overlap.sum())
    data = data.loc[~benchmark_overlap].copy()

    analysis_columns = ["Cohort_count", *RANK_COLUMNS.values()]
    for column in analysis_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    incomplete = data[analysis_columns].isna().any(axis=1)
    incomplete_excluded = int(incomplete.sum())
    data = data.loc[~incomplete].copy()
    return data, source_rows, benchmark_excluded, incomplete_excluded


def average_ranks(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average").to_numpy(dtype=float)


def pearson_columns(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x_centered = x - x.mean()
    y_centered = y - y.mean(axis=0)
    numerator = x_centered @ y_centered
    denominator = np.sqrt(
        np.dot(x_centered, x_centered) * np.sum(y_centered * y_centered, axis=0)
    )
    return numerator / denominator


def log_regularized_beta_small_x(a: float, b: float, x: float) -> float:
    """Log of I_x(a, b) when x < (a + 1) / (a + b + 2).

    The continued-fraction evaluation avoids underflow for the extremely small
    P values in this analysis. The implementation follows the standard
    modified-Lentz evaluation of the incomplete-beta continued fraction.
    """
    if not 0.0 < x < 1.0:
        raise ValueError("x must be strictly between 0 and 1")
    maximum_iterations = 10_000
    epsilon = 3.0e-14
    minimum = np.finfo(float).tiny / epsilon
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < minimum:
        d = minimum
    d = 1.0 / d
    fraction = d
    for iteration in range(1, maximum_iterations + 1):
        twice = 2.0 * iteration
        coefficient = (
            iteration * (b - iteration) * x /
            ((qam + twice) * (a + twice))
        )
        d = 1.0 + coefficient * d
        if abs(d) < minimum:
            d = minimum
        c = 1.0 + coefficient / c
        if abs(c) < minimum:
            c = minimum
        d = 1.0 / d
        fraction *= d * c

        coefficient = -(
            (a + iteration) * (qab + iteration) * x /
            ((a + twice) * (qap + twice))
        )
        d = 1.0 + coefficient * d
        if abs(d) < minimum:
            d = minimum
        c = 1.0 + coefficient / c
        if abs(c) < minimum:
            c = minimum
        d = 1.0 / d
        change = d * c
        fraction *= change
        if abs(change - 1.0) < epsilon:
            break
    else:
        raise RuntimeError("Incomplete-beta continued fraction did not converge")

    log_front = (
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    return log_front + math.log(fraction) - math.log(a)


def student_t_two_sided_log_p(t_statistic: float, df: int) -> float:
    """Two-sided Student-t log P without floating-point underflow."""
    beta_x = df / (df + t_statistic * t_statistic)
    return log_regularized_beta_small_x(df / 2.0, 0.5, beta_x)


def weighted_correlations(
    weights: np.ndarray,
    x: np.ndarray,
    x_squared: np.ndarray,
    y: np.ndarray,
    y_squared: np.ndarray,
    xy: np.ndarray,
) -> np.ndarray:
    total = float(weights.sum())
    sum_x = np.dot(weights, x)
    sum_x2 = np.dot(weights, x_squared)
    sum_y = weights @ y
    sum_y2 = weights @ y_squared
    sum_xy = weights @ xy
    covariance = sum_xy - sum_x * sum_y / total
    variance_x = sum_x2 - sum_x * sum_x / total
    variance_y = sum_y2 - sum_y * sum_y / total
    return covariance / np.sqrt(variance_x * variance_y)


def bootstrap_spearman(
    x_rank: np.ndarray,
    y_ranks: np.ndarray,
    n_bootstraps: int,
    seed: int,
) -> np.ndarray:
    """Paired bootstrap using multinomial case weights.

    Each multinomial vector represents the multiplicities of the N BSJ pairs
    in one nonparametric bootstrap sample. The same weights are used for all
    four percentile-rank variables within an iteration.
    """
    n = len(x_rank)
    probabilities = np.full(n, 1.0 / n, dtype=float)
    x_squared = x_rank * x_rank
    y_squared = y_ranks * y_ranks
    xy = x_rank[:, None] * y_ranks
    result = np.empty((n_bootstraps, y_ranks.shape[1]), dtype=float)
    rng = np.random.default_rng(seed)

    for iteration in range(n_bootstraps):
        weights = rng.multinomial(n, probabilities)
        result[iteration, :] = weighted_correlations(
            weights, x_rank, x_squared, y_ranks, y_squared, xy
        )
        if (iteration + 1) % 100 == 0 or iteration + 1 == n_bootstraps:
            print(
                f"Bootstrap {iteration + 1:,}/{n_bootstraps:,}",
                end="\r" if iteration + 1 < n_bootstraps else "\n",
                flush=True,
            )
    return result


def cohort_group(value: float) -> str:
    value = int(value)
    if 1 <= value <= 10:
        return str(value)
    if 11 <= value <= 15:
        return "11-15"
    if value >= 16:
        return ">=16"
    raise ValueError(f"Unexpected Cohort_count: {value}")


def median_iqr_table(data: pd.DataFrame) -> pd.DataFrame:
    order = [str(i) for i in range(1, 11)] + ["11-15", ">=16"]
    groups = data["Cohort_count"].map(cohort_group)
    rows = []
    for group_name in order:
        keep = groups == group_name
        for rank_name, column in RANK_COLUMNS.items():
            values = data.loc[keep, column].to_numpy(dtype=float)
            q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
            rows.append(
                {
                    "Cohort_group": group_name,
                    "Rank": rank_name,
                    "N_BSJ": len(values),
                    "Median": median,
                    "Q1": q1,
                    "Q3": q3,
                }
            )
    return pd.DataFrame(rows)


def strict_qc(
    data: pd.DataFrame,
    source_rows: int,
    benchmark_excluded: int,
    incomplete_excluded: int,
    observed_rho: np.ndarray,
    grouped: pd.DataFrame,
) -> None:
    if (source_rows, benchmark_excluded, incomplete_excluded, len(data)) != (
        500_535,
        1_813,
        0,
        498_722,
    ):
        raise RuntimeError("Fixed TCCIA population controls do not match.")
    expected_rho = np.array([EXPECTED_RHO[name] for name in RANK_COLUMNS])
    if not np.allclose(observed_rho, expected_rho, atol=5e-12, rtol=0):
        raise RuntimeError("Observed Spearman correlations do not match retained results.")
    counts = (
        grouped[["Cohort_group", "N_BSJ"]]
        .drop_duplicates()
        .set_index("Cohort_group")["N_BSJ"]
        .to_dict()
    )
    if counts != EXPECTED_GROUP_COUNTS:
        raise RuntimeError(f"Cohort-group counts do not match: {counts}")


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, source_rows, benchmark_excluded, incomplete_excluded = load_data(args.input)

    x_rank = average_ranks(data["Cohort_count"].to_numpy(dtype=float))
    y_ranks = np.column_stack(
        [average_ranks(data[column].to_numpy(dtype=float)) for column in RANK_COLUMNS.values()]
    )
    observed_rho = pearson_columns(x_rank, y_ranks)
    n = len(data)
    degrees_of_freedom = n - 2
    t_statistic = observed_rho * np.sqrt(
        degrees_of_freedom / (1.0 - observed_rho * observed_rho)
    )
    log_p = np.array(
        [student_t_two_sided_log_p(abs(value), degrees_of_freedom)
         for value in t_statistic],
        dtype=float,
    )
    log10_p = log_p / math.log(10.0)

    bootstrap = bootstrap_spearman(
        x_rank, y_ranks, args.bootstraps, args.seed
    )
    ci_low = np.quantile(bootstrap, 0.025, axis=0)
    ci_high = np.quantile(bootstrap, 0.975, axis=0)

    grouped = median_iqr_table(data)
    if not args.relaxed_qc:
        strict_qc(
            data,
            source_rows,
            benchmark_excluded,
            incomplete_excluded,
            observed_rho,
            grouped,
        )

    summary_rows = []
    for index, rank_name in enumerate(RANK_COLUMNS):
        summary_rows.append(
            {
                "Rank": rank_name,
                "N_BSJ": n,
                "Spearman_rho": observed_rho[index],
                "Bootstrap_95pct_CI_low": ci_low[index],
                "Bootstrap_95pct_CI_high": ci_high[index],
                "Bootstrap_resamples": args.bootstraps,
                "Random_seed": args.seed,
                "Source_rows": source_rows,
                "Benchmark_overlap_rows_excluded": benchmark_excluded,
                "Incomplete_rows_excluded_after_benchmark_filter": incomplete_excluded,
                "T_statistic": t_statistic[index],
                "Degrees_of_freedom": degrees_of_freedom,
                "P_value_two_sided": (
                    "<1e-300" if log10_p[index] < -300 else math.exp(log_p[index])
                ),
                "Log10_P_value_two_sided": log10_p[index],
            }
        )
    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        args.output_directory / "TCCIA_Spearman_bootstrap_summary.tsv",
        sep="\t",
        index=False,
    )
    grouped.to_csv(
        args.output_directory / "TCCIA_cohort_rank_median_IQR.tsv",
        sep="\t",
        index=False,
    )
    if not args.no_bootstrap_output:
        distribution = pd.DataFrame(bootstrap, columns=list(RANK_COLUMNS))
        distribution.insert(0, "Bootstrap_iteration", np.arange(1, args.bootstraps + 1))
        distribution.to_csv(
            args.output_directory /
            f"TCCIA_Spearman_bootstrap{args.bootstraps}_distribution.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
