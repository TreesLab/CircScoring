#!/usr/bin/env python3
"""CircScoring validation using two orthogonal experimental evidence sets.

Numerical workflow only: group counts, ECDF source data, observed RBC and
median differences, and 10,000 size-matched background-resampling tests.
No figure is generated.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RESAMPLES = 10_000
DEFAULT_SEED = 20260815

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

RT_INDEPENDENT_COLUMN = "RT_independent_NCL_event"
RT_NONRT_COLUMN = "Both RT/nonRT"

EVIDENCE_DEFINITIONS = {
    "RT-/non-RT evidence": RT_NONRT_COLUMN,
    "RT-independent evidence": RT_INDEPENDENT_COLUMN,
}

EXPECTED_COUNTS = {
    "Source rows": 4_455_197,
    "Benchmark-overlapping rows excluded": 2_282,
    "Incomplete-rank rows excluded": 0,
    "Background": 4_449_931,
    "RT-/non-RT evidence": 599,
    "RT-independent evidence": 2_493,
    "Overlap of two positive sets": 108,
}

EXPECTED_RBC = {
    ("RT-/non-RT evidence", "XGBoost CS-R"): 0.6392442256,
    ("RT-/non-RT evidence", "XGBoost CS-C"): 0.7377390465,
    ("RT-/non-RT evidence", "Elastic Net CS-R"): 0.6270580214,
    ("RT-/non-RT evidence", "Elastic Net CS-C"): 0.6944889884,
    ("RT-independent evidence", "XGBoost CS-R"): 0.9688858750,
    ("RT-independent evidence", "XGBoost CS-C"): 0.8644165523,
    ("RT-independent evidence", "Elastic Net CS-R"): 0.9717084467,
    ("RT-independent evidence", "Elastic Net CS-C"): 0.8751511594,
}

EXPECTED_POSITIVE_MEDIANS = {
    ("RT-/non-RT evidence", "XGBoost CS-R"): 98.52615945,
    ("RT-/non-RT evidence", "XGBoost CS-C"): 91.48703189,
    ("RT-/non-RT evidence", "Elastic Net CS-R"): 98.52615945,
    ("RT-/non-RT evidence", "Elastic Net CS-C"): 91.83335776,
    ("RT-independent evidence", "XGBoost CS-R"): 98.52615945,
    ("RT-independent evidence", "XGBoost CS-C"): 93.31259650,
    ("RT-independent evidence", "Elastic Net CS-R"): 98.52615945,
    ("RT-independent evidence", "Elastic Net CS-C"): 94.44525124,
}

EXPECTED_BACKGROUND_MEDIANS = {
    "XGBoost CS-R": 38.58723868,
    "XGBoost CS-C": 49.94904827,
    "Elastic Net CS-R": 29.70659883,
    "Elastic Net CS-C": 49.94367252,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "RT-independent and RT-/non-RT evidence ECDF, RBC, median, and "
            "size-matched background-resampling analyses."
        )
    )
    parser.add_argument("input", type=Path, help="One- or two-row-header circIG TSV[.gz].")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--resamples", type=int, default=DEFAULT_RESAMPLES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--relaxed-qc",
        action="store_true",
        help="Disable fixed-population controls; intended only for code testing.",
    )
    parser.add_argument(
        "--no-null-output",
        action="store_true",
        help="Do not save iteration-level null distributions.",
    )
    args = parser.parse_args()
    if args.resamples < 1:
        parser.error("--resamples must be a positive integer")
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
    first_header = first.split(delimiter)
    second_header = second.split(delimiter)
    required = {
        RT_INDEPENDENT_COLUMN,
        RT_NONRT_COLUMN,
        *RANK_COLUMNS.values(),
        *BENCHMARK_COLUMNS,
    }
    if required.issubset(first_header):
        return 0, delimiter
    if required.issubset(second_header):
        return 1, delimiter
    missing = sorted(required.difference(set(first_header) | set(second_header)))
    raise ValueError(f"Required columns were not found: {', '.join(missing)}")


def numeric_blank_zero(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float)


def load_analysis_data(
    path: Path,
) -> tuple[pd.DataFrame, dict[str, int]]:
    skiprows, delimiter = detect_header(path)
    usecols = [
        RT_INDEPENDENT_COLUMN,
        RT_NONRT_COLUMN,
        *RANK_COLUMNS.values(),
        *BENCHMARK_COLUMNS,
    ]
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

    for column in RANK_COLUMNS.values():
        data[column] = pd.to_numeric(data[column], errors="coerce")
    incomplete = data[list(RANK_COLUMNS.values())].isna().any(axis=1)
    incomplete_excluded = int(incomplete.sum())
    data = data.loc[~incomplete].copy()

    data[RT_INDEPENDENT_COLUMN] = numeric_blank_zero(data[RT_INDEPENDENT_COLUMN])
    data[RT_NONRT_COLUMN] = numeric_blank_zero(data[RT_NONRT_COLUMN])

    qc = {
        "Source rows": source_rows,
        "Benchmark-overlapping rows excluded": benchmark_excluded,
        "Incomplete-rank rows excluded": incomplete_excluded,
        "Rows retained after exclusions": len(data),
    }
    return data, qc


def define_groups(data: pd.DataFrame) -> dict[str, np.ndarray]:
    groups = {
        "Background": (
            (data[RT_INDEPENDENT_COLUMN].to_numpy() == 0)
            & (data[RT_NONRT_COLUMN].to_numpy() == 0)
        )
    }
    for evidence, column in EVIDENCE_DEFINITIONS.items():
        groups[evidence] = data[column].to_numpy() == 1
    return groups


def make_ecdf_table(
    data: pd.DataFrame,
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    frames = []
    for group_name, mask in groups.items():
        for rank_name, column in RANK_COLUMNS.items():
            values = np.sort(data.loc[mask, column].to_numpy(dtype=float))
            cumulative = np.searchsorted(values, thresholds, side="right")
            frames.append(
                pd.DataFrame(
                    {
                        "Group": group_name,
                        "Rank": rank_name,
                        "N_BSJ": len(values),
                        "Percentile_rank_threshold": thresholds,
                        "ECDF_percent": cumulative / len(values) * 100.0,
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


def observed_rbc(positive: np.ndarray, background: np.ndarray) -> np.ndarray:
    """RBC using pairwise wins plus half of tied pairs."""
    result = np.empty(positive.shape[1], dtype=float)
    n_positive = len(positive)
    n_background = len(background)
    for column in range(positive.shape[1]):
        sorted_background = np.sort(background[:, column])
        lower = np.searchsorted(sorted_background, positive[:, column], side="left")
        upper = np.searchsorted(sorted_background, positive[:, column], side="right")
        u_statistic = np.sum(lower + 0.5 * (upper - lower), dtype=float)
        result[column] = 2.0 * u_statistic / (n_positive * n_background) - 1.0
    return result


def background_average_ranks(background: np.ndarray) -> np.ndarray:
    return pd.DataFrame(background).rank(method="average", axis=0).to_numpy(dtype=float)


def rbc_for_sample_vs_complement(
    rank_matrix: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    n_sample = len(selected)
    n_complement = len(rank_matrix) - n_sample
    rank_sum = rank_matrix[selected, :].sum(axis=0)
    u_statistic = rank_sum - n_sample * (n_sample + 1) / 2.0
    return 2.0 * u_statistic / (n_sample * n_complement) - 1.0


def kth_after_removal(
    sorted_population: np.ndarray,
    sorted_removed: np.ndarray,
    zero_based_rank: int,
) -> float:
    """Exact order statistic after removing a small sampled multiset."""
    target_count = zero_based_rank + 1
    low = 0
    high = len(sorted_population) - 1
    while low < high:
        middle = (low + high) // 2
        value = sorted_population[middle]
        population_le = np.searchsorted(sorted_population, value, side="right")
        removed_le = np.searchsorted(sorted_removed, value, side="right")
        if population_le - removed_le >= target_count:
            high = middle
        else:
            low = middle + 1
    return float(sorted_population[low])


def complement_median(
    sorted_background: np.ndarray,
    sampled_values: np.ndarray,
) -> np.ndarray:
    n_remaining = len(sorted_background) - len(sampled_values)
    result = np.empty(sorted_background.shape[1], dtype=float)
    for column in range(sorted_background.shape[1]):
        removed = np.sort(sampled_values[:, column])
        if n_remaining % 2:
            result[column] = kth_after_removal(
                sorted_background[:, column], removed, n_remaining // 2
            )
        else:
            lower = kth_after_removal(
                sorted_background[:, column], removed, n_remaining // 2 - 1
            )
            upper = kth_after_removal(
                sorted_background[:, column], removed, n_remaining // 2
            )
            result[column] = (lower + upper) / 2.0
    return result


def summarize_null(
    evidence: str,
    statistic_name: str,
    observed: np.ndarray,
    null: np.ndarray,
    n_positive: int,
    n_background: int,
    n_resamples: int,
    seed: int,
    extra: dict[str, np.ndarray] | None = None,
) -> pd.DataFrame:
    rows = []
    extra = extra or {}
    for index, rank_name in enumerate(RANK_COLUMNS):
        values = null[:, index]
        extreme = int(np.count_nonzero(values >= observed[index]))
        row = {
            "Evidence": evidence,
            "Rank": rank_name,
            "N_positive": n_positive,
            "N_background": n_background,
            "N_random_sample": n_positive,
            "N_background_complement": n_background - n_positive,
            f"Observed_{statistic_name}": observed[index],
            "Null_mean": values.mean(),
            "Null_SD": values.std(ddof=1),
            "Null_2.5pct": np.quantile(values, 0.025),
            "Null_97.5pct": np.quantile(values, 0.975),
            "Null_max": values.max(),
            "Extreme_resamples_ge_observed": extreme,
            "Empirical_one_sided_P_plus1": (extreme + 1) / (n_resamples + 1),
            "N_resamples": n_resamples,
            "Random_seed": seed,
        }
        for name, extra_values in extra.items():
            row[name] = extra_values[index]
        rows.append(row)
    return pd.DataFrame(rows)


def strict_qc(
    qc: dict[str, int],
    groups: dict[str, np.ndarray],
    observed_rbcs: dict[str, np.ndarray],
    positive_medians: dict[str, np.ndarray],
    background_median: np.ndarray,
) -> None:
    observed_counts = {
        **qc,
        "Background": int(groups["Background"].sum()),
        "RT-/non-RT evidence": int(groups["RT-/non-RT evidence"].sum()),
        "RT-independent evidence": int(groups["RT-independent evidence"].sum()),
        "Overlap of two positive sets": int(
            np.count_nonzero(
                groups["RT-/non-RT evidence"]
                & groups["RT-independent evidence"]
            )
        ),
    }
    for name, expected in EXPECTED_COUNTS.items():
        if observed_counts[name] != expected:
            raise RuntimeError(
                f"Fixed-population control failed for {name}: "
                f"{observed_counts[name]} != {expected}"
            )

    rank_names = list(RANK_COLUMNS)
    for evidence, values in observed_rbcs.items():
        expected = np.array([EXPECTED_RBC[(evidence, rank)] for rank in rank_names])
        if not np.allclose(values, expected, atol=5e-10, rtol=0):
            raise RuntimeError(f"Observed RBC controls failed for {evidence}.")
    for evidence, values in positive_medians.items():
        expected = np.array(
            [EXPECTED_POSITIVE_MEDIANS[(evidence, rank)] for rank in rank_names]
        )
        if not np.allclose(values, expected, atol=5e-8, rtol=0):
            raise RuntimeError(f"Positive-median controls failed for {evidence}.")
    expected_background = np.array(
        [EXPECTED_BACKGROUND_MEDIANS[rank] for rank in rank_names]
    )
    if not np.allclose(background_median, expected_background, atol=5e-8, rtol=0):
        raise RuntimeError("Background-median controls failed.")


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)

    data, qc = load_analysis_data(args.input)
    groups = define_groups(data)
    rank_columns = list(RANK_COLUMNS.values())
    background = data.loc[groups["Background"], rank_columns].to_numpy(dtype=float)
    positives = {
        evidence: data.loc[groups[evidence], rank_columns].to_numpy(dtype=float)
        for evidence in EVIDENCE_DEFINITIONS
    }

    ecdf = make_ecdf_table(data, groups)
    ecdf.to_csv(
        args.output_directory / "RT_evidence_ECDF_source_data_0.1.tsv",
        sep="\t",
        index=False,
    )

    group_count_rows = [
        {"Group": name, "N_BSJ": int(mask.sum())} for name, mask in groups.items()
    ]
    group_count_rows.append(
        {
            "Group": "Overlap of two positive sets",
            "N_BSJ": int(
                np.count_nonzero(
                    groups["RT-/non-RT evidence"]
                    & groups["RT-independent evidence"]
                )
            ),
        }
    )
    pd.DataFrame(group_count_rows).to_csv(
        args.output_directory / "RT_evidence_group_counts.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [{"QC_item": name, "N_rows": value} for name, value in qc.items()]
    ).to_csv(
        args.output_directory / "RT_evidence_QC.tsv", sep="\t", index=False
    )

    observed_rbcs = {
        evidence: observed_rbc(values, background)
        for evidence, values in positives.items()
    }
    background_median = np.median(background, axis=0)
    positive_medians = {
        evidence: np.median(values, axis=0)
        for evidence, values in positives.items()
    }
    observed_median_differences = {
        evidence: values - background_median
        for evidence, values in positive_medians.items()
    }

    if not args.relaxed_qc:
        strict_qc(
            qc,
            groups,
            observed_rbcs,
            positive_medians,
            background_median,
        )

    # RBC randomization. Each evidence set is tested separately. A single
    # sampled index set is applied to all four ranks within each iteration.
    rank_matrix = background_average_ranks(background)
    rbc_rng = np.random.default_rng(args.seed)
    rbc_null: dict[str, np.ndarray] = {}
    rbc_summary = []
    for evidence, positive in positives.items():
        n_positive = len(positive)
        null = np.empty((args.resamples, len(RANK_COLUMNS)), dtype=float)
        for iteration in range(args.resamples):
            selected = rbc_rng.choice(len(background), size=n_positive, replace=False)
            null[iteration, :] = rbc_for_sample_vs_complement(rank_matrix, selected)
        rbc_null[evidence] = null
        rbc_summary.append(
            summarize_null(
                evidence,
                "RBC",
                observed_rbcs[evidence],
                null,
                n_positive,
                len(background),
                args.resamples,
                args.seed,
            )
        )
        print(f"RBC resampling completed: {evidence}", flush=True)
    del rank_matrix

    # Median randomization is a distinct Monte Carlo run initialized with the
    # same documented seed. The comparator is the background after removal of
    # the sampled BSJs, not the full unmodified background.
    sorted_background = np.sort(background, axis=0)
    median_rng = np.random.default_rng(args.seed)
    median_null: dict[str, np.ndarray] = {}
    median_summary = []
    for evidence, positive in positives.items():
        n_positive = len(positive)
        null = np.empty((args.resamples, len(RANK_COLUMNS)), dtype=float)
        for iteration in range(args.resamples):
            selected = median_rng.choice(len(background), size=n_positive, replace=False)
            sampled = background[selected, :]
            random_median = np.median(sampled, axis=0)
            remaining_median = complement_median(sorted_background, sampled)
            null[iteration, :] = random_median - remaining_median
        median_null[evidence] = null
        median_summary.append(
            summarize_null(
                evidence,
                "median_difference",
                observed_median_differences[evidence],
                null,
                n_positive,
                len(background),
                args.resamples,
                args.seed,
                extra={
                    "Positive_median": positive_medians[evidence],
                    "Background_median": background_median,
                },
            )
        )
        print(f"Median resampling completed: {evidence}", flush=True)

    pd.concat(rbc_summary, ignore_index=True).to_csv(
        args.output_directory / "RT_evidence_RBC_resampling_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.concat(median_summary, ignore_index=True).to_csv(
        args.output_directory / "RT_evidence_median_resampling_summary.tsv",
        sep="\t",
        index=False,
    )

    if not args.no_null_output:
        rbc_distribution = {}
        median_distribution = {}
        for evidence in EVIDENCE_DEFINITIONS:
            for index, rank_name in enumerate(RANK_COLUMNS):
                key = f"{evidence}__{rank_name}"
                rbc_distribution[key] = rbc_null[evidence][:, index]
                median_distribution[key] = median_null[evidence][:, index]
        pd.DataFrame(rbc_distribution).to_csv(
            args.output_directory /
            f"RT_evidence_null_RBC_{args.resamples}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        pd.DataFrame(median_distribution).to_csv(
            args.output_directory /
            f"RT_evidence_null_median_difference_{args.resamples}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
