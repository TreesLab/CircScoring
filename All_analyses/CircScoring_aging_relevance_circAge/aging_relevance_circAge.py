#!/usr/bin/env python3
"""CircScoring aging-relevance validation using circAge.

Numerical workflow only: five-tissue ECDF source data, descriptive summaries,
tie-corrected Wilcoxon rank-sum tests, observed RBCs, and 10,000 pooled-label
randomizations. No figure is generated.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RANDOMIZATIONS = 10_000

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

TISSUES = ["Aortic", "Lung", "Skin", "Umbilical", "white_blood_cell"]

TISSUE_SEEDS = {
    "Aortic": 2026091801,
    "Lung": 2026091802,
    "Skin": 2026091803,
    "Umbilical": 2026091804,
    "white_blood_cell": 2026091805,
}

EXPECTED_QC = {
    "Source rows": 171_686,
    "Benchmark-overlapping rows excluded": 1_504,
    "Incomplete-rank rows excluded": 0,
    "Rows retained after exclusions": 170_182,
}

EXPECTED_COUNTS = {
    "Aortic": (20_844, 142, 149_196),
    "Lung": (73_406, 3_914, 92_862),
    "Skin": (53_851, 702, 115_629),
    "Umbilical": (25_121, 130, 144_931),
    "white_blood_cell": (61_936, 411, 107_835),
}

EXPECTED_RBC = {
    ("Aortic", "XGBoost CS-R"): 0.219278,
    ("Aortic", "XGBoost CS-C"): 0.250992,
    ("Aortic", "Elastic Net CS-R"): 0.219732,
    ("Aortic", "Elastic Net CS-C"): 0.171713,
    ("Lung", "XGBoost CS-R"): 0.359007,
    ("Lung", "XGBoost CS-C"): 0.254858,
    ("Lung", "Elastic Net CS-R"): 0.360000,
    ("Lung", "Elastic Net CS-C"): 0.280342,
    ("Skin", "XGBoost CS-R"): 0.404993,
    ("Skin", "XGBoost CS-C"): 0.219870,
    ("Skin", "Elastic Net CS-R"): 0.404272,
    ("Skin", "Elastic Net CS-C"): 0.233211,
    ("Umbilical", "XGBoost CS-R"): 0.240372,
    ("Umbilical", "XGBoost CS-C"): 0.134132,
    ("Umbilical", "Elastic Net CS-R"): 0.241823,
    ("Umbilical", "Elastic Net CS-C"): 0.152149,
    ("white_blood_cell", "XGBoost CS-R"): 0.510627,
    ("white_blood_cell", "XGBoost CS-C"): 0.448555,
    ("white_blood_cell", "Elastic Net CS-R"): 0.509416,
    ("white_blood_cell", "Elastic Net CS-C"): 0.468315,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Five-tissue circAge ECDF, Wilcoxon, RBC, and pooled-label "
            "randomization analysis."
        )
    )
    parser.add_argument("input", type=Path, help="One- or two-row-header circAge TSV[.gz].")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument(
        "--randomizations", type=int, default=DEFAULT_RANDOMIZATIONS
    )
    parser.add_argument(
        "--relaxed-qc",
        action="store_true",
        help="Disable fixed-population controls; intended only for code testing.",
    )
    parser.add_argument(
        "--no-null-output",
        action="store_true",
        help="Do not save the iteration-level null distribution.",
    )
    args = parser.parse_args()
    if args.randomizations < 1:
        parser.error("--randomizations must be a positive integer")
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
        *RANK_COLUMNS.values(),
        *BENCHMARK_COLUMNS,
        *TISSUES,
    }
    if required.issubset(first_header):
        return 0, delimiter
    if required.issubset(second_header):
        return 1, delimiter
    missing = sorted(required.difference(set(first_header) | set(second_header)))
    raise ValueError(f"Required columns were not found: {', '.join(missing)}")


def numeric_blank_zero(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float)


def load_analysis_data(path: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    skiprows, delimiter = detect_header(path)
    usecols = [*RANK_COLUMNS.values(), *BENCHMARK_COLUMNS, *TISSUES]
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
    for tissue in TISSUES:
        data[tissue] = pd.to_numeric(data[tissue], errors="coerce")

    qc = {
        "Source rows": source_rows,
        "Benchmark-overlapping rows excluded": benchmark_excluded,
        "Incomplete-rank rows excluded": incomplete_excluded,
        "Rows retained after exclusions": len(data),
    }
    return data, qc


def average_rank_matrix(values: np.ndarray) -> np.ndarray:
    return pd.DataFrame(values).rank(method="average", axis=0).to_numpy(dtype=float)


def rbc_from_rank_sum(
    rank_sum_group_1: np.ndarray,
    n_group_1: int,
    n_group_0: int,
) -> np.ndarray:
    u_statistic = rank_sum_group_1 - n_group_1 * (n_group_1 + 1) / 2.0
    return 2.0 * u_statistic / (n_group_1 * n_group_0) - 1.0


def normal_log_sf(z: float) -> float:
    """Log upper-tail probability for a standard normal variate."""
    z = abs(float(z))
    if z < 8.0:
        return math.log(0.5 * math.erfc(z / math.sqrt(2.0)))
    inverse_square = 1.0 / (z * z)
    correction = (
        1.0
        - inverse_square
        + 3.0 * inverse_square**2
        - 15.0 * inverse_square**3
        + 105.0 * inverse_square**4
    )
    return (
        -0.5 * z * z
        - math.log(z)
        - 0.5 * math.log(2.0 * math.pi)
        + math.log(correction)
    )


def tie_corrected_wilcoxon(
    pooled_values: np.ndarray,
    pooled_ranks: np.ndarray,
    n_zero: int,
    n_one: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Two-sided Mann-Whitney/Wilcoxon rank-sum with tie correction.

    Group 1 is status=1. A 0.5 continuity correction is applied toward the null
    mean, matching the conventional asymptotic implementation.
    """
    rank_sum_one = pooled_ranks[n_zero:, :].sum(axis=0)
    u_one = rank_sum_one - n_one * (n_one + 1) / 2.0
    null_mean = n_zero * n_one / 2.0
    z_values = np.empty(pooled_values.shape[1], dtype=float)
    log10_p = np.empty(pooled_values.shape[1], dtype=float)
    p_numeric = np.empty(pooled_values.shape[1], dtype=float)

    total = n_zero + n_one
    for column in range(pooled_values.shape[1]):
        _, tie_counts = np.unique(pooled_values[:, column], return_counts=True)
        tie_sum = np.sum(tie_counts.astype(float) ** 3 - tie_counts)
        variance = (
            n_zero
            * n_one
            / 12.0
            * ((total + 1.0) - tie_sum / (total * (total - 1.0)))
        )
        difference = u_one[column] - null_mean
        corrected = difference - 0.5 * np.sign(difference)
        z_values[column] = corrected / math.sqrt(variance)
        log_p = math.log(2.0) + normal_log_sf(abs(z_values[column]))
        log10_p[column] = log_p / math.log(10.0)
        p_numeric[column] = 0.0 if log_p < math.log(np.finfo(float).tiny) else math.exp(log_p)
    return u_one, z_values, p_numeric, log10_p


def ecdf_source_data(data: pd.DataFrame) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    frames = []
    for tissue in TISSUES:
        for status in (0, 1):
            subset = data.loc[data[tissue] == status]
            for rank_name, column in RANK_COLUMNS.items():
                values = np.sort(subset[column].to_numpy(dtype=float))
                cumulative = np.searchsorted(values, thresholds, side="right")
                frames.append(
                    pd.DataFrame(
                        {
                            "Tissue": tissue,
                            "Age_related_status": status,
                            "Rank": rank_name,
                            "N_BSJ": len(values),
                            "Percentile_rank_threshold": thresholds,
                            "Cumulative_percent": cumulative / len(values) * 100.0,
                        }
                    )
                )
    return pd.concat(frames, ignore_index=True)


def strict_qc(
    qc: dict[str, int],
    counts: dict[str, tuple[int, int, int]],
    observed: dict[str, np.ndarray],
) -> None:
    if qc != EXPECTED_QC:
        raise RuntimeError(f"Fixed circAge input QC does not match: {qc}")
    if counts != EXPECTED_COUNTS:
        raise RuntimeError(f"circAge group counts do not match: {counts}")
    rank_names = list(RANK_COLUMNS)
    for tissue, values in observed.items():
        expected = np.array([EXPECTED_RBC[(tissue, rank)] for rank in rank_names])
        if not np.allclose(values, expected, atol=5e-7, rtol=0):
            raise RuntimeError(f"Observed RBC controls failed for {tissue}.")


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, qc = load_analysis_data(args.input)
    rank_columns = list(RANK_COLUMNS.values())

    counts = {
        tissue: (
            int(data[tissue].eq(0).sum()),
            int(data[tissue].eq(1).sum()),
            int(data[tissue].isna().sum()),
        )
        for tissue in TISSUES
    }

    tissue_objects = {}
    observed_rbcs = {}
    for tissue in TISSUES:
        zero = data.loc[data[tissue] == 0, rank_columns].to_numpy(dtype=float)
        one = data.loc[data[tissue] == 1, rank_columns].to_numpy(dtype=float)
        pooled = np.vstack([zero, one])
        pooled_ranks = average_rank_matrix(pooled)
        observed_rbcs[tissue] = rbc_from_rank_sum(
            pooled_ranks[len(zero):, :].sum(axis=0), len(one), len(zero)
        )
        tissue_objects[tissue] = (zero, one, pooled, pooled_ranks)

    if not args.relaxed_qc:
        strict_qc(qc, counts, observed_rbcs)

    ecdf = ecdf_source_data(data)
    descriptive_rows = []
    wilcoxon_rows = []
    rbc_summary_rows = []
    null_frames = []

    for tissue in TISSUES:
        zero, one, pooled, pooled_ranks = tissue_objects[tissue]
        n_zero, n_one = len(zero), len(one)

        for status, values in ((0, zero), (1, one)):
            for rank_index, rank_name in enumerate(RANK_COLUMNS):
                q1, median, q3 = np.quantile(
                    values[:, rank_index], [0.25, 0.50, 0.75]
                )
                descriptive_rows.append(
                    {
                        "Tissue": tissue,
                        "Age_related_status": status,
                        "Rank": rank_name,
                        "N_BSJ": len(values),
                        "Mean": values[:, rank_index].mean(),
                        "Q1": q1,
                        "Median": median,
                        "Q3": q3,
                    }
                )

        u_values, z_values, p_values, log10_p = tie_corrected_wilcoxon(
            pooled, pooled_ranks, n_zero, n_one
        )
        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            wilcoxon_rows.append(
                {
                    "Tissue": tissue,
                    "Rank": rank_name,
                    "N_status_0": n_zero,
                    "N_status_1": n_one,
                    "Mann_Whitney_U_status_1": u_values[rank_index],
                    "Tie_corrected_Z": z_values[rank_index],
                    "Two_sided_asymptotic_P": (
                        "<1e-300" if log10_p[rank_index] < -300
                        else p_values[rank_index]
                    ),
                    "Log10_two_sided_P": log10_p[rank_index],
                    "Continuity_correction": 0.5,
                }
            )

        rng = np.random.default_rng(TISSUE_SEEDS[tissue])
        null = np.empty((args.randomizations, len(RANK_COLUMNS)), dtype=float)
        for iteration in range(args.randomizations):
            randomized_one = rng.choice(
                len(pooled), size=n_one, replace=False, shuffle=False
            )
            null[iteration, :] = rbc_from_rank_sum(
                pooled_ranks[randomized_one, :].sum(axis=0), n_one, n_zero
            )

        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            null_values = null[:, rank_index]
            observed = observed_rbcs[tissue][rank_index]
            extreme = int(np.count_nonzero(null_values >= observed))
            rbc_summary_rows.append(
                {
                    "Tissue": tissue,
                    "Rank": rank_name,
                    "N_0": n_zero,
                    "N_1": n_one,
                    "Observed_RBC_1_vs_0": observed,
                    "Null_95pct_low": np.quantile(null_values, 0.025),
                    "Null_95pct_high": np.quantile(null_values, 0.975),
                    "Null_mean": null_values.mean(),
                    "Null_median": np.median(null_values),
                    "Extreme_randomizations_ge_observed": extreme,
                    "Empirical_one_sided_P_plus1": (
                        extreme + 1
                    ) / (args.randomizations + 1),
                    "N_randomizations": args.randomizations,
                    "Random_seed": TISSUE_SEEDS[tissue],
                }
            )
            null_frames.append(
                pd.DataFrame(
                    {
                        "Tissue": tissue,
                        "Rank": rank_name,
                        "Iteration": np.arange(1, args.randomizations + 1),
                        "Null_RBC": null_values,
                        "Random_seed": TISSUE_SEEDS[tissue],
                    }
                )
            )
        print(f"Randomization completed: {tissue}", flush=True)

    pd.DataFrame(
        [{"QC_item": item, "N_rows": value} for item, value in qc.items()]
    ).to_csv(args.output_directory / "circAge_QC.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "Tissue": tissue,
                "N_0_after_benchmark_exclusion": values[0],
                "N_1_after_benchmark_exclusion": values[1],
                "N_blank_after_benchmark_exclusion": values[2],
            }
            for tissue, values in counts.items()
        ]
    ).to_csv(
        args.output_directory / "circAge_five_tissue_group_counts.tsv",
        sep="\t",
        index=False,
    )
    ecdf.to_csv(
        args.output_directory / "circAge_five_tissue_ECDF_source_data_0.1.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(descriptive_rows).to_csv(
        args.output_directory / "circAge_five_tissue_rank_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(wilcoxon_rows).to_csv(
        args.output_directory / "circAge_five_tissue_Wilcoxon_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(rbc_summary_rows).to_csv(
        args.output_directory / "circAge_five_tissue_RBC_randomization_summary.tsv",
        sep="\t",
        index=False,
    )
    if not args.no_null_output:
        pd.concat(null_frames, ignore_index=True).to_csv(
            args.output_directory /
            f"circAge_five_tissue_RBC_pooled_label_{args.randomizations}_null.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
