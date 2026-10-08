#!/usr/bin/env python3
"""CircScoring validation of library preparation and reproducibility.

The script performs the numerical analyses only. It does not create figures.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260917
DEFAULT_RESAMPLES = 10_000

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

EXPECTED_COUNTS = {
    "Rows after benchmark exclusion": 4_452_915,
    "CircAI only": 299_005,
    "CircRic only": 30_720,
    "CircAI-CircRic shared": 3_892,
    "CSCD2 common": 291_809,
    "CSCD2 normal": 1_207_723,
    "CSCD2 cancer": 806_147,
}

EXPECTED_REPORTED_RBC = {
    ("CircAI only vs CircRic only", "XGBoost CS-R"): 0.83196412,
    ("CircAI only vs CircRic only", "XGBoost CS-C"): 0.34511617,
    ("CircAI only vs CircRic only", "Elastic Net CS-R"): 0.88579810,
    ("CircAI only vs CircRic only", "Elastic Net CS-C"): 0.30710553,
    ("Shared vs CircRic only", "XGBoost CS-R"): 0.9395284330571256,
    ("Shared vs CircRic only", "XGBoost CS-C"): 0.7391762566509292,
    ("Shared vs CircRic only", "Elastic Net CS-R"): 0.9577871157315219,
    ("Shared vs CircRic only", "Elastic Net CS-C"): 0.7513776924888660,
    ("Common vs Normal", "XGBoost CS-R"): 0.5288970798202597,
    ("Common vs Normal", "XGBoost CS-C"): 0.5138003102934152,
    ("Common vs Normal", "Elastic Net CS-R"): 0.5789637370208487,
    ("Common vs Normal", "Elastic Net CS-C"): 0.4849935857204195,
    ("Common vs Cancer", "XGBoost CS-R"): 0.5595793524643036,
    ("Common vs Cancer", "XGBoost CS-C"): 0.6332440619448263,
    ("Common vs Cancer", "Elastic Net CS-R"): 0.6478093406138741,
    ("Common vs Cancer", "Elastic Net CS-C"): 0.5705971083815806,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "CircAI-CircRic and CSCD2 ECDF, RBC, and pooled-label "
            "randomization analyses."
        )
    )
    parser.add_argument("input", type=Path, help="One- or two-row-header circIG TSV[.gz].")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--resamples", type=int, default=DEFAULT_RESAMPLES)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--all-pairs",
        action="store_true",
        help="Also calculate Shared versus CircAI only (not in the retained summary table).",
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
    if args.resamples < 1:
        parser.error("--resamples must be a positive integer")
    return args


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if path.suffix == ".gz" else path.open("rt", encoding="utf-8-sig", newline="")


def detect_header(path: Path) -> tuple[int, str]:
    with open_text(path) as handle:
        first = handle.readline().rstrip("\r\n")
        second = handle.readline().rstrip("\r\n")
    delimiter = "\t" if "\t" in first else ","
    h1 = first.split(delimiter)
    h2 = second.split(delimiter)
    required = {
        "CircAI", "CircRic", "CSCD2_group",
        *RANK_COLUMNS.values(), *BENCHMARK_COLUMNS,
    }
    if required.issubset(h1):
        return 0, delimiter
    if required.issubset(h2):
        return 1, delimiter
    missing = sorted(required.difference(set(h1) | set(h2)))
    raise ValueError(f"Required columns were not found: {', '.join(missing)}")


def numeric_blank_zero(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float)


def load_analysis_data(path: Path) -> tuple[pd.DataFrame, int, int]:
    skiprows, delimiter = detect_header(path)
    usecols = [
        "CircAI", "CircRic", "CSCD2_group",
        *RANK_COLUMNS.values(), *BENCHMARK_COLUMNS,
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
    excluded_benchmark_rows = int(benchmark_overlap.sum())
    data = data.loc[~benchmark_overlap].copy()

    for column in RANK_COLUMNS.values():
        data[column] = pd.to_numeric(data[column], errors="coerce")
    incomplete = data[list(RANK_COLUMNS.values())].isna().any(axis=1)
    incomplete_rank_rows = int(incomplete.sum())
    data = data.loc[~incomplete].copy()

    data["CircAI"] = numeric_blank_zero(data["CircAI"])
    data["CircRic"] = numeric_blank_zero(data["CircRic"])
    data["Library_group"] = "Neither"
    data.loc[
        (data["CircAI"] == 1) & (data["CircRic"] == 0), "Library_group"
    ] = "CircAI only"
    data.loc[
        (data["CircAI"] == 0) & (data["CircRic"] == 1), "Library_group"
    ] = "CircRic only"
    data.loc[
        (data["CircAI"] == 1) & (data["CircRic"] == 1), "Library_group"
    ] = "CircAI-CircRic shared"
    data["CSCD2_group_clean"] = (
        data["CSCD2_group"].fillna("").astype(str).str.strip().str.lower()
    )
    return data, excluded_benchmark_rows, incomplete_rank_rows


def pooled_average_ranks(group_1: np.ndarray, group_2: np.ndarray) -> np.ndarray:
    pooled = np.vstack([group_1, group_2])
    return pd.DataFrame(pooled).rank(method="average", axis=0).to_numpy(dtype=float)


def rbc_from_selected_positions(
    pooled_ranks: np.ndarray, selected_group_1: np.ndarray, n_group_1: int
) -> np.ndarray:
    n_group_2 = len(pooled_ranks) - n_group_1
    rank_sum = pooled_ranks[selected_group_1, :].sum(axis=0)
    u_group_1 = rank_sum - n_group_1 * (n_group_1 + 1) / 2
    return 2 * u_group_1 / (n_group_1 * n_group_2) - 1


def analyze_comparison(
    data: pd.DataFrame,
    analysis: str,
    group_column: str,
    comparison: str,
    group_1_name: str,
    group_2_name: str,
    rng: np.random.Generator,
    n_resamples: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = list(RANK_COLUMNS.values())
    group_1 = data.loc[data[group_column] == group_1_name, columns].to_numpy(dtype=float)
    group_2 = data.loc[data[group_column] == group_2_name, columns].to_numpy(dtype=float)
    n1, n2 = len(group_1), len(group_2)
    pooled_ranks = pooled_average_ranks(group_1, group_2)
    observed = rbc_from_selected_positions(pooled_ranks, np.arange(n1), n1)

    null = np.empty((n_resamples, len(columns)), dtype=float)
    for iteration in range(n_resamples):
        selected = rng.choice(n1 + n2, size=n1, replace=False)
        null[iteration, :] = rbc_from_selected_positions(pooled_ranks, selected, n1)

    summary_rows = []
    null_rows = []
    for rank_index, rank_name in enumerate(RANK_COLUMNS):
        null_values = null[:, rank_index]
        extreme = int(np.count_nonzero(null_values >= observed[rank_index]))
        summary_rows.append(
            {
                "Analysis": analysis,
                "Comparison": comparison,
                "Group_1": group_1_name,
                "Group_2": group_2_name,
                "Rank": rank_name,
                "N_group_1": n1,
                "N_group_2": n2,
                "Observed_RBC": observed[rank_index],
                "Null_95pct_low": np.quantile(null_values, 0.025),
                "Null_95pct_high": np.quantile(null_values, 0.975),
                "Null_mean": null_values.mean(),
                "Null_median": np.median(null_values),
                "Extreme_permutations_ge_observed": extreme,
                "Empirical_one_sided_P_plus1": (extreme + 1) / (n_resamples + 1),
                "N_resampling": n_resamples,
                "Random_seed": int(rng_seed),
            }
        )
        null_rows.append(
            pd.DataFrame(
                {
                    "Analysis": analysis,
                    "Comparison": comparison,
                    "Rank": rank_name,
                    "Iteration": np.arange(1, n_resamples + 1),
                    "Null_RBC": null_values,
                }
            )
        )
    return pd.DataFrame(summary_rows), pd.concat(null_rows, ignore_index=True)


def make_ecdf_table(
    data: pd.DataFrame,
    analysis: str,
    group_column: str,
    group_names: list[str],
    output_group_names: dict[str, str] | None = None,
) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    frames = []
    output_group_names = output_group_names or {}
    for group_name in group_names:
        group = data.loc[data[group_column] == group_name]
        for rank_name, column in RANK_COLUMNS.items():
            values = np.sort(group[column].to_numpy(dtype=float))
            cumulative = np.searchsorted(values, thresholds, side="right")
            frames.append(
                pd.DataFrame(
                    {
                        "Analysis": analysis,
                        "Group": output_group_names.get(group_name, group_name),
                        "Rank": rank_name,
                        "Threshold": thresholds,
                        "Cumulative_count": cumulative,
                        "ECDF_percent": 100 * cumulative / len(values),
                        "N_group": len(values),
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


def strict_qc(data: pd.DataFrame, summary: pd.DataFrame) -> None:
    observed_counts = {
        "Rows after benchmark exclusion": len(data),
        "CircAI only": int((data["Library_group"] == "CircAI only").sum()),
        "CircRic only": int((data["Library_group"] == "CircRic only").sum()),
        "CircAI-CircRic shared": int(
            (data["Library_group"] == "CircAI-CircRic shared").sum()
        ),
        "CSCD2 common": int((data["CSCD2_group_clean"] == "common").sum()),
        "CSCD2 normal": int((data["CSCD2_group_clean"] == "normal").sum()),
        "CSCD2 cancer": int((data["CSCD2_group_clean"] == "cancer").sum()),
    }
    if observed_counts != EXPECTED_COUNTS:
        raise RuntimeError(f"Fixed-population group counts do not match: {observed_counts}")
    for _, row in summary.iterrows():
        key = (row["Comparison"], row["Rank"])
        if key in EXPECTED_REPORTED_RBC:
            if not np.isclose(row["Observed_RBC"], EXPECTED_REPORTED_RBC[key], atol=5e-8):
                raise RuntimeError(f"Observed RBC does not match retained result for {key}.")


def main() -> None:
    global rng_seed
    args = parse_args()
    rng_seed = args.seed
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, benchmark_excluded, incomplete_ranks = load_analysis_data(args.input)

    comparisons = [
        (
            "CircAI_CircRic", "Library_group",
            "CircAI only vs CircRic only", "CircAI only", "CircRic only",
        ),
        (
            "CircAI_CircRic", "Library_group",
            "Shared vs CircRic only", "CircAI-CircRic shared", "CircRic only",
        ),
        (
            "CSCD2_group", "CSCD2_group_clean",
            "Common vs Normal", "common", "normal",
        ),
        (
            "CSCD2_group", "CSCD2_group_clean",
            "Common vs Cancer", "common", "cancer",
        ),
    ]
    if args.all_pairs:
        comparisons.insert(
            1,
            (
                "CircAI_CircRic", "Library_group",
                "Shared vs CircAI only", "CircAI-CircRic shared", "CircAI only",
            ),
        )

    rng = np.random.default_rng(args.seed)
    summaries, null_tables = [], []
    for analysis, group_column, comparison, group_1, group_2 in comparisons:
        summary, null = analyze_comparison(
            data, analysis, group_column, comparison,
            group_1, group_2, rng, args.resamples
        )
        summaries.append(summary)
        null_tables.append(null)
    summary = pd.concat(summaries, ignore_index=True)

    if not args.relaxed_qc:
        strict_qc(data, summary)

    group_counts = pd.DataFrame(
        {
            "Group": [
                "CircAI only", "CircRic only", "CircAI-CircRic shared",
                "CSCD2 common", "CSCD2 normal-specific", "CSCD2 cancer-specific",
            ],
            "Definition": [
                "CircAI = 1 and CircRic = 0",
                "CircAI = 0 and CircRic = 1",
                "CircAI = 1 and CircRic = 1",
                "CSCD2_group = common",
                "CSCD2_group = normal",
                "CSCD2_group = cancer",
            ],
            "N_after_benchmark_exclusion_and_complete_ranks": [
                int((data["Library_group"] == "CircAI only").sum()),
                int((data["Library_group"] == "CircRic only").sum()),
                int((data["Library_group"] == "CircAI-CircRic shared").sum()),
                int((data["CSCD2_group_clean"] == "common").sum()),
                int((data["CSCD2_group_clean"] == "normal").sum()),
                int((data["CSCD2_group_clean"] == "cancer").sum()),
            ],
        }
    )
    qc = pd.DataFrame(
        {
            "Item": [
                "Source rows",
                "Benchmark-overlap rows excluded",
                "Incomplete-rank rows excluded after benchmark exclusion",
                "Rows after benchmark and rank exclusions",
            ],
            "Value": [
                len(data) + benchmark_excluded + incomplete_ranks,
                benchmark_excluded,
                incomplete_ranks,
                len(data),
            ],
        }
    )

    summary.to_csv(
        args.output_directory / "Library_preparation_reproducibility_RBC_summary.tsv",
        sep="\t", index=False,
    )
    group_counts.to_csv(
        args.output_directory / "Library_preparation_reproducibility_group_counts.tsv",
        sep="\t", index=False,
    )
    qc.to_csv(
        args.output_directory / "Library_preparation_reproducibility_QC.tsv",
        sep="\t", index=False,
    )
    make_ecdf_table(
        data,
        "CircAI_CircRic",
        "Library_group",
        ["CircRic only", "CircAI only", "CircAI-CircRic shared"],
    ).to_csv(
        args.output_directory / "CircAI_CircRic_ECDF_data_0.1.tsv", sep="\t", index=False
    )
    make_ecdf_table(
        data,
        "CSCD2_group",
        "CSCD2_group_clean",
        ["cancer", "normal", "common"],
        {"cancer": "Cancer-specific", "normal": "Normal-specific", "common": "Common"},
    ).to_csv(
        args.output_directory / "CSCD2_group_ECDF_data_0.1.tsv", sep="\t", index=False
    )
    if not args.no_null_output:
        pd.concat(null_tables, ignore_index=True).to_csv(
            args.output_directory /
            f"Library_preparation_reproducibility_RBC_null_{args.resamples}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
