#!/usr/bin/env python3
"""CircScoring coding-potential validation using riboCIRC and TransCirc.

The two resources are analyzed independently and written to separate ECDF,
RBC-summary, and null-distribution files. No figure is generated.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


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

COLUMN_ALIASES = {
    "riboCIRC": ["riboCIRC"],
    "annotation-guided": [
        "riboCIRC-annotation-guided",
        "riboCIRC_annotation-guided",
        "annotation-guided",
    ],
    "context-specific": [
        "riboCIRC-context-specific",
        "riboCIRC_context-specific",
        "context-specific",
    ],
    "literature-reported": [
        "riboCIRC-literature-reported",
        "riboCIRC_literature-reported",
        "literature-reported",
    ],
    "TransCirc_evidences_num": [
        "TransCirc-evidences_type_number",
        "TransCirc_evidences_num",
        "TransCirc_evidences_type_number",
    ],
}

RIBOCIRC_EVIDENCE = {
    "Annotation-guided": "annotation-guided",
    "Context-specific": "context-specific",
    "Literature-reported": "literature-reported",
}

RIBOCIRC_SEEDS = {
    "Annotation-guided": 2026091701,
    "Context-specific": 2026091702,
    "Literature-reported": 2026091703,
}

TRANSCIRC_SEEDS = {
    "=1": 2026091601,
    "=2": 2026091602,
    "=3": 2026091603,
    "=4": 2026091604,
    ">=5": 2026091605,
}

EXPECTED_QC = {
    "Source rows": 4_455_197,
    "Benchmark-overlapping rows excluded": 2_282,
    "Incomplete-rank rows excluded": 0,
    "Rows retained after exclusions": 4_452_915,
}

EXPECTED_RIBOCIRC_COUNTS = {
    "Background": 4_449_774,
    "Annotation-guided": 3_085,
    "Context-specific": 56,
    "Literature-reported": 52,
}

EXPECTED_TRANSCIRC_COUNTS = {
    "blank": 4_162_531,
    "0": 929,
    "1": 14_122,
    "2": 16_485,
    "3": 219_878,
    "4": 37_309,
    ">=5": 1_661,
}

EXPECTED_RIBOCIRC_RBC = {
    ("Annotation-guided", "XGBoost CS-R"): 0.1476657406465527,
    ("Annotation-guided", "XGBoost CS-C"): 0.4425679143208743,
    ("Annotation-guided", "Elastic Net CS-R"): 0.1388716171165421,
    ("Annotation-guided", "Elastic Net CS-C"): 0.3597114427851347,
    ("Context-specific", "XGBoost CS-R"): 0.9072466015770047,
    ("Context-specific", "XGBoost CS-C"): 0.8544823247524160,
    ("Context-specific", "Elastic Net CS-R"): 0.9286341002936329,
    ("Context-specific", "Elastic Net CS-C"): 0.8763250672955525,
    ("Literature-reported", "XGBoost CS-R"): 0.7763111288175708,
    ("Literature-reported", "XGBoost CS-C"): 0.8331165505000064,
    ("Literature-reported", "Elastic Net CS-R"): 0.7703180802855640,
    ("Literature-reported", "Elastic Net CS-C"): 0.8413519860351766,
}

EXPECTED_TRANSCIRC_RBC = {
    ("=1", "XGBoost CS-R"): 0.129063,
    ("=1", "XGBoost CS-C"): 0.182079,
    ("=1", "Elastic Net CS-R"): 0.152978,
    ("=1", "Elastic Net CS-C"): 0.193771,
    ("=2", "XGBoost CS-R"): 0.158804,
    ("=2", "XGBoost CS-C"): 0.319843,
    ("=2", "Elastic Net CS-R"): 0.186931,
    ("=2", "Elastic Net CS-C"): 0.304103,
    ("=3", "XGBoost CS-R"): 0.386538,
    ("=3", "XGBoost CS-C"): 0.626283,
    ("=3", "Elastic Net CS-R"): 0.427214,
    ("=3", "Elastic Net CS-C"): 0.608373,
    ("=4", "XGBoost CS-R"): 0.428980,
    ("=4", "XGBoost CS-C"): 0.676888,
    ("=4", "Elastic Net CS-R"): 0.469458,
    ("=4", "Elastic Net CS-C"): 0.660239,
    (">=5", "XGBoost CS-R"): 0.617339,
    (">=5", "XGBoost CS-C"): 0.726241,
    (">=5", "Elastic Net CS-R"): 0.635315,
    (">=5", "Elastic Net CS-C"): 0.714539,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Independent riboCIRC and TransCirc ECDF/RBC coding-potential "
            "validation workflows."
        )
    )
    parser.add_argument("input", type=Path, help="One- or two-row-header circIG TSV[.gz].")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--resamples", type=int, default=DEFAULT_RESAMPLES)
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


def select_header(path: Path) -> tuple[int, str, list[str]]:
    with open_text(path) as handle:
        first = handle.readline().rstrip("\r\n")
        second = handle.readline().rstrip("\r\n")
    delimiter = "\t" if "\t" in first else ","
    first_header = first.split(delimiter)
    second_header = second.split(delimiter)
    fixed = {*RANK_COLUMNS.values(), *BENCHMARK_COLUMNS}
    if fixed.issubset(first_header):
        return 0, delimiter, first_header
    if fixed.issubset(second_header):
        return 1, delimiter, second_header
    missing = sorted(fixed.difference(set(first_header) | set(second_header)))
    raise ValueError(f"Required columns were not found: {', '.join(missing)}")


def resolve_aliases(header: list[str]) -> dict[str, str]:
    resolved = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        matches = [name for name in aliases if name in header]
        if not matches:
            raise ValueError(
                f"No input column found for {canonical}; accepted names: {aliases}"
            )
        resolved[canonical] = matches[0]
    return resolved


def numeric_blank_zero(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").fillna(0).to_numpy(dtype=float)


def load_analysis_data(path: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    skiprows, delimiter, header = select_header(path)
    aliases = resolve_aliases(header)
    usecols = [
        *RANK_COLUMNS.values(),
        *BENCHMARK_COLUMNS,
        *dict.fromkeys(aliases.values()),
    ]
    data = pd.read_csv(
        path,
        sep=delimiter,
        skiprows=skiprows,
        usecols=usecols,
        low_memory=False,
    )
    rename = {actual: canonical for canonical, actual in aliases.items()}
    data = data.rename(columns=rename)
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

    for column in COLUMN_ALIASES:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    qc = {
        "Source rows": source_rows,
        "Benchmark-overlapping rows excluded": benchmark_excluded,
        "Incomplete-rank rows excluded": incomplete_excluded,
        "Rows retained after exclusions": len(data),
    }
    return data, qc


def average_rank_matrix(values: np.ndarray) -> np.ndarray:
    return pd.DataFrame(values).rank(method="average", axis=0).to_numpy(dtype=float)


def rbc_from_group_1_ranks(
    rank_sum_group_1: np.ndarray,
    n_group_1: int,
    n_group_2: int,
) -> np.ndarray:
    u_statistic = rank_sum_group_1 - n_group_1 * (n_group_1 + 1) / 2.0
    return 2.0 * u_statistic / (n_group_1 * n_group_2) - 1.0


def observed_rbc_vs_background(
    positive: np.ndarray,
    background: np.ndarray,
) -> np.ndarray:
    result = np.empty(positive.shape[1], dtype=float)
    for column in range(positive.shape[1]):
        sorted_background = np.sort(background[:, column])
        lower = np.searchsorted(sorted_background, positive[:, column], side="left")
        upper = np.searchsorted(sorted_background, positive[:, column], side="right")
        u_statistic = np.sum(lower + 0.5 * (upper - lower), dtype=float)
        result[column] = (
            2.0 * u_statistic / (len(positive) * len(background)) - 1.0
        )
    return result


def summarize_group(values: np.ndarray, group_name: str) -> pd.DataFrame:
    rows = []
    for index, rank_name in enumerate(RANK_COLUMNS):
        q1, median, q3 = np.quantile(values[:, index], [0.25, 0.50, 0.75])
        rows.append(
            {
                "Group": group_name,
                "Rank": rank_name,
                "N": len(values),
                "Mean": values[:, index].mean(),
                "Q1": q1,
                "Median": median,
                "Q3": q3,
            }
        )
    return pd.DataFrame(rows)


def make_ribocirc_ecdf(
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    frames = []
    for group_name, values in groups.items():
        for index, rank_name in enumerate(RANK_COLUMNS):
            sorted_values = np.sort(values[:, index])
            cumulative = np.searchsorted(sorted_values, thresholds, side="right")
            frames.append(
                pd.DataFrame(
                    {
                        "Group": group_name,
                        "Rank": rank_name,
                        "Percentile_rank_threshold": thresholds,
                        "Cumulative_percent": cumulative / len(values) * 100.0,
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


def make_transcirc_ecdf(
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    rows = []
    labels = ["0", "1", "2", "3", "4", ">=5"]
    for index, rank_name in enumerate(RANK_COLUMNS):
        cumulative_by_group = {}
        for label in labels:
            values = np.sort(groups[label][:, index])
            cumulative_by_group[label] = (
                np.searchsorted(values, thresholds, side="right") / len(values) * 100.0
            )
        for threshold_index, threshold in enumerate(thresholds):
            row = {"Rank": rank_name, "Threshold": threshold}
            for label in labels:
                display = "TransCirc >= 5" if label == ">=5" else f"TransCirc = {label}"
                row[display] = cumulative_by_group[label][threshold_index]
            rows.append(row)
    return pd.DataFrame(rows)


def analyze_ribocirc(
    data: pd.DataFrame,
    n_resamples: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rank_columns = list(RANK_COLUMNS.values())
    positive_masks = {
        evidence: data[column].fillna(0).eq(1).to_numpy()
        for evidence, column in RIBOCIRC_EVIDENCE.items()
    }
    any_positive = np.logical_or.reduce(list(positive_masks.values()))
    background_mask = data["riboCIRC"].fillna(0).ne(1).to_numpy() & ~any_positive
    groups = {
        "Background": data.loc[background_mask, rank_columns].to_numpy(dtype=float)
    }
    groups.update(
        {
            evidence: data.loc[mask, rank_columns].to_numpy(dtype=float)
            for evidence, mask in positive_masks.items()
        }
    )

    background = groups["Background"]
    background_ranks = average_rank_matrix(background)
    summary_frames = []
    null_frames = []
    for evidence in RIBOCIRC_EVIDENCE:
        positive = groups[evidence]
        observed = observed_rbc_vs_background(positive, background)
        rng = np.random.default_rng(RIBOCIRC_SEEDS[evidence])
        null = np.empty((n_resamples, len(RANK_COLUMNS)), dtype=float)
        for iteration in range(n_resamples):
            selected = rng.choice(
                len(background), size=len(positive), replace=False, shuffle=False
            )
            rank_sum = background_ranks[selected, :].sum(axis=0)
            null[iteration, :] = rbc_from_group_1_ranks(
                rank_sum, len(positive), len(background) - len(positive)
            )

        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            null_values = null[:, rank_index]
            extreme = int(np.count_nonzero(null_values >= observed[rank_index]))
            summary_frames.append(
                {
                    "Evidence": evidence,
                    "Rank": rank_name,
                    "N_positive": len(positive),
                    "N_background": len(background),
                    "Observed_RBC_positive_vs_background": observed[rank_index],
                    "Null_95pct_low": np.quantile(null_values, 0.025),
                    "Null_95pct_high": np.quantile(null_values, 0.975),
                    "Null_mean": null_values.mean(),
                    "Null_median": np.median(null_values),
                    "Extreme_resamples_ge_observed": extreme,
                    "Empirical_one_sided_P_plus1": (extreme + 1) / (n_resamples + 1),
                    "N_resamples": n_resamples,
                    "Random_seed": RIBOCIRC_SEEDS[evidence],
                }
            )
            null_frames.append(
                pd.DataFrame(
                    {
                        "Evidence": evidence,
                        "Rank": rank_name,
                        "Iteration": np.arange(1, n_resamples + 1),
                        "Null_RBC": null_values,
                        "Random_seed": RIBOCIRC_SEEDS[evidence],
                    }
                )
            )
        print(f"riboCIRC resampling completed: {evidence}", flush=True)

    counts = pd.DataFrame(
        [{"Group": group, "N": len(values)} for group, values in groups.items()]
    )
    group_summary = pd.concat(
        [summarize_group(values, group) for group, values in groups.items()],
        ignore_index=True,
    )
    ecdf = make_ribocirc_ecdf(groups)
    return (
        pd.DataFrame(summary_frames),
        pd.concat(null_frames, ignore_index=True),
        counts,
        group_summary,
        ecdf,
    )


def transcirc_masks(data: pd.DataFrame) -> dict[str, np.ndarray]:
    values = data["TransCirc_evidences_num"].to_numpy(dtype=float)
    return {
        "blank": np.isnan(values),
        "0": values == 0,
        "1": values == 1,
        "2": values == 2,
        "3": values == 3,
        "4": values == 4,
        ">=5": values >= 5,
    }


def analyze_transcirc(
    data: pd.DataFrame,
    n_resamples: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rank_columns = list(RANK_COLUMNS.values())
    masks = transcirc_masks(data)
    groups = {
        label: data.loc[mask, rank_columns].to_numpy(dtype=float)
        for label, mask in masks.items()
    }
    n0 = groups["0"]
    summary_rows = []
    null_frames = []
    for positive_label in ["1", "2", "3", "4", ">=5"]:
        positive = groups[positive_label]
        # Historical ordering is N0 followed by the positive set. Under the
        # null, N0 labels are reassigned and the complementary rows receive the
        # positive label. This exactly preserves both observed group sizes.
        pooled = np.vstack([n0, positive])
        pooled_ranks = average_rank_matrix(pooled)
        n_zero = len(n0)
        n_positive = len(positive)
        total_rank_sum = pooled_ranks.sum(axis=0)
        observed_positive_rank_sum = pooled_ranks[n_zero:, :].sum(axis=0)
        observed = rbc_from_group_1_ranks(
            observed_positive_rank_sum, n_positive, n_zero
        )

        rng = np.random.default_rng(TRANSCIRC_SEEDS[f"={positive_label}" if positive_label != ">=5" else ">=5"])
        null = np.empty((n_resamples, len(RANK_COLUMNS)), dtype=float)
        for iteration in range(n_resamples):
            randomized_zero = rng.choice(
                len(pooled), size=n_zero, replace=False, shuffle=False
            )
            randomized_positive_rank_sum = (
                total_rank_sum - pooled_ranks[randomized_zero, :].sum(axis=0)
            )
            null[iteration, :] = rbc_from_group_1_ranks(
                randomized_positive_rank_sum, n_positive, n_zero
            )

        seed_key = f"={positive_label}" if positive_label != ">=5" else ">=5"
        seed = TRANSCIRC_SEEDS[seed_key]
        display_label = f"={positive_label}" if positive_label != ">=5" else ">=5"
        comparison = f"TransCirc_evidences_num {display_label} vs 0"
        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            null_values = null[:, rank_index]
            extreme = int(np.count_nonzero(null_values >= observed[rank_index]))
            summary_rows.append(
                {
                    "Comparison": comparison,
                    "Positive_group": display_label,
                    "Rank": rank_name,
                    "N_0": n_zero,
                    "N_positive": n_positive,
                    "Observed_RBC_positive_vs_0": observed[rank_index],
                    "Null_95pct_low": np.quantile(null_values, 0.025),
                    "Null_95pct_high": np.quantile(null_values, 0.975),
                    "Null_mean": null_values.mean(),
                    "Null_median": np.median(null_values),
                    "Extreme_permutations_ge_observed": extreme,
                    "Empirical_one_sided_P_plus1": (extreme + 1) / (n_resamples + 1),
                    "N_permutations": n_resamples,
                    "Random_seed": seed,
                }
            )
            null_frames.append(
                pd.DataFrame(
                    {
                        "Positive_group": display_label,
                        "Rank": rank_name,
                        "Iteration": np.arange(1, n_resamples + 1),
                        "Null_RBC": null_values,
                        "Random_seed": seed,
                    }
                )
            )
        print(f"TransCirc randomization completed: {display_label} versus 0", flush=True)

    counts = pd.DataFrame(
        [
            {"TransCirc_evidences_num": label, "N_after_benchmark_exclusion": len(values)}
            for label, values in groups.items()
        ]
    )
    group_summary = pd.concat(
        [
            summarize_group(groups[label], f"TransCirc {label}")
            for label in ["0", "1", "2", "3", "4", ">=5"]
        ],
        ignore_index=True,
    )
    ecdf = make_transcirc_ecdf(groups)
    return (
        pd.DataFrame(summary_rows),
        pd.concat(null_frames, ignore_index=True),
        counts,
        group_summary,
        ecdf,
    )


def strict_qc(
    qc: dict[str, int],
    ribo_summary: pd.DataFrame,
    ribo_counts: pd.DataFrame,
    trans_summary: pd.DataFrame,
    trans_counts: pd.DataFrame,
) -> None:
    if qc != EXPECTED_QC:
        raise RuntimeError(f"Fixed input QC does not match: {qc}")

    observed_ribo_counts = ribo_counts.set_index("Group")["N"].to_dict()
    if observed_ribo_counts != EXPECTED_RIBOCIRC_COUNTS:
        raise RuntimeError(f"riboCIRC counts do not match: {observed_ribo_counts}")
    observed_trans_counts = (
        trans_counts.set_index("TransCirc_evidences_num")[
            "N_after_benchmark_exclusion"
        ].to_dict()
    )
    if observed_trans_counts != EXPECTED_TRANSCIRC_COUNTS:
        raise RuntimeError(f"TransCirc counts do not match: {observed_trans_counts}")

    for row in ribo_summary.itertuples(index=False):
        expected = EXPECTED_RIBOCIRC_RBC[(row.Evidence, row.Rank)]
        if not np.isclose(
            row.Observed_RBC_positive_vs_background, expected, atol=5e-12, rtol=0
        ):
            raise RuntimeError(f"riboCIRC RBC control failed: {row.Evidence}, {row.Rank}")
    for row in trans_summary.itertuples(index=False):
        expected = EXPECTED_TRANSCIRC_RBC[(row.Positive_group, row.Rank)]
        if not np.isclose(
            row.Observed_RBC_positive_vs_0, expected, atol=5e-7, rtol=0
        ):
            raise RuntimeError(
                f"TransCirc RBC control failed: {row.Positive_group}, {row.Rank}"
            )


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, qc = load_analysis_data(args.input)

    (
        ribo_summary,
        ribo_null,
        ribo_counts,
        ribo_group_summary,
        ribo_ecdf,
    ) = analyze_ribocirc(data, args.resamples)
    (
        trans_summary,
        trans_null,
        trans_counts,
        trans_group_summary,
        trans_ecdf,
    ) = analyze_transcirc(data, args.resamples)

    if not args.relaxed_qc:
        strict_qc(qc, ribo_summary, ribo_counts, trans_summary, trans_counts)

    pd.DataFrame(
        [{"QC_item": item, "N_rows": value} for item, value in qc.items()]
    ).to_csv(args.output_directory / "coding_potential_QC.tsv", sep="\t", index=False)

    # riboCIRC outputs remain separate from TransCirc outputs.
    ribo_counts.to_csv(
        args.output_directory / "riboCIRC_group_counts.tsv", sep="\t", index=False
    )
    ribo_group_summary.to_csv(
        args.output_directory / "riboCIRC_group_rank_summary.tsv", sep="\t", index=False
    )
    ribo_ecdf.to_csv(
        args.output_directory / "riboCIRC_three_evidence_ECDF_data.tsv",
        sep="\t",
        index=False,
    )
    ribo_summary.to_csv(
        args.output_directory / "riboCIRC_three_evidence_RBC_summary.tsv",
        sep="\t",
        index=False,
    )

    trans_counts.to_csv(
        args.output_directory / "TransCirc_group_counts.tsv", sep="\t", index=False
    )
    trans_group_summary.to_csv(
        args.output_directory / "TransCirc_group_rank_summary.tsv", sep="\t", index=False
    )
    trans_ecdf.to_csv(
        args.output_directory / "TransCirc_evidence_levels_ECDF_data.tsv",
        sep="\t",
        index=False,
    )
    trans_summary.to_csv(
        args.output_directory / "TransCirc_RBC_summary.tsv", sep="\t", index=False
    )

    if not args.no_null_output:
        ribo_null.to_csv(
            args.output_directory /
            f"riboCIRC_three_evidence_RBC_{args.resamples}_null.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        trans_null.to_csv(
            args.output_directory /
            f"TransCirc_RBC_pooled_label_{args.resamples}_null.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
