#!/usr/bin/env python3
"""CircScoring validation using m6A2Circ and circTarget evidence.

Numerical workflow only: group counts, ECDF source data, RBC and median
summaries, and 10,000 size-matched background-resampling iterations.
No figure is generated.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RESAMPLES = 10_000
DEFAULT_SEED = 20260825

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
    "m6A_confidence": [
        "m6A2Circ-confidence_type_number",
        "m6A2Circ#confidence",
        "m6A2Circ_confidence",
    ],
    "circTarget_type_count": [
        "circTarget-#type of Seq (supported by >=2 chimeric reads)",
        ">=2#type of Seq",
        "circTarget->=2#type of Seq",
    ],
}

EVIDENCE_ORDER = [
    "m6A confidence = 1",
    "m6A confidence = 2",
    "m6A confidence >= 3",
    "circTarget",
]

EXPECTED_QC = {
    "Source rows": 4_455_197,
    "Benchmark-overlapping rows excluded": 2_282,
    "Incomplete-rank rows excluded": 0,
    "Rows retained after exclusions": 4_452_915,
}

EXPECTED_COUNTS = {
    "Background": 4_409_523,
    "m6A confidence = 1": 35_543,
    "m6A confidence = 2": 5_272,
    "m6A confidence >= 3": 1_204,
    "circTarget": 2_215,
}

EXPECTED_OVERLAPS = {
    "m6A confidence = 1 AND circTarget": 386,
    "m6A confidence = 2 AND circTarget": 294,
    "m6A confidence >= 3 AND circTarget": 162,
    "Any m6A confidence AND circTarget": 842,
}

EXPECTED_RBC = {
    ("m6A confidence = 1", "XGBoost CS-R"): 0.55377384,
    ("m6A confidence = 1", "XGBoost CS-C"): 0.58580469,
    ("m6A confidence = 1", "Elastic Net CS-R"): 0.57209807,
    ("m6A confidence = 1", "Elastic Net CS-C"): 0.56112549,
    ("m6A confidence = 2", "XGBoost CS-R"): 0.90179763,
    ("m6A confidence = 2", "XGBoost CS-C"): 0.79570450,
    ("m6A confidence = 2", "Elastic Net CS-R"): 0.91164068,
    ("m6A confidence = 2", "Elastic Net CS-C"): 0.79757420,
    ("m6A confidence >= 3", "XGBoost CS-R"): 0.95621953,
    ("m6A confidence >= 3", "XGBoost CS-C"): 0.83597174,
    ("m6A confidence >= 3", "Elastic Net CS-R"): 0.95565328,
    ("m6A confidence >= 3", "Elastic Net CS-C"): 0.84413472,
    ("circTarget", "XGBoost CS-R"): 0.69279395,
    ("circTarget", "XGBoost CS-C"): 0.69315450,
    ("circTarget", "Elastic Net CS-R"): 0.73584055,
    ("circTarget", "Elastic Net CS-C"): 0.70135579,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "m6A2Circ/circTarget ECDF, RBC, median, and size-matched "
            "background-resampling analysis."
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
    data = data.rename(
        columns={actual: canonical for canonical, actual in aliases.items()}
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

    data["m6A_confidence"] = pd.to_numeric(
        data["m6A_confidence"], errors="coerce"
    )
    data["circTarget_type_count"] = pd.to_numeric(
        data["circTarget_type_count"], errors="coerce"
    )

    qc = {
        "Source rows": source_rows,
        "Benchmark-overlapping rows excluded": benchmark_excluded,
        "Incomplete-rank rows excluded": incomplete_excluded,
        "Rows retained after exclusions": len(data),
    }
    return data, qc


def define_groups(data: pd.DataFrame) -> dict[str, np.ndarray]:
    m6a = data["m6A_confidence"].fillna(0).to_numpy(dtype=float)
    circ_target = data["circTarget_type_count"].fillna(0).to_numpy(dtype=float)
    return {
        "Background": (m6a == 0) & (circ_target == 0),
        "m6A confidence = 1": m6a == 1,
        "m6A confidence = 2": m6a == 2,
        "m6A confidence >= 3": m6a >= 3,
        "circTarget": circ_target >= 1,
    }


def average_rank_matrix(values: np.ndarray) -> np.ndarray:
    return pd.DataFrame(values).rank(method="average", axis=0).to_numpy(dtype=float)


def observed_rbc(positive: np.ndarray, background: np.ndarray) -> np.ndarray:
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


def random_sample_rbc(
    background_ranks: np.ndarray,
    selected: np.ndarray,
) -> np.ndarray:
    n_sample = len(selected)
    n_remaining = len(background_ranks) - n_sample
    rank_sum = background_ranks[selected, :].sum(axis=0)
    u_statistic = rank_sum - n_sample * (n_sample + 1) / 2.0
    return 2.0 * u_statistic / (n_sample * n_remaining) - 1.0


def group_rank_summary(groups: dict[str, np.ndarray]) -> pd.DataFrame:
    frames = []
    for group_name, values in groups.items():
        for index, rank_name in enumerate(RANK_COLUMNS):
            q1, median, q3 = np.quantile(values[:, index], [0.25, 0.50, 0.75])
            frames.append(
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
    return pd.DataFrame(frames)


def ecdf_source_data(groups: dict[str, np.ndarray]) -> pd.DataFrame:
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


def strict_qc(
    qc: dict[str, int],
    groups: dict[str, np.ndarray],
    masks: dict[str, np.ndarray],
    observed_rbcs: dict[str, np.ndarray],
) -> None:
    if qc != EXPECTED_QC:
        raise RuntimeError(f"Fixed input QC does not match: {qc}")
    counts = {name: len(values) for name, values in groups.items()}
    if counts != EXPECTED_COUNTS:
        raise RuntimeError(f"Regulatory-evidence group counts do not match: {counts}")
    overlap_counts = {
        "m6A confidence = 1 AND circTarget": int(
            np.count_nonzero(
                masks["m6A confidence = 1"] & masks["circTarget"]
            )
        ),
        "m6A confidence = 2 AND circTarget": int(
            np.count_nonzero(
                masks["m6A confidence = 2"] & masks["circTarget"]
            )
        ),
        "m6A confidence >= 3 AND circTarget": int(
            np.count_nonzero(
                masks["m6A confidence >= 3"] & masks["circTarget"]
            )
        ),
        "Any m6A confidence AND circTarget": int(
            np.count_nonzero(
                (
                    masks["m6A confidence = 1"]
                    | masks["m6A confidence = 2"]
                    | masks["m6A confidence >= 3"]
                )
                & masks["circTarget"]
            )
        ),
    }
    if overlap_counts != EXPECTED_OVERLAPS:
        raise RuntimeError(f"Evidence-overlap counts do not match: {overlap_counts}")
    rank_names = list(RANK_COLUMNS)
    for evidence, values in observed_rbcs.items():
        expected = np.array([EXPECTED_RBC[(evidence, rank)] for rank in rank_names])
        if not np.allclose(values, expected, atol=5e-8, rtol=0):
            raise RuntimeError(f"Observed RBC controls failed for {evidence}.")


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, qc = load_analysis_data(args.input)
    masks = define_groups(data)
    rank_columns = list(RANK_COLUMNS.values())
    groups = {
        name: data.loc[mask, rank_columns].to_numpy(dtype=float)
        for name, mask in masks.items()
    }
    background = groups["Background"]

    observed_rbcs = {
        evidence: observed_rbc(groups[evidence], background)
        for evidence in EVIDENCE_ORDER
    }
    if not args.relaxed_qc:
        strict_qc(qc, groups, masks, observed_rbcs)

    counts = pd.DataFrame(
        [{"Group": name, "N_BSJ": len(values)} for name, values in groups.items()]
    )
    overlap_counts = pd.DataFrame(
        [
            {
                "Overlap": name,
                "N_BSJ": value,
            }
            for name, value in {
                "m6A confidence = 1 AND circTarget": int(
                    np.count_nonzero(masks["m6A confidence = 1"] & masks["circTarget"])
                ),
                "m6A confidence = 2 AND circTarget": int(
                    np.count_nonzero(masks["m6A confidence = 2"] & masks["circTarget"])
                ),
                "m6A confidence >= 3 AND circTarget": int(
                    np.count_nonzero(masks["m6A confidence >= 3"] & masks["circTarget"])
                ),
                "Any m6A confidence AND circTarget": int(
                    np.count_nonzero(
                        (
                            masks["m6A confidence = 1"]
                            | masks["m6A confidence = 2"]
                            | masks["m6A confidence >= 3"]
                        )
                        & masks["circTarget"]
                    )
                ),
            }.items()
        ]
    )
    ranks_summary = group_rank_summary(groups)
    ecdf = ecdf_source_data(groups)

    background_ranks = average_rank_matrix(background)
    rng = np.random.default_rng(args.seed)
    rbc_summary_rows = []
    median_summary_rows = []
    rbc_null_frames = []
    median_null_frames = []

    # A single historical RNG stream is shared across evidence groups in this
    # fixed order. The same sampled rows yield both RBC and median statistics.
    for evidence in EVIDENCE_ORDER:
        positive = groups[evidence]
        n_positive = len(positive)
        observed_rbc_values = observed_rbcs[evidence]
        observed_medians = np.median(positive, axis=0)
        null_rbc = np.empty((args.resamples, len(RANK_COLUMNS)), dtype=float)
        null_median = np.empty((args.resamples, len(RANK_COLUMNS)), dtype=float)

        for iteration in range(args.resamples):
            selected = rng.choice(
                len(background),
                size=n_positive,
                replace=False,
                shuffle=False,
            )
            null_rbc[iteration, :] = random_sample_rbc(background_ranks, selected)
            null_median[iteration, :] = np.median(background[selected, :], axis=0)

        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            rbc_values = null_rbc[:, rank_index]
            rbc_extreme = int(
                np.count_nonzero(rbc_values >= observed_rbc_values[rank_index])
            )
            rbc_summary_rows.append(
                {
                    "Evidence": evidence,
                    "Rank": rank_name,
                    "N_positive": n_positive,
                    "N_background": len(background),
                    "N_background_complement": len(background) - n_positive,
                    "Observed_RBC": observed_rbc_values[rank_index],
                    "Null_95pct_low": np.quantile(rbc_values, 0.025),
                    "Null_95pct_high": np.quantile(rbc_values, 0.975),
                    "Null_mean": rbc_values.mean(),
                    "Null_median": np.median(rbc_values),
                    "Extreme_resamples_ge_observed": rbc_extreme,
                    "Empirical_one_sided_P_plus1": (
                        rbc_extreme + 1
                    ) / (args.resamples + 1),
                    "N_resamples": args.resamples,
                    "Random_seed": args.seed,
                }
            )

            median_values = null_median[:, rank_index]
            median_extreme = int(
                np.count_nonzero(median_values >= observed_medians[rank_index])
            )
            median_summary_rows.append(
                {
                    "Evidence": evidence,
                    "Rank": rank_name,
                    "N_positive": n_positive,
                    "N_background": len(background),
                    "Observed_positive_median": observed_medians[rank_index],
                    "Null_median_2.5pct": np.quantile(median_values, 0.025),
                    "Null_median_97.5pct": np.quantile(median_values, 0.975),
                    "Null_mean_of_random_set_medians": median_values.mean(),
                    "Null_median_of_random_set_medians": np.median(median_values),
                    "Extreme_resamples_ge_observed": median_extreme,
                    "Empirical_one_sided_P_plus1": (
                        median_extreme + 1
                    ) / (args.resamples + 1),
                    "N_resamples": args.resamples,
                    "Random_seed": args.seed,
                }
            )

            rbc_null_frames.append(
                pd.DataFrame(
                    {
                        "Evidence": evidence,
                        "Rank": rank_name,
                        "Permutation": np.arange(1, args.resamples + 1),
                        "RBC_perm": rbc_values,
                    }
                )
            )
            median_null_frames.append(
                pd.DataFrame(
                    {
                        "Evidence": evidence,
                        "Rank": rank_name,
                        "Permutation": np.arange(1, args.resamples + 1),
                        "Random_set_median": median_values,
                    }
                )
            )
        print(f"Resampling completed: {evidence}", flush=True)

    pd.DataFrame(
        [{"QC_item": item, "N_rows": value} for item, value in qc.items()]
    ).to_csv(args.output_directory / "regulatory_evidence_QC.tsv", sep="\t", index=False)
    counts.to_csv(
        args.output_directory / "m6A2Circ_circTarget_group_counts.tsv",
        sep="\t",
        index=False,
    )
    overlap_counts.to_csv(
        args.output_directory / "m6A2Circ_circTarget_overlap_counts.tsv",
        sep="\t",
        index=False,
    )
    ranks_summary.to_csv(
        args.output_directory / "m6A2Circ_circTarget_group_rank_summary.tsv",
        sep="\t",
        index=False,
    )
    ecdf.to_csv(
        args.output_directory / "m6A2Circ_circTarget_ECDF_source_data_0.1.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(rbc_summary_rows).to_csv(
        args.output_directory / "m6A2Circ_circTarget_RBC_resampling_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(median_summary_rows).to_csv(
        args.output_directory / "m6A2Circ_circTarget_median_resampling_summary.tsv",
        sep="\t",
        index=False,
    )

    if not args.no_null_output:
        pd.concat(rbc_null_frames, ignore_index=True).to_csv(
            args.output_directory /
            f"m6A2Circ_circTarget_RBC_{args.resamples}_null.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        pd.concat(median_null_frames, ignore_index=True).to_csv(
            args.output_directory /
            f"m6A2Circ_circTarget_median_{args.resamples}_null.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    print(f"Analysis completed: {args.output_directory.resolve()}")


if __name__ == "__main__":
    main()
