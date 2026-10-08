#!/usr/bin/env python3
"""CircScoring validation using disease-association evidence.

Numerical workflow only: benchmark exclusion, group definitions, ECDF source
data, descriptive summaries, observed rank-biserial correlations (RBCs), and
10,000 background-based size-matched resamples for RBC and median rank.
No figure is generated.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_RESAMPLES = 10_000
DEFAULT_SEED = 20260827

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

EVIDENCE_ALIASES = {
    "circad": ["circad"],
    "LncRNADisease": ["LncRNADisease"],
    "CircR2Disease": ["CircR2Disease"],
    "Circ2Disease": ["Circ2Disease"],
    "unique_disease_count": [
        "unique_disease_count",
        "circRNADisease-unique_disease_count",
        "circRNADisease_unique_disease_count",
    ],
}

GROUP_ORDER = [
    "CircR2Disease",
    "Circ2Disease",
    "LncRNADisease",
    "circad",
    "Unique diseases =1",
    "Unique diseases =2",
    "Unique diseases >=3",
]

EXPECTED_QC = {
    "Source rows": 4_455_197,
    "Benchmark-overlapping rows excluded": 2_282,
    "Incomplete-rank rows excluded after benchmark exclusion": 0,
    "Benchmark-free complete-rank rows": 4_452_915,
    "Background rows": 4_449_090,
}

EXPECTED_GROUP_COUNTS = {
    "CircR2Disease": 1_275,
    "Circ2Disease": 166,
    "LncRNADisease": 196,
    "circad": 391,
    "Unique diseases =1": 2_572,
    "Unique diseases =2": 505,
    "Unique diseases >=3": 387,
}

EXPECTED_RBC = {
    ("CircR2Disease", "XGBoost CS-R"): 0.613396,
    ("CircR2Disease", "XGBoost CS-C"): 0.710855,
    ("CircR2Disease", "Elastic Net CS-R"): 0.596895,
    ("CircR2Disease", "Elastic Net CS-C"): 0.666302,
    ("Circ2Disease", "XGBoost CS-R"): 0.759566,
    ("Circ2Disease", "XGBoost CS-C"): 0.761455,
    ("Circ2Disease", "Elastic Net CS-R"): 0.742787,
    ("Circ2Disease", "Elastic Net CS-C"): 0.751667,
    ("LncRNADisease", "XGBoost CS-R"): 0.617018,
    ("LncRNADisease", "XGBoost CS-C"): 0.754385,
    ("LncRNADisease", "Elastic Net CS-R"): 0.593060,
    ("LncRNADisease", "Elastic Net CS-C"): 0.712468,
    ("circad", "XGBoost CS-R"): 0.667732,
    ("circad", "XGBoost CS-C"): 0.746976,
    ("circad", "Elastic Net CS-R"): 0.645847,
    ("circad", "Elastic Net CS-C"): 0.723924,
    ("Unique diseases =1", "XGBoost CS-R"): 0.667589,
    ("Unique diseases =1", "XGBoost CS-C"): 0.767055,
    ("Unique diseases =1", "Elastic Net CS-R"): 0.653562,
    ("Unique diseases =1", "Elastic Net CS-C"): 0.725763,
    ("Unique diseases =2", "XGBoost CS-R"): 0.764787,
    ("Unique diseases =2", "XGBoost CS-C"): 0.772310,
    ("Unique diseases =2", "Elastic Net CS-R"): 0.759434,
    ("Unique diseases =2", "Elastic Net CS-C"): 0.760513,
    ("Unique diseases >=3", "XGBoost CS-R"): 0.821341,
    ("Unique diseases >=3", "XGBoost CS-C"): 0.795855,
    ("Unique diseases >=3", "Elastic Net CS-R"): 0.817166,
    ("Unique diseases >=3", "Elastic Net CS-C"): 0.780512,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Disease-association ECDF, RBC and median-rank analysis with "
            "background-based size-matched resampling."
        )
    )
    parser.add_argument("input", type=Path, help="Annotated circIG TSV/CSV file")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument(
        "--resamples", type=int, default=DEFAULT_RESAMPLES,
        help="Number of background resamples (default: 10000)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--relaxed-qc", action="store_true",
        help="Do not stop if retained-data QC values differ",
    )
    parser.add_argument(
        "--no-null-output", action="store_true",
        help="Do not write iteration-level gzip-compressed null tables",
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
    first_fields = first.split(delimiter)
    second_fields = second.split(delimiter)
    required = set(RANK_COLUMNS.values())
    if required.issubset(second_fields):
        return 1, delimiter
    if required.issubset(first_fields):
        return 0, delimiter
    raise ValueError("Could not identify the header containing the four rank columns")


def resolve_columns(path: Path, skiprows: int, delimiter: str) -> dict[str, str]:
    header = pd.read_csv(path, sep=delimiter, skiprows=skiprows, nrows=0).columns.tolist()
    missing = [c for c in [*RANK_COLUMNS.values(), *BENCHMARK_COLUMNS] if c not in header]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    resolved = {}
    for logical, aliases in EVIDENCE_ALIASES.items():
        match = next((name for name in aliases if name in header), None)
        if match is None:
            raise ValueError(f"Missing evidence column for {logical}; accepted names: {aliases}")
        resolved[logical] = match
    return resolved


def load_data(path: Path):
    skiprows, delimiter = detect_header(path)
    evidence_columns = resolve_columns(path, skiprows, delimiter)
    usecols = [
        *RANK_COLUMNS.values(), *BENCHMARK_COLUMNS, *evidence_columns.values()
    ]
    data = pd.read_csv(
        path, sep=delimiter, skiprows=skiprows, usecols=usecols, low_memory=False
    )
    data = data.apply(pd.to_numeric, errors="coerce")

    benchmark_overlap = data[BENCHMARK_COLUMNS].fillna(0).ne(0).any(axis=1)
    complete_ranks = data[list(RANK_COLUMNS.values())].notna().all(axis=1)
    eligible = (~benchmark_overlap) & complete_ranks

    evidence = pd.DataFrame(
        {logical: data[column].fillna(0) for logical, column in evidence_columns.items()}
    )
    background = eligible & evidence.eq(0).all(axis=1)
    groups = {
        "CircR2Disease": eligible & evidence["CircR2Disease"].eq(1),
        "Circ2Disease": eligible & evidence["Circ2Disease"].eq(1),
        "LncRNADisease": eligible & evidence["LncRNADisease"].eq(1),
        "circad": eligible & evidence["circad"].eq(1),
        "Unique diseases =1": eligible & evidence["unique_disease_count"].eq(1),
        "Unique diseases =2": eligible & evidence["unique_disease_count"].eq(2),
        "Unique diseases >=3": eligible & evidence["unique_disease_count"].ge(3),
    }
    qc = {
        "Source rows": len(data),
        "Benchmark-overlapping rows excluded": int(benchmark_overlap.sum()),
        "Incomplete-rank rows excluded after benchmark exclusion": int(
            ((~benchmark_overlap) & (~complete_ranks)).sum()
        ),
        "Benchmark-free complete-rank rows": int(eligible.sum()),
        "Background rows": int(background.sum()),
    }
    return data, background, groups, qc


def observed_rbc(positive: np.ndarray, background: np.ndarray) -> np.ndarray:
    """RBC for positive vs background; positive values favor the positive set."""
    result = np.empty(positive.shape[1], dtype=float)
    n_background = len(background)
    for column in range(positive.shape[1]):
        reference = np.sort(background[:, column])
        values = positive[:, column]
        lower = np.searchsorted(reference, values, side="left")
        upper = np.searchsorted(reference, values, side="right")
        u = np.sum(lower + 0.5 * (upper - lower), dtype=float)
        result[column] = 2.0 * u / (len(values) * n_background) - 1.0
    return result


def background_average_ranks(background: np.ndarray) -> np.ndarray:
    return pd.DataFrame(background).rank(method="average", axis=0).to_numpy(float)


def ecdf_table(group_values: dict[str, np.ndarray]) -> pd.DataFrame:
    thresholds = np.round(np.arange(0.0, 100.0 + 0.05, 0.1), 1)
    rows = []
    for group, matrix in group_values.items():
        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            values = np.sort(matrix[:, rank_index])
            percentages = np.searchsorted(values, thresholds, side="right") / len(values) * 100
            rows.append(pd.DataFrame({
                "Group": group,
                "Rank": rank_name,
                "Threshold": thresholds,
                "ECDF_percent": percentages,
                "N": len(values),
            }))
    return pd.concat(rows, ignore_index=True)


def strict_qc(qc, groups, observed):
    errors = []
    for key, expected in EXPECTED_QC.items():
        if qc[key] != expected:
            errors.append(f"{key}: observed {qc[key]}, expected {expected}")
    for group, expected in EXPECTED_GROUP_COUNTS.items():
        observed_count = int(groups[group].sum())
        if observed_count != expected:
            errors.append(f"{group}: observed N={observed_count}, expected N={expected}")
    for group in GROUP_ORDER:
        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            expected = EXPECTED_RBC[(group, rank_name)]
            value = observed[group][rank_index]
            if abs(value - expected) > 5e-7:
                errors.append(
                    f"{group} / {rank_name}: observed RBC={value}, expected {expected}"
                )
    if errors:
        raise RuntimeError("Strict retained-data QC failed:\n- " + "\n- ".join(errors))


def main() -> None:
    args = parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    data, background_mask, group_masks, qc = load_data(args.input)
    rank_columns = list(RANK_COLUMNS.values())
    background = data.loc[background_mask, rank_columns].to_numpy(float)
    positive_values = {
        group: data.loc[group_masks[group], rank_columns].to_numpy(float)
        for group in GROUP_ORDER
    }
    observed = {
        group: observed_rbc(positive_values[group], background)
        for group in GROUP_ORDER
    }
    if not args.relaxed_qc:
        strict_qc(qc, group_masks, observed)

    group_values = {"Background": background, **positive_values}
    count_rows = [
        {"Group": group, "N": len(matrix)} for group, matrix in group_values.items()
    ]
    descriptive_rows = []
    for group, matrix in group_values.items():
        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            values = matrix[:, rank_index]
            descriptive_rows.append({
                "Group": group,
                "Rank": rank_name,
                "N": len(values),
                "Mean": float(np.mean(values)),
                "Q1": float(np.quantile(values, 0.25)),
                "Median": float(np.median(values)),
                "Q3": float(np.quantile(values, 0.75)),
            })

    overlap_rows = []
    for left_index, left in enumerate(GROUP_ORDER):
        for right in GROUP_ORDER[left_index + 1:]:
            overlap_rows.append({
                "Group_1": left,
                "Group_2": right,
                "N_overlap": int((group_masks[left] & group_masks[right]).sum()),
            })

    # One RNG is used in the preserved group order. shuffle=False is required
    # for exact reproduction of the retained NumPy Generator stream.
    rng = np.random.default_rng(args.seed)
    background_ranks = background_average_ranks(background)
    n_background = len(background)
    rbc_null_wide = {"resample_id": np.arange(1, args.resamples + 1)}
    median_null_wide = {"resample_id": np.arange(1, args.resamples + 1)}
    rbc_summary_rows = []
    median_summary_rows = []

    for group in GROUP_ORDER:
        n_positive = len(positive_values[group])
        n_remaining = n_background - n_positive
        null_rbc = np.empty((args.resamples, len(RANK_COLUMNS)), dtype=float)
        null_median = np.empty_like(null_rbc)
        for iteration in range(args.resamples):
            sampled = rng.choice(
                n_background, size=n_positive, replace=False, shuffle=False
            )
            rank_sum = background_ranks[sampled, :].sum(axis=0)
            u = rank_sum - n_positive * (n_positive + 1) / 2.0
            null_rbc[iteration, :] = 2.0 * u / (n_positive * n_remaining) - 1.0
            null_median[iteration, :] = np.median(background[sampled, :], axis=0)

        for rank_index, rank_name in enumerate(RANK_COLUMNS):
            label = f"{group} | {rank_name}"
            rbc_null_wide[label] = null_rbc[:, rank_index]
            median_null_wide[label] = null_median[:, rank_index]

            observed_rbc_value = observed[group][rank_index]
            observed_median = float(np.median(positive_values[group][:, rank_index]))
            rbc_extreme = int(np.sum(null_rbc[:, rank_index] >= observed_rbc_value))
            median_extreme = int(np.sum(null_median[:, rank_index] >= observed_median))
            rbc_summary_rows.append({
                "Evidence": group,
                "Rank": rank_name,
                "N_positive": n_positive,
                "N_background": n_background,
                "N_background_remaining_per_resample": n_remaining,
                "Observed_RBC": observed_rbc_value,
                "Null_95pct_low": float(np.quantile(null_rbc[:, rank_index], 0.025)),
                "Null_95pct_high": float(np.quantile(null_rbc[:, rank_index], 0.975)),
                "Null_mean": float(np.mean(null_rbc[:, rank_index])),
                "Null_median": float(np.median(null_rbc[:, rank_index])),
                "Extreme_resamples_ge_observed": rbc_extreme,
                "Empirical_one_sided_P_plus1": (rbc_extreme + 1) / (args.resamples + 1),
                "N_resamples": args.resamples,
                "Random_seed": args.seed,
            })
            median_summary_rows.append({
                "Evidence": group,
                "Rank": rank_name,
                "N_positive": n_positive,
                "N_background": n_background,
                "Observed_positive_median": observed_median,
                "Null_95pct_low": float(np.quantile(null_median[:, rank_index], 0.025)),
                "Null_95pct_high": float(np.quantile(null_median[:, rank_index], 0.975)),
                "Null_mean": float(np.mean(null_median[:, rank_index])),
                "Null_median": float(np.median(null_median[:, rank_index])),
                "Extreme_resamples_ge_observed": median_extreme,
                "Empirical_one_sided_P_plus1": (median_extreme + 1) / (args.resamples + 1),
                "N_resamples": args.resamples,
                "Random_seed": args.seed,
            })
        print(f"Resampling completed: {group}", flush=True)

    pd.DataFrame([{"Metric": key, "Value": value} for key, value in qc.items()]).to_csv(
        args.output_directory / "disease_association_QC.tsv", sep="\t", index=False
    )
    pd.DataFrame(count_rows).to_csv(
        args.output_directory / "disease_association_group_counts.tsv", sep="\t", index=False
    )
    pd.DataFrame(overlap_rows).to_csv(
        args.output_directory / "disease_association_positive_group_overlap_counts.tsv",
        sep="\t", index=False,
    )
    pd.DataFrame(descriptive_rows).to_csv(
        args.output_directory / "disease_association_rank_summary.tsv", sep="\t", index=False
    )
    ecdf_table(group_values).to_csv(
        args.output_directory / "disease_association_ECDF_thresholds_0.1.tsv",
        sep="\t", index=False,
    )
    pd.DataFrame(rbc_summary_rows).to_csv(
        args.output_directory / "disease_association_RBC_resampling_summary.tsv",
        sep="\t", index=False,
    )
    pd.DataFrame(median_summary_rows).to_csv(
        args.output_directory / "disease_association_median_resampling_summary.tsv",
        sep="\t", index=False,
    )
    if not args.no_null_output:
        pd.DataFrame(rbc_null_wide).to_csv(
            args.output_directory /
            f"disease_association_RBC_{args.resamples}_null_distributions.tsv.gz",
            sep="\t", index=False,
        )
        pd.DataFrame(median_null_wide).to_csv(
            args.output_directory /
            f"disease_association_median_{args.resamples}_null_distributions.tsv.gz",
            sep="\t", index=False,
        )
    print(f"Analysis completed: {args.output_directory}")


if __name__ == "__main__":
    main()
