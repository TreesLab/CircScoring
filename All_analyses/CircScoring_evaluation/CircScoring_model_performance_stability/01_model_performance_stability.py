#!/usr/bin/env python3
"""Reproduce the CircScoring model-performance stability analysis.

The analysis evaluates four *fixed* CircScoring score vectors in 1,000
independent, unstratified, random 80:20 partitions of the P1N1 dataset.
Models are not refitted and hyperparameters are not re-optimized.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


MASTER_SEED = 20260801
N_ITERATIONS = 1000
TRAIN_FRACTION = 0.80
SCORE_COLUMNS = {
    "XGBoost_CS-R": "XGBoost_CS-R_score",
    "XGBoost_CS-C": "XGBoost_CS-C_score",
    "ElasticNet_CS-R": "ElasticNet_CS-R_score",
    "ElasticNet_CS-C": "ElasticNet_CS-C_score",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def prepare_p1n1(data: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    required = {"P1", "N1", *SCORE_COLUMNS.values()}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    p1 = pd.to_numeric(data["P1"], errors="raise").to_numpy()
    n1 = pd.to_numeric(data["N1"], errors="raise").to_numpy()
    keep = (p1 == 1) | (n1 == 1)
    selected = data.loc[keep].reset_index(drop=True).copy()
    if ((selected["P1"] == 1) & (selected["N1"] == 1)).any():
        raise ValueError("P1 and N1 must be mutually exclusive.")

    y = (pd.to_numeric(selected["P1"], errors="raise") == 1).astype(int).to_numpy()
    for column in SCORE_COLUMNS.values():
        selected[column] = pd.to_numeric(selected[column], errors="raise")
        if not np.isfinite(selected[column].to_numpy()).all():
            raise ValueError(f"Non-finite score detected in {column}.")

    if len(selected) != 1301 or int(y.sum()) != 797 or int((1 - y).sum()) != 504:
        raise ValueError(
            "Unexpected P1N1 composition; expected n=1,301 (797 P1 and 504 N1), "
            f"observed n={len(selected)} ({int(y.sum())} P1 and {int((1-y).sum())} N1)."
        )
    return selected, y


def auroc(y: np.ndarray, score: np.ndarray) -> float:
    """AUROC by the Mann-Whitney formulation with average ranks for ties."""
    ranks = pd.Series(score).rank(method="average").to_numpy(dtype=float)
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC requires both outcome classes.")
    return float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def auprc(y: np.ndarray, score: np.ndarray) -> float:
    """Trapezoidal area under the threshold-based precision-recall curve.

    Observations are sorted by decreasing score. Curve points are evaluated at
    the end of each tied-score block, with the conventional origin (recall=0,
    precision=1) prepended. This is deliberately not average precision.
    """
    order = np.argsort(-score, kind="stable")
    ys = y[order]
    scores = score[order]
    tp = np.cumsum(ys)
    fp = np.cumsum(1 - ys)
    threshold_ends = np.r_[np.flatnonzero(np.diff(scores) != 0), len(scores) - 1]
    recall = tp[threshold_ends] / tp[-1]
    precision = tp[threshold_ends] / (tp[threshold_ends] + fp[threshold_ends])
    recall = np.r_[0.0, recall]
    precision = np.r_[1.0, precision]
    return float(np.trapezoid(precision, recall))


def metric_pair(y: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    return auroc(y, score), auprc(y, score)


def run_analysis(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    p1n1, y = prepare_p1n1(data)
    scores = {name: p1n1[col].to_numpy(dtype=float) for name, col in SCORE_COLUMNS.items()}
    n = len(p1n1)
    n80 = int(np.floor(TRAIN_FRACTION * n))

    master_rng = np.random.default_rng(MASTER_SEED)
    seeds = master_rng.integers(1, np.iinfo(np.int32).max, size=N_ITERATIONS)
    rows: list[dict[str, float | int]] = []

    for iteration, seed_value in enumerate(seeds, start=1):
        seed = int(seed_value)
        order = np.random.default_rng(seed).permutation(n)
        idx80, idx20 = order[:n80], order[n80:]
        y80, y20 = y[idx80], y[idx20]
        row: dict[str, float | int] = {
            "Iteration": iteration,
            "Random_seed": seed,
            "Subset80_n": len(idx80),
            "Subset80_positive_n": int(y80.sum()),
            "Subset80_negative_n": int(len(y80) - y80.sum()),
            "Subset20_n": len(idx20),
            "Subset20_positive_n": int(y20.sum()),
            "Subset20_negative_n": int(len(y20) - y20.sum()),
        }
        for name, score in scores.items():
            auc80, pr80 = metric_pair(y80, score[idx80])
            auc20, pr20 = metric_pair(y20, score[idx20])
            row[f"{name}_AUROC_80"] = auc80
            row[f"{name}_AUROC_20"] = auc20
            row[f"{name}_AUROC_diff_80minus20"] = auc80 - auc20
            row[f"{name}_AUPRC_80"] = pr80
            row[f"{name}_AUPRC_20"] = pr20
            row[f"{name}_AUPRC_diff_80minus20"] = pr80 - pr20
        rows.append(row)

    results = pd.DataFrame(rows)
    reference_rows = []
    for name, score in scores.items():
        model, score_name = name.split("_", maxsplit=1)
        roc, pr = metric_pair(y, score)
        reference_rows.append(
            {
                "Model": model,
                "Score": score_name,
                "n": n,
                "Positive_n": int(y.sum()),
                "Negative_n": int((1 - y).sum()),
                "AUROC": roc,
                "AUPRC": pr,
            }
        )
    return results, pd.DataFrame(reference_rows)


def summary_statistics(results: pd.DataFrame) -> pd.DataFrame:
    metric_columns = list(results.columns[8:])
    records = []
    for column in metric_columns:
        x = results[column]
        records.append(
            {
                "Metric": column,
                "Mean": x.mean(),
                "SD": x.std(ddof=1),
                "Median": x.median(),
                "2.5th percentile": x.quantile(0.025),
                "97.5th percentile": x.quantile(0.975),
                "Minimum": x.min(),
                "Maximum": x.max(),
            }
        )
    return pd.DataFrame(records)


def verify_against_expected(observed: pd.DataFrame, expected_path: Path) -> None:
    if expected_path.suffix.lower() in {".xlsx", ".xlsm"}:
        expected = pd.read_excel(expected_path, sheet_name="1000_partitions")
    else:
        expected = pd.read_csv(expected_path)
    if list(observed.columns) != list(expected.columns):
        raise AssertionError("Observed and expected column names/order differ.")
    if observed.shape != expected.shape:
        raise AssertionError(f"Shape differs: observed {observed.shape}, expected {expected.shape}.")

    exact_columns = list(observed.columns[:8])
    for column in exact_columns:
        if not np.array_equal(observed[column].to_numpy(), expected[column].to_numpy()):
            raise AssertionError(f"Exact verification failed for {column}.")
    numeric_columns = list(observed.columns[8:])
    max_abs = float(
        np.nanmax(np.abs(observed[numeric_columns].to_numpy() - expected[numeric_columns].to_numpy()))
    )
    if max_abs > 5e-14:
        raise AssertionError(f"Numerical verification failed; maximum absolute difference={max_abs:.3g}.")
    print(f"Verification passed: maximum absolute metric difference={max_abs:.3g}")


def write_workbook(
    output_path: Path,
    results: pd.DataFrame,
    summary: pd.DataFrame,
    reference: pd.DataFrame,
    input_path: Path,
) -> None:
    readme = pd.DataFrame(
        [
            ("Analysis", "Model performance stability using 1,000 independent random 80:20 partitions of P1N1"),
            ("Input file", input_path.name),
            ("Input SHA-256", sha256(input_path)),
            ("P1N1 sample size", 1301),
            ("P1 positive BSJs", 797),
            ("N1 negative BSJs", 504),
            ("Partition method", "Unstratified simple random partition without replacement within each iteration"),
            ("Subset sizes", "80% subset: n=1,040; 20% subset: n=261"),
            ("Iterations", N_ITERATIONS),
            ("Master random seed", MASTER_SEED),
            ("Model fitting", "No refitting or hyperparameter optimization; four fixed scores were evaluated"),
            ("AUROC", "Mann-Whitney/rank-based AUROC; ties assigned average ranks"),
            ("AUPRC", "Trapezoidal integration of the threshold-based precision-recall curve"),
        ],
        columns=["Item", "Description"],
    )
    metadata = pd.DataFrame(
        [
            ("source_path", str(input_path.resolve())),
            ("source_sha256", sha256(input_path)),
            ("master_seed", MASTER_SEED),
            ("iterations", N_ITERATIONS),
            ("partition_type", "unstratified simple random 80:20 without replacement"),
            ("implementation", "01_model_performance_stability.py"),
        ],
        columns=["Key", "Value"],
    )
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="README", index=False)
        results.to_excel(writer, sheet_name="1000_partitions", index=False)
        summary.to_excel(writer, sheet_name="Summary_statistics", index=False)
        reference.to_excel(writer, sheet_name="Full_P1N1_reference", index=False)
        metadata.to_excel(writer, sheet_name="_metadata", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True, help="CSV/TSV containing P1, N1 and four fixed score columns")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--expected", type=Path, help="Optional retained 1000-split CSV for exact verification")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = read_table(args.input)
    results, reference = run_analysis(data)
    summary = summary_statistics(results)
    if args.expected:
        verify_against_expected(results, args.expected)
    results.to_csv(args.outdir / "model_performance_stability_1000_partitions.csv", index=False)
    summary.to_csv(args.outdir / "model_performance_stability_summary.csv", index=False)
    reference.to_csv(args.outdir / "model_performance_stability_full_P1N1_reference.csv", index=False)
    write_workbook(
        args.outdir / "model_performance_stability_1000_complete.xlsx",
        results,
        summary,
        reference,
        args.input,
    )
    print(f"Completed {N_ITERATIONS} partitions; outputs written to {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
