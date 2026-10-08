#!/usr/bin/env python3
"""Reproduce the joint CircScoring confidence landscape analysis.

The analysis constructs 10 x 10 CS-C-by-CS-R matrices for XGBoost and Elastic
Net, compares each of 11 benchmark-free evidence-positive sets with the fixed
circIG reference, performs 10,000 size-matched reference resamples,
tests positive cell enrichment, applies BH correction across 100 cells, and
tests the global distributional difference using total variation distance.

Example
-------
python Joint_confidence_enrichment.py --input circIG.tsv --output results
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


N_RESAMPLES_DEFAULT = 10_000
RANDOM_SEED_BASE = 2_026_091_500

BENCHMARKS = [
    "positive1_P1N1", "positive2_P4N4_L", "positive3_P4N4_H",
    "negative1_P1N1", "negative2_P4N4_L", "negative3_P4N4_H",
]

ALIASES: Mapping[str, Sequence[str]] = {
    "riboCIRC": ["riboCIRC"],
    "TransCirc": ["TransCirc_evidences_num", "TransCirc-evidences_type_number"],
    "m6A": ["m6A2Circ#confidence", "m6A2Circ-confidence_type_number"],
    "circTarget": [
        ">=2#type of Seq",
        "circTarget-#type of Seq (supported by >=2 chimeric reads)",
    ],
    "Disease": ["Disease"],
    "Aortic": ["Aortic", "circAge-Aortic"],
    "Lung": ["Lung", "circAge-Lung"],
    "Skin": ["Skin", "circAge-Skin"],
    "Umbilical": ["Umbilical", "circAge-Umbilical"],
    "White blood cell": ["white_blood_cell", "circAge-White_blood_cell"],
    "RT-independent": ["RT_independent_NCL_event"],
    "Both RT/nonRT": ["Both RT/nonRT"],
    "XGB_R": ["XGBoost_CS-R_percentile_rank"],
    "XGB_C": ["XGBoost_CS-C_percentile_rank"],
    "EN_R": ["ElasticNet_CS-R_percentile_rank"],
    "EN_C": ["ElasticNet_CS-C_percentile_rank"],
    **{name: [name] for name in BENCHMARKS},
}

EVIDENCE_ORDER = [
    "RT-independent", "Both RT/nonRT", "Translation", "m6A2Circ >=1",
    "circTarget >=1", "Aortic", "Lung", "Skin", "Umbilical",
    "White blood cell", "Disease",
]

EVIDENCE_DEFINITIONS = {
    "RT-independent": "RT_independent_NCL_event = 1",
    "Both RT/nonRT": "Both RT/nonRT = 1",
    "Translation": "riboCIRC = 1 OR TransCirc evidence count >= 1",
    "m6A2Circ >=1": "m6A2Circ confidence count >= 1",
    "circTarget >=1": "circTarget #type of Seq supported by >=2 chimeric reads >= 1",
    "Aortic": "Aortic = 1", "Lung": "Lung = 1", "Skin": "Skin = 1",
    "Umbilical": "Umbilical = 1", "White blood cell": "white_blood_cell = 1",
    "Disease": "Disease = 1",
}

MODELS = {
    "XGBoost": ("XGB_R", "XGB_C"),
    "Elastic Net": ("EN_R", "EN_C"),
}

EXPECTED_COUNTS = {
    "N_total": 4_455_197, "N_benchmark_overlap": 2_282,
    "N_reference": 4_142_933, "RT-independent": 2_493,
    "Both RT/nonRT": 599, "Translation": 291_730,
    "m6A2Circ >=1": 42_019, "circTarget >=1": 2_215,
    "Aortic": 142, "Lung": 3_914, "Skin": 702, "Umbilical": 130,
    "White blood cell": 411, "Disease": 3_825,
}


def detect_header_row(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for i in range(10):
            fields = handle.readline().rstrip("\r\n").split("\t")
            if len(fields) >= 4 and fields[:4] == ["chrom", "start", "end", "strand"]:
                return i
    raise ValueError("Could not locate the chrom/start/end/strand header row.")


def resolve_columns(path: Path, header_row: int) -> Dict[str, str]:
    columns = list(pd.read_csv(path, sep="\t", header=header_row, nrows=0).columns)
    resolved: Dict[str, str] = {}
    for logical, alternatives in ALIASES.items():
        match = next((x for x in alternatives if x in columns), None)
        if match is None:
            raise KeyError(f"Missing required column {logical!r}; tried {alternatives}")
        resolved[logical] = match
    return resolved


def load_arrays(path: Path, chunksize: int) -> Dict[str, np.ndarray]:
    header_row = detect_header_row(path)
    resolved = resolve_columns(path, header_row)
    reverse = {actual: logical for logical, actual in resolved.items()}
    pieces: Dict[str, List[np.ndarray]] = {logical: [] for logical in resolved}
    for chunk in pd.read_csv(
        path, sep="\t", header=header_row, usecols=list(reverse),
        chunksize=chunksize, low_memory=False,
    ):
        for actual, logical in reverse.items():
            pieces[logical].append(
                pd.to_numeric(chunk[actual], errors="coerce").to_numpy(np.float64)
            )
    return {name: np.concatenate(parts) for name, parts in pieces.items()}


def eq1(x: np.ndarray) -> np.ndarray:
    return np.isfinite(x) & (x == 1.0)


def ge1(x: np.ndarray) -> np.ndarray:
    return np.isfinite(x) & (x >= 1.0)


def define_populations(data: Mapping[str, np.ndarray]) -> Tuple[np.ndarray, Dict[str, np.ndarray], np.ndarray]:
    benchmark_overlap = np.zeros(len(data["riboCIRC"]), dtype=bool)
    for name in BENCHMARKS:
        benchmark_overlap |= eq1(data[name])
    benchmark_free = ~benchmark_overlap

    raw = {
        "RT-independent": eq1(data["RT-independent"]),
        "Both RT/nonRT": eq1(data["Both RT/nonRT"]),
        "Translation": eq1(data["riboCIRC"]) | ge1(data["TransCirc"]),
        "m6A2Circ >=1": ge1(data["m6A"]),
        "circTarget >=1": ge1(data["circTarget"]),
        "Aortic": eq1(data["Aortic"]), "Lung": eq1(data["Lung"]),
        "Skin": eq1(data["Skin"]), "Umbilical": eq1(data["Umbilical"]),
        "White blood cell": eq1(data["White blood cell"]),
        "Disease": eq1(data["Disease"]),
    }
    positives = {name: benchmark_free & raw[name] for name in EVIDENCE_ORDER}
    any_evidence = np.zeros_like(benchmark_overlap)
    for mask in raw.values():
        any_evidence |= mask
    reference = benchmark_free & ~any_evidence
    return reference, positives, benchmark_overlap


def decile_index(values: np.ndarray) -> np.ndarray:
    """Return 0-based deciles: [0,10), ..., [90,100], with 100 in decile 10."""
    if np.any(~np.isfinite(values)):
        raise ValueError("Missing percentile ranks encountered in an analysis population.")
    return np.clip(np.floor(values / 10.0).astype(np.int16), 0, 9)


def joint_cells(ranks_r: np.ndarray, ranks_c: np.ndarray) -> np.ndarray:
    return decile_index(ranks_c) * 10 + decile_index(ranks_r)


def cell_counts(cells: np.ndarray) -> np.ndarray:
    return np.bincount(cells, minlength=100).astype(np.int64)


def interval(decile_zero_based: int) -> str:
    lo = decile_zero_based * 10
    return f"[{lo},{lo + 10}]" if decile_zero_based == 9 else f"[{lo},{lo + 10})"


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def write_tsv(path: Path, rows: List[Mapping[str, object]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def write_gzip_tsv(path: Path, rows: List[Mapping[str, object]], fields: Sequence[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def resample_cells(reference_counts: np.ndarray, sample_size: int, n_resamples: int,
                   rng: np.random.Generator) -> np.ndarray:
    """Exact sampling without replacement from cell counts."""
    return rng.multivariate_hypergeometric(
        reference_counts, sample_size, size=n_resamples, method="marginals"
    ).astype(np.int64)


def analyze_pair(evidence: str, model: str, positive_counts: np.ndarray,
                 reference_counts: np.ndarray, n_resamples: int, seed: int):
    n_pos = int(positive_counts.sum())
    n_ref = int(reference_counts.sum())
    pos_prop = positive_counts / n_pos
    ref_prop = reference_counts / n_ref
    expected = n_pos * ref_prop
    observed_log2 = np.log2((positive_counts + 0.5) / (expected + 0.5))
    observed_tvd = 0.5 * float(np.abs(pos_prop - ref_prop).sum())

    rng = np.random.default_rng(seed)
    sampled = resample_cells(reference_counts, n_pos, n_resamples, rng)
    remaining = reference_counts[None, :] - sampled
    remaining_n = n_ref - n_pos
    sampled_prop = sampled / n_pos
    remaining_prop = remaining / remaining_n
    null_expected = n_pos * remaining_prop
    null_log2 = np.log2((sampled + 0.5) / (null_expected + 0.5))
    null_tvd = 0.5 * np.abs(sampled_prop - remaining_prop).sum(axis=1)

    extreme = (null_log2 >= observed_log2[None, :]).sum(axis=0)
    raw_upper = (extreme + 1.0) / (n_resamples + 1.0)
    enriched = observed_log2 > 0
    empirical_p = np.where(enriched, raw_upper, 1.0)
    fdr = bh_adjust(empirical_p)

    cell_rows = []
    for cell in range(100):
        c = cell // 10; r = cell % 10
        cell_rows.append({
            "Evidence": evidence, "Model": model,
            "CS_C_decile": c + 1, "CS_C_interval": interval(c),
            "CS_R_decile": r + 1, "CS_R_interval": interval(r),
            "Positive_count": int(positive_counts[cell]),
            "Positive_proportion": float(pos_prop[cell]),
            "Reference_count": int(reference_counts[cell]),
            "Reference_proportion": float(ref_prop[cell]),
            "N_positive": n_pos, "N_reference_pool": n_ref,
            "Expected_positive_count_from_reference": float(expected[cell]),
            "Observed_log2_enrichment": float(observed_log2[cell]),
            "Null_log2_mean": float(null_log2[:, cell].mean()),
            "Null_log2_median": float(np.median(null_log2[:, cell])),
            "Null_log2_95pct_low": float(np.quantile(null_log2[:, cell], 0.025)),
            "Null_log2_95pct_high": float(np.quantile(null_log2[:, cell], 0.975)),
            "Extreme_resamples_ge_observed": int(extreme[cell]),
            "Raw_upper_tail_P_before_direction_filter": float(raw_upper[cell]),
            "Empirical_enrichment_one_sided_P_plus1": float(empirical_p[cell]),
            "BH_FDR_across_100_cells": float(fdr[cell]),
            "Enriched_log2_gt_0": bool(enriched[cell]),
            "FDR_lt_0.05": bool(enriched[cell] and fdr[cell] < 0.05),
            "P_lt_0.05_FDR_ge_0.05": bool(enriched[cell] and empirical_p[cell] < 0.05 and fdr[cell] >= 0.05),
            "N_resamples": n_resamples, "Random_seed": seed,
        })

    tvd_extreme = int((null_tvd >= observed_tvd).sum())
    tvd_row = {
        "Evidence": evidence, "Model": model, "N_positive": n_pos,
        "N_reference_pool": n_ref, "Observed_TVD": observed_tvd,
        "Null_TVD_mean": float(null_tvd.mean()),
        "Null_TVD_median": float(np.median(null_tvd)),
        "Null_TVD_95pct_low": float(np.quantile(null_tvd, 0.025)),
        "Null_TVD_95pct_high": float(np.quantile(null_tvd, 0.975)),
        "Extreme_resamples_ge_observed": tvd_extreme,
        "Empirical_one_sided_P_plus1": (tvd_extreme + 1.0) / (n_resamples + 1.0),
        "N_resamples": n_resamples, "Random_seed": seed,
    }
    null_rows = [
        {"Evidence": evidence, "Model": model, "Iteration": i + 1,
         "Null_TVD": float(value), "Random_seed": seed}
        for i, value in enumerate(null_tvd)
    ]
    return cell_rows, tvd_row, null_rows


def validation_rows(reference: np.ndarray, positives: Mapping[str, np.ndarray],
                    benchmark_overlap: np.ndarray) -> List[Dict[str, object]]:
    observed = {
        "N_total": len(reference), "N_benchmark_overlap": int(benchmark_overlap.sum()),
        "N_reference": int(reference.sum()),
        **{name: int(mask.sum()) for name, mask in positives.items()},
    }
    rows = []
    for key, expected in EXPECTED_COUNTS.items():
        rows.append({"Check": key, "Expected": expected, "Observed": observed[key],
                     "Pass": observed[key] == expected})
    return rows


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--n-resamples", type=int, default=N_RESAMPLES_DEFAULT)
    p.add_argument("--chunksize", type=int, default=250_000)
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--evidence", action="append", choices=EVIDENCE_ORDER)
    p.add_argument("--model", action="append", choices=list(MODELS))
    return p.parse_args()


def main() -> None:
    args = parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    data = load_arrays(args.input, args.chunksize)
    reference, positives, benchmark_overlap = define_populations(data)
    checks = validation_rows(reference, positives, benchmark_overlap)
    write_tsv(args.output / "Validation_checks.tsv", checks,
              ["Check", "Expected", "Observed", "Pass"])
    if not all(row["Pass"] for row in checks):
        raise RuntimeError(f"Population validation failed: {[x for x in checks if not x['Pass']]}")

    summary_rows = [
        {"Evidence": e, "Definition": EVIDENCE_DEFINITIONS[e],
         "N_positive": int(positives[e].sum()), "N_reference": int(reference.sum())}
        for e in EVIDENCE_ORDER
    ]
    write_tsv(args.output / "Summary.tsv", summary_rows,
              ["Evidence", "Definition", "N_positive", "N_reference"])
    if args.validate_only:
        return

    evidence_names = [e for e in EVIDENCE_ORDER if not args.evidence or e in args.evidence]
    model_names = [m for m in MODELS if not args.model or m in args.model]
    reference_tables: Dict[str, np.ndarray] = {}
    reference_rows = []; rank_rows = []

    for model, (r_name, c_name) in MODELS.items():
        ref_cells = joint_cells(data[r_name][reference], data[c_name][reference])
        counts = cell_counts(ref_cells); reference_tables[model] = counts
        for cell in range(100):
            c = cell // 10; r = cell % 10
            reference_rows.append({
                "Model": model, "CS_C_decile": c + 1, "CS_C_interval": interval(c),
                "CS_R_decile": r + 1, "CS_R_interval": interval(r),
                "Reference_count": int(counts[cell]),
                "Population_proportion_percent": 100.0 * counts[cell] / counts.sum(),
            })
        for rank_label, values in [(r_name, data[r_name][reference]), (c_name, data[c_name][reference])]:
            d = decile_index(values); dc = np.bincount(d, minlength=10)
            for j in range(10):
                rank_rows.append({"Group": "reference", "Rank": rank_label,
                                  "Rank_decile": j + 1, "Rank_interval": interval(j),
                                  "Count": int(dc[j]), "Percent": 100.0 * dc[j] / len(d),
                                  "N_group": len(d)})

    write_tsv(args.output / "Reference_10x10.tsv", reference_rows,
              ["Model", "CS_C_decile", "CS_C_interval", "CS_R_decile",
               "CS_R_interval", "Reference_count", "Population_proportion_percent"])

    cell_output = []; tvd_output = []; tvd_null_output = []
    for ei, evidence in enumerate(EVIDENCE_ORDER, start=1):
        if evidence not in evidence_names:
            continue
        for mi, model in enumerate(MODELS, start=1):
            if model not in model_names:
                continue
            r_name, c_name = MODELS[model]
            mask = positives[evidence]
            pos_cells = joint_cells(data[r_name][mask], data[c_name][mask])
            pos_counts = cell_counts(pos_cells)
            seed = RANDOM_SEED_BASE + ei * 10 + mi
            rows, tvd, null_rows = analyze_pair(
                evidence, model, pos_counts, reference_tables[model],
                args.n_resamples, seed,
            )
            cell_output.extend(rows); tvd_output.append(tvd); tvd_null_output.extend(null_rows)
            for rank_label, values in [(r_name, data[r_name][mask]), (c_name, data[c_name][mask])]:
                d = decile_index(values); dc = np.bincount(d, minlength=10)
                for j in range(10):
                    rank_rows.append({"Group": evidence, "Rank": rank_label,
                                      "Rank_decile": j + 1, "Rank_interval": interval(j),
                                      "Count": int(dc[j]), "Percent": 100.0 * dc[j] / len(d),
                                      "N_group": len(d)})

    cell_fields = [
        "Evidence", "Model", "CS_C_decile", "CS_C_interval", "CS_R_decile",
        "CS_R_interval", "Positive_count", "Positive_proportion", "Reference_count",
        "Reference_proportion", "N_positive", "N_reference_pool",
        "Expected_positive_count_from_reference", "Observed_log2_enrichment",
        "Null_log2_mean", "Null_log2_median", "Null_log2_95pct_low",
        "Null_log2_95pct_high", "Extreme_resamples_ge_observed",
        "Raw_upper_tail_P_before_direction_filter", "Empirical_enrichment_one_sided_P_plus1",
        "BH_FDR_across_100_cells", "Enriched_log2_gt_0", "FDR_lt_0.05",
        "P_lt_0.05_FDR_ge_0.05", "N_resamples", "Random_seed",
    ]
    write_gzip_tsv(args.output / "Cell_results.tsv.gz", cell_output, cell_fields)
    write_tsv(args.output / "TVD_summary.tsv", tvd_output, [
        "Evidence", "Model", "N_positive", "N_reference_pool", "Observed_TVD",
        "Null_TVD_mean", "Null_TVD_median", "Null_TVD_95pct_low",
        "Null_TVD_95pct_high", "Extreme_resamples_ge_observed",
        "Empirical_one_sided_P_plus1", "N_resamples", "Random_seed",
    ])
    write_gzip_tsv(args.output / "TVD_10000_null.tsv.gz", tvd_null_output,
                   ["Evidence", "Model", "Iteration", "Null_TVD", "Random_seed"])
    write_tsv(args.output / "Rank_distribution.tsv", rank_rows,
              ["Group", "Rank", "Rank_decile", "Rank_interval", "Count", "Percent", "N_group"])

    figure_rows = []
    for evidence in evidence_names:
        for model in model_names:
            rows = [x for x in cell_output if x["Evidence"] == evidence and x["Model"] == model]
            figure_rows.append({
                "Evidence": evidence, "Model": model,
                "Positive_enrichment_cells": sum(x["Enriched_log2_gt_0"] for x in rows),
                "Displayed_significant_positive_cells": sum(x["FDR_lt_0.05"] or x["P_lt_0.05_FDR_ge_0.05"] for x in rows),
                "FDR_lt_0_05_cells": sum(x["FDR_lt_0.05"] for x in rows),
                "Raw_P_lt_0_05_FDR_ge_0_05_cells": sum(x["P_lt_0.05_FDR_ge_0.05"] for x in rows),
            })
    write_tsv(args.output / "Figure_points.tsv", figure_rows,
              ["Evidence", "Model", "Positive_enrichment_cells",
               "Displayed_significant_positive_cells", "FDR_lt_0_05_cells",
               "Raw_P_lt_0_05_FDR_ge_0_05_cells"])

    methods = [
        {"Item": "Analysis population", "Value": "Each evidence-positive set was analyzed separately after excluding overlap with any of the six benchmark sets."},
        {"Item": "Translation potential", "Value": "riboCIRC=1 or TransCirc evidence count>=1."},
        {"Item": "Reference", "Value": "Benchmark-free and negative for all 11 evidence definitions."},
        {"Item": "Joint landscape", "Value": "Atlas-wide CS-R and CS-C percentile ranks assigned to 10-point intervals [0,10), ..., [90,100]."},
        {"Item": "Randomization", "Value": "Size-matched samples drawn without replacement from the reference; each sample compared with the reference remaining after its removal."},
        {"Item": "Cell enrichment", "Value": "log2[(observed+0.5)/(expected+0.5)]; upper-tail empirical P with +1 correction only for observed log2 enrichment>0."},
        {"Item": "Multiple testing", "Value": "Benjamini-Hochberg FDR across 100 cells separately for each evidence set and model."},
        {"Item": "TVD", "Value": "0.5 times the sum of absolute differences between unsmoothed 10x10 proportions; upper-tail empirical P with +1 correction."},
        {"Item": "Resamples", "Value": args.n_resamples},
    ]
    write_tsv(args.output / "Methods.tsv", methods, ["Item", "Value"])
    manifest = {
        "input": str(args.input.resolve()), "n_resamples": args.n_resamples,
        "seed_base": RANDOM_SEED_BASE,
        "software": {"python": sys.version, "platform": platform.platform(),
                     "numpy": np.__version__, "pandas": pd.__version__},
        "counts": {x["Check"]: x["Observed"] for x in checks},
    }
    (args.output / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
