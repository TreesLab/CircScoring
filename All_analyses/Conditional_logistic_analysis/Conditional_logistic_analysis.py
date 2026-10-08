#!/usr/bin/env python3
"""Conditional logistic-regression analysis for CircScoring.

For each of 11 evidence sets and each model (XGBoost and Elastic Net), this
script combines evidence-positive BSJs with the fixed circIG
reference, fits CS-R-only, CS-C-only, and mutually adjusted logistic models,
checks linearity against restricted cubic splines (RCS), and quantifies the
incremental contribution of each score dimension.

The script accepts both the original and the later prefixed column names used
in the circIG TSV.  The TSV has two header rows; the row beginning with
"chrom, start, end, strand" is detected automatically.

Run:
  python Sec6_conditional_logistic_analysis.py --input circIG.tsv --output results

Use --validate-only to reproduce population counts and spline knots without
fitting the logistic models.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


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

EXPECTED_COUNTS = {
    "N_total": 4_455_197,
    "N_reference": 4_142_933,
    "RT-independent": 2_493,
    "Both RT/nonRT": 599,
    "Translation": 291_730,
    "m6A2Circ >=1": 42_019,
    "circTarget >=1": 2_215,
    "Aortic": 142,
    "Lung": 3_914,
    "Skin": 702,
    "Umbilical": 130,
    "White blood cell": 411,
    "Disease": 3_825,
}

MODEL_SCORES = {
    "XGBoost": ("XGB_R", "XGB_C"),
    "Elastic Net": ("EN_R", "EN_C"),
}


def detect_header_row(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for i in range(10):
            line = handle.readline()
            if not line:
                break
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 4 and fields[:4] == ["chrom", "start", "end", "strand"]:
                return i
    raise ValueError("Could not locate the TSV header row beginning with chrom/start/end/strand.")


def resolve_columns(path: Path, header_row: int) -> Dict[str, str]:
    columns = list(pd.read_csv(path, sep="\t", header=header_row, nrows=0).columns)
    resolved: Dict[str, str] = {}
    for logical, alternatives in ALIASES.items():
        match = next((name for name in alternatives if name in columns), None)
        if match is None:
            raise KeyError(f"Required column for {logical!r} not found. Tried: {alternatives}")
        resolved[logical] = match
    return resolved


def load_arrays(path: Path, chunksize: int) -> Dict[str, np.ndarray]:
    header_row = detect_header_row(path)
    resolved = resolve_columns(path, header_row)
    reverse = {actual: logical for logical, actual in resolved.items()}
    usecols = list(reverse)
    pieces: Dict[str, List[np.ndarray]] = {logical: [] for logical in resolved}

    for chunk in pd.read_csv(
        path, sep="\t", header=header_row, usecols=usecols,
        chunksize=chunksize, low_memory=False,
    ):
        for actual, logical in reverse.items():
            values = pd.to_numeric(chunk[actual], errors="coerce").to_numpy(np.float64)
            pieces[logical].append(values)

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
        "Aortic": eq1(data["Aortic"]),
        "Lung": eq1(data["Lung"]),
        "Skin": eq1(data["Skin"]),
        "Umbilical": eq1(data["Umbilical"]),
        "White blood cell": eq1(data["White blood cell"]),
        "Disease": eq1(data["Disease"]),
    }
    positives = {name: benchmark_free & raw[name] for name in EVIDENCE_ORDER}

    any_evidence = np.zeros_like(benchmark_overlap)
    for mask in raw.values():
        any_evidence |= mask
    reference = benchmark_free & ~any_evidence
    return reference, positives, benchmark_overlap


def choose_knots(x_reference: np.ndarray) -> Tuple[np.ndarray, str]:
    x = x_reference[np.isfinite(x_reference)]
    four = np.quantile(x, [0.05, 0.35, 0.65, 0.95], method="linear")
    if len(np.unique(four)) == 4 and np.all(np.diff(four) > 0):
        return four, "reference percentiles (5,35,65,95)"
    distinct = np.unique(x)
    three = np.quantile(distinct, [0.10, 0.50, 0.90], method="linear")
    if len(np.unique(three)) != 3 or not np.all(np.diff(three) > 0):
        raise ValueError("Unable to construct distinct RCS knots.")
    return three, "three-knot distinct-value fallback (10,50,90)"


def rcs_basis(x: np.ndarray, knots_original: np.ndarray) -> np.ndarray:
    """Harrell restricted cubic spline basis: linear term + K-2 nonlinear terms."""
    x10 = np.asarray(x, dtype=np.float64) / 10.0
    knots = np.asarray(knots_original, dtype=np.float64) / 10.0
    k_last2, k_last = knots[-2], knots[-1]
    scale = (knots[-1] - knots[0]) ** 2
    columns = [x10]
    for kj in knots[:-2]:
        term = np.maximum(x10 - kj, 0.0) ** 3
        term -= np.maximum(x10 - k_last2, 0.0) ** 3 * (k_last - kj) / (k_last - k_last2)
        term += np.maximum(x10 - k_last, 0.0) ** 3 * (k_last2 - kj) / (k_last - k_last2)
        columns.append(term / scale)
    return np.column_stack(columns)


def expit_stable(eta: np.ndarray) -> np.ndarray:
    out = np.empty_like(eta)
    pos = eta >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-eta[pos]))
    e = np.exp(eta[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def log_likelihood(y: np.ndarray, eta: np.ndarray) -> float:
    return float(np.sum(y * eta - np.logaddexp(0.0, eta)))


@dataclass
class Fit:
    structure: str
    names: List[str]
    beta: np.ndarray
    covariance: np.ndarray
    ll: float
    aic: float
    mcfadden_r2: float
    converged: bool
    iterations: int
    last_max_abs_step: float

    @property
    def df_model(self) -> int:
        return len(self.beta) - 1


def make_design(r: np.ndarray, c: np.ndarray, structure: str,
                knots_r: np.ndarray, knots_c: np.ndarray) -> Tuple[np.ndarray, List[str]]:
    columns = [np.ones(len(r), dtype=np.float64)]
    names = ["(Intercept)"]
    if "linR" in structure:
        columns.append(r / 10.0); names.append("R_lin")
    elif "splR" in structure:
        br = rcs_basis(r, knots_r)
        for j in range(br.shape[1]):
            columns.append(br[:, j]); names.append(f"R_rcs{j + 1}")
    if "linC" in structure:
        columns.append(c / 10.0); names.append("C_lin")
    elif "splC" in structure:
        bc = rcs_basis(c, knots_c)
        for j in range(bc.shape[1]):
            columns.append(bc[:, j]); names.append(f"C_rcs{j + 1}")
    return np.column_stack(columns), names


def fit_logistic(y: np.ndarray, design: np.ndarray, names: List[str], structure: str,
                 max_iter: int = 100, tolerance: float = 1e-8,
                 max_step: float = 8.0, initial: np.ndarray | None = None) -> Fit:
    beta = np.zeros(design.shape[1], dtype=np.float64) if initial is None else np.asarray(initial, dtype=np.float64).copy()
    prevalence = min(max(float(y.mean()), 1e-12), 1.0 - 1e-12)
    if initial is None:
        beta[0] = math.log(prevalence / (1.0 - prevalence))
    converged = False
    last_step = math.inf
    iterations = 0

    for iterations in range(1, max_iter + 1):
        eta = design @ beta
        prob = expit_stable(eta)
        weight = np.maximum(prob * (1.0 - prob), 1e-12)
        gradient = design.T @ (y - prob)
        information = design.T @ (design * weight[:, None])
        try:
            step = np.linalg.solve(information, gradient)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(information) @ gradient
        step = np.clip(step, -max_step, max_step)
        old_ll = log_likelihood(y, eta)
        factor = 1.0
        while factor >= 2.0 ** -20:
            candidate = beta + factor * step
            if log_likelihood(y, design @ candidate) >= old_ll - 1e-10:
                break
            factor *= 0.5
        step *= factor
        beta += step
        last_step = float(np.max(np.abs(step)))
        if last_step < tolerance:
            converged = True
            break

    eta = design @ beta
    ll = log_likelihood(y, eta)
    prob = expit_stable(eta)
    weight = np.maximum(prob * (1.0 - prob), 1e-12)
    information = design.T @ (design * weight[:, None])
    covariance = np.linalg.pinv(information)
    null_eta = np.full(len(y), math.log(prevalence / (1.0 - prevalence)))
    null_ll = log_likelihood(y, null_eta)
    return Fit(
        structure=structure, names=names, beta=beta, covariance=covariance,
        ll=ll, aic=-2.0 * ll + 2.0 * len(beta),
        mcfadden_r2=1.0 - ll / null_ll, converged=converged,
        iterations=iterations, last_max_abs_step=last_step,
    )


def starting_values(target_names: Sequence[str], base_fit: Fit | None) -> np.ndarray | None:
    if base_fit is None:
        return None
    initial = np.zeros(len(target_names), dtype=np.float64)
    base = {name: value for name, value in zip(base_fit.names, base_fit.beta)}
    for i, name in enumerate(target_names):
        if name in base:
            initial[i] = base[name]
        elif name == "R_rcs1" and "R_lin" in base:
            initial[i] = base["R_lin"]
        elif name == "C_rcs1" and "C_lin" in base:
            initial[i] = base["C_lin"]
    return initial


def chi_square_sf(x: float, df: int) -> float:
    if not np.isfinite(x) or x < 0:
        return math.nan
    z = math.sqrt(x / 2.0)
    if df == 1:
        return math.erfc(z)
    if df == 2:
        return math.exp(-x / 2.0)
    if df == 3:
        return math.erfc(z) + math.sqrt(2.0 * x / math.pi) * math.exp(-x / 2.0)
    raise ValueError(f"chi_square_sf currently supports df 1-3, received {df}.")


def lr_compare(reduced: Fit, full: Fit) -> Dict[str, float]:
    lr = max(0.0, 2.0 * (full.ll - reduced.ll))
    df = full.df_model - reduced.df_model
    return {
        "LR_chisq": lr, "df": df, "P": chi_square_sf(lr, df),
        "Delta_AIC": full.aic - reduced.aic,
        "Delta_R2": full.mcfadden_r2 - reduced.mcfadden_r2,
    }


def substantial(test: Mapping[str, float], spline_fit: Fit) -> bool:
    return bool(
        spline_fit.converged and test["P"] < 0.001
        and test["Delta_AIC"] <= -10.0 and test["Delta_R2"] >= 0.001
    )


def normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def safe_exp(x: float) -> float:
    if x > 709:
        return math.inf
    if x < -745:
        return 0.0
    return math.exp(x)


def coefficient_row(fit: Fit, term: str) -> Tuple[float, float, float, float, float, float]:
    i = fit.names.index(term)
    beta = float(fit.beta[i])
    se = math.sqrt(max(float(fit.covariance[i, i]), 0.0))
    z = beta / se if se > 0 else math.inf
    return safe_exp(beta), safe_exp(beta - 1.96 * se), safe_exp(beta + 1.96 * se), normal_two_sided_p(z), beta, se


def contrast_rows(fit: Fit, predictor: str, form: str,
                  score_values: np.ndarray, knots: np.ndarray) -> List[Dict[str, object]]:
    if form == "Linear":
        term = f"{predictor}_lin"
        or_, low, high, _, _, _ = coefficient_row(fit, term)
        return [{"10_rank_contrast": "constant per 10-rank increase", "OR": or_,
                 "CI95_low": low, "CI95_high": high}]

    finite = score_values[np.isfinite(score_values)]
    start = int(math.ceil(float(finite.min()) / 10.0) * 10)
    end = int(math.floor(float(finite.max()) / 10.0) * 10)
    indices = [i for i, name in enumerate(fit.names) if name.startswith(f"{predictor}_rcs")]
    rows = []
    for a in range(start, end, 10):
        b = a + 10
        if b > end:
            continue
        basis = rcs_basis(np.asarray([a, b], dtype=float), knots)
        delta_basis = basis[1] - basis[0]
        delta = np.zeros(len(fit.beta), dtype=float)
        delta[indices] = delta_basis
        log_or = float(delta @ fit.beta)
        variance = max(float(delta @ fit.covariance @ delta), 0.0)
        se = math.sqrt(variance)
        rows.append({"10_rank_contrast": f"{a}-{b}", "OR": safe_exp(log_or),
                     "CI95_low": safe_exp(log_or - 1.96 * se),
                     "CI95_high": safe_exp(log_or + 1.96 * se)})
    return rows


def write_tsv(path: Path, rows: List[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_analysis(data: Mapping[str, np.ndarray], reference: np.ndarray,
                 positives: Mapping[str, np.ndarray], output: Path,
                 evidence_filter: Sequence[str] | None = None,
                 model_filter: Sequence[str] | None = None,
                 model_specification: Mapping[str, object] | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    selected_models: List[Dict[str, object]] = []
    nonlinear_rows: List[Dict[str, object]] = []
    incremental_rows: List[Dict[str, object]] = []
    linear_rows: List[Dict[str, object]] = []
    contrast_output: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    saved_json: Dict[str, object] = {}

    evidence_names = [e for e in EVIDENCE_ORDER if not evidence_filter or e in evidence_filter]
    model_names = [m for m in MODEL_SCORES if not model_filter or m in model_filter]

    for model_name in model_names:
        r_name, c_name = MODEL_SCORES[model_name]
        knots_r, method_r = choose_knots(data[r_name][reference])
        knots_c, method_c = choose_knots(data[c_name][reference])
        for evidence in evidence_names:
            positive = positives[evidence]
            analysis_mask = reference | positive
            complete = analysis_mask & np.isfinite(data[r_name]) & np.isfinite(data[c_name])
            y = positive[complete].astype(np.float64)
            r = data[r_name][complete]
            c = data[c_name][complete]
            n_positive = int(y.sum())
            n_reference = int(len(y) - y.sum())

            structures = ["linR", "linC", "linR_linC", "splR", "splC",
                          "splR_linC", "linR_splC", "splR_splC"]
            fits: Dict[str, Fit] = {}
            for structure in structures:
                design, names = make_design(r, c, structure, knots_r, knots_c)
                if structure == "splR":
                    base_fit = fits["linR"]
                elif structure == "splC":
                    base_fit = fits["linC"]
                elif structure in ("splR_linC", "linR_splC", "splR_splC"):
                    base_fit = fits["linR_linC"]
                else:
                    base_fit = None
                fits[structure] = fit_logistic(
                    y, design, names, structure,
                    initial=starting_values(names, base_fit),
                )
                del design

            comparisons = {
                ("Univariable", "CS-R"): (fits["linR"], fits["splR"]),
                ("Univariable", "CS-C"): (fits["linC"], fits["splC"]),
                ("Mutually adjusted", "CS-R"): (fits["linR_linC"], fits["splR_linC"]),
                ("Mutually adjusted", "CS-C"): (fits["linR_linC"], fits["linR_splC"]),
            }
            decisions: Dict[Tuple[str, str], bool] = {}
            for (context, predictor), (linear_fit, spline_fit) in comparisons.items():
                test = lr_compare(linear_fit, spline_fit)
                decision = substantial(test, spline_fit)
                decisions[(context, predictor)] = decision
                nonlinear_rows.append({
                    "Evidence": evidence, "Model": model_name, "Context": context,
                    "Predictor": predictor, "LR_chisq_nonlinearity": test["LR_chisq"] if spline_fit.converged else "",
                    "df": test["df"], "P_nonlinearity": test["P"] if spline_fit.converged else "",
                    "Delta_AIC_spline_minus_linear": test["Delta_AIC"] if spline_fit.converged else "",
                    "Delta_McFadden_R2": test["Delta_R2"] if spline_fit.converged else "",
                    "Substantial_nonlinearity": decision,
                    "Linear_converged": linear_fit.converged, "Spline_converged": spline_fit.converged,
                    "Spline_last_max_abs_step": spline_fit.last_max_abs_step,
                })

            specification_key = f"{evidence}|{model_name}"
            if model_specification and specification_key in model_specification:
                locked = model_specification[specification_key].get("decisions", {})
                for context, predictor in list(decisions):
                    saved_key = f"{context}|{predictor}"
                    if saved_key in locked:
                        decisions[(context, predictor)] = bool(locked[saved_key])
                for row in nonlinear_rows[-4:]:
                    row["Substantial_nonlinearity"] = decisions[(row["Context"], row["Predictor"])]

            form_r = "RCS" if decisions[("Mutually adjusted", "CS-R")] else "Linear"
            form_c = "RCS" if decisions[("Mutually adjusted", "CS-C")] else "Linear"
            full_structure = ("splR" if form_r == "RCS" else "linR") + "_" + ("splC" if form_c == "RCS" else "linC")
            full = fits[full_structure]
            reduced_r = fits["splR" if form_r == "RCS" else "linR"]
            reduced_c = fits["splC" if form_c == "RCS" else "linC"]
            add_c = lr_compare(reduced_r, full)
            add_r = lr_compare(reduced_c, full)

            for added, reduced, test in [
                ("CS-C added to CS-R", reduced_r, add_c),
                ("CS-R added to CS-C", reduced_c, add_r),
            ]:
                incremental_rows.append({
                    "Evidence": evidence, "Model": model_name, "Added_dimension": added,
                    "Reduced_model": reduced.structure, "Full_model": full.structure,
                    "LR_chisq": test["LR_chisq"], "df": test["df"], "LR_P": test["P"],
                    "Delta_AIC_full_minus_reduced": test["Delta_AIC"],
                    "Incremental_McFadden_R2": test["Delta_R2"],
                })

            selected_models.append({
                "Evidence": evidence, "Model": model_name, "N_positive": n_positive,
                "N_reference": n_reference, "CS-R_form": form_r, "CS-C_form": form_c,
                "Selected_structure": full.structure, "Log_likelihood": full.ll,
                "AIC": full.aic, "McFadden_R2": full.mcfadden_r2,
                "Converged": full.converged, "Last_max_abs_step": full.last_max_abs_step,
            })
            summary_rows.append({
                "Evidence": evidence, "Model": model_name, "N_positive": n_positive,
                "N_reference": n_reference, "CS-R form": form_r, "CS-C form": form_c,
                "McFadden_R2": full.mcfadden_r2,
                "P: add CS-C to CS-R": add_c["P"], "Delta R2: add CS-C": add_c["Delta_R2"],
                "P: add CS-R to CS-C": add_r["P"], "Delta R2: add CS-R": add_r["Delta_R2"],
                "Converged": full.converged,
            })

            # Always report conventional linear ORs.
            for analysis, fit_key in [("Univariable", "linR"), ("Univariable", "linC"),
                                      ("Mutually adjusted", "linR_linC")]:
                fit = fits[fit_key]
                predictors = ["CS-R"] if fit_key == "linR" else (["CS-C"] if fit_key == "linC" else ["CS-R", "CS-C"])
                for predictor in predictors:
                    term = "R_lin" if predictor == "CS-R" else "C_lin"
                    or_, low, high, p, beta, se = coefficient_row(fit, term)
                    linear_rows.append({
                        "Evidence": evidence, "Model": model_name, "N_positive": n_positive,
                        "N_reference": n_reference, "Analysis": analysis, "Predictor": predictor,
                        "OR_per_10_rank": or_, "CI95_low": low, "CI95_high": high,
                        "Wald_two_sided_P": p, "Beta_per_10_rank": beta, "SE": se,
                        "Log_likelihood": fit.ll, "AIC": fit.aic, "Converged": fit.converged,
                        "Iterations": fit.iterations, "Last_max_abs_step": fit.last_max_abs_step,
                    })

            # Report selected-form 10-rank contrasts.
            selected_contexts = [
                ("Univariable", "CS-R", "splR" if decisions[("Univariable", "CS-R")] else "linR"),
                ("Univariable", "CS-C", "splC" if decisions[("Univariable", "CS-C")] else "linC"),
                ("Mutually adjusted", "CS-R", full_structure),
                ("Mutually adjusted", "CS-C", full_structure),
            ]
            for analysis, predictor, fit_key in selected_contexts:
                fit = fits[fit_key]
                form = (form_r if predictor == "CS-R" else form_c) if analysis == "Mutually adjusted" else ("RCS" if decisions[(analysis, predictor)] else "Linear")
                values = r if predictor == "CS-R" else c
                knots = knots_r if predictor == "CS-R" else knots_c
                for row in contrast_rows(fit, "R" if predictor == "CS-R" else "C", form, values, knots):
                    contrast_output.append({
                        "Evidence": evidence, "Model": model_name, "Analysis": analysis,
                        "Predictor": predictor, "Functional_form": form, **row,
                    })

            saved_json[f"{evidence}|{model_name}"] = {
                "knots_r": knots_r.tolist(), "knots_c": knots_c.tolist(),
                "knot_method_r": method_r, "knot_method_c": method_c,
                "decisions": {f"{ctx}|{pred}": value for (ctx, pred), value in decisions.items()},
                "selected": full_structure,
            }
            del r, c, y, fits

    write_tsv(output / "Summary.tsv", summary_rows, [
        "Evidence", "Model", "N_positive", "N_reference", "CS-R form", "CS-C form",
        "McFadden_R2", "P: add CS-C to CS-R", "Delta R2: add CS-C",
        "P: add CS-R to CS-C", "Delta R2: add CS-R", "Converged",
    ])
    write_tsv(output / "Linear_ORs.tsv", linear_rows, [
        "Evidence", "Model", "N_positive", "N_reference", "Analysis", "Predictor",
        "OR_per_10_rank", "CI95_low", "CI95_high", "Wald_two_sided_P",
        "Beta_per_10_rank", "SE", "Log_likelihood", "AIC", "Converged",
        "Iterations", "Last_max_abs_step",
    ])
    write_tsv(output / "Nonlinearity.tsv", nonlinear_rows, [
        "Evidence", "Model", "Context", "Predictor", "LR_chisq_nonlinearity", "df",
        "P_nonlinearity", "Delta_AIC_spline_minus_linear", "Delta_McFadden_R2",
        "Substantial_nonlinearity", "Linear_converged", "Spline_converged",
        "Spline_last_max_abs_step",
    ])
    write_tsv(output / "Selected_models.tsv", selected_models, [
        "Evidence", "Model", "N_positive", "N_reference", "CS-R_form", "CS-C_form",
        "Selected_structure", "Log_likelihood", "AIC", "McFadden_R2", "Converged",
        "Last_max_abs_step",
    ])
    write_tsv(output / "Selected_form_10rank_ORs.tsv", contrast_output, [
        "Evidence", "Model", "Analysis", "Predictor", "Functional_form",
        "10_rank_contrast", "OR", "CI95_low", "CI95_high",
    ])
    write_tsv(output / "Incremental_tests.tsv", incremental_rows, [
        "Evidence", "Model", "Added_dimension", "Reduced_model", "Full_model",
        "LR_chisq", "df", "LR_P", "Delta_AIC_full_minus_reduced",
        "Incremental_McFadden_R2",
    ])
    (output / "conditional_results.json").write_text(
        json.dumps(saved_json, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def write_validation(output: Path, data: Mapping[str, np.ndarray], reference: np.ndarray,
                     positives: Mapping[str, np.ndarray], benchmark_overlap: np.ndarray) -> None:
    output.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, object]] = [
        {"Item": "N_total", "Value": len(reference)},
        {"Item": "N_benchmark_overlap", "Value": int(benchmark_overlap.sum())},
        {"Item": "N_Sec6_reference", "Value": int(reference.sum())},
        {"Item": "New_translation_definition", "Value": "riboCIRC=1 OR TransCirc evidence count>=1"},
    ]
    for model, (r_name, c_name) in MODEL_SCORES.items():
        for dimension, score_name in [("CS-R", r_name), ("CS-C", c_name)]:
            x = data[score_name][reference]
            knots, method = choose_knots(x)
            rows.extend([
                {"Item": f"{model}_{dimension}_unique_values", "Value": len(np.unique(x[np.isfinite(x)]))},
                {"Item": f"{model}_{dimension}_knots", "Value": ",".join(f"{v:.12g}" for v in knots)},
                {"Item": f"{model}_{dimension}_knot_method", "Value": method},
            ])
    for evidence in EVIDENCE_ORDER:
        rows.append({"Item": f"N_positive_{evidence}", "Value": int(positives[evidence].sum())})
    write_tsv(output / "QC.tsv", rows, ["Item", "Value"])

    observed = {"N_total": len(reference), "N_reference": int(reference.sum()),
                **{name: int(mask.sum()) for name, mask in positives.items()}}
    checks = []
    for key, expected in EXPECTED_COUNTS.items():
        actual = observed[key]
        checks.append({"Check": key, "Expected": expected, "Observed": actual, "Pass": actual == expected})
    write_tsv(output / "Validation_checks.tsv", checks, ["Check", "Expected", "Observed", "Pass"])
    if not all(row["Pass"] for row in checks):
        failed = [row for row in checks if not row["Pass"]]
        raise RuntimeError(f"Population validation failed: {failed}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Two-header-row circIG TSV")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument("--validate-only", action="store_true", help="Only reproduce counts and knots")
    parser.add_argument("--evidence", action="append", choices=EVIDENCE_ORDER,
                        help="Optional evidence subset; may be supplied repeatedly")
    parser.add_argument("--model", action="append", choices=list(MODEL_SCORES),
                        help="Optional model subset; may be supplied repeatedly")
    parser.add_argument(
        "--model-specification", type=Path,
        default=Path(__file__).with_name("model_specification.json"),
        help="Frozen diagnostic decisions used for the published analysis",
    )
    parser.add_argument(
        "--reselect-models", action="store_true",
        help="Ignore the frozen specification and select forms from the current run",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_arrays(args.input, args.chunksize)
    reference, positives, benchmark_overlap = define_populations(data)
    write_validation(args.output, data, reference, positives, benchmark_overlap)
    metadata = {
        "input": str(args.input.resolve()), "python": sys.version,
        "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__,
        "reference_definition": "Benchmark-free and negative for all 11 evidence definitions",
        "positive_definition": "Evidence-specific positive after exclusion of all six benchmark sets",
        "nonlinearity_thresholds": {"LR_P": 0.001, "Delta_AIC_max": -10.0, "Delta_McFadden_R2_min": 0.001},
    }
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not args.validate_only:
        model_specification = None
        if not args.reselect_models and args.model_specification.exists():
            model_specification = json.loads(args.model_specification.read_text(encoding="utf-8"))
        run_analysis(
            data, reference, positives, args.output, args.evidence, args.model,
            model_specification=model_specification,
        )


if __name__ == "__main__":
    main()
