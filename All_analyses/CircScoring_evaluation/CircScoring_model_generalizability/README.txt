# CircScoring model-generalizability analysis

This directory contains the submission-ready source code for the CircScoring
**model generalizability** analysis.

## Historical fidelity

The historical source script itself was not retained. This is a faithful
reconstruction from the preserved fixed-score table, bootstrap result CSV,
figure assets, sample counts, iteration count, and random seed. The numerical
implementation was verified against all retained estimates and confidence
limits. It should not be described as a byte-for-byte copy of the lost script.

## Exact analysis reconstructed from the retained outputs

- The four finalized scores were evaluated without model refitting.
- Generalizability was examined in P1N1 (797/504), P2N2 (805/484), and
  P2′N2′ (203/162) positive/negative benchmark pairs.
- AUROC used the Mann-Whitney formulation with average ranks for ties.
- The result column labeled `AUPRC` used non-interpolated average precision,
  `AP = sum((R_n - R_(n-1)) * P_n)`, evaluated at distinct score thresholds.
- A single NumPy generator was initialized with seed `20260803`.
- For each metric, 1,000 bootstrap samples were drawn separately within the
  positive and negative classes, preserving the original class sizes.
- The 95% confidence interval was defined by the 2.5th and 97.5th percentiles
  of the bootstrap distribution.
- For exact historical reproducibility, the benchmark/score/model/metric loop
  order in the source code must not be changed because all analyses consume one
  sequential random-number stream.

## Files

- `01_model_generalizability_bootstrap.py`: numerical analysis and exact
  regression test against the retained result table.
- `02_plot_model_generalizability.R`: transparent 12 × 7.59 cm SVG/PNG figure.
- `data/generalizability_fixed_scores.tsv`: labels and the four fixed score
  vectors used in the analysis.
- `expected/CircScoring_model_generalizability_bootstrap_results_retained.csv`:
  retained reference results.
- `expected/figure_retained.svg` and `expected/figure_retained.png`: retained
  final figure assets for visual comparison with regenerated output.
- `METHODS_text.md`: concise manuscript-ready description.
- `requirements.txt` and `SOFTWARE_ENVIRONMENT.tsv`: verified software setup.

## Run

```bash
python 01_model_generalizability_bootstrap.py \
  --input data/generalizability_fixed_scores.tsv \
  --expected expected/CircScoring_model_generalizability_bootstrap_results_retained.csv \
  --output results/CircScoring_model_generalizability_bootstrap_results.csv

Rscript 02_plot_model_generalizability.R \
  results/CircScoring_model_generalizability_bootstrap_results.csv \
  results/figures
```

The numerical script stops with an error if any categorical field, sample
count, estimate, or confidence limit differs from the retained output beyond
floating-point precision.

## Requirements

Component	Version	Role
Python	3.12.14	Bootstrap analysis and CSV export
NumPy	2.3.5	Deterministic random-number generation and numerical operations
pandas	3.0.1	Input/output and average-rank calculation
R	4.6.1	Transparent SVG/PNG figure generation using base R

numpy==2.3.5
pandas==3.0.1
