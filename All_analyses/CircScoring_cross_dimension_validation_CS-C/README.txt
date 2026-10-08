# Cross-dimension validation of CS-C

This compact source bundle reproduces the numerical analysis underlying the
cross-dimension validation of CS-C. It intentionally contains no plotting code
because the final figure has already been prepared.

## Historical fidelity

The historical execution script was not retained. The included R script is a
faithful reconstruction from the preserved analysis outputs, event counts,
figure, and fixed circIG field definitions. The analysis specification and all
retained numerical controls are preserved; the script should not be described
as a byte-for-byte copy of the lost historical file.

## Analysis population and required input

The analysis uses the complete analysis-eligible circIG population
(`N = 4,455,197` BSJs). The input must be a tab-delimited text file, optionally
gzip-compressed, containing these fields:

- `XGBoost_CS-C_score`
- `ElasticNet_CS-C_score`
- `FL-circAS or circFL_seq`
- `mouse conserved` or `mouse_conserved`
- `#7_DB` or `#7_DBs`

The reader accepts either a conventional one-row header or the two-row circIG
export header in which the first row contains column groups. Blank evidence
cells are treated as absence of evidence. CS-C predictor values must be present
for every analyzed BSJ.

The script checks the fixed atlas size and the following retained event counts
before fitting models:

| Outcome | Event count |
|---|---:|
| Full-length support | 800,687 |
| Mouse conservation | 94,840 |
| Database recurrence ≥1 | 1,291,639 |
| Database recurrence ≥2 | 348,111 |
| Database recurrence ≥3 | 214,436 |
| Database recurrence ≥4 | 146,985 |
| Database recurrence ≥5 | 101,435 |
| Database recurrence ≥6 | 45,429 |
| Database recurrence ≥7 | 10,287 |

These checks prevent silent use of a different circIG release. The optional
`--relaxed-qc` argument is intended only for testing the code on another input
and should not be used to reproduce the reported analysis.

The retained atlas-wide means and sample standard deviations of both CS-C
predictors are also checked internally before model fitting.

## Statistical procedure

The XGBoost and Elastic Net CS-C scores are standardized independently across
all 4,455,197 BSJs using the sample standard deviation (denominator `N - 1`).
Separate univariable logistic regression models are then fitted for each model
and each binary outcome:

1. full-length circRNA support;
2. mouse conservation; and
3. database recurrence at cumulative thresholds `#7_DB ≥ 1` through
   `#7_DB ≥ 7`.

For each regression,

`logit[Pr(evidence = 1)] = β0 + β1 × standardized CS-C`.

The reported odds ratio is `exp(β1)` and therefore represents the change in
odds per 1-SD increase in CS-C. The 95% confidence interval is the model-based
Wald interval `exp[β1 ± 1.96 × SE(β1)]`; the P value is the two-sided Wald test
of `H0: β1 = 0`. Models are fitted separately, without covariate adjustment,
resampling, or multiple-testing correction.

## Software environment and requirements

The reconstruction was prepared and syntax-checked with:

- R 4.6.1
- base R packages only: `stats` and `utils` (distributed with R)

No external R packages are required. Memory sufficient to hold the selected
atlas-wide columns is required; approximately 2–4 GB of available RAM is
recommended for the fixed input.

## Run

```bash
Rscript cross_dimension_validation_cs_c.R \
  cross_dimension_validation_input.tsv.gz \
  results
```

The script writes two tab-delimited files:

- `Cross_dimension_validation_CS-C_logistic_regression.tsv`
- `Cross_dimension_validation_CS-C_event_counts.tsv`

## Manuscript-ready methods text

Cross-dimension validation of CS-C was performed using orthogonal measures of
circRNA reliability, comprising full-length circRNA support, mouse
conservation, and recurrence across seven circRNA databases. XGBoost and
Elastic Net CS-C scores were standardized independently across all
analysis-eligible circIG BSJs using the sample standard deviation. Separate
univariable logistic regression models were fitted for each evidence outcome
and model. Database recurrence was evaluated at cumulative thresholds ranging
from at least one to at least seven supporting databases. Odds ratios and
model-based Wald 95% confidence intervals were reported per 1-SD increase in
CS-C, together with two-sided Wald P values.
