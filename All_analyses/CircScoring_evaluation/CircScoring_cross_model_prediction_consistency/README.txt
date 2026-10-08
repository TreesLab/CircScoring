# CircScoring cross-model prediction consistency

This compact submission bundle reproduces the **cross-model prediction
consistency** analysis comparing the finalized XGBoost and Elastic Net CS-R and
CS-C scores.

## Historical fidelity

The historical execution script was not retained. The included script is a
faithful reconstruction from the preserved fixed-score table and
`Cross_model_consistency_summary.xlsx`. It was verified against all 5,910
retained BSJ-level decile assignments and all six reported Spearman
correlations, P values, and quadratic-weighted kappa estimates. It should not be
described as a byte-for-byte copy of the lost historical script.

## Analysis

Consistency was evaluated separately in three benchmark pairs:

- P1N1: 1,301 BSJs
- P2N2: 1,289 BSJs
- P2′N2′: 365 BSJs

For CS-R and CS-C separately, the script:

1. calculates Spearman's rank correlation between the continuous XGBoost and
   Elastic Net scores, with a two-sided asymptotic P value (`exact = FALSE`);
2. orders BSJs by ascending score, resolving tied scores by ascending
   `circRNA_id`;
3. assigns the BSJ at ordered position `r` among `n` BSJs to decile
   `ceil(10 × r / n)`;
4. calculates quadratic-weighted Cohen's kappa between the XGBoost and Elastic
   Net deciles; and
5. exports the corresponding 10 × 10 count and overall-percentage matrix.

No model refitting, bootstrap procedure, or permutation test is used in this
analysis.

## Included files

- `cross_model_prediction_consistency.R`: complete numerical analysis.
- `cross_model_fixed_scores.tsv`: fixed model scores and benchmark labels used
  as input.
- `README.md`: methods, software environment, requirements, and run command.

The script produces three machine-readable outputs: the six-row summary,
BSJ-level decile assignments, and the six 10 × 10 matrices. Plotting code is
intentionally omitted because the final figure is already available.

## Software environment and requirements

The reconstruction was verified with:

- R 4.6.1
- base R packages only (`stats` and `utils`, distributed with R)

No external R packages are required.

## Run

```bash
Rscript cross_model_prediction_consistency.R \
  cross_model_fixed_scores.tsv \
  results
```

## Manuscript-ready methods text

Cross-model prediction consistency was assessed separately for CS-R and CS-C
within the P1N1, P2N2, and P2′N2′ benchmark sets. Spearman's rank correlation
was calculated between the continuous XGBoost and Elastic Net scores, with
two-sided asymptotic P values. For categorical agreement assessment, BSJs were
sorted in ascending order of score for each model; tied scores were resolved
deterministically by ascending BSJ identifier. Equal-count deciles were assigned
according to the ordered position of each BSJ, and agreement between XGBoost and
Elastic Net deciles was quantified using quadratic-weighted Cohen's kappa.
