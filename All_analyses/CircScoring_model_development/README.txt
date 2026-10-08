# FEATURE DEFINITIONS

feature_id	feature_name	CircScoring_dimension
f1	donor_site_at_the_annotated_boundary	CS-C
f2	acceptor_site_at_the_annotated_boundary	CS-C
f3	donor_acceptor_sites_at_the_same_transcript_isoform	CS-C
f4	has_AS_event(donor)	CS-C
f5	has_AS_event(acceptor)	CS-C
f6	phyloP(acceptor_in)	CS-C
f7	phyloP(acceptor_out)	CS-C
f8	phyloP(donor_in)	CS-C
f9	phyloP(donor_out)	CS-C
f10	phastCons(acceptor_in)	CS-C
f11	phastCons(acceptor_out)	CS-C
f12	phastCons(donor_in)	CS-C
f13	phastCons(donor_out)	CS-C
f14	MAXENT(donor)	CS-C
f15	MAXENT(acceptor)	CS-C
f16	7 DB	CS-R
f17	FL-circAS or circFL_seq	CS-R
f18	mouse_conserved	CS-R

# CircScoring model-development source code

This directory reproduces the development of the four finalized CircScoring
models used by `predictCS`:

1. XGBoost CS-R, using f16-f18;
2. XGBoost CS-C, using f1-f15;
3. Elastic Net CS-R, using f16-f18; and
4. Elastic Net CS-C, using f1-f15.

The authoritative model specification is the one recorded in
`tunning_P1N1 1.xlsx` and embedded in the released `predictCS.R` and
`predictCS.py` programs. The alternative archive named
`P1N1_retraining_corrected_AUROC_AUPRC1.zip` used a different candidate grid
and a different model-selection rule; it is not the source of the deployed
models and is therefore not used here.

## Analysis population

The development set contains 1,301 unique BSJs from the `6cols` worksheet of
`Benchmark_features_lite.xlsx`:

- 797 P1-positive BSJs;
- 504 N1-negative BSJs;
- no P1/N1 overlap; and
- no missing values in f1-f18.

`data/P1N1_training_data.tsv` is a frozen, machine-readable copy of the exact
analysis population. `prepare_training_data.py` regenerates it from the source
workbook and verifies its class counts, identifiers and feature completeness.

## Reproducible workflow

Use R 4.4.1 with `xgboost` 3.2.1.1 and `glmnet` 4.1-10.

```bash
Rscript install_exact_packages.R
Rscript CircScoring_model_development.R \
  --input data/P1N1_training_data.tsv \
  --output results \
  --threads 1
```

On Windows, the same command can be launched with:

```text
run_model_development.bat
```

`run_refit_only.bat` regenerates the four final models from the archived
selected hyperparameters without repeating the full candidate search.

The default run performs the complete search:

- stratified five-fold cross-validation with seed 123;
- 2,592 XGBoost candidates for each score dimension;
- 550 Elastic Net candidates for each score dimension;
- candidate ranking by descending mean fold AUROC, with descending mean fold
  AUPRC as the prespecified tie-breaker;
- out-of-fold predictions for the selected models; and
- refitting of all four selected models using all 1,301 development BSJs.

The complete search fits many models and can take several hours. For a quick
installation and data-contract check, use `--quick-test`. To refit only the
published models without repeating tuning, use `--refit-only`.

## Main outputs

The program writes the following files under the selected output directory:

- `fold_assignments.tsv`
- `XGBoost_tuning.tsv`
- `ElasticNet_tuning.tsv`
- `best_parameters.tsv`
- `out_of_fold_predictions.tsv`
- `out_of_fold_performance.tsv`
- `XGBoost_feature_importance.tsv`
- `ElasticNet_coefficients.tsv`
- `training_fitted_scores.tsv`
- `model_manifest.tsv`
- `reproducibility_checks.tsv`
- `sessionInfo.txt`
- four fitted-model files under `models/`

The expected final hyperparameters and Elastic Net coefficients are stored in
`expected/`. `reproducibility_checks.tsv` compares a completed full run with
these archived values.

## Statistical definitions

AUROC is calculated from the Mann-Whitney/rank-sum identity with average ranks
for tied predictions. AUPRC is calculated as average precision over distinct
prediction thresholds. Fold-level metrics are summarized by their arithmetic
mean and sample standard deviation.

Elastic Net models use binomial logistic regression with an intercept and
`glmnet`'s default predictor standardization during fitting. Reported final
coefficients are returned on the original predictor scale. XGBoost models use
the `binary:logistic` objective and output positive-class probabilities.

## Repository contents

- `CircScoring_model_development.R`: complete tuning, evaluation and refitting
  workflow.
- `prepare_training_data.py`: deterministic extraction of the P1/N1
  development population from the source workbook.
- `install_exact_packages.R`: installs the archived R package versions.
- `Dockerfile`: container recipe for R 4.4.1 and the model workflow.
- `METHODS_model_development.md`: concise manuscript-ready methods text.
- `PACKAGE_VERSIONS_AND_CITATIONS.md`: software versions and citations.
- `PROVENANCE.md`: authoritative sources and excluded alternative artifacts.
- `FEATURE_DEFINITIONS.tsv`: f1-f18 names and CircScoring dimensions.
- `SOFTWARE_ENVIRONMENT.tsv`: software versions used by the workflow.
- `source_artifacts/`: the exact source and archived tuning workbooks.
- `deployment_reference/`: the released prediction programs used to confirm
  that the reconstructed workflow targets the deployed models.

## Reproducibility notes

- The input row order is retained and recorded.
- The archived seed-123 fold membership is supplied under `expected/`; the
  same deterministic, class-stratified algorithm is used if another seed is
  requested.
- The same fold assignments are used for all four models and all candidates.
- XGBoost should be run with one thread for the strictest numerical
  reproducibility across repeated runs.
- Minor floating-point differences can occur across operating systems or CPU
  architectures even when software versions and seeds are fixed.

# SOFTWARE ENVIRONMENT
software	version	purpose
R	4.4.1	Model tuning, evaluation and final fitting
xgboost	3.2.1.1	XGBoost model development
glmnet	4.1-10	Elastic Net binomial logistic regression
Python	>=3.10	Optional extraction of the frozen TSV from Excel
openpyxl	3.1.5	Optional reading of Benchmark_features_lite.xlsx

