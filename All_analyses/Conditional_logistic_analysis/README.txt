# Conditional logistic-regression analysis

This folder contains the reproducible source code for the CircScoring conditional analysis of 11 evidence datasets against the fixed circIG reference.

## Analysis implemented

- Six benchmark sets are excluded from both the reference and evidence-positive populations.
- Translation potential is defined as `riboCIRC = 1` or `TransCirc evidence count >= 1`.
- circTarget support is defined as `circTarget #type>= 1 (Seq supported by >=2 chimeric reads)`.
- The reference is benchmark-free and negative for all 11 evidence definitions.
- Positive evidence sets are analyzed separately and may overlap other evidence sets.
- Original atlas-wide CS-R and CS-C percentile ranks are modeled as continuous predictors.
- Logistic-regression ORs are scaled per 10-percentile-rank increase.
- Restricted cubic splines use reference-derived 5th, 35th, 65th and 95th percentiles. When duplicated knots occur for discrete CS-R distributions, the 10th, 50th and 90th percentiles of the distinct reference values are used.
- Substantial nonlinearity requires all three criteria: likelihood-ratio P < 0.001, delta AIC <= -10 and delta McFadden R2 >= 0.001.
- Incremental contributions are tested by likelihood-ratio comparison of the mutually adjusted model with the corresponding single-score reduced model, using the functional forms selected for the full model.

## Run

```text
python Conditional_logistic_analysis.py --input path/to/circIG.tsv --output path/to/results
```

By default, the script reads `model_specification.json`, which preserves the
functional-form decisions from the original diagnostic run. This ensures that
the submitted analysis reproduces the reported models even where an original
spline fit was explicitly recorded as nonconverged. To perform model selection
again using the current run, add `--reselect-models`.

Validate population definitions and spline knots without fitting models:

```text
python Conditional_logistic_analysis.py --input path/to/circIG.tsv --output path/to/validation --validate-only
```

For a short computational check, one evidence/model combination can be selected:

```text
python Conditional_logistic_analysis.py --input path/to/circIG.tsv --output path/to/test --evidence Aortic --model XGBoost
```

## Required software

- Python 3.12 or later
- NumPy
- pandas

The script writes the package and Python versions used for each run to
`run_metadata.json`. Logistic regression, Wald intervals, likelihood-ratio
tests, restricted cubic splines and McFadden R2 are implemented explicitly in
the script to make the numerical procedure auditable and minimize dependencies.

## Principal outputs

- `Summary.tsv`
- `Linear_ORs.tsv`
- `Nonlinearity.tsv`
- `Selected_models.tsv`
- `Selected_form_10rank_ORs.tsv`
- `Incremental_tests.tsv`
- `conditional_results.json`
- `QC.tsv`
- `Validation_checks.tsv`
- `run_metadata.json`

The complete analysis is computationally intensive because each evidence set is
compared with more than four million reference BSJs.
