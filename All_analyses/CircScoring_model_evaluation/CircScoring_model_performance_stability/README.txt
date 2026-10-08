# CircScoring model-performance stability analysis

This directory contains the submission-ready, reproducible source code for the
CircScoring **model performance stability** analysis.

## Historical fidelity

The historical script itself was not retained. The code here is a faithful
reconstruction from the preserved source table, complete 1,000-partition result
table, workbook metadata, random seeds, and published figure assets. The
reconstruction was verified against all 1,000 retained partitions and all 24,000
reported AUROC/AUPRC/difference values. It should not be described as a
byte-for-byte copy of the lost historical script.

The retained records establish the following exact procedure:

- P1N1 contained 1,301 BSJs: 797 P1 positives and 504 N1 negatives.
- The four finalized CircScoring score vectors were fixed before this analysis.
- No model refitting or hyperparameter optimization occurred during stability
  assessment.
- A master NumPy generator initialized with seed `20260801` generated 1,000
  unique iteration seeds.
- In each iteration, all 1,301 records were randomly permuted without
  stratification or replacement. The first 1,040 records formed the 80% subset;
  the remaining 261 records formed the mutually exclusive 20% subset.
- AUROC used the Mann-Whitney/rank formulation with average ranks for ties.
- AUPRC used trapezoidal integration of the threshold-based precision-recall
  curve; it was not calculated as average precision.
- Stability was summarized using the distributions of the 80% and 20% metrics
  and their within-iteration differences (`80% - 20%`).

## Files

- `01_model_performance_stability.py`: exact numerical analysis, workbook/CSV
  export, and optional verification against the retained result table.
- `02_plot_model_performance_stability.R`: base-R plotting script for transparent
  SVG and PNG violin/box plots.
- `data/P1N1_fixed_scores.tsv`: the 1,301 P1N1 labels and four fixed scores used
  by the stability analysis.
- `expected/model_performance_stability_1000_retained.xlsx`: retained full-
  precision reference output used for exact regression testing.
- `METHODS_text.md`: concise manuscript-ready description.
- `requirements.txt`: Python dependency versions for the verified environment.

## Run

From this directory:

```bash
python 01_model_performance_stability.py \
  --input data/P1N1_fixed_scores.tsv \
  --expected expected/model_performance_stability_1000_retained.xlsx \
  --outdir results

Rscript 02_plot_model_performance_stability.R \
  results/model_performance_stability_1000_partitions.csv \
  results/figures
```

The first command stops with an error if the regenerated partitions, seeds,
class counts, or metrics differ from the retained results.

## Requirements

numpy==2.3.5
pandas==3.0.1
openpyxl==3.1.5
