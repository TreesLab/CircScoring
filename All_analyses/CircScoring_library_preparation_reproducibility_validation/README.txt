# CircScoring validation: library preparation and reproducibility

This compact source bundle reproduces the joint numerical workflow used for:

1. RNA-seq library-preparation validation using the RNase R-treated CircAI and
   poly(A)-selected CircRic resources; and
2. cross-condition reproducibility using normal/cancer classifications from
   CSCD2.

The two analyses share one RBC and pooled-label randomization implementation.
Their ECDF source data are exported separately, matching the separate ECDF
panels used in the study. Plotting code is intentionally omitted.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved source table, two ECDF datasets,
complete 10,000-iteration null distributions, summary workbook, final figures,
and manuscript methods. All six group counts and all 16 reported observed RBCs
were verified against the retained outputs. Both ECDF output tables were also
verified row-by-row against their retained counterparts. The script should not
be described as a byte-for-byte copy of the lost historical file.

The retained null-distribution file was generated within an older combined run
that also evaluated a subsequently discarded Type_Length comparison before the
reported analyses, thereby advancing the shared random-number stream. The
present submission script starts the documented seed at the beginning of the
reported CircAI/CircRic and CSCD2 workflow. It implements the same
randomization test, but individual Monte Carlo draws and the last decimal places
of the null limits need not be identical to the older combined-run file.
Observed RBCs and ECDF values are unchanged.

## Required input

The fixed input is the analysis-eligible circIG export used for the reported
analysis. It may have either a conventional header or the two-row grouped header
produced by circIG. Required columns are:

- `CircAI`, `CircRic`, and `CSCD2_group`;
- `XGBoost_CS-R_percentile_rank` and `XGBoost_CS-C_percentile_rank`;
- `ElasticNet_CS-R_percentile_rank` and
  `ElasticNet_CS-C_percentile_rank`;
- `positive1_P1N1`, `positive2_P4N4_L`, `positive3_P4N4_H`;
- `negative1_P1N1`, `negative2_P4N4_L`, `negative3_P4N4_H`.

The retained input was `filtered-20260917-120220.tsv`. This 473-MB atlas-wide
file is not duplicated in the compact source-code bundle.

## Exclusions and group definitions

BSJs overlapping any of the six benchmark sets are excluded before either
analysis. Rows lacking any of the four percentile ranks are then excluded.
Blank CircAI or CircRic fields are treated as zero.

Library-preparation groups:

- CircAI only: `CircAI = 1` and `CircRic = 0` (`n = 299,005`).
- CircRic only: `CircAI = 0` and `CircRic = 1` (`n = 30,720`).
- CircAI–CircRic shared: `CircAI = 1` and `CircRic = 1` (`n = 3,892`).

CSCD2 reproducibility groups:

- common to normal and cancer samples (`n = 291,809`);
- normal-specific (`n = 1,207,723`); and
- cancer-specific (`n = 806,147`).

The source contains 4,455,197 BSJs. After excluding 2,282 benchmark-overlapping
BSJs, 4,452,915 remain; no additional rows are excluded for incomplete ranks.

## ECDF analyses

ECDF values are calculated independently for the three CircAI/CircRic groups
and the three CSCD2 groups at fixed percentile-rank thresholds from 0 to 100 in
increments of 0.1. Two separate source-data files are written:

- `CircAI_CircRic_ECDF_data_0.1.tsv`; and
- `CSCD2_group_ECDF_data_0.1.tsv`.

The script does not combine or plot these ECDFs.

## RBC and randomization analyses

For each of the four CircScoring percentile ranks, rank-biserial correlation is
calculated as

`RBC = 2U / (n1 × n2) - 1`,

where `U` is the Mann–Whitney statistic for Group 1 calculated from pooled
average ranks. Positive RBC indicates higher percentile ranks in Group 1.

The four reported comparisons are:

1. CircAI only versus CircRic only;
2. CircAI–CircRic shared versus CircRic only;
3. CSCD2 common versus normal-specific; and
4. CSCD2 common versus cancer-specific.

For each comparison, both groups are pooled and Group 1 labels are reassigned
without replacement while preserving both observed group sizes. The same
randomized labels are applied to all four percentile-rank variables within each
iteration. Comparisons are processed by the same script and the same NumPy
random-number generator. The procedure is repeated 10,000 times using
`default_rng` with seed `20260917`.

Null 95% intervals are the 2.5th and 97.5th percentiles of the randomized RBCs.
The one-sided empirical P value uses the +1 correction:

`P = [1 + count(RBCnull ≥ RBCobserved)] / (10,000 + 1)`.

The optional `--all-pairs` switch additionally evaluates CircAI–CircRic shared
versus CircAI only. This comparison appears in one Methods draft but is absent
from the retained summary table and final RBC panel; it is not part of the
default reproduction.

## Software environment and requirements

The reconstruction was executed and checked with:

- Python 3.12.14
- NumPy 2.3.5
- pandas 3.0.1

Install the exact package versions with:

```bash
python -m pip install numpy==2.3.5 pandas==3.0.1
```

Approximately 4–8 GB of available RAM is recommended for the atlas-wide TSV.

## Run

```bash
python library_preparation_and_reproducibility_validation.py \
  filtered-20260917-120220.tsv \
  results
```

The script writes separate ECDF source tables, combined group counts and QC,
the combined RBC summary, and the gzip-compressed iteration-level null
distributions. No figure is generated.

For a quick code test, reduce the resampling count, for example:

```bash
python library_preparation_and_reproducibility_validation.py \
  filtered-20260917-120220.tsv test_results \
  --resamples 20 --no-null-output
```

## Manuscript-ready methods text

CircScoring was evaluated against expected library-preparation and
cross-condition reproducibility patterns using BSJs from CircAI, CircRic and
CSCD2. BSJs overlapping any of the six benchmark sets were excluded before
analysis. For library-preparation validation, BSJs were classified as
CircAI-only, shared between CircAI and CircRic, or CircRic-only. For
cross-condition reproducibility, CSCD2 BSJs were classified as common to normal
and cancer samples, normal-specific, or cancer-specific. ECDFs were calculated
separately for the CircAI/CircRic and CSCD2 groupings. Pairwise differences were
quantified using rank-biserial correlation for CircAI-only versus CircRic-only,
shared versus CircRic-only, common versus normal-specific and common versus
cancer-specific BSJs. Null intervals and one-sided empirical P values were
estimated using 10,000 pooled-label randomization iterations that preserved the
observed group sizes, with the same randomized labels applied to all four
percentile-rank variables within each iteration.
