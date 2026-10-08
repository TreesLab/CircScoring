# Orthogonal experimental evidence: RT-independent and RT-/non-RT validation

This compact source bundle reproduces the numerical CircScoring validation
using two orthogonal experimental evidence sets:

1. RT-independent NCL events; and
2. BSJs validated by both reverse-transcription-dependent and
   reverse-transcription-independent approaches (RT-/non-RT evidence).

The script exports group counts and ECDF source data, calculates observed
rank-biserial correlations (RBCs) and median-rank differences, and performs the
corresponding 10,000-iteration size-matched background-resampling tests.
Plotting code is intentionally omitted because the final figures have already
been prepared.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved atlas-wide source table, complete
10,000-iteration null distributions, summary tables, final figures, and
manuscript methods. All population counts, all eight observed RBCs, and all
positive/background medians were independently verified against the retained
outputs. The script should not be described as a byte-for-byte copy of the
lost historical file.

The retained historical Monte Carlo files were generated from an earlier
atlas-wide export whose row order is no longer available. Because random
sampling operates on row positions, a rerun from the current equivalent table
can produce a different—but statistically equivalent—sequence of null draws
and minor differences in Monte Carlo interval endpoints. The observed RBCs,
medians, group sizes, test definition, and empirical-P-value calculation are
unchanged.

## Required input

The fixed input is the analysis-eligible circIG table containing:

- `RT_independent_NCL_event`;
- `Both RT/nonRT`;
- `XGBoost_CS-R_percentile_rank` and
  `XGBoost_CS-C_percentile_rank`;
- `ElasticNet_CS-R_percentile_rank` and
  `ElasticNet_CS-C_percentile_rank`;
- `positive1_P1N1`, `positive2_P4N4_L`, `positive3_P4N4_H`;
- `negative1_P1N1`, `negative2_P4N4_L`, `negative3_P4N4_H`.

The reader accepts either a conventional header or the two-row grouped circIG
header. The retained equivalent input is
`used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv`.
This atlas-wide file is not duplicated in the compact source-code bundle.

## Exclusions and group definitions

BSJs overlapping any of the six benchmark sets are excluded before analysis.
Rows lacking any of the four percentile ranks are subsequently excluded. Blank
evidence fields are treated as zero.

- Background: `RT_independent_NCL_event = 0` and `Both RT/nonRT = 0`.
- RT-independent evidence: `RT_independent_NCL_event = 1`.
- RT-/non-RT evidence: `Both RT/nonRT = 1`.

The positive sets are analyzed independently; a BSJ satisfying both evidence
definitions is retained in both applicable positive sets. In the fixed input:

- 4,455,197 BSJs are present;
- 2,282 benchmark-overlapping BSJs are excluded;
- 4,449,931 BSJs form the background;
- 2,493 BSJs have RT-independent evidence; and
- 599 BSJs have RT-/non-RT evidence.

After benchmark exclusion, 108 BSJs belong to both positive sets; consistent
with the prespecified analysis, they are retained in each separate comparison.

## ECDF source data

For each group and each of the four CircScoring percentile ranks, ECDF values
are evaluated at thresholds from 0 to 100 in 0.1-rank increments. The resulting
source table is exported but is not plotted by this script.

## Rank-biserial correlation

For each positive set and rank, RBC against the background is calculated as

`RBC = 2U / (Npositive × Nbackground) - 1`,

where `U` is the Mann–Whitney statistic for the positive set, with half credit
assigned to tied positive-background pairs. Positive RBC indicates higher
percentile ranks in the evidence-positive set.

For the null distribution, a set equal in size to the relevant positive set is
sampled without replacement from the background. Its RBC is calculated against
the remaining background after removal of the sampled BSJs, ensuring mutually
exclusive groups. The same sampled BSJs are used for all four ranks within an
iteration.

## Median-rank analysis

The observed statistic is

`median(positive set) - median(full background)`.

For each null iteration, the median of the size-matched random background set
is compared with the median of the remaining background after removal of that
set:

`median(random set) - median(background minus random set)`.

The complement median is calculated exactly; it is not approximated by the
median of the unmodified background.

## Monte Carlo inference

The two evidence sets are tested separately. RBC and median analyses are
implemented as separate Monte Carlo randomization runs, each initialized using
NumPy `default_rng` with seed `20260815`. Each run uses 10,000 size-matched
background resamples.

Null 95% intervals are the 2.5th and 97.5th percentiles of the null estimates.
The one-sided empirical P value uses the +1 correction:

`P = [1 + count(null statistic >= observed statistic)] / (10,000 + 1)`.

This procedure is most precisely described as a Monte Carlo randomization test
based on size-matched background resampling, rather than a conventional
pooled-label permutation test.

## Software environment and requirements

The reconstruction was prepared and checked with:

- Python 3.12.14
- NumPy 2.3.5
- pandas 3.0.1

Install the exact package versions with:

```bash
python -m pip install numpy==2.3.5 pandas==3.0.1
```

Approximately 4–8 GB of available RAM is recommended. The complete
10,000-iteration RBC and median analyses are computationally intensive.

## Run

```bash
python orthogonal_experimental_evidence.py \
  used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv \
  results
```

The script writes:

- group-count and QC tables;
- ECDF source data at 0.1-rank intervals;
- RBC resampling summary;
- median-rank resampling summary; and
- optional gzip-compressed iteration-level null distributions.

No figure is generated.

For a quick execution check:

```bash
python orthogonal_experimental_evidence.py input.tsv test_results \
  --resamples 2 --no-null-output
```

## Manuscript-ready methods text

To evaluate CircScoring against orthogonal experimental evidence, BSJs with
RT-independent NCL support or validation by both RT-dependent and
RT-independent approaches were compared separately with a background lacking
either form of evidence. BSJs overlapping any of the six benchmark sets were
excluded. Differences in percentile-rank distributions were quantified using
rank-biserial correlation and the difference in median percentile rank.
Statistical significance was assessed using a Monte Carlo randomization test
based on 10,000 size-matched background resamples. In each iteration, a random
set matching the corresponding positive-set size was sampled without
replacement from the background and compared with the remaining background
after removal of the sampled BSJs. Null 95% intervals were defined by the 2.5th
and 97.5th percentiles of the resampled statistics, and one-sided empirical P
values were calculated using the +1 correction.
