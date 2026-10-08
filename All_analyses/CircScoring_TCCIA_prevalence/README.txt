# Reproducibility: circRNA prevalence in TCCIA

This compact source bundle reproduces the numerical analysis relating the
number of TCCIA cohorts supporting each BSJ to its CircScoring percentile
ranks. It performs the BSJ-level Spearman analysis, paired nonparametric
bootstrap confidence intervals, and cohort-group median/IQR summaries. Plotting
code is intentionally omitted because the final figure has already been
prepared.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved fixed input, complete
10,000-iteration bootstrap distribution, summary table, median/IQR table,
figure, and manuscript methods. The retained population counts, four observed
Spearman correlations, the bootstrap algorithm and the cohort-group summaries
were independently verified. The script should not be described as a
byte-for-byte copy of the lost historical file.

## Required input

The retained input is `TCCIA_with_Cohort_count.tsv`, containing 500,535 BSJs.
The reader accepts either a conventional header or the two-row grouped circIG
header. Required columns are:

- `Cohort_count`;
- `XGBoost_CS-R_percentile_rank` and
  `XGBoost_CS-C_percentile_rank`;
- `ElasticNet_CS-R_percentile_rank` and
  `ElasticNet_CS-C_percentile_rank`;
- `positive1_P1N1`, `positive2_P4N4_L`, `positive3_P4N4_H`;
- `negative1_P1N1`, `negative2_P4N4_L`, `negative3_P4N4_H`.

The 61-MB fixed input is not duplicated in this compact source-code bundle.

## Exclusions and analysis population

BSJs overlapping any of the six benchmark sets are excluded. Rows missing
`Cohort_count` or any of the four percentile ranks are subsequently excluded.
The retained analysis contains:

- 500,535 source BSJs;
- 1,813 benchmark-overlapping BSJs excluded;
- no additional incomplete rows; and
- 498,722 BSJs in the final analysis.

## Spearman correlation and P value

For each of the four CircScoring percentile ranks, Spearman's correlation is
calculated at the individual-BSJ level between `Cohort_count` and percentile
rank. Average ranks are assigned to tied values, and Spearman's rho is the
Pearson correlation between these two average-rank variables. It is not
calculated from cohort-group medians.

The two-sided asymptotic P value uses

`t = rho × sqrt[(N - 2) / (1 - rho²)]`, with `df = N - 2`.

The logarithmic Student-t survival probability is retained so that extremely
small P values can be reported without numerical underflow.

## Paired nonparametric bootstrap

Confidence intervals are estimated using 10,000 paired BSJ-level bootstrap
resamples with seed `20260917`. In each iteration, a multinomial vector of BSJ
multiplicities is sampled from

`Multinomial[N; (1/N, ..., 1/N)]`.

This is mathematically equivalent to sampling N BSJ rows with replacement. A
single weight vector is applied simultaneously to `Cohort_count` and all four
percentile-rank variables, thereby preserving the pairing among measurements
from the same BSJ. Bootstrap rho is calculated as the weighted Pearson
correlation of the original average-rank variables. The 2.5th and 97.5th
percentiles of the 10,000 bootstrap estimates define the 95% confidence
interval. This implementation exactly reproduces the retained bootstrap
algorithm.

## Cohort-recurrence summaries

For presentation, BSJs are grouped as cohort counts 1 through 10, 11–15, and
at least 16. Within each group and for each rank, the script reports the number
of BSJs, median, first quartile and third quartile. These grouped summaries are
used only for presentation and are not used to calculate Spearman's rho.

## Software environment and requirements

The reconstruction was prepared and checked with:

- Python 3.12.14
- NumPy 2.3.5
- pandas 3.0.1

Install the requirements with:

```bash
python -m pip install numpy==2.3.5 pandas==3.0.1
```

Approximately 2–4 GB of available RAM is recommended. The complete 10,000-run
bootstrap is computationally intensive; `--bootstraps` may be reduced only for
code testing, not for reproducing the reported results.

## Run

```bash
python tccia_prevalence_validation.py \
  TCCIA_with_Cohort_count.tsv \
  results
```

The script writes:

- the four-row Spearman/bootstrap summary;
- the cohort-group median/IQR table; and
- the gzip-compressed iteration-level bootstrap distribution.

No figure is generated.

For a quick execution check:

```bash
python tccia_prevalence_validation.py \
  TCCIA_with_Cohort_count.tsv test_results \
  --bootstraps 2 --no-bootstrap-output
```

## Manuscript-ready methods text

To assess the relationship between circRNA prevalence and CircScoring
confidence, the number of TCCIA cohorts in which each BSJ was detected was
compared with its CS-R and CS-C percentile ranks after excluding BSJs
overlapping any of the six benchmark sets. Spearman's rank correlation was
calculated at the individual-BSJ level between cohort count and the
corresponding percentile rank. Confidence intervals for rho were estimated
using 10,000 paired nonparametric bootstrap resamples of BSJs, with the 2.5th
and 97.5th percentiles of the bootstrap distribution defining the 95%
confidence interval. For visualization, percentile ranks were summarized by
cohort-recurrence category using the median and interquartile range; these
grouped summaries were not used in the correlation analysis.
