# Aging relevance: circAge

This compact source bundle reproduces the final five-tissue CircScoring aging-
relevance analysis based on circAge. It exports ECDF source data and
descriptive rank summaries, performs two-sided tie-corrected Wilcoxon rank-sum
tests, calculates observed rank-biserial correlations (RBCs), and generates
10,000-iteration pooled-label null distributions. Plotting code is
intentionally omitted because the final figures have already been prepared.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved fixed input, complete
iteration-level null distribution, summary workbook, final figures and
manuscript methods. All tissue-specific sample sizes and all 20 observed RBCs
were verified against the retained results. The historical group ordering,
five random-number seeds, status-label direction and `shuffle=False`
implementation were recovered; the first three retained null draws for all
four ranks were reproduced exactly. The script should not be described as a
byte-for-byte copy of the lost historical file.

## Required input

The retained input is
`合併circAge_circIG座標rank_tissue_age_status.tsv`. It contains:

- four CircScoring percentile ranks;
- the six benchmark-set indicators; and
- explicit 0/1/blank circAge status for Aortic, Lung, Skin, Umbilical and
  white-blood-cell tissues.

The reader accepts either a conventional header or the retained two-row grouped
header. The 19-MB fixed input is not duplicated in this compact source bundle.

## Exclusions and group definitions

BSJs overlapping any of the six benchmark sets are excluded before every
tissue-specific comparison. Rows missing any of the four percentile ranks are
then excluded.

Within each tissue, explicit status 1 (age-related) is compared with explicit
status 0 (not age-related). Blank tissue status is excluded from that tissue's
comparison rather than treated as zero. Because each tissue is analyzed
independently, a BSJ may contribute to more than one tissue-specific analysis.

After benchmark exclusion, the group sizes are:

- Aortic: N0 = 20,844; N1 = 142;
- Lung: N0 = 73,406; N1 = 3,914;
- Skin: N0 = 53,851; N1 = 702;
- Umbilical: N0 = 25,121; N1 = 130; and
- white blood cell: N0 = 61,936; N1 = 411.

The input contains 171,686 rows; 1,504 benchmark-overlapping rows are excluded,
leaving 170,182 rows before tissue-specific removal of blank statuses.

## ECDF and Wilcoxon analyses

For each tissue, status and rank, ECDF values are evaluated at thresholds from
0 to 100 in increments of 0.1. Mean, first quartile, median and third quartile
are also exported.

Status 1 and status 0 distributions are compared using a two-sided
Mann–Whitney/Wilcoxon rank-sum test based on pooled average ranks. The
asymptotic variance is corrected for ties, and a 0.5 continuity correction is
applied toward the null mean. The output includes U, the tie-corrected Z
statistic, the two-sided P value and log10(P). This test compares the full rank
distributions; it is not specifically a test of medians.

## RBC and pooled-label randomization

RBC is calculated as

`RBC = 2U / (N1 × N0) - 1`,

where status 1 is Group 1 and average ranks are assigned to ties. Positive RBC
indicates higher percentile ranks among age-related BSJs.

For each tissue, status-0 and status-1 BSJs are pooled and status-1 labels are
randomly reassigned without replacement while preserving N0 and N1 exactly.
The same randomized labels are used for all four ranks within an iteration.
The procedure is repeated 10,000 times using independent fixed seeds:

- Aortic: `2026091801`;
- Lung: `2026091802`;
- Skin: `2026091803`;
- Umbilical: `2026091804`; and
- white blood cell: `2026091805`.

Null 95% intervals are the 2.5th and 97.5th percentiles of the randomized RBCs.
The one-sided empirical P value uses the +1 correction:

`P = [1 + count(RBCnull >= RBCobserved)] / (10,000 + 1)`.

The earlier size-matched background sampling of median ranks was exploratory
and is not part of this final pooled-label circAge workflow.

## Software environment and requirements

The reconstruction was prepared and checked with:

- Python 3.12.14
- NumPy 2.3.5
- pandas 3.0.1

Install the exact package versions with:

```bash
python -m pip install numpy==2.3.5 pandas==3.0.1
```

Approximately 1–2 GB of available RAM is recommended.

## Run

```bash
python aging_relevance_circAge.py \
  合併circAge_circIG座標rank_tissue_age_status.tsv \
  results
```

The script writes:

- input QC and five-tissue group counts;
- ECDF source data;
- tissue/status rank summaries;
- tie-corrected Wilcoxon results;
- observed RBC and randomization summaries; and
- an optional gzip-compressed iteration-level RBC null distribution.

No figure is generated.

For a quick execution check:

```bash
python aging_relevance_circAge.py input.tsv test_results \
  --randomizations 3 --no-null-output
```

## Manuscript-ready methods text

Aging relevance was assessed using circAge annotations for aortic, lung, skin,
umbilical and white-blood-cell tissues. BSJs overlapping any of the six
benchmark sets were excluded. Within each tissue, BSJs explicitly annotated as
age-related were compared with those explicitly annotated as not age-related;
blank tissue annotations were excluded. Percentile-rank distributions were
summarized using ECDFs and compared using two-sided Wilcoxon rank-sum tests with
tie-corrected asymptotic P values. Effect sizes were quantified using
rank-biserial correlation. Statistical significance of RBC was assessed by
10,000 tissue-specific pooled-label randomizations preserving the observed
group sizes. Null 95% intervals were defined by the 2.5th and 97.5th
percentiles, and one-sided empirical P values were calculated using the +1
correction.
