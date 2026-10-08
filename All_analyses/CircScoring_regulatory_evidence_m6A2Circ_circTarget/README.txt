# Regulatory evidence: m6A2Circ and circTarget

This compact source bundle reproduces the numerical CircScoring validation
using m6A2Circ confidence levels and circTarget interaction evidence. It
exports group counts, rank summaries and ECDF source data, calculates observed
rank-biserial correlations (RBCs) and median percentile ranks, and performs the
corresponding 10,000-iteration size-matched background-resampling analyses.
Plotting code is intentionally omitted because the final figures have already
been prepared.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved atlas-wide table, complete
10,000-iteration RBC and median null distributions, final figures and
manuscript methods. All group counts and all 16 observed RBCs were independently
verified. The historical seed, group order, shared random-number stream,
sampling direction and `shuffle=False` implementation were recovered from the
retained iteration-level results. The first three retained RBC draws for all
four ranks in the first evidence group were reproduced exactly, and the median
draws were confirmed to use the same sampled BSJs. The script should not be
described as a byte-for-byte copy of the lost historical file.

The retained null-median file was generated from an earlier equivalent circIG
export whose stored percentile ranks differ from the current fixed table only
in inconsequential final decimal places. Consequently, a rerun from the current
table reproduces the analysis and RBC results but may differ in the last shown
decimal of individual random-set medians.

## Required input

The fixed input is an analysis-eligible circIG table containing:

- the m6A2Circ confidence/evidence-count field;
- the number of circTarget sequence types supported by at least two chimeric
  reads;
- four CircScoring percentile ranks; and
- the six benchmark-set indicators.

The reader accepts either a conventional header or the two-row grouped circIG
header and recognizes the column-name variants used in the retained project
files. The retained equivalent input is
`used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv`.
The atlas-wide input is not duplicated in this compact source bundle.

## Exclusions and evidence groups

BSJs overlapping any of the six benchmark sets are excluded before analysis.
Rows lacking any of the four percentile ranks are subsequently excluded. Blank
m6A2Circ and circTarget fields are treated as absence of the corresponding
evidence.

The four evidence-positive groups are:

- m6A2Circ confidence = 1 (`n = 35,543`);
- m6A2Circ confidence = 2 (`n = 5,272`);
- m6A2Circ confidence >= 3 (`n = 1,204`); and
- circTarget with at least one sequence type supported by at least two chimeric
  reads (`n = 2,215`).

The common background consists of BSJs with neither m6A2Circ evidence nor the
specified circTarget evidence (`n = 4,409,523`). The evidence groups are tested
separately; BSJs satisfying both an m6A2Circ and circTarget definition remain
in each applicable positive group.

The m6A2Circ–circTarget overlaps comprise 386 confidence-1, 294 confidence-2
and 162 confidence >=3 BSJs (842 unique overlapping BSJs in total); these are
reported separately by the script and are not removed from either applicable
positive analysis.

The fixed source contains 4,455,197 BSJs. After excluding 2,282
benchmark-overlapping rows, 4,452,915 remain; no additional rows are excluded
for incomplete ranks.

## ECDF and descriptive summaries

For the background and each evidence-positive group, ECDF values are calculated
for all four ranks at thresholds from 0 to 100 in 0.1-rank increments. Mean,
first quartile, median and third quartile are also exported. No figure is
generated.

## Rank-biserial correlation

For each positive set and rank, RBC against the common background is calculated
as

`RBC = 2U / (Npositive × Nbackground) - 1`,

where `U` is the Mann–Whitney statistic with average ranks for ties. Positive
RBC indicates higher ranks in the evidence-positive group.

For the null distribution, a random set matching the relevant positive-set
size is sampled without replacement from the background and compared with the
remaining background after removal of the sampled BSJs. The same sampled BSJs
are used for all four ranks within an iteration.

## Median-rank analysis

The observed statistic is the median percentile rank of the evidence-positive
set. For each null iteration, the statistic is the median percentile rank of
the size-matched random background set. The median is calculated directly for
the sampled set; it is not expressed as a difference from another median.

The same random set is used to calculate both the null RBC and null median in a
given iteration.

## Monte Carlo inference

The evidence groups are processed in this fixed order: m6A confidence = 1,
m6A confidence = 2, m6A confidence >= 3, and circTarget. All groups share one
NumPy `default_rng` stream initialized with seed `20260825`. Each group uses
10,000 size-matched background resamples.

Null 95% intervals are the 2.5th and 97.5th percentiles of the corresponding
null distribution. Upper-tail empirical P values use the +1 correction:

`P = [1 + count(null statistic >= observed statistic)] / (10,000 + 1)`.

This procedure is most precisely described as a Monte Carlo randomization test
based on size-matched background resampling.

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
10,000-iteration analysis is computationally intensive.

## Run

```bash
python regulatory_evidence_m6A2Circ_circTarget.py \
  used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv \
  results
```

The script writes:

- input QC and group-count tables;
- group mean, Q1, median and Q3 summaries;
- ECDF source data;
- RBC resampling summary;
- median-rank resampling summary; and
- optional gzip-compressed iteration-level null distributions.

No figure is generated.

For a quick execution check:

```bash
python regulatory_evidence_m6A2Circ_circTarget.py input.tsv test_results \
  --resamples 3 --no-null-output
```

## Manuscript-ready methods text

Regulatory-evidence validation was performed using m6A2Circ confidence levels
and circTarget interaction evidence after excluding BSJs overlapping any of the
six benchmark sets. BSJs with m6A2Circ confidence levels 1, 2 or at least 3,
and BSJs having at least one circTarget sequence type supported by at least two
chimeric reads, were evaluated separately against a common background lacking
either form of evidence. Differences were quantified using rank-biserial
correlation and median percentile rank. Statistical significance was assessed
using a Monte Carlo randomization test based on 10,000 size-matched background
resamples. Each random set was compared with the remaining background for RBC,
whereas its median percentile rank was calculated directly for the median
analysis. Null 95% intervals were defined by the 2.5th and 97.5th percentiles,
and one-sided empirical P values were calculated using the +1 correction.
