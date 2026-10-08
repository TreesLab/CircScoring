# Coding-potential validation: riboCIRC v2.0 and TransCirc

This compact source bundle reproduces the numerical CircScoring analyses using
riboCIRC v2.0 and TransCirc coding-potential evidence. The two resources use
different prespecified comparison and randomization designs and therefore
remain separate throughout the script. Separate ECDF data, RBC summaries and
null distributions are written for riboCIRC and TransCirc. Plotting code is
intentionally omitted because their RBC results were displayed in separate
figures that have already been prepared.

## Historical fidelity

The original execution scripts were not retained. The included Python script
is a faithful reconstruction from the preserved atlas-wide input, ECDF source
tables, complete iteration-level null distributions, RBC summary tables,
workbooks, figures and analysis notes. All group counts and all 32 observed
RBCs were checked against the retained outputs. The retained null draws were
also reproduced exactly, confirming the group ordering, random-number seeds,
sampling implementation and label-randomization direction. The script should not be
described as a byte-for-byte copy of the lost historical files.

## Required input and common exclusions

The fixed input is the analysis-eligible circIG table. The reader accepts a
conventional header or the two-row grouped circIG header. Required information
includes:

- `riboCIRC`;
- riboCIRC annotation-guided, context-specific and literature-reported flags;
- the TransCirc evidence-count field;
- four CircScoring percentile ranks; and
- the six benchmark-set indicators.

Common column-name variants used in the project are recognized automatically.
The retained input is
`used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv`.
The atlas-wide input is not duplicated in this compact source bundle.

Before either analysis, BSJs overlapping any of the six benchmark sets are
excluded. Rows lacking any of the four percentile ranks are also excluded. The
fixed input contains 4,455,197 BSJs; 2,282 benchmark-overlapping BSJs are
excluded, leaving 4,452,915 analysis-eligible rows.

## riboCIRC v2.0 analysis

Three evidence-positive sets are analyzed separately:

- annotation-guided (`n = 3,085`);
- context-specific (`n = 56`); and
- literature-reported (`n = 52`).

The common riboCIRC background comprises all eligible circIG BSJs with
`riboCIRC != 1`, additionally excluding every BSJ positive for any of the three
analyzed evidence categories. This gives 4,449,774 background BSJs and ensures
that each positive-background comparison is mutually exclusive. Two
literature-reported BSJs have `riboCIRC = 0`; they remain literature positives
and are excluded from the background.

For each evidence set and rank, RBC against the background is calculated as

`RBC = 2U / (Npositive × Nbackground) - 1`,

with average ranks for ties. Positive RBC indicates higher ranks in the
evidence-positive set.

For null inference, a set matching the evidence-positive sample size is drawn
without replacement from the background and compared with the remaining
background after removal of the sampled BSJs. The same sampled rows are used
for all four ranks within an iteration. This is repeated 10,000 times.
Independent fixed seeds are used for the three evidence sets:

- annotation-guided: `2026091701`;
- context-specific: `2026091702`;
- literature-reported: `2026091703`.

## TransCirc analysis

Only the 932 coordinate-verified BSJs explicitly assigned
`TransCirc_evidences_num = 0` are eligible for N0. Three overlap a benchmark
set, leaving 929 N0 BSJs. Blank TransCirc values are not treated as zero and
are excluded from all TransCirc comparisons.

N0 is compared separately with:

- evidence count `= 1` (`n = 14,122`);
- evidence count `= 2` (`n = 16,485`);
- evidence count `= 3` (`n = 219,878`);
- evidence count `= 4` (`n = 37,309`); and
- evidence count `>= 5` (`n = 1,661`).

For each comparison, N0 and the positive group are pooled. Labels are randomly
reassigned while preserving both observed group sizes exactly, and RBC is
recalculated with the positive group as Group 1. The same randomized labels are
used for all four ranks within an iteration. Each comparison uses 10,000
iterations and its own fixed seed (`2026091601` through `2026091605`). This is a
pooled-label randomization test, unlike the background-resampling design used
for riboCIRC.

## ECDFs and Monte Carlo inference

ECDF source values are evaluated independently for each resource, evidence
group and rank at thresholds from 0 to 100 in increments of 0.1. Separate
source-data files are exported; no figure is generated.

For both analyses, null 95% intervals are the 2.5th and 97.5th percentiles of
the 10,000 null RBCs. One-sided empirical P values use the +1 correction:

`P = [1 + count(RBCnull >= RBCobserved)] / (10,000 + 1)`.

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
10,000-iteration workflow is computationally intensive.

## Run

```bash
python coding_potential_riboCIRC_TransCirc.py \
  used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv \
  results
```

The script writes separate riboCIRC and TransCirc files for:

- group counts and rank summaries;
- ECDF source data;
- observed RBC and resampling summaries; and
- optional gzip-compressed iteration-level null distributions.

No figure is generated, and riboCIRC and TransCirc RBC outputs are not merged.

For a quick execution check:

```bash
python coding_potential_riboCIRC_TransCirc.py input.tsv test_results \
  --resamples 2 --no-null-output
```

## Manuscript-ready methods text

Coding-potential validation was performed separately using riboCIRC v2.0 and
TransCirc after excluding BSJs overlapping any of the six benchmark sets. For
riboCIRC, annotation-guided, context-specific and literature-reported BSJs were
each compared with a mutually exclusive circIG background lacking riboCIRC
membership and any of the three positive evidence flags. Null RBC distributions
were generated using 10,000 size-matched background resamples, in which each
random set was compared with the remaining background after removal of the
sampled BSJs. For TransCirc, BSJs with one, two, three, four or at least five
evidence types were compared separately with coordinate-verified BSJs having an
explicit evidence count of zero; blank evidence-count fields were excluded.
For each TransCirc comparison, group labels were randomly reassigned 10,000
times within the pooled N0 and evidence-positive sets while preserving the
observed group sizes. Null 95% intervals were defined by the 2.5th and 97.5th
percentiles, and one-sided empirical P values were calculated using the +1
correction.
