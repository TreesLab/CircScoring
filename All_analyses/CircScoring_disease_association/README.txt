# Disease association

This compact source bundle reproduces the final CircScoring disease-
association analysis. It exports ECDF source data and descriptive rank
summaries, calculates observed rank-biserial correlations (RBCs), and performs
10,000 background-based size-matched resamples for RBC and median percentile
rank. Plotting code is intentionally omitted because the final figures have
already been prepared.

## Historical fidelity

The original execution script was not retained. The included Python script is
a faithful reconstruction from the preserved annotated circIG data, complete
iteration-level null distributions, method record and final figures. The
retained sample sizes, all 28 observed RBCs, random-number seed, group order and
`shuffle=False` sampling implementation were recovered. The retained null RBC
and median values are reproduced when the script is run with 10,000 resamples.
The script should not be described as a byte-for-byte copy of the lost
historical file.

## Required input

The retained analysis can be run from
`used_circIG_with_Translation_Disease_circAge_tissue_and_experiment_evidence_TransCirc_N0_corrected_annotation_guided.tsv`.
The earlier input name was
`used_circIG_with_CircR2Disease_and_Circ2Disease.tsv`. The reader supports the
retained two-row grouped header and a conventional one-row header. It also
accepts either `unique_disease_count` or
`circRNADisease-unique_disease_count` for the disease-count field.

The large fixed input is not duplicated in this source bundle.

## Exclusions and group definitions

BSJs positive in any of the six benchmark-set columns are excluded from all
groups. Rows missing any of the four CircScoring percentile ranks are also
excluded.

The common background consists of benchmark-free BSJs with no evidence in any
of the following fields:

- `CircR2Disease`;
- `Circ2Disease`;
- `LncRNADisease`;
- `circad`; and
- `unique_disease_count`.

Seven evidence groups are compared independently with this background:

- `CircR2Disease = 1`;
- `Circ2Disease = 1`;
- `LncRNADisease = 1`;
- `circad = 1`;
- `unique_disease_count = 1`;
- `unique_disease_count = 2`; and
- `unique_disease_count >= 3`.

Evidence groups may overlap one another. They are not made mutually exclusive.
After benchmark exclusion, the shared background contains 4,449,090 BSJs. The
seven positive-set sizes are 1,275, 166, 196, 391, 2,572, 505 and 387,
respectively. The identical `N_background` reported for every comparison is
therefore intentional.

## ECDF and descriptive summaries

ECDF values are evaluated separately for all four CircScoring percentile ranks
at thresholds from 0 to 100 in increments of 0.1. The script also exports N,
mean, first quartile, median and third quartile for every group and rank.

## RBC and background-based size-matched resampling

Observed RBC compares an evidence-positive set with the complete shared
background. RBC is oriented as evidence positive relative to background, so a
positive value indicates higher percentile ranks in the evidence-positive set.
For a positive set of size N1 and background of size N0,

`RBC = 2U / (N1 × N0) - 1`,

with average ranks assigned to ties.

For each evidence group, N1 BSJs are sampled without replacement from the
shared background. The sampled BSJs form the random positive set and are
compared with `background minus the sampled BSJs`, keeping the two randomized
comparison groups mutually exclusive. The same sampled BSJs are used for all
four ranks within a resample. This procedure is repeated 10,000 times.

Null 95% intervals are the 2.5th and 97.5th percentiles of the resampled RBCs.
The one-sided empirical P value uses the +1 correction:

`P = [1 + count(RBCnull >= RBCobserved)] / (10,000 + 1)`.

## Median-rank resampling

For the median analysis, the median percentile rank is calculated directly for
each size-matched random set drawn from the shared background. It is not a
difference in medians. The null 95% interval and upper-tail empirical P value
are calculated in the same way as above, using the observed positive-set
median as the test statistic.

## Reproducibility details

The final analysis used `numpy.random.default_rng(20260827)`. Evidence groups
were processed in the order shown above. Sampling used
`Generator.choice(..., replace=False, shuffle=False)`. One continuous generator
stream was used across all seven groups. These details are required to
reproduce the retained iteration-level null tables exactly.

## Software environment and requirements

The reconstruction was prepared and checked with:

- Python 3.12.14
- NumPy 2.3.5
- pandas 3.0.1

Install the exact package versions with:

```bash
python -m pip install numpy==2.3.5 pandas==3.0.1
```

At least 3–4 GB of available RAM is recommended for the 4.46-million-row
input. The full 10,000-resample run may take several minutes depending on CPU
and storage speed.

## Run

```bash
python disease_association.py annotated_circIG.tsv results
```

The script writes:

- input QC, group counts and positive-group overlap counts;
- ECDF source data at 0.1-rank intervals;
- group/rank descriptive summaries;
- observed RBC, null intervals and empirical P values;
- observed median ranks, null intervals and empirical P values; and
- optional gzip-compressed iteration-level RBC and median null distributions.

No figure is generated. For a short execution check:

```bash
python disease_association.py annotated_circIG.tsv test_results \
  --resamples 3 --no-null-output
```

Because one continuous random-number stream is used across groups, a run with
fewer than 10,000 resamples is only a functional test; its later-group random
draws are not prefixes of the retained 10,000-resample tables.

## Manuscript-ready methods text

Disease association was evaluated after excluding BSJs overlapping any of the
six benchmark sets. A shared background was defined as benchmark-free BSJs
without evidence in CircR2Disease, Circ2Disease, LncRNADisease, circad or
circRNADisease v3.0. Seven evidence groups were evaluated independently:
CircR2Disease, Circ2Disease, LncRNADisease and circad annotations, and
circRNADisease v3.0 support for one, two or at least three unique diseases.
CircScoring percentile-rank distributions were summarized using ECDFs. Effect
sizes were quantified using rank-biserial correlation relative to the shared
background. Statistical significance was assessed using 10,000
background-based size-matched resamples. In each iteration, a random set equal
in size to the corresponding positive set was sampled without replacement from
the background and compared with the remaining background after removal of the
sampled BSJs. Null 95% intervals were defined by the 2.5th and 97.5th
percentiles, and one-sided empirical P values were calculated using the +1
correction. Median-rank null distributions were generated by directly
calculating the median percentile rank of each size-matched random set.
