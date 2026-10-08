# Joint CircScoring confidence landscape

This folder contains the reproducible analysis and figure source code for the
joint CS-R/CS-C confidence landscape and 10 x 10 enrichment analyses.

## Files

- `Joint_confidence_enrichment.py`: population construction, 10 x 10 matrices, 10,000 size-matched reference resamples, positive-enrichment tests, BH-FDR and TVD.
- `plot_Joint_confidence_landscape.py`: transparent 6 x 6 cm reference and evidence bubble matrices.
- `requirements.txt`: Python dependencies.

## Analysis definitions

- All evidence-positive sets exclude BSJs overlapping any of the six benchmark sets.
- Translation potential is `riboCIRC = 1 OR TransCirc evidence count >= 1`.
- circTarget is `#type of Seq supported by >=2 chimeric reads >= 1`.
- The reference is benchmark-free and negative for all 11 evidence definitions.
- Atlas-wide CS-R and CS-C percentile ranks are divided into `[0,10)`, ..., `[90,100]`.
- Each null set is sampled without replacement from the reference and compared with the remaining reference after removal of the sampled BSJs.
- Cell enrichment is `log2[(observed + 0.5)/(expected + 0.5)]`.
- Only cells with observed log2 enrichment greater than zero receive an upper-tail empirical P value; other cells are assigned P = 1.
- Empirical P values use `(extreme + 1)/(10,000 + 1)`.
- BH-FDR is calculated across all 100 cells separately within each evidence–model analysis.
- TVD is `0.5 * sum(abs(p_positive - p_reference))` using unsmoothed proportions.

## Run analysis

```text
python Joint_confidence_enrichment.py --input path/to/circIG.tsv --output results
```

Validate population definitions without running resampling:

```text
python Joint_confidence_enrichment.py --input path/to/circIG.tsv --output validation --validate-only
```

Run one evidence/model combination for testing:

```text
python Joint_confidence_enrichment.py --input path/to/circIG.tsv --output test --evidence RT-independent --model XGBoost
```

## Draw figures

```text
python plot_Joint_confidence_landscape.py --results results --output figures
```

Each run records software versions, population counts, random seeds and analysis settings.

## Requirements:
numpy==2.3.5
pandas==3.0.1
matplotlib>=3.9