# circRNA integration Snakemake workflow

This workflow starts from the clean tables produced for each database. It builds
the human and mouse tier 1 through tier 4 base tables, then generates the
mouse-to-human conservation table. Raw-data downloads and database-specific
cleaning remain the responsibility of `download.py` and `clean.py` in each
database directory.

## Main targets

```bash
snakemake -s workflow/Snakefile --cores 2 all
snakemake -s workflow/Snakefile --cores 2 human_v4
snakemake -s workflow/Snakefile --cores 2 mouse_v4
snakemake -s workflow/Snakefile --cores 2 mouse_conservation
snakemake -s workflow/Snakefile --cores 2 coordinate_qc
```

Production results are written to `results/` under the repository root by
default. The coordinate-system manifests are locked decisions for production
builds. The `coordinate_qc` target only regenerates QC reports and does not
modify those decisions.

The clean stage for `21_CircTarget` and `43_circRNADisease` produces circRNA IDs
only. The workflow resolves their coordinates using the tier 1 circBase and
circAtlas final tables, then generates six-column clean tables. Before running
the workflow, run the latest `clean.py` for both databases.
`43_circRNADisease/clean.py` requires `openpyxl`.

## Locked coordinate assignments

The mouse data from `30_CircR2Disease` and `06_CircFunBase` contain mixed
coordinate systems, so production conversion requires row-level coordinate
classification. The workflow stores the verified classification snapshots in
`workflow/config/locked_coordinate_assignments/` so production builds do not
depend on untracked legacy QC outputs. These snapshots are identical to the
following original files:

- `30_CircR2Disease/mouse/coordinate_check/30_CircR2Disease.coordinate_check.tsv`
- `06_CircFunBase/mouse/coordinate_check/06_CircFunBase.coordinate_check.tsv`
