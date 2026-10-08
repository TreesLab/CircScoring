#!/usr/bin/env python3
"""Clean CIRIonco human circRNA raw data into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW = BASE_DIR / "raw" / "human_circrnas_api_merged.tsv"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "37_CIRIonco.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "37_CIRIonco.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "37_CIRIonco.circRNAs.clean.summary.tsv"

RAW_COLUMNS = [
    "CIRC_ID",
    "regulation_status",
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Gene_Name",
    "CIRC_type",
    "Tissue_number",
    "Tissue_cancer_number",
    "Tissue_BSJ_number",
    "Tissue_cancer_BSJ_number",
    "Tissue_names",
    "cancer_tissue_names",
    "start_exon",
    "end_exon",
    "start_exon_boundary",
    "end_exon_boundary",
    "exon_composition",
    "system",
    "tissue",
    "disease",
    "system_regulation",
    "tissue_regulation",
    "disease_regulation",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["line_number", *RAW_COLUMNS, "invalid_reason"]
SUMMARY_KEYS = [
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "unique_gene_values",
    "missing_gene_values",
    "non_chr_chromosome_rows",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CIRIonco circRNA API TSV into unified six-column TSV.")
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW, help=f"Raw merged TSV. Default: {DEFAULT_RAW}")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"Clean output TSV. Default: {DEFAULT_OUTPUT}")
    parser.add_argument(
        "--invalid-output",
        type=Path,
        default=DEFAULT_INVALID_OUTPUT,
        help=f"Invalid/audit output TSV. Default: {DEFAULT_INVALID_OUTPUT}",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Summary output TSV. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def parse_int(value: str, field: str) -> str:
    normalized = normalize_missing(value)
    if normalized == "NA":
        raise ValueError(f"missing {field}")
    try:
        parsed = int(normalized)
    except ValueError as exc:
        raise ValueError(f"invalid integer {field}: {value!r}") from exc
    if parsed < 0:
        raise ValueError(f"negative {field}: {parsed}")
    return str(parsed)


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def write_summary(path: Path, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in SUMMARY_KEYS:
            writer.writerow([key, counts[key]])


def main() -> int:
    args = parse_args()
    raw_path = args.raw.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([output, invalid_output, summary_output], args.force)

    if not raw_path.exists():
        raise FileNotFoundError(f"Missing raw merged TSV: {raw_path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    gene_values: set[str] = set()

    with (
        raw_path.open(newline="") as raw_handle,
        output.open("w", newline="") as output_handle,
        invalid_output.open("w", newline="") as invalid_handle,
    ):
        reader = csv.DictReader(raw_handle, delimiter="\t")
        missing_columns = [column for column in RAW_COLUMNS if column not in (reader.fieldnames or [])]
        if missing_columns:
            raise ValueError(f"Missing raw column(s): {', '.join(missing_columns)}")

        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for line_number, raw in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                circRNA_id = normalize_missing(raw.get("CIRC_ID"))
                chrom = normalize_missing(raw.get("Chromosome"))
                start = parse_int(raw.get("Start", ""), "Start")
                end = parse_int(raw.get("End", ""), "End")
                strand = normalize_missing(raw.get("Strand"))
                gene_symbol = normalize_missing(raw.get("Gene_Name"))

                if circRNA_id == "NA":
                    raise ValueError("missing CIRC_ID")
                if chrom == "NA":
                    raise ValueError("missing Chromosome")
                if int(start) > int(end):
                    raise ValueError(f"Start is greater than End: {start}>{end}")
                if strand not in {"+", "-"}:
                    raise ValueError(f"invalid strand: {strand!r}")

                row = {
                    "circRNA_id": circRNA_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_symbol": gene_symbol,
                }
            except ValueError as exc:
                invalid = {"line_number": str(line_number), "invalid_reason": str(exc)}
                invalid.update({column: raw.get(column, "") for column in RAW_COLUMNS})
                invalid_writer.writerow(invalid)
                counts["invalid_rows"] += 1
                continue

            counts["clean_rows_before_uniq"] += 1
            key = tuple(row[column] for column in OUTPUT_COLUMNS)
            if key in seen_rows:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen_rows.add(key)
            writer.writerow(row)

            circ_ids.add(row["circRNA_id"])
            if row["gene_symbol"] == "NA":
                counts["missing_gene_values"] += 1
            else:
                gene_values.add(row["gene_symbol"])
            if not row["chrom"].startswith("chr"):
                counts["non_chr_chromosome_rows"] += 1
            counts["stranded_rows"] += 1

    counts["clean_rows"] = len(seen_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_gene_values"] = len(gene_values)
    write_summary(summary_output, counts)

    print(f"raw_rows\t{counts['raw_rows']}")
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    print(f"duplicate_full_rows_removed\t{counts['duplicate_full_rows_removed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
