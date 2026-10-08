#!/usr/bin/env python3
"""Clean TCCIA Ensemble circRNA annotations into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_FILE = BASE_DIR / "raw" / "TCCIA_Ensemble_circRNAs.tsv.gz"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "27_TCCIA.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "27_TCCIA.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "27_TCCIA.circRNAs.clean.summary.tsv"

EXPECTED_HEADER = ["ID", "Gene", "Mean", "Median", "Max", "SD", "N", "Cohort"]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = EXPECTED_HEADER + ["line_number", "invalid_reason"]
SUMMARY_KEYS = [
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "unique_genes",
    "unique_cohorts",
]

ID_RE = re.compile(r"^(?P<gene>.+):(?P<strand>[+-]):(?P<chrom>chr[^:]+):(?P<start>\d+):(?P<end>\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean TCCIA Ensemble circRNA table into unified six-column TSV.")
    parser.add_argument("--raw-file", type=Path, default=DEFAULT_RAW_FILE, help=f"Raw input file. Default: {DEFAULT_RAW_FILE}")
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


def parse_id(raw_id: str) -> tuple[str, str, str, str, str]:
    match = ID_RE.match(raw_id.strip())
    if not match:
        raise ValueError(f"invalid ID format: {raw_id!r}")
    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"start is greater than end: {raw_id!r}")
    return match.group("gene"), match.group("strand"), chrom, str(start), str(end)


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
    raw_file = args.raw_file.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([output, invalid_output, summary_output], args.force)

    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    genes: set[str] = set()
    cohorts: set[str] = set()

    with (
        gzip.open(raw_file, "rt", newline="") as input_handle,
        output.open("w", newline="") as output_handle,
        invalid_output.open("w", newline="") as invalid_handle,
    ):
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_HEADER:
            raise ValueError(f"Unexpected header in {raw_file}: {reader.fieldnames!r}")

        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for line_number, raw in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            cohorts.add(normalize_missing(raw.get("Cohort")))
            try:
                raw_id = normalize_missing(raw.get("ID"))
                if raw_id == "NA":
                    raise ValueError("missing ID")
                id_gene, strand, chrom, start, end = parse_id(raw_id)
                gene_symbol = normalize_missing(raw.get("Gene"))
                if gene_symbol == "NA":
                    gene_symbol = normalize_missing(id_gene)
                row = {
                    "circRNA_id": raw_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_symbol": gene_symbol,
                }
            except ValueError as exc:
                invalid = {column: raw.get(column, "") for column in EXPECTED_HEADER}
                invalid["line_number"] = str(line_number)
                invalid["invalid_reason"] = str(exc)
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
            genes.add(row["gene_symbol"])
            if row["strand"] in {"+", "-"}:
                counts["stranded_rows"] += 1
            else:
                counts["unstranded_rows"] += 1

    counts["clean_rows"] = len(seen_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_genes"] = len(genes)
    counts["unique_cohorts"] = len(cohorts)
    write_summary(summary_output, counts)

    print(f"raw_rows\t{counts['raw_rows']}")
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
