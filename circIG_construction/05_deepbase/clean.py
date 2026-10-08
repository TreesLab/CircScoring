#!/usr/bin/env python3
"""Clean deepBase human circRNA browser table into unified pipeline columns."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "deepbase_human_circRNA_browser_api.tsv"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "05_deepbase.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "05_deepbase.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "05_deepbase.circRNAs.clean.summary.tsv"
REQUIRED_COLUMNS = ["Name", "Chromosome", "Start", "End", "Strand", "Source Gene Symbol"]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
SUMMARY_KEYS = [
    "input_file",
    "output_file",
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
]


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def parse_coordinate(row: dict[str, str]) -> dict[str, str]:
    circ_id = normalize_missing(row.get("Name"))
    chrom = normalize_missing(row.get("Chromosome"))
    start_text = normalize_missing(row.get("Start"))
    end_text = normalize_missing(row.get("End"))
    strand = normalize_missing(row.get("Strand"))
    gene_symbol = normalize_missing(row.get("Source Gene Symbol"))

    if circ_id == "NA":
        raise ValueError("missing Name")
    if not re.fullmatch(r"chr[\w.]+", chrom):
        raise ValueError(f"invalid Chromosome: {chrom!r}")
    if not start_text.isdigit() or not end_text.isdigit():
        raise ValueError(f"invalid coordinates: {start_text}-{end_text}")
    start = int(start_text)
    end = int(end_text)
    if start > end:
        raise ValueError(f"start is greater than end: {start}>{end}")
    if strand not in {"+", "-"}:
        raise ValueError(f"invalid Strand: {strand!r}")

    return {
        "circRNA_id": circ_id,
        "chrom": chrom,
        "start": str(start),
        "end": str(end),
        "strand": strand,
        "gene_symbol": gene_symbol,
    }


def validate_header(path: Path, fieldnames: list[str] | None) -> None:
    if fieldnames is None:
        raise ValueError(f"Input table has no header: {path}")
    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"Output exists; use --force to overwrite: {', '.join(existing)}")


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "input_file": str(args.input),
        "output_file": str(args.output),
        "raw_rows": str(counts["raw_rows"]),
        "clean_rows_before_uniq": str(counts["clean_rows_before_uniq"]),
        "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
        "clean_rows": str(counts["clean_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "unique_circRNA_ids": str(counts["unique_circRNA_ids"]),
        "stranded_rows": str(counts["stranded_rows"]),
        "unstranded_rows": str(counts["unstranded_rows"]),
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_KEYS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw browser table TSV. Default: {DEFAULT_INPUT}")
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


def main() -> int:
    args = parse_args()
    refuse_overwrite([args.output, args.invalid_output, args.summary_output], args.force)

    rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    counts: Counter[str] = Counter()

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(args.input, reader.fieldnames)
        for line_number, raw_row in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                clean_row = parse_coordinate(raw_row)
            except ValueError as exc:
                invalid = dict(raw_row)
                invalid["line_number"] = str(line_number)
                invalid["invalid_reason"] = str(exc)
                invalid_rows.append(invalid)
                counts["invalid_rows"] += 1
                continue

            counts["clean_rows_before_uniq"] += 1
            key = tuple(clean_row[column] for column in OUTPUT_COLUMNS)
            if key in seen_rows:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen_rows.add(key)
            rows.append(clean_row)
            circ_ids.add(clean_row["circRNA_id"])
            if clean_row["strand"] in {"+", "-"}:
                counts["stranded_rows"] += 1
            else:
                counts["unstranded_rows"] += 1

    counts["clean_rows"] = len(rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    invalid_columns = list(invalid_rows[0].keys()) if invalid_rows else ["line_number", "invalid_reason"]
    write_rows(args.output, rows, OUTPUT_COLUMNS)
    write_rows(args.invalid_output, invalid_rows, invalid_columns)
    write_summary(args.summary_output, args, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
