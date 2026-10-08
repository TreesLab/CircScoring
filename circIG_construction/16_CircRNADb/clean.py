#!/usr/bin/env python3
"""Clean CircRNADb circRNA dataset into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "circRNA_dataset" / "circRNA_dataset.txt"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "16_CircRNADb.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "16_CircRNADb.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "16_CircRNADb.circRNAs.clean.summary.tsv"

RAW_COLUMNS = [
    "circ_id",
    "chrom",
    "start",
    "end",
    "strand",
    "gene_symbol",
    "genomic_length",
    "best_transcript",
    "spliced_length",
    "exon_count",
    "exon_lengths",
    "exon_offsets",
    "transcript_exon_annotation",
    "samples",
    "pmids",
]
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


def normalize_chrom(value: str | None) -> str:
    chrom = normalize_missing(value)
    if chrom == "NA":
        return chrom
    if chrom.lower().startswith("chr"):
        return "chr" + chrom[3:]
    return f"chr{chrom}"


def normalize_strand(value: str | None) -> str:
    strand = normalize_missing(value)
    return strand if strand in {"+", "-"} else "NA"


def parse_coordinate(value: str | None, column: str) -> str:
    value = normalize_missing(value)
    if value == "NA":
        raise ValueError(f"missing {column}")
    try:
        coord = int(value)
    except ValueError as exc:
        raise ValueError(f"{column} is not an integer: {value!r}") from exc
    if coord < 0:
        raise ValueError(f"{column} is negative: {coord}")
    return str(coord)


def parse_raw_row(fields: list[str], line_number: int) -> dict[str, str]:
    if len(fields) != len(RAW_COLUMNS):
        raise ValueError(f"expected {len(RAW_COLUMNS)} fields, found {len(fields)}")

    raw = dict(zip(RAW_COLUMNS, fields))
    circ_id = normalize_missing(raw["circ_id"])
    if circ_id == "NA":
        raise ValueError("missing circ_id")
    chrom = normalize_chrom(raw["chrom"])
    if chrom == "NA":
        raise ValueError("missing chrom")
    start = parse_coordinate(raw["start"], "start")
    end = parse_coordinate(raw["end"], "end")
    if int(start) > int(end):
        raise ValueError(f"start is greater than end: {chrom}:{start}-{end}")
    strand = normalize_strand(raw["strand"])

    return {
        "circRNA_id": circ_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "gene_symbol": normalize_missing(raw["gene_symbol"]),
    }


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {key: str(counts[key]) for key in SUMMARY_KEYS}
    row["input_file"] = str(args.input)
    row["output_file"] = str(args.output)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_KEYS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CircRNADb raw coordinate table into unified six-column TSV.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw coordinate table. Default: {DEFAULT_INPUT}")
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
    args.input = args.input.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.invalid_output = args.invalid_output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([args.output, args.invalid_output, args.summary_output], args.force)

    rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    counts: Counter[str] = Counter()

    with args.input.open(errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, start=1):
            if not fields:
                continue
            counts["raw_rows"] += 1
            try:
                clean_row = parse_raw_row(fields, line_number)
            except ValueError as exc:
                invalid = {f"raw_col_{index + 1}": value for index, value in enumerate(fields)}
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
