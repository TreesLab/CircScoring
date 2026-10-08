#!/usr/bin/env python3
"""Clean ExoRBase circRNA annotation into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "circRNAs_anno.csv"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "10_ExoRBase.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "10_ExoRBase.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "10_ExoRBase.circRNAs.clean.summary.tsv"

REQUIRED_COLUMNS = ["circID", "Genomic position", "Strand", "Gene symbol"]
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

POSITION_RE = re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def normalize_strand(value: str | None) -> str:
    value = normalize_missing(value)
    return value if value in {"+", "-"} else "NA"


def parse_position(position: str) -> tuple[str, str, str]:
    match = POSITION_RE.match(position.strip())
    if not match:
        raise ValueError(f"invalid Genomic position format: {position!r}")
    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"start is greater than end: {chrom}:{start}-{end}")
    return chrom, str(start), str(end)


def validate_header(path: Path, fieldnames: list[str] | None) -> None:
    if fieldnames is None:
        raise ValueError(f"Input table has no header: {path}")
    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")


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
    parser = argparse.ArgumentParser(description="Clean ExoRBase circRNA annotation into unified six-column TSV.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw annotation CSV. Default: {DEFAULT_INPUT}")
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

    with args.input.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle)
        validate_header(args.input, reader.fieldnames)

        for line_number, raw_row in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                circ_id = normalize_missing(raw_row.get("circID"))
                if circ_id == "NA":
                    raise ValueError("missing circID")
                chrom, start, end = parse_position(raw_row.get("Genomic position", ""))
                clean_row = {
                    "circRNA_id": circ_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": normalize_strand(raw_row.get("Strand")),
                    "gene_symbol": normalize_missing(raw_row.get("Gene symbol")),
                }
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
