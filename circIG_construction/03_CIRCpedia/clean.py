#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "Homo-sapiens.txt"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "03_CIRCpedia.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "03_CIRCpedia.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "03_CIRCpedia.circRNAs.clean.summary.tsv"
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
LOCATION_RE = re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)$")


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.lower() in {"none", "null", "na", "nan"}:
        return "NA"
    return value


def parse_location(location: str) -> tuple[str, str, str, str]:
    match = LOCATION_RE.match(location.strip())
    if not match:
        raise ValueError(f"invalid Location format: {location!r}")

    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    strand = match.group("strand")
    if start >= end:
        raise ValueError(f"invalid interval: {chrom}:{start}-{end}")
    return chrom, str(start), str(end), strand


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in [
            "raw_rows",
            "clean_rows_before_uniq",
            "duplicate_full_rows_removed",
            "clean_rows",
            "invalid_rows",
            "unique_circRNA_ids",
            "stranded_rows",
            "unstranded_rows",
        ]:
            writer.writerow([key, counts[key]])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CIRCpedia raw Homo sapiens circRNA annotation into unified columns.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw Homo-sapiens.txt. Default: {DEFAULT_INPUT}")
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
    for output in [args.output, args.invalid_output, args.summary_output]:
        if output.exists() and not args.force:
            raise FileExistsError(f"Output exists; use --force to overwrite: {output}")

    counts: Counter[str] = Counter()
    rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {args.input}")

        for line_number, row in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                chrom, start, end, strand = parse_location(row.get("Location", ""))
                clean_row = {
                    "circRNA_id": normalize_missing(row.get("circID")),
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_symbol": normalize_missing(row.get("gene_Refseq")),
                }
                if clean_row["circRNA_id"] == "NA":
                    raise ValueError("missing circID")
            except ValueError as exc:
                invalid = dict(row)
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

    invalid_columns = []
    if invalid_rows:
        invalid_columns = list(invalid_rows[0].keys())

    write_rows(args.output, rows, OUTPUT_COLUMNS)
    write_rows(args.invalid_output, invalid_rows, invalid_columns or ["line_number", "invalid_reason"])
    write_summary(args.summary_output, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
