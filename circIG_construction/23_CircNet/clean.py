#!/usr/bin/env python3
"""Clean CircNet raw cancer CSV files into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "23_CircNet.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "23_CircNet.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "23_CircNet.circRNAs.clean.summary.tsv"

CIRCID_RE = re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)\|(?P<end>\d+)$")
EXPECTED_HEADER = ["", "CircID", "Strand", "CircType", "HostGene", "Algorithm", "Sequence"]
EXCLUDED_NAMES = {"samples.csv", "sample_0913.csv", "mir_update_0914.csv"}
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["source_file", "line_number", "invalid_reason", "CircID", "Strand", "HostGene"]
SUMMARY_COLUMNS = [
    "raw_dir",
    "output_file",
    "input_file_count",
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


def normalize_hostgene(value: str | None) -> str:
    value = normalize_missing(value)
    if value == "NA":
        return value
    parts = [part.strip() for part in re.split(r"[;|]+", value) if part.strip()]
    return ";".join(parts) if parts else "NA"


def parse_circid(value: str | None) -> tuple[str, str, str]:
    circid = normalize_missing(value)
    if circid == "NA":
        raise ValueError("missing CircID")
    match = CIRCID_RE.match(circid)
    if not match:
        raise ValueError(f"invalid CircID format: {circid!r}")
    chrom = "chr" + match.group("chrom")[3:]
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"start is greater than end: {circid!r}")
    return chrom, str(start), str(end)


def normalize_strand(value: str | None) -> str:
    strand = normalize_missing(value)
    if strand not in {"+", "-"}:
        raise ValueError(f"invalid strand: {strand!r}")
    return strand


def iter_input_files(raw_dir: Path) -> list[Path]:
    return [
        path
        for path in sorted(raw_dir.glob("*.csv"))
        if path.name.lower() not in EXCLUDED_NAMES
    ]


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {key: str(counts[key]) for key in SUMMARY_COLUMNS}
    row["raw_dir"] = str(args.raw_dir)
    row["output_file"] = str(args.output)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CircNet raw CSV files into unified six-column TSV.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw input directory. Default: {DEFAULT_RAW_DIR}")
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
    args.raw_dir = args.raw_dir.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.invalid_output = args.invalid_output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([args.output, args.invalid_output, args.summary_output], args.force)

    input_files = iter_input_files(args.raw_dir)
    if not input_files:
        raise FileNotFoundError(f"No CircNet cancer CSV files found in {args.raw_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter(input_file_count=len(input_files))
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()

    with args.output.open("w", newline="") as output_handle, args.invalid_output.open("w", newline="") as invalid_handle:
        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        invalid_writer.writeheader()

        for path in input_files:
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames != EXPECTED_HEADER:
                    raise ValueError(f"Unexpected header in {path}: {reader.fieldnames!r}")

                for line_number, row in enumerate(reader, start=2):
                    counts["raw_rows"] += 1
                    try:
                        circid = normalize_missing(row.get("CircID"))
                        chrom, start, end = parse_circid(circid)
                        strand = normalize_strand(row.get("Strand"))
                        clean_row = {
                            "circRNA_id": circid,
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "gene_symbol": normalize_hostgene(row.get("HostGene")),
                        }
                    except ValueError as exc:
                        invalid_writer.writerow(
                            {
                                "source_file": str(path.relative_to(args.raw_dir.parent)),
                                "line_number": str(line_number),
                                "invalid_reason": str(exc),
                                "CircID": normalize_missing(row.get("CircID")),
                                "Strand": normalize_missing(row.get("Strand")),
                                "HostGene": normalize_missing(row.get("HostGene")),
                            }
                        )
                        counts["invalid_rows"] += 1
                        continue

                    counts["clean_rows_before_uniq"] += 1
                    key = tuple(clean_row[column] for column in OUTPUT_COLUMNS)
                    if key in seen_rows:
                        counts["duplicate_full_rows_removed"] += 1
                        continue
                    seen_rows.add(key)
                    writer.writerow(clean_row)
                    circ_ids.add(clean_row["circRNA_id"])
                    if clean_row["strand"] in {"+", "-"}:
                        counts["stranded_rows"] += 1
                    else:
                        counts["unstranded_rows"] += 1

    counts["clean_rows"] = len(seen_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    write_summary(args.summary_output, args, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
