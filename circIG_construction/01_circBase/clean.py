#!/usr/bin/env python3
"""Clean circBase human circRNA raw tables into the unified pipeline format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "raw"
DEFAULT_OUTPUT = SCRIPT_DIR / "clean" / "01_circBase.circRNAs.clean.tsv"
DEFAULT_SUMMARY_OUTPUT = SCRIPT_DIR / "clean" / "01_circBase.circRNAs.clean.summary.tsv"

RAW_FILES = [
    "hsa_hg19_circRNA.txt",
    "hsa_hg19_Rybak2015.txt",
    "hsa_hg19_Jeck2013.txt",
    "hsa_hg19_Memczak2013.txt",
    "hsa_hg19_Salzman2013.txt",
    "hsa_hg19_Zhang2013_H9.txt",
]

REQUIRED_COLUMNS = ["# chrom", "start", "end", "strand", "circRNA ID", "gene symbol"]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
SUMMARY_COLUMNS = [
    "raw_dir",
    "output",
    "raw_files",
    "raw_data_rows",
    "duplicate_full_raw_rows",
    "unique_full_raw_rows",
    "duplicate_clean_rows",
    "invalid_rows",
    "clean_output_rows",
    "unique_circRNA_ids",
    "unique_coord_strand_id_keys",
]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge and clean circBase Homo sapiens hg19 raw TXT tables into "
            "the unified six-column circRNA format."
        )
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Directory containing circBase raw TXT files. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Clean TSV output path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Summary TSV output path. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE"}:
        return "NA"
    return value


def validate_raw_files(raw_dir: Path) -> list[Path]:
    paths = [raw_dir / filename for filename in RAW_FILES]
    missing = [path for path in paths if not path.exists() or path.stat().st_size == 0]
    if missing:
        names = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing or empty raw file(s): {names}")
    return paths


def validate_header(path: Path, fieldnames: list[str] | None) -> list[str]:
    if fieldnames is None:
        raise ValueError(f"Raw table has no header: {path}")
    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")
    return fieldnames


def raw_row_key(fieldnames: list[str], row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(column, "") for column in fieldnames)


def is_valid_row(row: dict[str, str]) -> bool:
    chrom = normalize_missing(row.get("# chrom"))
    start = normalize_missing(row.get("start"))
    end = normalize_missing(row.get("end"))
    strand = normalize_missing(row.get("strand"))
    circRNA_id = normalize_missing(row.get("circRNA ID"))

    if chrom == "NA" or circRNA_id == "NA":
        return False
    if strand not in {"+", "-"}:
        return False
    try:
        int(start)
        int(end)
    except ValueError:
        return False
    return True


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {
        "circRNA_id": normalize_missing(row.get("circRNA ID")),
        "chrom": normalize_missing(row.get("# chrom")),
        "start": normalize_missing(row.get("start")),
        "end": normalize_missing(row.get("end")),
        "strand": normalize_missing(row.get("strand")),
        "gene_symbol": normalize_missing(row.get("gene symbol")),
    }


def read_clean_rows(raw_paths: list[Path]) -> tuple[list[dict[str, str]], Counter[str]]:
    counts: Counter[str] = Counter(raw_files=len(raw_paths))
    seen_full_rows: set[tuple[str, ...]] = set()
    seen_clean_rows: set[tuple[str, ...]] = set()
    clean_rows: list[dict[str, str]] = []

    for path in raw_paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = validate_header(path, reader.fieldnames)
            for row in reader:
                counts["raw_data_rows"] += 1
                key = raw_row_key(fieldnames, row)
                if key in seen_full_rows:
                    counts["duplicate_full_raw_rows"] += 1
                    continue
                seen_full_rows.add(key)
                counts["unique_full_raw_rows"] += 1

                if not is_valid_row(row):
                    counts["invalid_rows"] += 1
                    continue
                cleaned = clean_row(row)
                clean_key = tuple(cleaned[column] for column in OUTPUT_COLUMNS)
                if clean_key in seen_clean_rows:
                    counts["duplicate_clean_rows"] += 1
                    continue
                seen_clean_rows.add(clean_key)
                clean_rows.append(cleaned)

    counts["clean_output_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len({row["circRNA_id"] for row in clean_rows})
    counts["unique_coord_strand_id_keys"] = len(
        {
            (row["chrom"], row["start"], row["end"], row["strand"], row["circRNA_id"])
            for row in clean_rows
        }
    )
    return clean_rows, counts


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "raw_dir": str(args.raw_dir),
        "output": str(args.output),
        "raw_files": str(counts["raw_files"]),
        "raw_data_rows": str(counts["raw_data_rows"]),
        "duplicate_full_raw_rows": str(counts["duplicate_full_raw_rows"]),
        "unique_full_raw_rows": str(counts["unique_full_raw_rows"]),
        "duplicate_clean_rows": str(counts["duplicate_clean_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "clean_output_rows": str(counts["clean_output_rows"]),
        "unique_circRNA_ids": str(counts["unique_circRNA_ids"]),
        "unique_coord_strand_id_keys": str(counts["unique_coord_strand_id_keys"]),
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def refuse_overwrite(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Output already exists; use --force to overwrite: {path}")


def main() -> int:
    args = parse_args()
    args.raw_dir = args.raw_dir.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()

    try:
        refuse_overwrite(args.output, args.force)
        refuse_overwrite(args.summary_output, args.force)
        raw_paths = validate_raw_files(args.raw_dir)
        clean_rows, counts = read_clean_rows(raw_paths)
        write_table(args.output, clean_rows)
        write_summary(args.summary_output, args, counts)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote clean table: {args.output}")
    print(f"Wrote summary: {args.summary_output}")
    print(
        "Summary: "
        f"{counts['raw_data_rows']} raw rows, "
        f"{counts['clean_output_rows']} clean rows, "
        f"{counts['invalid_rows']} invalid rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
