#!/usr/bin/env python3
"""Clean BloodCircR CircRNA Atlas into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "CircRNA_Atlas"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "11_BloodCircR.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "11_BloodCircR.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "11_BloodCircR.circRNAs.clean.summary.tsv"

REQUIRED_COLUMNS = ["BSJ_ID", "chr", "circRNA_start", "circRNA_end", "strand", "host_gene"]
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
    value = normalize_missing(value)
    if value == "NA":
        return value
    if value.lower().startswith("chr"):
        return "chr" + value[3:]
    return "chr" + value


def normalize_strand(value: str | None) -> str:
    value = normalize_missing(value)
    return value if value in {"+", "-"} else "NA"


def parse_coord(value: str | None, column: str) -> str:
    value = normalize_missing(value)
    if value == "NA":
        raise ValueError(f"missing {column}")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"non-integer {column}: {value!r}") from exc
    if parsed < 0:
        raise ValueError(f"negative {column}: {value!r}")
    return str(parsed)


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
    parser = argparse.ArgumentParser(description="Clean BloodCircR CircRNA Atlas into unified six-column TSV.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw atlas. Default: {DEFAULT_INPUT}")
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen_rows: set[tuple[str, str, str, str, str, str]] = set()
    unique_ids: set[str] = set()

    with args.input.open(newline="", errors="replace") as input_handle, args.output.open("w", newline="") as output_handle, args.invalid_output.open("w", newline="") as invalid_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        validate_header(args.input, reader.fieldnames)
        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_columns = list(reader.fieldnames or []) + ["failure_reason"]
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=invalid_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        invalid_writer.writeheader()

        for row in reader:
            counts["raw_rows"] += 1
            try:
                chrom = normalize_chrom(row.get("chr"))
                start = parse_coord(row.get("circRNA_start"), "circRNA_start")
                end = parse_coord(row.get("circRNA_end"), "circRNA_end")
                if int(start) > int(end):
                    raise ValueError(f"start is greater than end: {chrom}:{start}-{end}")
                out = {
                    "circRNA_id": normalize_missing(row.get("BSJ_ID")),
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": normalize_strand(row.get("strand")),
                    "gene_symbol": normalize_missing(row.get("host_gene")),
                }
                if out["circRNA_id"] == "NA" or out["chrom"] == "NA":
                    raise ValueError("missing circRNA_id or chrom")
            except ValueError as exc:
                counts["invalid_rows"] += 1
                invalid_row = dict(row)
                invalid_row["failure_reason"] = str(exc)
                invalid_writer.writerow(invalid_row)
                continue

            counts["clean_rows_before_uniq"] += 1
            key = tuple(out[column] for column in OUTPUT_COLUMNS)
            if key in seen_rows:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen_rows.add(key)
            writer.writerow(out)
            unique_ids.add(out["circRNA_id"])
            counts["clean_rows"] += 1
            if out["strand"] in {"+", "-"}:
                counts["stranded_rows"] += 1
            else:
                counts["unstranded_rows"] += 1

    counts["unique_circRNA_ids"] = len(unique_ids)
    write_summary(args.summary_output, args, counts)

    print(f"wrote {args.output}")
    print(f"clean_rows={counts['clean_rows']}")
    print(f"invalid_rows={counts['invalid_rows']}")
    print(f"duplicate_full_rows_removed={counts['duplicate_full_rows_removed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
