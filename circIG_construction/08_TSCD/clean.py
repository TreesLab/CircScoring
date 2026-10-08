#!/usr/bin/env python3
"""Clean TSCD human hg38 circRNA raw tables into the unified pipeline format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_ADULT = BASE_DIR / "raw" / "hg38_adult_TS_circRNAs.txt"
DEFAULT_FETAL = BASE_DIR / "raw" / "hg38_fetal_TS_circRNAs.txt"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "08_TSCD.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "08_TSCD.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "08_TSCD.circRNAs.clean.summary.tsv"

RAW_COLUMNS = [
    "tissue_or_sample",
    "coordinates",
    "chromosome",
    "donor_site",
    "acceptor_site",
    "junction_reads",
    "strand",
    "algorithms",
    "srptm",
    "gene_annotation",
    "gene_body",
    "gene_type",
    "gene_strand",
    "microRNAs",
    "RBPs",
    "flanking_length",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
SUMMARY_COLUMNS = [
    "adult_input",
    "fetal_input",
    "output_file",
    "raw_rows",
    "adult_rows",
    "fetal_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "raw_plus_rows",
    "raw_minus_rows",
    "raw_mixed_rows",
    "raw_other_strand_rows",
]
NULL_VALUES = {"", "-", ".", "n/a", "N/A", "NA", "None", "none", "NULL", "null", "Unknown"}


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if value in NULL_VALUES:
        return "NA"
    return value


def normalize_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return value if value in {"+", "-"} else "NA"


def normalize_chrom(value: str | None) -> str:
    value = normalize_missing(value)
    if value == "NA":
        return value
    if value.lower().startswith("chr"):
        return "chr" + value[3:]
    return f"chr{value}"


def parse_int(value: str | None, column: str) -> int:
    value = normalize_missing(value)
    if value == "NA":
        raise ValueError(f"missing {column}")
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"non-integer {column}: {value!r}") from exc


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def iter_raw_rows(path: Path, stage: str):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) != len(RAW_COLUMNS):
                raise ValueError(f"{path}:{line_number}: expected {len(RAW_COLUMNS)} columns, got {len(row)}")
            record = dict(zip(RAW_COLUMNS, row))
            record["developmental_stage"] = stage
            record["source_file"] = str(path)
            record["line_number"] = str(line_number)
            yield record


def raw_strand_bucket(value: str | None) -> str:
    value = (value or "").strip()
    if value == "+":
        return "raw_plus_rows"
    if value == "-":
        return "raw_minus_rows"
    if "," in value:
        return "raw_mixed_rows"
    return "raw_other_strand_rows"


def to_clean_row(record: dict[str, str]) -> dict[str, str]:
    chrom = normalize_chrom(record.get("chromosome"))
    if chrom == "NA":
        raise ValueError("missing chromosome")
    start = parse_int(record.get("donor_site"), "donor_site")
    end = parse_int(record.get("acceptor_site"), "acceptor_site")
    if start > end:
        start, end = end, start
    if start < 1:
        raise ValueError(f"start must be >= 1: {start}")
    coordinates = normalize_missing(record.get("coordinates"))
    if coordinates == "NA":
        raise ValueError("missing coordinates")
    return {
        "circRNA_id": coordinates,
        "chrom": chrom,
        "start": str(start),
        "end": str(end),
        "strand": normalize_strand(record.get("strand")),
        "gene_symbol": normalize_missing(record.get("gene_annotation")),
    }


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "adult_input": str(args.adult),
        "fetal_input": str(args.fetal),
        "output_file": str(args.output),
    }
    for column in SUMMARY_COLUMNS[3:]:
        row[column] = str(counts[column])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean TSCD human hg38 adult/fetal circRNA tables into unified six-column TSV.")
    parser.add_argument("--adult", type=Path, default=DEFAULT_ADULT, help=f"Adult raw input. Default: {DEFAULT_ADULT}")
    parser.add_argument("--fetal", type=Path, default=DEFAULT_FETAL, help=f"Fetal raw input. Default: {DEFAULT_FETAL}")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"Clean output TSV. Default: {DEFAULT_OUTPUT}")
    parser.add_argument("--invalid-output", type=Path, default=DEFAULT_INVALID_OUTPUT, help=f"Invalid/audit output TSV. Default: {DEFAULT_INVALID_OUTPUT}")
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT, help=f"Summary output TSV. Default: {DEFAULT_SUMMARY_OUTPUT}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.adult = args.adult.expanduser().resolve()
    args.fetal = args.fetal.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.invalid_output = args.invalid_output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([args.output, args.invalid_output, args.summary_output], args.force)

    rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    counts: Counter[str] = Counter()

    for path, stage in ((args.adult, "adult"), (args.fetal, "fetal")):
        for record in iter_raw_rows(path, stage):
            counts["raw_rows"] += 1
            counts[f"{stage}_rows"] += 1
            counts[raw_strand_bucket(record.get("strand"))] += 1
            try:
                clean_row = to_clean_row(record)
            except ValueError as exc:
                invalid = dict(record)
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

    invalid_columns = list(invalid_rows[0].keys()) if invalid_rows else ["source_file", "line_number", "invalid_reason"]
    write_rows(args.output, rows, OUTPUT_COLUMNS)
    write_rows(args.invalid_output, invalid_rows, invalid_columns)
    write_summary(args.summary_output, args, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
