#!/usr/bin/env python3
"""Clean circExp web annotation tables into per-GSE unified circRNA tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = BASE_DIR / "raw" / "web_annotation_tables"
DEFAULT_OUTPUT_DIR = BASE_DIR / "clean" / "web_annotation_tables"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "28_circExp.circRNAs.clean.summary.tsv"
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
ASCRP_COLUMN = "circExp_ASCRP"
ASCRP_OUTPUT_COLUMNS = [*OUTPUT_COLUMNS, ASCRP_COLUMN]
RAW_REQUIRED_COLUMNS = [
    "source_file",
    "geo_dataset",
    "original_id",
    "circbase_id",
    "genomic_location",
    "strand",
    "feature",
    "parental_gene",
]
SUMMARY_COLUMNS = [
    "geo_dataset",
    "input_file",
    "output_file",
    "ascrp_output_file",
    "invalid_output_file",
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "coordinate_from_genomic_location",
    "coordinate_from_original_id",
    "rows_with_ascrp",
    "unique_ascrp_ids",
]
NULL_VALUES = {"", "-", ".", "*", "NA", "N/A", "na", "n/a", "NULL", "null", "None", "none"}
COORD_RE = re.compile(r"^(chr[^:\s]+):(\d+)-(\d+)(?::[+-])?$")
ASCRP_RE = re.compile(r"^ASCRP\d+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR, help=f"Input web table directory. Default: {DEFAULT_INPUT_DIR}")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Per-GSE clean output directory. Default: {DEFAULT_OUTPUT_DIR}")
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Clean summary TSV. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return "NA" if value in NULL_VALUES else value


def normalize_strand(value: str | None) -> str:
    value = normalize_missing(value)
    return value if value in {"+", "-"} else "NA"


def normalize_chrom(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        return "chr" + value[3:]
    return f"chr{value}"


def parse_coordinate(value: str | None) -> tuple[str, str, str]:
    text = normalize_missing(value)
    if text == "NA":
        raise ValueError("missing coordinate")
    match = COORD_RE.match(text)
    if not match:
        raise ValueError(f"invalid coordinate format: {text!r}")
    chrom, start, end = match.groups()
    start_int = int(start)
    end_int = int(end)
    if start_int < 0 or end_int < 0:
        raise ValueError(f"negative coordinate: {text!r}")
    if start_int > end_int:
        raise ValueError(f"start greater than end: {text!r}")
    return normalize_chrom(chrom), str(start_int), str(end_int)


def choose_coordinate(row: dict[str, str], counts: Counter[str]) -> tuple[str, str, str]:
    if row.get("geo_dataset") == "GSE136113":
        try:
            coord = parse_coordinate(row.get("original_id"))
            counts["coordinate_from_original_id"] += 1
            return coord
        except ValueError as first_error:
            try:
                coord = parse_coordinate(row.get("genomic_location"))
                counts["coordinate_from_genomic_location"] += 1
                return coord
            except ValueError as second_error:
                raise ValueError(f"cannot parse original_id or genomic_location: {first_error}; {second_error}") from second_error

    try:
        coord = parse_coordinate(row.get("genomic_location"))
        counts["coordinate_from_genomic_location"] += 1
        return coord
    except ValueError as first_error:
        try:
            coord = parse_coordinate(row.get("original_id"))
            counts["coordinate_from_original_id"] += 1
            return coord
        except ValueError as second_error:
            raise ValueError(f"cannot parse genomic_location or original_id: {first_error}; {second_error}") from second_error


def clean_row(row: dict[str, str], counts: Counter[str]) -> dict[str, str]:
    circ_id = normalize_missing(row.get("circbase_id"))
    if circ_id == "NA":
        circ_id = normalize_missing(row.get("original_id"))
    if circ_id == "NA":
        raise ValueError("missing circRNA ID")

    chrom, start, end = choose_coordinate(row, counts)
    original_id = (row.get("original_id") or "").strip()
    return {
        "circRNA_id": circ_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": normalize_strand(row.get("strand")),
        "gene_symbol": normalize_missing(row.get("parental_gene")),
        ASCRP_COLUMN: original_id if ASCRP_RE.fullmatch(original_id) else "",
    }


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"Output exists; use --force to overwrite: {', '.join(existing)}")


def write_dict_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def process_file(input_path: Path, output_dir: Path) -> tuple[dict[str, str], list[dict[str, str]], list[dict[str, str]]]:
    gse = input_path.name.removesuffix(".web_annotation.tsv")
    output_path = output_dir / f"{gse}.circRNAs.clean.tsv"
    ascrp_output_path = output_dir / f"{gse}.circRNAs.clean.ASCRP.tsv"
    invalid_path = output_dir / f"{gse}.circRNAs.clean.invalid.tsv"

    counts: Counter[str] = Counter()
    rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    ascrp_by_key: dict[tuple[str, ...], set[str]] = {}
    circ_ids: set[str] = set()

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {input_path}")
        missing = [column for column in RAW_REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"{input_path} is missing required column(s): {', '.join(missing)}")

        invalid_columns = ["line_number", *reader.fieldnames, "invalid_reason"]
        for line_number, raw in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                row = clean_row(raw, counts)
            except ValueError as exc:
                counts["invalid_rows"] += 1
                invalid_rows.append({"line_number": str(line_number), **raw, "invalid_reason": str(exc)})
                continue

            counts["clean_rows_before_uniq"] += 1
            key = tuple(row[column] for column in OUTPUT_COLUMNS)
            if key in seen:
                counts["duplicate_full_rows_removed"] += 1
                if row[ASCRP_COLUMN]:
                    ascrp_by_key[key].add(row[ASCRP_COLUMN])
                continue
            seen.add(key)
            rows.append(row)
            ascrp_by_key[key] = {row[ASCRP_COLUMN]} if row[ASCRP_COLUMN] else set()
            circ_ids.add(row["circRNA_id"])
            if row["strand"] in {"+", "-"}:
                counts["stranded_rows"] += 1
            else:
                counts["unstranded_rows"] += 1

    ascrp_rows = []
    unique_ascrp_ids: set[str] = set()
    for row in rows:
        key = tuple(row[column] for column in OUTPUT_COLUMNS)
        ascrp_ids = sorted(ascrp_by_key[key])
        unique_ascrp_ids.update(ascrp_ids)
        ascrp_rows.append({**row, ASCRP_COLUMN: ",".join(ascrp_ids)})

    write_dict_rows(output_path, OUTPUT_COLUMNS, rows)
    write_dict_rows(ascrp_output_path, ASCRP_OUTPUT_COLUMNS, ascrp_rows)
    write_dict_rows(invalid_path, invalid_columns if invalid_rows else ["line_number", *RAW_REQUIRED_COLUMNS, "invalid_reason"], invalid_rows)
    counts["clean_rows"] = len(rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["rows_with_ascrp"] = sum(bool(row[ASCRP_COLUMN]) for row in ascrp_rows)
    counts["unique_ascrp_ids"] = len(unique_ascrp_ids)

    summary = {
        "geo_dataset": gse,
        "input_file": str(input_path),
        "output_file": str(output_path),
        "ascrp_output_file": str(ascrp_output_path),
        "invalid_output_file": str(invalid_path),
    }
    for column in SUMMARY_COLUMNS:
        if column not in summary:
            summary[column] = str(counts[column])
    return summary, rows, invalid_rows


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Missing input directory: {input_dir}")

    input_files = sorted(input_dir.glob("*.web_annotation.tsv"))
    if not input_files:
        raise FileNotFoundError(f"No *.web_annotation.tsv files found in {input_dir}")

    output_paths = [output_dir / f"{path.name.removesuffix('.web_annotation.tsv')}.circRNAs.clean.tsv" for path in input_files]
    output_paths += [output_dir / f"{path.name.removesuffix('.web_annotation.tsv')}.circRNAs.clean.ASCRP.tsv" for path in input_files]
    output_paths += [output_dir / f"{path.name.removesuffix('.web_annotation.tsv')}.circRNAs.clean.invalid.tsv" for path in input_files]
    output_paths.append(summary_output)
    refuse_overwrite(output_paths, args.force)

    summaries: list[dict[str, str]] = []
    for input_path in input_files:
        summary, _rows, _invalid = process_file(input_path, output_dir)
        summaries.append(summary)
        print(f"{summary['geo_dataset']}\tclean_rows={summary['clean_rows']}\tinvalid_rows={summary['invalid_rows']}")

    write_dict_rows(summary_output, SUMMARY_COLUMNS, summaries)
    print(f"summary\t{summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
