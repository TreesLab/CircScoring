#!/usr/bin/env python3
"""Clean CircRic official expression CSV into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = [
    "circ_location",
    "ccl_id",
    "ccl_name",
    "cancer_type",
    "number_bsReads",
    "circ_gene",
    "invalid_reason",
]
SUMMARY_COLUMNS = [
    "input",
    "output",
    "raw_rows",
    "invalid_rows",
    "duplicate_full_rows_removed",
    "clean_rows",
    "unique_circRNA_ids",
    "unique_gene_symbols",
    "raw_stranded_rows",
    "raw_unstranded_rows",
]
NULL_VALUES = {"", "-", ".", "*", "NA", "N/A", "None", "none", "null", "NULL"}
LOCATION_RE = re.compile(r"^(?P<chrom>chr)?(?P<chrom_body>[0-9XYM]+)_(?P<start>[0-9]+)_(?P<end>[0-9]+)\|(?P<gene>.+)$")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "raw" / "circRNA_expression.csv.zip")
    parser.add_argument("--member", default="circRNA_expression.csv")
    parser.add_argument("--output", type=Path, default=root / "clean" / "42_CircRic.circRNAs.clean.tsv")
    parser.add_argument("--invalid-output", type=Path, default=root / "clean" / "42_CircRic.circRNAs.clean.invalid.tsv")
    parser.add_argument("--summary-output", type=Path, default=root / "clean" / "42_CircRic.circRNAs.clean.summary.tsv")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def clean_text(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return "NA" if text in NULL_VALUES else text


def ensure_can_write(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists; use --force to overwrite")


def parse_location(value: str) -> tuple[str, str, str, str] | None:
    match = LOCATION_RE.match(value)
    if not match:
        return None
    chrom_body = match.group("chrom_body")
    chrom = chrom_body if chrom_body.startswith("chr") else f"chr{chrom_body}"
    start = int(match.group("start"))
    end = int(match.group("end"))
    gene = clean_text(match.group("gene"))
    if start > end or gene == "NA":
        return None
    return chrom, str(start), str(end), gene


def open_csv_from_zip(path: Path, member: str):
    archive = zipfile.ZipFile(path)
    try:
        handle = archive.open(member)
    except KeyError:
        archive.close()
        members = ", ".join(archive.namelist())
        raise KeyError(f"{path}: missing member {member!r}; available members: {members}") from None
    text_handle = io.TextIOWrapper(handle, encoding="utf-8-sig", newline="")
    return archive, text_handle


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    for path in (args.output, args.invalid_output, args.summary_output):
        ensure_can_write(path, args.force)

    counts: Counter[str] = Counter()
    clean_rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    gene_symbols: set[str] = set()

    archive, text_handle = open_csv_from_zip(input_path, args.member)
    with archive, text_handle:
        reader = csv.DictReader(text_handle)
        expected = ["circ_location", "ccl_id", "ccl_name", "cancer_type", "number_bsReads", "circ_gene"]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected header in {input_path}:{args.member}: {reader.fieldnames!r}; expected {expected!r}")

        for row in reader:
            counts["raw_rows"] += 1
            circ_location = clean_text(row.get("circ_location"))
            circ_gene = clean_text(row.get("circ_gene"))
            parsed = parse_location(circ_location)
            reasons: list[str] = []
            if parsed is None:
                reasons.append("invalid_circ_location")
            else:
                _chrom, _start, _end, gene_symbol = parsed
                if circ_gene != "NA" and circ_gene != gene_symbol:
                    reasons.append("circ_gene_mismatch")

            if reasons:
                invalid = {column: clean_text(row.get(column)) for column in INVALID_COLUMNS if column != "invalid_reason"}
                invalid["invalid_reason"] = ",".join(reasons)
                invalid_rows.append(invalid)
                counts["invalid_rows"] += 1
                continue

            assert parsed is not None
            chrom, start, end, gene_symbol = parsed
            clean_row = {
                "circRNA_id": circ_location,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": "NA",
                "gene_symbol": gene_symbol,
            }
            key = tuple(clean_row[column] for column in OUTPUT_COLUMNS)
            counts["raw_unstranded_rows"] += 1
            if key in seen:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen.add(key)
            clean_rows.append(clean_row)
            gene_symbols.add(gene_symbol)

    counts["clean_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len({row["circRNA_id"] for row in clean_rows})
    counts["unique_gene_symbols"] = len(gene_symbols)

    write_tsv(args.output, OUTPUT_COLUMNS, clean_rows)
    write_tsv(args.invalid_output, INVALID_COLUMNS, invalid_rows)
    summary = {
        "input": f"{input_path}:{args.member}",
        "output": str(args.output),
        "raw_rows": str(counts["raw_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
        "clean_rows": str(counts["clean_rows"]),
        "unique_circRNA_ids": str(counts["unique_circRNA_ids"]),
        "unique_gene_symbols": str(counts["unique_gene_symbols"]),
        "raw_stranded_rows": str(counts["raw_stranded_rows"]),
        "raw_unstranded_rows": str(counts["raw_unstranded_rows"]),
    }
    write_tsv(args.summary_output, SUMMARY_COLUMNS, [summary])
    for key in SUMMARY_COLUMNS[2:]:
        print(f"{key}\t{summary[key]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
