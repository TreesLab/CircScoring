#!/usr/bin/env python3
"""Clean MiOncoCirc raw release table into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["chrom", "start", "end", "reads", "symbol", "sample", "release", "invalid_reason"]
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


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=root / "raw" / "v0.1.release.txt")
    parser.add_argument("--output", type=Path, default=root / "clean" / "29_MiOncoCirc.circRNAs.clean.tsv")
    parser.add_argument("--invalid-output", type=Path, default=root / "clean" / "29_MiOncoCirc.circRNAs.clean.invalid.tsv")
    parser.add_argument("--summary-output", type=Path, default=root / "clean" / "29_MiOncoCirc.circRNAs.clean.summary.tsv")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def clean_text(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return "NA" if text in NULL_VALUES else text


def parse_int(value: str | None) -> int | None:
    text = clean_text(value)
    if text == "NA":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def ensure_can_write(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"{path} exists; use --force to overwrite")


def circ_id(chrom: str, start: str, end: str, gene_symbol: str) -> str:
    return f"{chrom}:{start}|{end}:{gene_symbol}"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    args.input = args.input.expanduser().resolve()
    for path in (args.output, args.invalid_output, args.summary_output):
        ensure_can_write(path, args.force)

    counts: Counter[str] = Counter()
    clean_rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    gene_symbols: set[str] = set()

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["chr", "start", "end", "reads", "symbol", "sample", "release"]
        if reader.fieldnames != expected:
            raise ValueError(f"Unexpected header in {args.input}: {reader.fieldnames!r}; expected {expected!r}")

        for row in reader:
            counts["raw_rows"] += 1
            chrom = clean_text(row.get("chr"))
            start = parse_int(row.get("start"))
            end = parse_int(row.get("end"))
            gene_symbol = clean_text(row.get("symbol"))
            reasons: list[str] = []
            if chrom == "NA":
                reasons.append("missing_chrom")
            if start is None:
                reasons.append("invalid_start")
            if end is None:
                reasons.append("invalid_end")
            if start is not None and end is not None and end < start:
                reasons.append("end_before_start")
            if gene_symbol == "NA":
                reasons.append("missing_gene_symbol")

            if reasons:
                invalid = {key: clean_text(row.get(key)) for key in ["chr", "start", "end", "reads", "symbol", "sample", "release"]}
                invalid["chrom"] = invalid.pop("chr")
                invalid["invalid_reason"] = ",".join(reasons)
                invalid_rows.append(invalid)
                counts["invalid_rows"] += 1
                continue

            assert start is not None and end is not None
            clean_row = {
                "circRNA_id": circ_id(chrom, str(start), str(end), gene_symbol),
                "chrom": chrom,
                "start": str(start),
                "end": str(end),
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
        "input": str(args.input),
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
