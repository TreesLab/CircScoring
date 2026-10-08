#!/usr/bin/env python3
"""Clean circVAR GRCh38 BED5 circRNA table into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW = BASE_DIR / "raw" / "circ_GRCH38_bedfile.txt"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "35_circVAR.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "35_circVAR.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "35_circVAR.circRNAs.clean.summary.tsv"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["line_number", "invalid_reason", "raw_row"]
VALID_RAW_MISSING_STRANDS = {"x", "not_available", ".", "", "NA", "N/A"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW, help=f"Raw GRCh38 BED5 file. Default: {DEFAULT_RAW}")
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


def chrom_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    if value.isdigit():
        return (int(value), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (special.get(value, 1000), value)


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        raise FileExistsError("Output exists; use --force to overwrite: " + ", ".join(str(path) for path in existing))


def normalize_strand(value: str) -> str:
    value = value.strip()
    if value in {"+", "-"}:
        return value
    if value in VALID_RAW_MISSING_STRANDS:
        return "NA"
    raise ValueError(f"invalid strand: {value!r}")


def convert_row(row: list[str]) -> tuple[str, str, str, str, str, str]:
    if len(row) != 5:
        raise ValueError(f"expected 5 columns, got {len(row)}")
    chrom, start_text, end_text, strand_text, circ_id = [value.strip() for value in row]
    if not chrom.startswith("chr"):
        raise ValueError(f"invalid chrom: {chrom!r}")
    if not circ_id:
        raise ValueError("missing circRNA_id")
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise ValueError(f"invalid coordinates: start={start_text!r}, end={end_text!r}") from exc
    if start < 0 or end < 0:
        raise ValueError(f"negative coordinate: start={start}, end={end}")
    if start >= end:
        raise ValueError(f"0-based interval must satisfy start < end: start={start}, end={end}")
    strand = normalize_strand(strand_text)
    return circ_id, chrom, str(start), str(end), strand, "NA"


def write_summary(path: Path, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = [
        "raw_rows",
        "clean_rows_before_uniq",
        "duplicate_full_rows_removed",
        "clean_rows",
        "invalid_rows",
        "unique_circRNA_ids",
        "stranded_rows",
        "unstranded_rows",
        "gene_symbol_na_rows",
        "strand:+",
        "strand:-",
        "raw_strand:x",
        "raw_strand:not_available",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in keys:
            writer.writerow([key, counts[key]])


def main() -> int:
    args = parse_args()
    raw = args.raw.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()

    refuse_overwrite([output, invalid_output, summary_output], args.force)
    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    rows: list[tuple[str, str, str, str, str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    circ_ids: set[str] = set()

    with raw.open(newline="") as raw_handle, invalid_output.open("w", newline="") as invalid_handle:
        reader = csv.reader(raw_handle, delimiter="\t")
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            counts["raw_rows"] += 1
            if len(row) >= 4:
                counts[f"raw_strand:{row[3].strip()}"] += 1
            try:
                clean_row = convert_row(row)
            except ValueError as exc:
                counts["invalid_rows"] += 1
                invalid_writer.writerow({"line_number": line_number, "invalid_reason": str(exc), "raw_row": repr(row)})
                continue

            counts["clean_rows_before_uniq"] += 1
            if clean_row in seen:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen.add(clean_row)
            rows.append(clean_row)
            circ_ids.add(clean_row[0])
            if clean_row[4] in {"+", "-"}:
                counts["stranded_rows"] += 1
                counts[f"strand:{clean_row[4]}"] += 1
            else:
                counts["unstranded_rows"] += 1
            counts["gene_symbol_na_rows"] += 1

    rows.sort(key=lambda row: (chrom_key(row[1]), int(row[2]), int(row[3]), row[4], row[0]))
    counts["clean_rows"] = len(rows)
    counts["unique_circRNA_ids"] = len(circ_ids)

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(rows)

    write_summary(summary_output, counts)
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    print(f"output\t{output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
