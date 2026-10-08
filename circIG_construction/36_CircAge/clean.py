#!/usr/bin/env python3
"""Clean CircAge human circRNA raw tables into the unified six-column format."""

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
DEFAULT_OUTPUT = BASE_DIR / "clean" / "36_CircAge.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "36_CircAge.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "36_CircAge.circRNAs.clean.summary.tsv"

RAW_FILES = [
    "human_Aortic_main.txt",
    "human_Lung_main.txt",
    "human_Skin_main.txt",
    "human_Umbilical_main.txt",
    "human_white_blood_cell_main.txt",
]

RAW_COLUMNS = [
    "circRNA",
    "species",
    "tissue",
    "strand",
    "hostGene",
    "circType",
    "geneType",
    "conservation",
    "algorithm",
    "ageRelated",
    "expression",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["source_file", "line_number", *RAW_COLUMNS, "invalid_reason"]
SUMMARY_KEYS = [
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "unique_genes",
    "unique_tissues",
]

COORD_RE = re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)\|(?P<end>\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CircAge human circRNA raw tables into unified six-column TSV.")
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


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def parse_coord(raw_coord: str) -> tuple[str, str, str]:
    match = COORD_RE.match(raw_coord.strip())
    if not match:
        raise ValueError(f"invalid circRNA coordinate format: {raw_coord!r}")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"start is greater than end: {raw_coord!r}")
    return match.group("chrom"), str(start), str(end)


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def write_summary(path: Path, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in SUMMARY_KEYS:
            writer.writerow([key, counts[key]])


def raw_rows(raw_dir: Path):
    for filename in RAW_FILES:
        path = raw_dir / filename
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for line_number, values in enumerate(reader, start=1):
                raw = {column: values[index] if index < len(values) else "" for index, column in enumerate(RAW_COLUMNS)}
                yield filename, line_number, values, raw


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([output, invalid_output, summary_output], args.force)

    missing = [str(raw_dir / filename) for filename in RAW_FILES if not (raw_dir / filename).exists()]
    if missing:
        raise FileNotFoundError(f"Missing raw file(s): {', '.join(missing)}")

    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()
    genes: set[str] = set()
    tissues: set[str] = set()

    with (
        output.open("w", newline="") as output_handle,
        invalid_output.open("w", newline="") as invalid_handle,
    ):
        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for source_file, line_number, values, raw in raw_rows(raw_dir):
            counts["raw_rows"] += 1
            tissues.add(normalize_missing(raw.get("tissue")))
            try:
                if len(values) != len(RAW_COLUMNS):
                    raise ValueError(f"expected {len(RAW_COLUMNS)} columns, observed {len(values)}")

                circRNA_id = normalize_missing(raw.get("circRNA"))
                if circRNA_id == "NA":
                    raise ValueError("missing circRNA")
                chrom, start, end = parse_coord(circRNA_id)

                strand = normalize_missing(raw.get("strand"))
                if strand not in {"+", "-"}:
                    raise ValueError(f"invalid strand: {strand!r}")

                row = {
                    "circRNA_id": circRNA_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_symbol": normalize_missing(raw.get("hostGene")),
                }
            except ValueError as exc:
                invalid = {"source_file": source_file, "line_number": str(line_number), "invalid_reason": str(exc)}
                invalid.update({column: raw.get(column, "") for column in RAW_COLUMNS})
                invalid_writer.writerow(invalid)
                counts["invalid_rows"] += 1
                continue

            counts["clean_rows_before_uniq"] += 1
            key = tuple(row[column] for column in OUTPUT_COLUMNS)
            if key in seen_rows:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen_rows.add(key)
            writer.writerow(row)
            circ_ids.add(row["circRNA_id"])
            genes.add(row["gene_symbol"])
            counts["stranded_rows"] += 1

    counts["clean_rows"] = len(seen_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_genes"] = len(genes)
    counts["unique_tissues"] = len(tissues)
    write_summary(summary_output, counts)

    print(f"raw_rows\t{counts['raw_rows']}")
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
