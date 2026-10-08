#!/usr/bin/env python3
"""Clean CSCD2 raw circRNA annotations into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "26_CSCD2.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "26_CSCD2.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "26_CSCD2.circRNAs.clean.summary.tsv"

GROUPS = {"cancer", "normal", "common"}
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
EXPECTED_HEADER = [
    "circRNA",
    "circbase_id",
    "type",
    "sample_source",
    "position",
    "host_gene",
    "reads_counts",
    "algorithm",
    "chain",
    "subcellular_location",
    "gene_type",
    "average_read_counts",
    "average_log2SRPTM",
    "log2SRPTM",
    "alternative_splicing",
    "ratios",
    "works_cited",
    "mioncocirc",
    "number_of_algorithm",
    "sequence",
]
SUMMARY_KEYS = [
    "input_file_count",
    "raw_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "cancer_rows",
    "normal_rows",
    "common_rows",
    "source_chromosome_mismatch_rows",
]

CIRCRNA_RE = re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)\|(?P<end>\d+)$")
SOURCE_RE = re.compile(r"^hg38_(?P<group>cancer|normal|common)_circrna_circ_(?P<chrom>chr[^.]+)\.txt\.gz$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CSCD2 raw annotation files into unified six-column TSV.")
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
    if not value or value == "---" or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def parse_circrna(value: str) -> tuple[str, str, str]:
    match = CIRCRNA_RE.match(value.strip())
    if not match:
        raise ValueError(f"invalid circRNA coordinate: {value!r}")
    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start > end:
        raise ValueError(f"start is greater than end: {value!r}")
    return chrom, str(start), str(end)


def parse_source(path: Path) -> tuple[str, str]:
    match = SOURCE_RE.match(path.name)
    if not match:
        raise ValueError(f"unexpected source filename: {path.name}")
    return match.group("group"), match.group("chrom")


def iter_input_files(raw_dir: Path) -> list[Path]:
    files: list[Path] = []
    for group in sorted(GROUPS):
        files.extend(sorted((raw_dir / group).glob("*.txt.gz")))
    if not files:
        raise FileNotFoundError(f"No CSCD2 .txt.gz files found under {raw_dir}")
    return files


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


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([output, invalid_output, summary_output], args.force)

    input_files = iter_input_files(raw_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    counts["input_file_count"] = len(input_files)
    seen_rows: set[tuple[str, ...]] = set()
    circ_ids: set[str] = set()

    with (
        output.open("w", newline="") as output_handle,
        invalid_output.open("w", newline="") as invalid_handle,
    ):
        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        invalid_columns = EXPECTED_HEADER + ["source_file", "line_number", "invalid_reason"]
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=invalid_columns, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for path in input_files:
            source_group, source_chrom = parse_source(path)
            with gzip.open(path, "rt", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != EXPECTED_HEADER:
                    raise ValueError(f"Unexpected header in {path}: {reader.fieldnames!r}")

                for line_number, raw in enumerate(reader, start=2):
                    counts["raw_rows"] += 1
                    counts[f"{source_group}_rows"] += 1
                    try:
                        circ_id = normalize_missing(raw.get("circRNA"))
                        if circ_id == "NA":
                            raise ValueError("missing circRNA")
                        chrom, start, end = parse_circrna(circ_id)
                        if chrom != source_chrom:
                            counts["source_chromosome_mismatch_rows"] += 1
                        strand = normalize_missing(raw.get("chain"))
                        if strand not in {"+", "-"}:
                            strand = "NA"
                        row = {
                            "circRNA_id": circ_id,
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "strand": strand,
                            "gene_symbol": normalize_missing(raw.get("host_gene")),
                        }
                    except ValueError as exc:
                        invalid = {column: raw.get(column, "") for column in EXPECTED_HEADER}
                        invalid["source_file"] = str(path.relative_to(raw_dir))
                        invalid["line_number"] = str(line_number)
                        invalid["invalid_reason"] = str(exc)
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
                    if row["strand"] in {"+", "-"}:
                        counts["stranded_rows"] += 1
                    else:
                        counts["unstranded_rows"] += 1

    counts["clean_rows"] = len(seen_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    write_summary(summary_output, counts)

    print(f"input_files\t{counts['input_file_count']}")
    print(f"raw_rows\t{counts['raw_rows']}")
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
