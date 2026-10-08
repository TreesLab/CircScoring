#!/usr/bin/env python3
"""Clean CircAI human circRNA raw zip into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_ZIP = BASE_DIR / "raw" / "human_circrna.zip"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "40_CircAI.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "40_CircAI.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "40_CircAI.circRNAs.clean.summary.tsv"
ZIP_MEMBER = "human_circrna.txt"

RAW_COLUMNS = [
    "CircAI_ID",
    "Chromosome",
    "Start_of_circRNA",
    "End_of_circRNA",
    "Name_of_circular_transcript",
    "Strand",
    "Parental_gene_name",
    "Exon_sizes",
    "Exon_offsets",
    "seq",
    "length",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["line_number", *RAW_COLUMNS, "invalid_reason"]
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
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean CircAI human circRNA raw zip into unified six-column TSV.")
    parser.add_argument("--raw-zip", type=Path, default=DEFAULT_RAW_ZIP, help=f"Raw zip input. Default: {DEFAULT_RAW_ZIP}")
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


def normalize_chrom(value: str | None) -> str:
    value = normalize_missing(value)
    if value == "NA":
        return value
    return value if value.startswith("chr") else f"chr{value}"


def normalize_strand(value: str | None) -> str:
    value = normalize_missing(value)
    if value in {"+", "-"}:
        return value
    return "NA"


def parse_int(value: str | None, field: str) -> str:
    value = normalize_missing(value)
    if value == "NA":
        raise ValueError(f"missing {field}")
    number = int(value)
    if number < 0:
        raise ValueError(f"negative {field}: {value}")
    return str(number)


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


def raw_rows(raw_zip: Path):
    list_result = subprocess.run(
        ["unzip", "-Z", "-1", str(raw_zip)],
        check=True,
        capture_output=True,
        text=True,
    )
    if ZIP_MEMBER not in list_result.stdout.splitlines():
        raise FileNotFoundError(f"{raw_zip} does not contain {ZIP_MEMBER}")

    process = subprocess.Popen(
        ["unzip", "-p", str(raw_zip), ZIP_MEMBER],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
    )
    assert process.stdout is not None
    reader = csv.DictReader(process.stdout, delimiter="\t")
    if reader.fieldnames != RAW_COLUMNS:
        process.kill()
        process.wait()
        raise ValueError(f"Unexpected header in {ZIP_MEMBER}: {reader.fieldnames}")

    try:
        for line_number, row in enumerate(reader, start=2):
            yield line_number, row
    finally:
        stderr = process.stderr.read() if process.stderr else ""
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"unzip -p failed with exit code {return_code}: {stderr.strip()}")


def main() -> int:
    args = parse_args()
    raw_zip = args.raw_zip.expanduser().resolve()
    output = args.output.expanduser().resolve()
    invalid_output = args.invalid_output.expanduser().resolve()
    summary_output = args.summary_output.expanduser().resolve()
    refuse_overwrite([output, invalid_output, summary_output], args.force)

    if not raw_zip.exists():
        raise FileNotFoundError(f"Missing raw zip: {raw_zip}")

    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    clean_rows: list[tuple[str, str, str, str, str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    circ_ids: set[str] = set()
    genes: set[str] = set()

    with invalid_output.open("w", newline="") as invalid_handle:
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for line_number, raw in raw_rows(raw_zip):
            counts["raw_rows"] += 1
            try:
                circ_id = normalize_missing(raw.get("CircAI_ID"))
                chrom = normalize_chrom(raw.get("Chromosome"))
                start = parse_int(raw.get("Start_of_circRNA"), "start")
                end = parse_int(raw.get("End_of_circRNA"), "end")
                strand = normalize_strand(raw.get("Strand"))
                gene = normalize_missing(raw.get("Parental_gene_name"))

                if circ_id == "NA":
                    raise ValueError("missing CircAI_ID")
                if chrom == "NA":
                    raise ValueError("missing chromosome")
                if int(start) > int(end):
                    raise ValueError(f"start is greater than end: {start} > {end}")
                if strand == "NA":
                    raise ValueError(f"invalid strand: {raw.get('Strand')!r}")

                row = (circ_id, chrom, start, end, strand, gene)
                counts["clean_rows_before_uniq"] += 1
                if row in seen:
                    counts["duplicate_full_rows_removed"] += 1
                    continue
                seen.add(row)
                clean_rows.append(row)
                circ_ids.add(circ_id)
                if strand in {"+", "-"}:
                    counts["stranded_rows"] += 1
                else:
                    counts["unstranded_rows"] += 1
                if gene != "NA":
                    genes.add(gene)
            except Exception as exc:
                counts["invalid_rows"] += 1
                invalid_writer.writerow({"line_number": line_number, **raw, "invalid_reason": str(exc)})

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(clean_rows)

    counts["clean_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_genes"] = len(genes)
    write_summary(summary_output, counts)

    print(f"raw rows: {counts['raw_rows']}")
    print(f"clean rows: {counts['clean_rows']}")
    print(f"invalid rows: {counts['invalid_rows']}")
    print(f"duplicate full rows removed: {counts['duplicate_full_rows_removed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
