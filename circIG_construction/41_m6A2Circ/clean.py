#!/usr/bin/env python3
"""Clean m6A2Circ human circRNA raw zip into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_ZIP = BASE_DIR / "raw" / "huamn_circRNAs_Information.zip"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "41_m6A2Circ.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "41_m6A2Circ.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "41_m6A2Circ.circRNAs.clean.summary.tsv"
ZIP_MEMBER = "huamn_circRNAs_Information.txt"

RAW_COLUMNS = [
    "m6A2CircID",
    "Organism",
    "Position",
    "CancerNormal",
    "Tissue",
    "Tools",
    "ConfidenceLevel",
    "OtherDatabaseID",
    "HostGene",
    "Strand",
    "gene_id",
    "transcript_id",
    "gene_biotype",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["line_number", *RAW_COLUMNS, "invalid_reason"]
SUMMARY_KEYS = [
    "raw_rows",
    "human_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "unique_genes",
]
POSITION_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)(?::([+-]))?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean m6A2Circ human circRNA raw zip into unified six-column TSV.")
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
    if not value or value in {"."} or value.upper() in {"NA", "N/A", "NULL", "NONE", "NAN"}:
        return "NA"
    return value


def normalize_strand(value: str | None) -> str:
    value = normalize_missing(value)
    if value in {"+", "-"}:
        return value
    return "NA"


def parse_position(value: str | None) -> tuple[str, str, str]:
    value = normalize_missing(value)
    match = POSITION_RE.match(value)
    if not match:
        raise ValueError(f"invalid Position: {value!r}")
    chrom, start, end, _position_strand = match.groups()
    start_int = int(start)
    end_int = int(end)
    if start_int < 0 or end_int < 0:
        raise ValueError(f"negative coordinate in Position: {value!r}")
    if start_int > end_int:
        raise ValueError(f"start is greater than end in Position: {value!r}")
    return chrom, str(start_int), str(end_int)


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
                organism = normalize_missing(raw.get("Organism"))
                if organism != "Homo sapiens":
                    raise ValueError(f"non-human organism: {organism}")
                counts["human_rows"] += 1

                circ_id = normalize_missing(raw.get("m6A2CircID"))
                chrom, start, end = parse_position(raw.get("Position"))
                strand = normalize_strand(raw.get("Strand"))
                gene = normalize_missing(raw.get("HostGene"))

                if circ_id == "NA":
                    raise ValueError("missing m6A2CircID")

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
