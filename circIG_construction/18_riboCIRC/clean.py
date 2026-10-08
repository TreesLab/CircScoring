#!/usr/bin/env python3
"""Clean riboCIRC raw archives into the unified six-column circRNA table."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Iterator

from coordinate_rules import normalize_hg38_0based_start
from species_rules import classify_literature_human


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "18_riboCIRC.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "18_riboCIRC.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "18_riboCIRC.circRNAs.clean.summary.tsv"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["source_table", "line_number", "invalid_reason", "raw_row"]
SUMMARY_KEYS = [
    "raw_rows",
    "human_rows",
    "non_human_rows",
    "clean_rows_before_uniq",
    "duplicate_full_rows_removed",
    "clean_rows",
    "invalid_rows",
    "unique_circRNA_ids",
    "stranded_rows",
    "unstranded_rows",
    "unique_gene_symbols",
    "annotation_guided_human_rows",
    "context_specific_human_rows",
    "literature_reported_human_rows",
    "literature_reported_human_by_species_rows",
    "literature_reported_human_by_hg38_hsa_id_rows",
    "coordinate_exception_rows",
]

POSITION_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")
MISSING_VALUES = {"", "NA", "N/A", "NULL", "NONE", "NAN"}

SOURCES = [
    {
        "name": "annotation_guided",
        "rar": "All_annotation-guided_ribo-circRNA.rar",
        "member": "All_annotation-guided_ribo-circRNA/basic_info_Annotation-guided.txt",
        "species_value": "human",
        "type": "basic_info",
    },
    {
        "name": "context_specific",
        "rar": "All_context-specific_ribo-circRNAs.rar",
        "member": "All_context-specific_ribo-circRNAs/basic_info_Context-specific.txt",
        "species_value": "human",
        "type": "basic_info",
    },
    {
        "name": "literature_reported",
        "rar": "literature-reported_trans-circRNAs.rar",
        "member": "literature-reported_trans-circRNAs/literature-reported_trans-circRNAs.txt",
        "species_value": "Homo sapiens",
        "type": "literature",
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean riboCIRC raw RAR archives into unified six-column TSV.")
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


def normalize_missing(value: str | None, *, dash_is_missing: bool = False) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if value.upper() in MISSING_VALUES or (dash_is_missing and value == "-"):
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
        raise ValueError(f"invalid Genome_position: {value!r}")
    chrom, start, end = match.groups()
    start_int = int(start)
    end_int = int(end)
    if start_int < 0 or end_int < 0:
        raise ValueError(f"negative coordinate in Genome_position: {value!r}")
    if start_int > end_int:
        raise ValueError(f"start is greater than end in Genome_position: {value!r}")
    return chrom, str(start_int), str(end_int)


def choose_literature_id(row: dict[str, str]) -> str:
    for field in ("riboCIRC_ID", "circBase_ID", "circRNA_name"):
        value = normalize_missing(row.get(field), dash_is_missing=True)
        if value != "NA":
            return value
    return "NA"


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def iter_rar_tsv(raw_dir: Path, rar_name: str, member: str) -> Iterator[tuple[int, dict[str, str]]]:
    rar_path = raw_dir / rar_name
    if not rar_path.exists():
        raise FileNotFoundError(f"Missing raw RAR: {rar_path}")

    process = subprocess.Popen(
        ["bsdtar", "-xOf", str(rar_path), member],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
    )
    assert process.stdout is not None
    reader = csv.DictReader(process.stdout, delimiter="\t")
    try:
        for line_number, row in enumerate(reader, start=2):
            yield line_number, row
    finally:
        stderr = process.stderr.read() if process.stderr else ""
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"bsdtar failed for {rar_path}:{member} with exit code {return_code}: {stderr.strip()}")


def convert_row(
    source: dict[str, str], row: dict[str, str]
) -> tuple[tuple[str, str, str, str, str, str], bool]:
    chrom, start, end = parse_position(row.get("Genome_position"))
    strand = normalize_strand(row.get("Strand"))
    if strand == "NA":
        raise ValueError(f"invalid strand: {row.get('Strand')!r}")

    if source["type"] == "basic_info":
        circ_id = normalize_missing(row.get("riboCIRC_ID"), dash_is_missing=True)
        gene_symbol = normalize_missing(row.get("Host_gene_Symbol"), dash_is_missing=True)
    elif source["type"] == "literature":
        circ_id = choose_literature_id(row)
        gene_symbol = "NA"
    else:
        raise ValueError(f"unsupported source type: {source['type']}")

    if circ_id == "NA":
        raise ValueError("missing circRNA_id")

    start_int, coordinate_corrected = normalize_hg38_0based_start(
        source_name=source["name"],
        circbase_id=row.get("circBase_ID"),
        assembly=row.get("Genome_assembly"),
        chrom=chrom,
        start=int(start),
        end=int(end),
        strand=strand,
    )
    clean_row = circ_id, chrom, str(start_int), end, strand, gene_symbol
    return clean_row, coordinate_corrected


def classify_human_row(source: dict[str, str], row: dict[str, str]) -> tuple[bool, str | None]:
    if source["type"] == "literature":
        return classify_literature_human(row)
    if row.get("Species") == source["species_value"]:
        return True, "species"
    return False, None


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

    output.parent.mkdir(parents=True, exist_ok=True)
    invalid_output.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    clean_rows: list[tuple[str, str, str, str, str, str]] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    circ_ids: set[str] = set()
    gene_symbols: set[str] = set()

    with invalid_output.open("w", newline="") as invalid_handle:
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for source in SOURCES:
            for line_number, raw in iter_rar_tsv(raw_dir, source["rar"], source["member"]):
                counts["raw_rows"] += 1
                is_human, classification_reason = classify_human_row(source, raw)
                if not is_human:
                    counts["non_human_rows"] += 1
                    continue
                counts["human_rows"] += 1
                counts[f"{source['name']}_human_rows"] += 1
                if source["type"] == "literature":
                    counts[f"literature_reported_human_by_{classification_reason}_rows"] += 1
                try:
                    row, coordinate_corrected = convert_row(source, raw)
                    if coordinate_corrected:
                        counts["coordinate_exception_rows"] += 1
                    counts["clean_rows_before_uniq"] += 1
                    if row in seen:
                        counts["duplicate_full_rows_removed"] += 1
                        continue
                    seen.add(row)
                    clean_rows.append(row)
                    circ_ids.add(row[0])
                    if row[4] in {"+", "-"}:
                        counts["stranded_rows"] += 1
                    else:
                        counts["unstranded_rows"] += 1
                    if row[5] != "NA":
                        gene_symbols.add(row[5])
                except Exception as exc:
                    counts["invalid_rows"] += 1
                    invalid_writer.writerow(
                        {
                            "source_table": source["name"],
                            "line_number": line_number,
                            "invalid_reason": str(exc),
                            "raw_row": json.dumps(raw, ensure_ascii=True, sort_keys=True),
                        }
                    )

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(clean_rows)

    counts["clean_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_gene_symbols"] = len(gene_symbols)
    write_summary(summary_output, counts)

    print(f"raw rows: {counts['raw_rows']}")
    print(f"human rows: {counts['human_rows']}")
    print(f"clean rows: {counts['clean_rows']}")
    print(f"invalid rows: {counts['invalid_rows']}")
    print(f"duplicate full rows removed: {counts['duplicate_full_rows_removed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
