#!/usr/bin/env python3
"""Clean CircMine raw circRNA records into the unified six-column table."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


POSITION_RE = re.compile(r"^(chr[0-9A-Za-z]+):([0-9]+)-([0-9]+)$")
NA_VALUES = {"", "NA", "N/A", "None", "none", "null", "NULL"}
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Clean CircMine raw merged TSV.")
    parser.add_argument("--input", type=Path, default=root / "raw" / "human_circRNAs_of_circMine.raw_merged.tsv")
    parser.add_argument("--out-dir", type=Path, default=root / "clean")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--invalid-output", type=Path, default=None)
    parser.add_argument("--dataset-audit-output", type=Path, default=None)
    parser.add_argument("--summary-output", type=Path, default=None)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def clean_text(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return "NA" if text in NA_VALUES else text


def parse_position(value: str) -> tuple[str, str, str] | None:
    match = POSITION_RE.match(value)
    if not match:
        return None
    chrom, start, end = match.groups()
    if int(start) > int(end):
        return None
    return chrom, start, end


def choose_circrna_id(row: dict[str, str]) -> str:
    circrna = clean_text(row.get("circrna"))
    if circrna != "NA":
        return circrna
    return clean_text(row.get("circbank"))


def ensure_can_write(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise SystemExit(f"{path} exists; use --force to overwrite")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.output or args.out_dir / "33_CircMine.circRNAs.clean.tsv"
    invalid_output = args.invalid_output or args.out_dir / "33_CircMine.circRNAs.clean.invalid.tsv"
    dataset_audit_output = args.dataset_audit_output or args.out_dir / "33_CircMine.circRNAs.dataset_ids.tsv"
    summary_output = args.summary_output or args.out_dir / "33_CircMine.clean.summary.tsv"

    for path in [output, invalid_output, dataset_audit_output, summary_output]:
        ensure_can_write(path, args.force)

    raw_rows = 0
    invalid_rows = 0
    duplicate_full_rows = 0
    strand_counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str, str, str]] = set()
    datasets_by_clean_key: dict[tuple[str, str, str, str, str, str], set[str]] = defaultdict(set)

    with args.input.open(newline="") as input_handle, output.open("w", newline="") as out_handle, invalid_output.open("w", newline="") as invalid_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        invalid_fields = (reader.fieldnames or []) + ["invalid_reason"]
        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=invalid_fields, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        invalid_writer.writeheader()

        for row in reader:
            raw_rows += 1
            circRNA_id = choose_circrna_id(row)
            parsed = parse_position(clean_text(row.get("hg38")))
            strand = clean_text(row.get("strand"))
            gene_symbol = clean_text(row.get("hostGene"))
            dataset_id = clean_text(row.get("circmine"))

            reasons = []
            if circRNA_id == "NA":
                reasons.append("missing_circrna_id")
            if parsed is None:
                reasons.append("invalid_hg38")
            if strand not in {"+", "-"}:
                reasons.append("invalid_strand")
            if reasons:
                invalid = dict(row)
                invalid["invalid_reason"] = ",".join(reasons)
                invalid_writer.writerow(invalid)
                invalid_rows += 1
                continue

            chrom, start, end = parsed
            clean_row = (circRNA_id, chrom, start, end, strand, gene_symbol)
            datasets_by_clean_key[clean_row].add(dataset_id)
            strand_counts[strand] += 1
            if clean_row in seen:
                duplicate_full_rows += 1
                continue
            seen.add(clean_row)
            writer.writerow(clean_row)

    with dataset_audit_output.open("w", newline="") as audit_handle:
        audit_writer = csv.writer(audit_handle, delimiter="\t", lineterminator="\n")
        audit_writer.writerow(OUTPUT_COLUMNS + ["dataset_ids", "dataset_count"])
        for clean_row in sorted(datasets_by_clean_key, key=lambda r: (r[1], int(r[2]), int(r[3]), r[4], r[0], r[5])):
            dataset_ids = sorted(x for x in datasets_by_clean_key[clean_row] if x != "NA")
            audit_writer.writerow(list(clean_row) + [",".join(dataset_ids) if dataset_ids else "NA", len(dataset_ids)])

    summary_rows = [
        ("raw_rows", raw_rows),
        ("clean_rows", len(seen)),
        ("invalid_rows", invalid_rows),
        ("duplicate_full_rows_removed", duplicate_full_rows),
        ("strand_plus_raw_valid_rows", strand_counts.get("+", 0)),
        ("strand_minus_raw_valid_rows", strand_counts.get("-", 0)),
        ("clean_rows_with_multiple_datasets", sum(1 for datasets in datasets_by_clean_key.values() if len(datasets) > 1)),
        ("max_dataset_count_for_clean_row", max((len(datasets) for datasets in datasets_by_clean_key.values()), default=0)),
    ]
    with summary_output.open("w", newline="") as summary_handle:
        summary_writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        summary_writer.writerow(["metric", "value"])
        summary_writer.writerows(summary_rows)

    for key, value in summary_rows:
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
