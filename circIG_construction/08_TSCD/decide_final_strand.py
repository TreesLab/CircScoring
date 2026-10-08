#!/usr/bin/env python3
"""Decide final TSCD strand and return the unified six-column table."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "processed" / "08_TSCD.circRNAs.with_strand_presence_master.tsv"
DEFAULT_OUTPUT = BASE_DIR / "processed" / "08_TSCD.circRNAs.final.tsv"
DEFAULT_AUDIT_OUTPUT = BASE_DIR / "processed" / "08_TSCD.circRNAs.final.audit.tsv"
DEFAULT_CONFLICT_OUTPUT = BASE_DIR / "processed" / "08_TSCD.circRNAs.strand_conflicts.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "processed" / "08_TSCD.circRNAs.final.summary.tsv"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
AUDIT_COLUMNS = [
    "circRNA_id",
    "chrom",
    "start",
    "end",
    "gene_symbol",
    "original_strand",
    "strand_presence_master",
    "final_strand",
    "strand_decision_method",
]
VALID_STRANDS = {"+", "-"}


def clean_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return value if value in VALID_STRANDS else "NA"


def decide(row: dict[str, str]) -> tuple[str, str]:
    original = clean_strand(row.get("strand"))
    presence_master = clean_strand(row.get("strand_presence_master"))

    if original in VALID_STRANDS:
        return original, "original"

    if presence_master in VALID_STRANDS:
        return presence_master, "presence_master"
    return "NA", "unresolved"


def build_audit_row(row: dict[str, str], final_strand: str, method: str) -> dict[str, str]:
    return {
        "circRNA_id": row.get("circRNA_id", "NA") or "NA",
        "chrom": row.get("chrom", "NA") or "NA",
        "start": row.get("start", "NA") or "NA",
        "end": row.get("end", "NA") or "NA",
        "gene_symbol": row.get("gene_symbol", "NA") or "NA",
        "original_strand": clean_strand(row.get("strand")),
        "strand_presence_master": clean_strand(row.get("strand_presence_master")),
        "final_strand": final_strand,
        "strand_decision_method": method,
    }


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    columns = [
        "input_file",
        "output_file",
        "audit_output",
        "total_rows",
        "final_rows",
        "final_stranded_rows",
        "final_unstranded_rows",
        "original_rows",
        "presence_master_rows",
        "conflict_rows",
        "unresolved_rows",
    ]
    row = {
        "input_file": str(args.input),
        "output_file": str(args.output),
        "audit_output": str(args.audit_output),
        "total_rows": str(counts["total_rows"]),
        "final_rows": str(counts["final_rows"]),
        "final_stranded_rows": str(counts["final_stranded_rows"]),
        "final_unstranded_rows": str(counts["final_unstranded_rows"]),
        "original_rows": str(counts["method_original"]),
        "presence_master_rows": str(counts["method_presence_master"]),
        "conflict_rows": str(counts["method_conflict"]),
        "unresolved_rows": str(counts["method_unresolved"]),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Choose final TSCD strand from original/presence master evidence.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Input TSV. Default: {DEFAULT_INPUT}")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"Final TSV. Default: {DEFAULT_OUTPUT}")
    parser.add_argument("--audit-output", type=Path, default=DEFAULT_AUDIT_OUTPUT, help=f"Final audit TSV. Default: {DEFAULT_AUDIT_OUTPUT}")
    parser.add_argument("--conflict-output", type=Path, default=DEFAULT_CONFLICT_OUTPUT, help=f"Conflict audit TSV. Default: {DEFAULT_CONFLICT_OUTPUT}")
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT, help=f"Summary TSV. Default: {DEFAULT_SUMMARY_OUTPUT}")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.input = args.input.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.audit_output = args.audit_output.expanduser().resolve()
    args.conflict_output = args.conflict_output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.conflict_output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()

    with args.input.open(newline="") as input_handle, args.output.open("w", newline="") as output_handle, args.audit_output.open("w", newline="") as audit_handle, args.conflict_output.open("w", newline="") as conflict_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {args.input}")
        missing = [column for column in OUTPUT_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"{args.input} is missing required column(s): {', '.join(missing)}")

        writer = csv.DictWriter(output_handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        audit_writer = csv.DictWriter(audit_handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        conflict_columns = reader.fieldnames + ["final_strand", "strand_decision_method"]
        conflict_writer = csv.DictWriter(conflict_handle, fieldnames=conflict_columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        audit_writer.writeheader()
        conflict_writer.writeheader()

        for row in reader:
            counts["total_rows"] += 1
            final_strand, method = decide(row)
            output_row = {column: row.get(column, "NA") for column in OUTPUT_COLUMNS}
            output_row["strand"] = final_strand
            writer.writerow(output_row)
            audit_writer.writerow(build_audit_row(row, final_strand, method))
            counts["final_rows"] += 1
            counts[f"method_{method}"] += 1
            if final_strand in VALID_STRANDS:
                counts["final_stranded_rows"] += 1
            else:
                counts["final_unstranded_rows"] += 1
            if method == "conflict":
                conflict_row = dict(row)
                conflict_row["final_strand"] = final_strand
                conflict_row["strand_decision_method"] = method
                conflict_writer.writerow(conflict_row)

    write_summary(args.summary_output, args, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
