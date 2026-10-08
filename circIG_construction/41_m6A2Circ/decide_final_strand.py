#!/usr/bin/env python3
"""Decide final m6A2Circ strand and write audit outputs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "processed" / "41_m6A2Circ.circRNAs.with_strand_presence_master.tsv"
DEFAULT_OUTPUT = BASE_DIR / "processed" / "41_m6A2Circ.circRNAs.final.tsv"
DEFAULT_AUDIT_OUTPUT = BASE_DIR / "processed" / "41_m6A2Circ.circRNAs.final.audit.tsv"
DEFAULT_CONFLICT_OUTPUT = BASE_DIR / "processed" / "41_m6A2Circ.circRNAs.strand_conflicts.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "processed" / "41_m6A2Circ.circRNAs.final.summary.tsv"

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--audit-output", type=Path, default=DEFAULT_AUDIT_OUTPUT)
    parser.add_argument("--conflict-output", type=Path, default=DEFAULT_CONFLICT_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    return parser.parse_args()


def clean_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return value if value in VALID_STRANDS else "NA"


def decide_strand(row: dict[str, str]) -> tuple[str, str]:
    raw = clean_strand(row.get("strand"))
    presence = clean_strand(row.get("strand_presence_master"))
    if raw in VALID_STRANDS:
        return raw, "raw"
    if presence in VALID_STRANDS:
        return presence, "presence_master"
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


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    for attr in ("input", "output", "audit_output", "conflict_output", "summary_output"):
        setattr(args, attr, getattr(args, attr).expanduser().resolve())

    counts: Counter[str] = Counter()
    final_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    conflict_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{args.input}: missing header")
        missing = [column for column in OUTPUT_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"{args.input}: missing required columns: {', '.join(missing)}")
        header = list(reader.fieldnames)

        for row in reader:
            counts["input_rows"] += 1
            final_strand, method = decide_strand(row)
            final_row = {column: row.get(column, "NA") or "NA" for column in OUTPUT_COLUMNS}
            final_row["strand"] = final_strand
            key = tuple(final_row[column] for column in OUTPUT_COLUMNS)
            if key in seen:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen.add(key)
            final_rows.append(final_row)
            audit_rows.append(build_audit_row(row, final_strand, method))
            counts["final_rows"] += 1
            counts[f"method_{method}"] += 1
            if final_strand in VALID_STRANDS:
                counts["final_stranded_rows"] += 1
            else:
                counts["final_unstranded_rows"] += 1
            if method == "conflict":
                conflict_rows.append({**row, "final_strand": final_strand, "strand_decision_method": method})

    write_tsv(args.output, OUTPUT_COLUMNS, final_rows)
    write_tsv(args.audit_output, AUDIT_COLUMNS, audit_rows)
    write_tsv(args.conflict_output, ["final_strand", "strand_decision_method", *header], conflict_rows)
    summary = {
        "input_file": str(args.input),
        "output_file": str(args.output),
        "audit_output": str(args.audit_output),
        "input_rows": str(counts["input_rows"]),
        "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
        "final_rows": str(counts["final_rows"]),
        "final_stranded_rows": str(counts["final_stranded_rows"]),
        "final_unstranded_rows": str(counts["final_unstranded_rows"]),
        "raw_rows": str(counts["method_raw"]),
        "presence_master_rows": str(counts["method_presence_master"]),
        "conflict_rows": str(counts["method_conflict"]),
        "unresolved_rows": str(counts["method_unresolved"]),
        "conflict_output": str(args.conflict_output),
    }
    write_tsv(args.summary_output, list(summary), [summary])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
