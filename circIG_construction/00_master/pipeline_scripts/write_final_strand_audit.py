#!/usr/bin/env python3
"""Write final strand decision audit from an evidence table.

This is a lightweight backfill/validation helper. It replays the same final
strand decision rule used by the current per-database pipelines and checks
that the reconstructed six-column final table matches the existing final.tsv.
Gene-annotation strand is retained for audit only and never participates in
the final decision.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
AUDIT_COLUMNS = [
    "circRNA_id",
    "chrom",
    "start",
    "end",
    "gene_symbol",
    "original_strand",
    "strand_gene_annotation",
    "strand_presence_master",
    "final_strand",
    "strand_decision_method",
]
VALID_STRANDS = {"+", "-"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--conflict-output", type=Path)
    parser.add_argument(
        "--raw-method-label",
        choices=["raw", "original"],
        default="raw",
        help="Decision-method label to use when the input strand is already valid.",
    )
    parser.add_argument(
        "--no-deduplicate",
        action="store_true",
        help="Do not remove duplicate final six-column rows while replaying decisions.",
    )
    return parser.parse_args()


def clean_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return value if value in VALID_STRANDS else "NA"


def decide_strand(row: dict[str, str], raw_method_label: str) -> tuple[str, str]:
    raw = clean_strand(row.get("strand"))
    presence = clean_strand(row.get("strand_presence_master"))
    if raw in VALID_STRANDS:
        return raw, raw_method_label
    if presence in VALID_STRANDS:
        return presence, "presence_master"
    return "NA", "unresolved"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return list(reader.fieldnames), list(reader)


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


def build_audit_row(row: dict[str, str], final_strand: str, method: str) -> dict[str, str]:
    return {
        "circRNA_id": row.get("circRNA_id", "NA") or "NA",
        "chrom": row.get("chrom", "NA") or "NA",
        "start": row.get("start", "NA") or "NA",
        "end": row.get("end", "NA") or "NA",
        "gene_symbol": row.get("gene_symbol", "NA") or "NA",
        "original_strand": clean_strand(row.get("strand")),
        "strand_gene_annotation": clean_strand(row.get("strand_gene_annotation")),
        "strand_presence_master": clean_strand(row.get("strand_presence_master")),
        "final_strand": final_strand,
        "strand_decision_method": method,
    }


def main() -> int:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    final_path = args.final.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    header, rows = read_tsv(input_path)
    missing = [column for column in OUTPUT_COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{input_path}: missing required columns: {', '.join(missing)}")

    counts: Counter[str] = Counter()
    reconstructed_final: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    conflict_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()

    for row in rows:
        counts["input_rows"] += 1
        final_strand, method = decide_strand(row, args.raw_method_label)
        final_row = {column: row.get(column, "NA") or "NA" for column in OUTPUT_COLUMNS}
        final_row["strand"] = final_strand
        key = tuple(final_row[column] for column in OUTPUT_COLUMNS)
        if not args.no_deduplicate and key in seen:
            counts["duplicate_full_rows_removed"] += 1
            continue
        seen.add(key)
        reconstructed_final.append(final_row)
        audit_rows.append(build_audit_row(row, final_strand, method))
        counts["final_rows"] += 1
        counts[f"method_{method}"] += 1
        if final_strand in VALID_STRANDS:
            counts["final_stranded_rows"] += 1
        else:
            counts["final_unstranded_rows"] += 1
        if method == "conflict":
            conflict_rows.append({**row, "final_strand": final_strand, "strand_decision_method": method})

    final_header, final_rows = read_tsv(final_path)
    final_projection = [
        {column: row.get(column, "NA") or "NA" for column in OUTPUT_COLUMNS}
        for row in final_rows
    ]
    if final_header[: len(OUTPUT_COLUMNS)] != OUTPUT_COLUMNS:
        raise ValueError(f"{final_path}: unexpected final table columns")
    if reconstructed_final != final_projection:
        raise ValueError(
            f"Reconstructed final table does not match {final_path}: "
            f"reconstructed={len(reconstructed_final)} existing={len(final_projection)}"
        )

    write_tsv(output_path, AUDIT_COLUMNS, audit_rows)

    if args.conflict_output:
        conflict_path = args.conflict_output.expanduser().resolve()
        conflict_header = ["final_strand", "strand_decision_method", *header]
        write_tsv(conflict_path, conflict_header, conflict_rows)

    if args.summary_output:
        summary_path = args.summary_output.expanduser().resolve()
        summary = {
            "input_file": str(input_path),
            "final_file": str(final_path),
            "audit_output": str(output_path),
            "input_rows": str(counts["input_rows"]),
            "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
            "final_rows": str(counts["final_rows"]),
            "final_stranded_rows": str(counts["final_stranded_rows"]),
            "final_unstranded_rows": str(counts["final_unstranded_rows"]),
            "raw_or_original_rows": str(counts["method_raw"] + counts["method_original"]),
            "presence_master_rows": str(counts["method_presence_master"]),
            "unresolved_rows": str(counts["method_unresolved"]),
        }
        write_tsv(summary_path, list(summary), [summary])

    print(f"audit_rows\t{len(audit_rows)}")
    print(f"conflict_rows\t{counts['method_conflict']}")
    print(f"unresolved_rows\t{counts['method_unresolved']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
