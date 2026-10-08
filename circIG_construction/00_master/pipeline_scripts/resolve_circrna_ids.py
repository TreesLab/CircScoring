#!/usr/bin/env python3
"""Resolve circRNA IDs against one or more standardized reference tables."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
AUDIT_COLUMNS = ["circRNA_id", "issue", "candidate_count", "coordinate_sources", "coordinates"]


@dataclass(frozen=True)
class Reference:
    name: str
    path: Path


def parse_reference(value: str) -> Reference:
    if "=" not in value:
        raise argparse.ArgumentTypeError("reference must use NAME=PATH")
    name, raw_path = value.split("=", 1)
    name = name.strip()
    raw_path = raw_path.strip()
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("reference must use non-empty NAME=PATH")
    return Reference(name, Path(raw_path))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="ID-only TSV containing circRNA_id.")
    parser.add_argument(
        "--reference",
        type=parse_reference,
        action="append",
        required=True,
        help="Standardized reference as NAME=PATH; may be repeated.",
    )
    parser.add_argument("--output", type=Path, required=True, help="Resolved six-column TSV.")
    parser.add_argument("--audit", type=Path, required=True, help="Unmatched and ambiguous ID audit TSV.")
    parser.add_argument("--summary-output", type=Path, required=True, help="Resolution summary TSV.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def normalize(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return text if text else "NA"


def read_ids(path: Path) -> tuple[set[str], Counter[str]]:
    counts: Counter[str] = Counter()
    circ_ids: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "circRNA_id" not in (reader.fieldnames or []):
            raise ValueError(f"{path} missing required column: circRNA_id")
        for row in reader:
            counts["input_rows"] += 1
            circ_id = normalize(row.get("circRNA_id"))
            if circ_id == "NA":
                counts["invalid_input_rows"] += 1
                continue
            if circ_id in circ_ids:
                counts["duplicate_input_ids"] += 1
            circ_ids.add(circ_id)
    counts["unique_input_ids"] = len(circ_ids)
    return circ_ids, counts


def load_references(
    references: list[Reference],
) -> tuple[dict[str, set[tuple[str, str, str, str, str]]], dict[str, set[str]]]:
    candidates: dict[str, set[tuple[str, str, str, str, str]]] = defaultdict(set)
    sources: dict[str, set[str]] = defaultdict(set)
    seen_names: set[str] = set()

    for reference in references:
        if reference.name in seen_names:
            raise ValueError(f"duplicate reference name: {reference.name}")
        seen_names.add(reference.name)
        with reference.path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = set(COLUMNS) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{reference.path} missing required columns: {', '.join(sorted(missing))}")
            for row in reader:
                circ_id = normalize(row.get("circRNA_id"))
                if circ_id == "NA":
                    continue
                candidate = tuple(normalize(row.get(column)) for column in COLUMNS[1:])
                if "NA" in candidate[:4]:
                    raise ValueError(f"{reference.path} contains incomplete coordinates for {circ_id}")
                try:
                    int(candidate[1])
                    int(candidate[2])
                except ValueError as exc:
                    raise ValueError(
                        f"{reference.path} contains non-integer coordinates for {circ_id}: {candidate[1:3]}"
                    ) from exc
                candidates[circ_id].add(candidate)
                sources[circ_id].add(reference.name)
    return candidates, sources


def chrom_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    if value.isdigit():
        return int(value), ""
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return special.get(value, 1000), value


def sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            chrom_key(row["chrom"]),
            int(row["start"]),
            int(row["end"]),
            row["strand"],
            row["circRNA_id"],
        ),
    )


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def resolve_ids(
    input_path: Path,
    references: list[Reference],
    output_path: Path,
    audit_path: Path,
    summary_path: Path,
    force: bool,
) -> Counter[str]:
    existing = [str(path) for path in (output_path, audit_path, summary_path) if path.exists()]
    if existing and not force:
        raise FileExistsError("Output exists; use --force to overwrite: " + ", ".join(existing))

    circ_ids, counts = read_ids(input_path)
    candidates, sources = load_references(references)
    resolved_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []

    for circ_id in sorted(circ_ids):
        matches = candidates.get(circ_id, set())
        source_text = ",".join(sorted(sources.get(circ_id, set()))) or "NA"
        if len(matches) == 1:
            chrom, start, end, strand, gene_symbol = next(iter(matches))
            resolved_rows.append(
                {
                    "circRNA_id": circ_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_symbol": gene_symbol,
                }
            )
            counts["resolved_ids"] += 1
            continue

        issue = "unmatched" if not matches else "ambiguous_reference"
        counts[f"{issue}_ids"] += 1
        coordinate_text = ";".join("|".join(candidate) for candidate in sorted(matches)) or "NA"
        audit_rows.append(
            {
                "circRNA_id": circ_id,
                "issue": issue,
                "candidate_count": str(len(matches)),
                "coordinate_sources": source_text,
                "coordinates": coordinate_text,
            }
        )

    resolved_rows = sort_rows(resolved_rows)
    counts["output_rows"] = len(resolved_rows)
    counts["audit_rows"] = len(audit_rows)
    write_tsv(output_path, COLUMNS, resolved_rows)
    write_tsv(audit_path, AUDIT_COLUMNS, audit_rows)
    write_tsv(
        summary_path,
        ["metric", "value"],
        [{"metric": key, "value": str(value)} for key, value in counts.items()],
    )
    return counts


def main() -> int:
    args = parse_args()
    references = [Reference(item.name, item.path.expanduser().resolve()) for item in args.reference]
    counts = resolve_ids(
        args.input.expanduser().resolve(),
        references,
        args.output.expanduser().resolve(),
        args.audit.expanduser().resolve(),
        args.summary_output.expanduser().resolve(),
        args.force,
    )
    for key, value in counts.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
