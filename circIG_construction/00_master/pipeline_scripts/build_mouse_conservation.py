#!/usr/bin/env python3
"""Annotate the human presence table with lifted mouse circRNAs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


NULL_STRINGS = {"", "nan", "na", "none", "null", "<na>"}
HUMAN_COLUMNS = ("chrom", "start", "end", "strand")
MOUSE_COLUMNS = (
    "chrom", "start", "end", "strand", "original_chrom", "original_start",
    "original_end", "original_strand",
)
OUTPUT_COLUMNS = (*HUMAN_COLUMNS, "circ_id", "mouse_circ_id")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--human", required=True, type=Path)
    parser.add_argument("--mouse-lifted", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args()


def normalize_strand(value: str | None) -> str:
    text = (value or "").strip()
    return "NA" if text.lower() in NULL_STRINGS else text


def require_columns(path: Path, columns: list[str] | None, required: tuple[str, ...]) -> None:
    missing = [column for column in required if column not in (columns or [])]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def coordinate_id(chrom: str, start: str, end: str, strand: str) -> str:
    return f"{chrom}:{start}|{end}({strand})"


def load_mouse(path: Path, counts: Counter[str]) -> dict[str, str]:
    grouped: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames, MOUSE_COLUMNS)
        for row in reader:
            counts["mouse_lifted_rows"] += 1
            strand = normalize_strand(row["strand"])
            original_strand = normalize_strand(row["original_strand"])
            counts["mouse_NA_strand_rows"] += strand == "NA"
            counts["mouse_NA_original_strand_rows"] += original_strand == "NA"
            human_id = coordinate_id(row["chrom"], row["start"], row["end"], strand)
            mouse_id = coordinate_id(
                row["original_chrom"], row["original_start"], row["original_end"],
                original_strand,
            )
            pair = (human_id, mouse_id)
            if pair not in seen:
                seen.add(pair)
                grouped[human_id].append(mouse_id)
    counts["mouse_grouped_circ_ids"] = len(grouped)
    return {key: ",".join(values) for key, values in grouped.items()}


def build(human: Path, mouse: Path, output: Path, summary: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    mouse_by_human = load_mouse(mouse, counts)
    output.parent.mkdir(parents=True, exist_ok=True)
    with human.open(newline="") as source, output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        require_columns(human, reader.fieldnames, HUMAN_COLUMNS)
        writer = csv.DictWriter(target, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            counts["human_rows"] += 1
            strand = normalize_strand(row["strand"])
            counts["human_NA_strand_rows"] += strand == "NA"
            circ_id = coordinate_id(row["chrom"], row["start"], row["end"], strand)
            mouse_ids = mouse_by_human.get(circ_id, "")
            counts["conserved_human_rows"] += bool(mouse_ids)
            writer.writerow(
                {
                    "chrom": row["chrom"], "start": row["start"], "end": row["end"],
                    "strand": strand, "circ_id": circ_id, "mouse_circ_id": mouse_ids,
                }
            )
    counts["output_rows"] = counts["human_rows"]
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows((key, value) for key, value in counts.items())
    return counts


def main() -> int:
    args = parse_args()
    counts = build(
        args.human.expanduser().resolve(), args.mouse_lifted.expanduser().resolve(),
        args.output.expanduser().resolve(), args.summary_output.expanduser().resolve(),
    )
    for key, value in counts.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
