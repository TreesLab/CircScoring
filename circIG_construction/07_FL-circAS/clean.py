#!/usr/bin/env python3
"""Clean FL-circAS human BSJ table into the unified pipeline format."""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "raw" / "BSJ_human.csv"
DEFAULT_OUTPUT = SCRIPT_DIR / "clean" / "07_FL-circAS.circRNAs.clean.tsv"
DEFAULT_SUMMARY_OUTPUT = SCRIPT_DIR / "clean" / "07_FL-circAS.circRNAs.clean.summary.tsv"
DEFAULT_DUPLICATE_COORD_OUTPUT = SCRIPT_DIR / "clean" / "07_FL-circAS.duplicated_coordinates.tsv"

REQUIRED_COLUMNS = [
    "species",
    "BSJ_ID",
    "gene_name",
    "chr",
    "BSJ_start",
    "BSJ_end",
    "strand",
]
OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
SUMMARY_COLUMNS = [
    "input",
    "output",
    "raw_data_rows",
    "non_human_rows",
    "invalid_rows",
    "duplicate_clean_rows",
    "clean_output_rows",
    "unique_circRNA_ids",
    "unique_coord_strand_keys",
    "duplicated_coord_strand_keys",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean the FL-circAS Homo sapiens BSJ-level CSV into the unified "
            "six-column circRNA format. Coordinates are kept in the original "
            "FL-circAS hg38_1based system."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"FL-circAS BSJ_human.csv. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Clean TSV output path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Summary TSV output path. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument(
        "--duplicate-coord-output",
        type=Path,
        default=DEFAULT_DUPLICATE_COORD_OUTPUT,
        help=(
            "Audit TSV for coordinate+strand groups that have multiple circRNA IDs. "
            f"Default: {DEFAULT_DUPLICATE_COORD_OUTPUT}"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )
    return parser.parse_args()


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if not value or value.upper() in {"NA", "N/A", "NULL", "NONE", "UNKNOWN"}:
        return "NA"
    return value


def normalize_int(value: str | None) -> str | None:
    value = normalize_missing(value)
    if value == "NA":
        return None
    try:
        parsed = int(value)
    except ValueError:
        return None
    if parsed < 1:
        return None
    return str(parsed)


def validate_input(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input file: {path}")


def validate_header(path: Path, fieldnames: list[str] | None) -> None:
    if fieldnames is None:
        raise ValueError(f"Input table has no header: {path}")
    missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")


def clean_row(row: dict[str, str]) -> dict[str, str] | None:
    species = normalize_missing(row.get("species")).lower()
    if species != "human":
        return None

    circRNA_id = normalize_missing(row.get("BSJ_ID"))
    chrom = normalize_missing(row.get("chr"))
    start = normalize_int(row.get("BSJ_start"))
    end = normalize_int(row.get("BSJ_end"))
    strand = normalize_missing(row.get("strand"))
    gene_symbol = normalize_missing(row.get("gene_name"))

    if (
        circRNA_id == "NA"
        or chrom == "NA"
        or start is None
        or end is None
        or strand not in {"+", "-"}
    ):
        return None

    if int(start) > int(end):
        return None

    return {
        "circRNA_id": circRNA_id,
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "gene_symbol": gene_symbol,
    }


def read_clean_rows(
    path: Path,
) -> tuple[list[dict[str, str]], Counter[str], dict[tuple[str, str, str, str], set[str]]]:
    counts: Counter[str] = Counter()
    seen_clean_rows: set[tuple[str, ...]] = set()
    coord_to_ids: dict[tuple[str, str, str, str], set[str]] = {}
    clean_rows: list[dict[str, str]] = []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        validate_header(path, reader.fieldnames)
        for row in reader:
            counts["raw_data_rows"] += 1
            species = normalize_missing(row.get("species")).lower()
            cleaned = clean_row(row)
            if cleaned is None:
                if species != "human":
                    counts["non_human_rows"] += 1
                else:
                    counts["invalid_rows"] += 1
                continue

            clean_key = tuple(cleaned[column] for column in OUTPUT_COLUMNS)
            if clean_key in seen_clean_rows:
                counts["duplicate_clean_rows"] += 1
                continue
            seen_clean_rows.add(clean_key)
            clean_rows.append(cleaned)

            coord_key = (
                cleaned["chrom"],
                cleaned["start"],
                cleaned["end"],
                cleaned["strand"],
            )
            coord_to_ids.setdefault(coord_key, set()).add(cleaned["circRNA_id"])

    counts["clean_output_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len({row["circRNA_id"] for row in clean_rows})
    counts["unique_coord_strand_keys"] = len(coord_to_ids)
    counts["duplicated_coord_strand_keys"] = sum(
        1 for circ_ids in coord_to_ids.values() if len(circ_ids) > 1
    )
    return clean_rows, counts, coord_to_ids


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_duplicate_coord_audit(
    path: Path, coord_to_ids: dict[tuple[str, str, str, str], set[str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fieldnames = ["chrom", "start", "end", "strand", "circRNA_ids", "circRNA_id_count"]
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for coord_key, circ_ids in sorted(coord_to_ids.items()):
            if len(circ_ids) <= 1:
                continue
            chrom, start, end, strand = coord_key
            writer.writerow(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "circRNA_ids": ",".join(sorted(circ_ids)),
                    "circRNA_id_count": str(len(circ_ids)),
                }
            )


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "input": str(args.input),
        "output": str(args.output),
        "raw_data_rows": str(counts["raw_data_rows"]),
        "non_human_rows": str(counts["non_human_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "duplicate_clean_rows": str(counts["duplicate_clean_rows"]),
        "clean_output_rows": str(counts["clean_output_rows"]),
        "unique_circRNA_ids": str(counts["unique_circRNA_ids"]),
        "unique_coord_strand_keys": str(counts["unique_coord_strand_keys"]),
        "duplicated_coord_strand_keys": str(counts["duplicated_coord_strand_keys"]),
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)


def refuse_overwrite(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Output already exists; use --force to overwrite: {path}")


def main() -> int:
    args = parse_args()
    args.input = args.input.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()
    args.duplicate_coord_output = args.duplicate_coord_output.expanduser().resolve()

    try:
        validate_input(args.input)
        refuse_overwrite(args.output, args.force)
        refuse_overwrite(args.summary_output, args.force)
        refuse_overwrite(args.duplicate_coord_output, args.force)
        clean_rows, counts, coord_to_ids = read_clean_rows(args.input)
        write_table(args.output, clean_rows)
        write_duplicate_coord_audit(args.duplicate_coord_output, coord_to_ids)
        write_summary(args.summary_output, args, counts)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {counts['clean_output_rows']} clean rows to {args.output}")
    print(f"Wrote summary to {args.summary_output}")
    print(f"Wrote duplicate-coordinate audit to {args.duplicate_coord_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
