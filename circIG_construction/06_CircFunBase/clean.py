#!/usr/bin/env python3
"""Clean CircFunBase human API results into the unified pipeline format."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any


csv.field_size_limit(sys.maxsize)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_API_RESULTS = SCRIPT_DIR / "raw" / "human_api_results.jsonl"
DEFAULT_OUTPUT = SCRIPT_DIR / "clean" / "06_CircFunBase.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = SCRIPT_DIR / "clean" / "06_CircFunBase.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = SCRIPT_DIR / "clean" / "06_CircFunBase.circRNAs.clean.summary.tsv"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["circRNA_id", "species", "position", "gene_symbol", "reason"]
SUMMARY_COLUMNS = [
    "input",
    "output",
    "raw_api_rows",
    "non_human_rows",
    "invalid_coordinate_rows",
    "duplicate_clean_rows",
    "clean_output_rows",
    "unique_circRNA_ids",
    "raw_stranded_rows",
    "raw_unstranded_rows",
]

POSITION_RE = re.compile(r"^(chr[^:\s]+):(\d+)-(\d+)$")
NULL_VALUES = {"", "-", ".", "*", "NA", "N/A", "None", "none", "null", "NULL"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse CircFunBase CircRNA.getinfo JSONL into the unified six-column "
            "circRNA format. Coordinates are kept in the original raw coordinate system."
        )
    )
    parser.add_argument(
        "--api-results",
        type=Path,
        default=DEFAULT_API_RESULTS,
        help=f"Raw API JSONL. Default: {DEFAULT_API_RESULTS}",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--invalid-output", type=Path, default=DEFAULT_INVALID_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def clean(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    if text in NULL_VALUES:
        return "NA"
    return text


def normalize_strand(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    return text if text in {"+", "-"} else "NA"


def get_info(record: dict[str, Any]) -> dict[str, Any]:
    response = record.get("response")
    if not isinstance(response, dict):
        return {}
    data = response.get("data")
    if not isinstance(data, dict):
        return {}
    data2 = data.get("data")
    if not isinstance(data2, dict):
        return {}
    info = data2.get("Info")
    return info if isinstance(info, dict) else {}


def parse_position(position: str) -> tuple[str, str, str] | None:
    match = POSITION_RE.match(position)
    if not match:
        return None
    chrom, start, end = match.groups()
    start_i = int(start)
    end_i = int(end)
    if start_i < 0 or end_i <= start_i:
        return None
    return chrom, str(start_i), str(end_i)


def load_clean_rows(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
    clean_rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    counts: Counter[str] = Counter()

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            counts["raw_api_rows"] += 1
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                counts["invalid_coordinate_rows"] += 1
                invalid_rows.append(
                    {
                        "circRNA_id": "NA",
                        "species": "NA",
                        "position": "NA",
                        "gene_symbol": "NA",
                        "reason": f"invalid_json_line_{line_number}: {exc}",
                    }
                )
                continue

            info = get_info(record)
            circRNA_id = clean(info.get("circRNA") or record.get("circRNA"))
            species = clean(info.get("Species"))
            position = clean(info.get("Position"))
            gene_symbol = clean(info.get("Gene"))
            strand = normalize_strand(info.get("Strand"))

            if strand in {"+", "-"}:
                counts["raw_stranded_rows"] += 1
            else:
                counts["raw_unstranded_rows"] += 1

            if species != "Homo sapiens":
                counts["non_human_rows"] += 1
                invalid_rows.append(
                    {
                        "circRNA_id": circRNA_id,
                        "species": species,
                        "position": position,
                        "gene_symbol": gene_symbol,
                        "reason": "non_human_or_missing_species",
                    }
                )
                continue

            parsed = parse_position(position)
            if parsed is None:
                counts["invalid_coordinate_rows"] += 1
                invalid_rows.append(
                    {
                        "circRNA_id": circRNA_id,
                        "species": species,
                        "position": position,
                        "gene_symbol": gene_symbol,
                        "reason": "unparseable_position",
                    }
                )
                continue

            chrom, start, end = parsed
            row = {
                "circRNA_id": circRNA_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_symbol": gene_symbol,
            }
            key = tuple(row[column] for column in OUTPUT_COLUMNS)
            if key in seen:
                counts["duplicate_clean_rows"] += 1
                continue
            seen.add(key)
            clean_rows.append(row)

    counts["clean_output_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len({row["circRNA_id"] for row in clean_rows})
    return clean_rows, invalid_rows, counts


def refuse_overwrite(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Output already exists; use --force to overwrite: {path}")


def write_table(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    row = {
        "input": str(args.api_results),
        "output": str(args.output),
        "raw_api_rows": str(counts["raw_api_rows"]),
        "non_human_rows": str(counts["non_human_rows"]),
        "invalid_coordinate_rows": str(counts["invalid_coordinate_rows"]),
        "duplicate_clean_rows": str(counts["duplicate_clean_rows"]),
        "clean_output_rows": str(counts["clean_output_rows"]),
        "unique_circRNA_ids": str(counts["unique_circRNA_ids"]),
        "raw_stranded_rows": str(counts["raw_stranded_rows"]),
        "raw_unstranded_rows": str(counts["raw_unstranded_rows"]),
    }
    write_table(path, [row], SUMMARY_COLUMNS)


def main() -> int:
    args = parse_args()
    args.api_results = args.api_results.expanduser().resolve()
    args.output = args.output.expanduser().resolve()
    args.invalid_output = args.invalid_output.expanduser().resolve()
    args.summary_output = args.summary_output.expanduser().resolve()

    try:
        if not args.api_results.exists() or args.api_results.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {args.api_results}")
        for path in (args.output, args.invalid_output, args.summary_output):
            refuse_overwrite(path, args.force)
        clean_rows, invalid_rows, counts = load_clean_rows(args.api_results)
        write_table(args.output, clean_rows, OUTPUT_COLUMNS)
        write_table(args.invalid_output, invalid_rows, INVALID_COLUMNS)
        write_summary(args.summary_output, args, counts)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote clean table: {args.output}")
    print(f"Wrote invalid/audit table: {args.invalid_output}")
    print(f"Wrote summary: {args.summary_output}")
    print(
        "Summary: "
        f"{counts['raw_api_rows']} raw API rows, "
        f"{counts['clean_output_rows']} clean rows, "
        f"{counts['invalid_coordinate_rows']} invalid coordinate rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
