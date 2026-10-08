#!/usr/bin/env python3
"""Clean Circ2Traits raw circRNA table into the unified six-column format."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW = BASE_DIR / "raw" / "circRNA_info_all.txt"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "38_Circ2Traits.circRNAs.clean.tsv"
DEFAULT_INVALID_OUTPUT = BASE_DIR / "clean" / "38_Circ2Traits.circRNAs.clean.invalid.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "38_Circ2Traits.circRNAs.clean.summary.tsv"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
INVALID_COLUMNS = ["line_number", "invalid_reason", "raw_row"]
MISSING_VALUES = {"", "NA", "N/A", "NULL", "NONE", "NAN"}
GENE_MISSING_VALUES = {"", "NA", "N/A", "NULL", "NAN"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Clean Circ2Traits circRNA coordinate table.")
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW, help=f"Raw input table. Default: {DEFAULT_RAW}")
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
    if value.upper() in MISSING_VALUES:
        return "NA"
    return value


def normalize_gene_symbol(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if value == "None":
        return "intergenic"
    if value.upper() in GENE_MISSING_VALUES:
        return "NA"
    return value


def chrom_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr")
    if value.isdigit():
        return (int(value), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (special.get(value, 1000), value)


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; use --force to overwrite: {names}")


def convert_row(row: dict[str, str]) -> tuple[str, str, str, str, str, str]:
    circ_id = normalize_missing(row.get("name"))
    chrom = normalize_missing(row.get("#chrom"))
    start = normalize_missing(row.get("start"))
    end = normalize_missing(row.get("end"))
    strand = normalize_missing(row.get("strand"))
    gene_symbol = normalize_gene_symbol(row.get("gene"))

    if circ_id == "NA":
        raise ValueError("missing circRNA_id")
    if chrom == "NA":
        raise ValueError("missing chrom")
    if strand not in {"+", "-"}:
        raise ValueError(f"invalid strand: {strand!r}")

    try:
        start_int = int(start)
        end_int = int(end)
    except ValueError as exc:
        raise ValueError(f"invalid coordinates: start={start!r}, end={end!r}") from exc
    if start_int < 0 or end_int < 0:
        raise ValueError(f"negative coordinate: start={start_int}, end={end_int}")
    if start_int > end_int:
        raise ValueError(f"start greater than end: start={start_int}, end={end_int}")

    return circ_id, chrom, str(start_int), str(end_int), strand, gene_symbol


def write_summary(path: Path, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in [
            "raw_rows",
            "clean_rows_before_uniq",
            "duplicate_full_rows_removed",
            "clean_rows",
            "invalid_rows",
            "unique_circRNA_ids",
            "stranded_rows",
            "unstranded_rows",
            "gene_symbol_na_rows",
            "unique_gene_symbols",
            "strand:+",
            "strand:-",
        ]:
            writer.writerow([key, counts[key]])


def main() -> int:
    args = parse_args()
    raw = args.raw.expanduser().resolve()
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

    with raw.open(newline="") as raw_handle, invalid_output.open("w", newline="") as invalid_handle:
        reader = csv.DictReader(raw_handle, delimiter="\t")
        required = {"#chrom", "start", "end", "name", "strand", "gene"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{raw} missing required columns: {', '.join(sorted(missing))}")

        invalid_writer = csv.DictWriter(invalid_handle, fieldnames=INVALID_COLUMNS, delimiter="\t", lineterminator="\n")
        invalid_writer.writeheader()

        for line_number, row in enumerate(reader, start=2):
            counts["raw_rows"] += 1
            try:
                clean_row = convert_row(row)
            except ValueError as exc:
                counts["invalid_rows"] += 1
                invalid_writer.writerow(
                    {
                        "line_number": line_number,
                        "invalid_reason": str(exc),
                        "raw_row": repr(row),
                    }
                )
                continue

            counts["clean_rows_before_uniq"] += 1
            if clean_row in seen:
                counts["duplicate_full_rows_removed"] += 1
                continue
            seen.add(clean_row)
            clean_rows.append(clean_row)
            circ_ids.add(clean_row[0])
            if clean_row[4] in {"+", "-"}:
                counts["stranded_rows"] += 1
                counts[f"strand:{clean_row[4]}"] += 1
            else:
                counts["unstranded_rows"] += 1
            if clean_row[5] == "NA":
                counts["gene_symbol_na_rows"] += 1
            else:
                gene_symbols.add(clean_row[5])

    clean_rows.sort(key=lambda row: (chrom_key(row[1]), int(row[2]), int(row[3]), row[4], row[0]))
    counts["clean_rows"] = len(clean_rows)
    counts["unique_circRNA_ids"] = len(circ_ids)
    counts["unique_gene_symbols"] = len(gene_symbols)

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        writer.writerows(clean_rows)

    write_summary(summary_output, counts)
    print(f"clean_rows\t{counts['clean_rows']}")
    print(f"invalid_rows\t{counts['invalid_rows']}")
    print(f"output\t{output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
