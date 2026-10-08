#!/usr/bin/env python3
"""Write a standard six-column final table from a postprocessed circRNA table."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


FINAL_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract the standard six columns from a headered TSV and write a circRNA final table.",
    )
    parser.add_argument("--input", required=True, type=Path, help="Input headered TSV.")
    parser.add_argument("--output", required=True, type=Path, help="Output six-column final TSV.")
    return parser.parse_args()


def write_final(input_path: Path, output_path: Path) -> int:
    input_path = input_path.expanduser().resolve()
    output_path = output_path.expanduser().resolve()

    if input_path == output_path:
        raise ValueError("--input and --output must not refer to the same file")
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None

    try:
        with input_path.open("r", encoding="utf-8", newline="") as input_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Input file has no header: {input_path}")

            missing = [column for column in FINAL_COLUMNS if column not in reader.fieldnames]
            if missing:
                raise ValueError(f"Input file is missing required columns: {', '.join(missing)}")

            with tempfile.NamedTemporaryFile(
                "w",
                encoding="utf-8",
                newline="",
                dir=output_path.parent,
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                delete=False,
            ) as output_handle:
                temporary_path = Path(output_handle.name)
                writer = csv.DictWriter(
                    output_handle,
                    fieldnames=FINAL_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="ignore",
                )
                writer.writeheader()

                row_count = 0
                for line_number, row in enumerate(reader, start=2):
                    missing_values = [column for column in FINAL_COLUMNS if row.get(column) is None]
                    if missing_values:
                        raise ValueError(
                            f"Row {line_number} is missing values for: {', '.join(missing_values)}"
                        )
                    writer.writerow({column: row[column] for column in FINAL_COLUMNS})
                    row_count += 1

        os.replace(temporary_path, output_path)
        temporary_path = None
        return row_count
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()
    row_count = write_final(args.input, args.output)
    print(f"output_rows\t{row_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
