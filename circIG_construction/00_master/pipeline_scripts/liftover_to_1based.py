#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


csv.field_size_limit(sys.maxsize)

PROVENANCE_COLUMNS = [
    "original_chrom",
    "original_start",
    "original_end",
    "original_strand",
    "original_coordinate_system",
    "conversion_status",
    "conversion_note",
]
UNMAPPED_COLUMNS = ["failure_reason", "original_strand", "original_coordinate_system"]
SUMMARY_COLUMNS = [
    "input_file",
    "output_file",
    "unmapped_output_file",
    "source_coordinate_system",
    "source_assembly",
    "target_assembly",
    "total_rows",
    "converted_rows",
    "unmapped_rows",
    "invalid_rows",
    "failed_liftover_rows",
    "chain",
    "liftover_bin",
]


@dataclass(frozen=True)
class CoordinateSystem:
    assembly: str
    base: str

    @classmethod
    def parse(cls, value: str) -> "CoordinateSystem":
        if value.endswith("_1based"):
            return cls(value[: -len("_1based")], "1based")
        if value.endswith("_0based"):
            return cls(value[: -len("_0based")], "0based")
        raise ValueError(f"Coordinate system must end with _1based or _0based: {value}")


@dataclass(frozen=True)
class Coordinate:
    chrom: str
    start: int
    end: int
    strand: str


@dataclass(frozen=True)
class RowCoordinate:
    row_index: int
    row: dict[str, str]
    coord: Coordinate


def normalize_chrom(value: str | None) -> str:
    if value is None:
        return ""
    chrom = value.strip()
    if not chrom:
        return ""
    if chrom.lower().startswith("chr"):
        return "chr" + chrom[3:]
    return f"chr{chrom}"


def normalize_strand(value: str | None) -> str:
    if value is None:
        return "."
    strand = value.strip()
    return strand if strand in {"+", "-", "."} else "."


def parse_int(value: str | None, column: str) -> int:
    if value is None:
        raise ValueError(f"missing required column value: {column}")
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"column {column!r} contains a non-integer coordinate: {value!r}") from exc


def read_table(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {path}")
        return list(reader.fieldnames), rows


def validate_required_columns(header: list[str], columns: list[str]) -> None:
    missing = [column for column in columns if column and column not in header]
    if missing:
        raise ValueError(f"Input table is missing required column(s): {', '.join(missing)}")


def parse_coordinate(
    row: dict[str, str],
    chrom_col: str,
    start_col: str,
    end_col: str,
    strand_col: str | None,
    coordinate_system: CoordinateSystem,
) -> Coordinate:
    chrom = normalize_chrom(row.get(chrom_col))
    start = parse_int(row.get(start_col), start_col)
    end = parse_int(row.get(end_col), end_col)
    strand = normalize_strand(row.get(strand_col) if strand_col else None)

    if not chrom:
        raise ValueError("missing chromosome")
    if coordinate_system.base == "1based" and start > end:
        raise ValueError(f"start is greater than end: {chrom}:{start}-{end}")
    if coordinate_system.base == "0based" and start >= end:
        raise ValueError(f"0-based interval must satisfy start < end: {chrom}:{start}-{end}")
    if coordinate_system.base == "1based" and start < 1:
        raise ValueError(f"1-based start must be >= 1: {chrom}:{start}-{end}")
    if coordinate_system.base == "0based" and start < 0:
        raise ValueError(f"0-based start must be >= 0: {chrom}:{start}-{end}")
    return Coordinate(chrom, start, end, strand)


def provenance(row: dict[str, str], coord: Coordinate, source_coordinate_system: str, status: str, note: str = "NA") -> dict[str, str]:
    out = dict(row)
    out.update(
        {
            "original_chrom": coord.chrom,
            "original_start": str(coord.start),
            "original_end": str(coord.end),
            "original_strand": row.get("strand", coord.strand) or "NA",
            "original_coordinate_system": source_coordinate_system,
            "conversion_status": status,
            "conversion_note": note,
        }
    )
    return out


def failure_row(row: dict[str, str], source_coordinate_system: str, reason: str) -> dict[str, str]:
    out = dict(row)
    out.update(
        {
            "failure_reason": reason,
            "original_strand": row.get("strand", "NA") or "NA",
            "original_coordinate_system": source_coordinate_system,
        }
    )
    return out


def write_table(path: Path, header: list[str], rows: list[dict[str, str]], delimiter: str, extra_columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(header) + [column for column in extra_columns if column not in header]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, source_system: CoordinateSystem, counts: Counter[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "input_file": str(args.input),
        "output_file": str(args.output),
        "unmapped_output_file": str(args.unmapped_output),
        "source_coordinate_system": args.source_coordinate_system,
        "source_assembly": source_system.assembly,
        "target_assembly": args.target_assembly,
        "total_rows": str(counts["total_rows"]),
        "converted_rows": str(counts["converted_rows"]),
        "unmapped_rows": str(counts["unmapped_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "failed_liftover_rows": str(counts["failed_liftover_rows"]),
        "chain": str(args.chain) if args.chain else "NA",
        "liftover_bin": args.liftover_bin,
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def boundary_intervals(coord: Coordinate, coordinate_system: CoordinateSystem) -> tuple[tuple[int, int], tuple[int, int]]:
    if coordinate_system.base == "1based":
        start_interval = (coord.start - 1, coord.start)
        end_interval = (coord.end - 1, coord.end)
    elif coordinate_system.base == "0based":
        start_interval = (coord.start, coord.start + 1)
        end_interval = (coord.end - 1, coord.end)
    else:
        raise ValueError(f"Unsupported coordinate base: {coordinate_system.base}")

    if start_interval[0] < 0 or start_interval[0] >= start_interval[1]:
        raise ValueError(f"invalid start boundary interval: {coord.chrom}:{coord.start}-{coord.end}")
    if end_interval[0] < 0 or end_interval[0] >= end_interval[1]:
        raise ValueError(f"invalid end boundary interval: {coord.chrom}:{coord.start}-{coord.end}")
    return start_interval, end_interval


def write_boundary_bed(row_coords: list[RowCoordinate], coordinate_system: CoordinateSystem, path: Path) -> None:
    with path.open("w", newline="") as handle:
        for item in row_coords:
            start_interval, end_interval = boundary_intervals(item.coord, coordinate_system)
            for boundary, interval in (("start", start_interval), ("end", end_interval)):
                handle.write(
                    "\t".join(
                        [
                            item.coord.chrom,
                            str(interval[0]),
                            str(interval[1]),
                            f"{item.row_index}|{boundary}",
                            "0",
                            item.coord.strand,
                        ]
                    )
                    + "\n"
                )


def run_liftover(liftover_bin: str, input_bed: Path, chain: Path, mapped_bed: Path, unmapped_bed: Path) -> None:
    if shutil.which(liftover_bin) is None:
        raise FileNotFoundError(f"liftOver executable not found: {liftover_bin}")
    if not chain.exists():
        raise FileNotFoundError(f"liftOver chain file not found: {chain}")
    subprocess.run(
        [liftover_bin, str(input_bed), str(chain), str(mapped_bed), str(unmapped_bed)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def parse_mapped_bed(path: Path) -> dict[tuple[int, str], tuple[str, int, int, str]]:
    mapped: dict[tuple[int, str], tuple[str, int, int, str]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start, end, name, _score, strand = fields[:6]
            row_index_text, boundary = name.split("|", 1)
            mapped[(int(row_index_text), boundary)] = (normalize_chrom(chrom), int(start), int(end), normalize_strand(strand))
    return mapped


def parse_unmapped_bed(path: Path) -> set[tuple[int, str]]:
    unmapped: set[tuple[int, str]] = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or "|" not in fields[3]:
                continue
            row_index_text, boundary = fields[3].split("|", 1)
            unmapped.add((int(row_index_text), boundary))
    return unmapped


def convert_direct(
    header: list[str],
    rows: list[dict[str, str]],
    args: argparse.Namespace,
    source_system: CoordinateSystem,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
    converted: list[dict[str, str]] = []
    unmapped: list[dict[str, str]] = []
    counts: Counter[str] = Counter(total_rows=len(rows))

    for row in rows:
        try:
            coord = parse_coordinate(row, args.chrom_col, args.start_col, args.end_col, args.strand_col, source_system)
        except ValueError as exc:
            counts["invalid_rows"] += 1
            unmapped.append(failure_row(row, args.source_coordinate_system, str(exc)))
            continue

        status = f"direct_{args.target_assembly}_1based" if source_system.base == "1based" else f"converted_from_{args.target_assembly}_0based"
        out = provenance(row, coord, args.source_coordinate_system, status)
        out[args.chrom_col] = coord.chrom
        out[args.start_col] = str(coord.start + 1 if source_system.base == "0based" else coord.start)
        out[args.end_col] = str(coord.end)
        converted.append(out)
        counts["converted_rows"] += 1

    counts["unmapped_rows"] = len(unmapped)
    return converted, unmapped, counts


def rebuild_lifted_rows(
    row_coords: list[RowCoordinate],
    mapped: dict[tuple[int, str], tuple[str, int, int, str]],
    unmapped_boundaries: set[tuple[int, str]],
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
    converted: list[dict[str, str]] = []
    failed: list[dict[str, str]] = []
    counts: Counter[str] = Counter()

    for item in row_coords:
        start_key = (item.row_index, "start")
        end_key = (item.row_index, "end")
        if start_key in unmapped_boundaries or end_key in unmapped_boundaries:
            counts["one_or_both_boundaries_unmapped"] += 1
            failed.append(failure_row(item.row, args.source_coordinate_system, "one_or_both_boundaries_unmapped"))
            continue
        if start_key not in mapped or end_key not in mapped:
            counts["missing_boundary_after_liftover"] += 1
            failed.append(failure_row(item.row, args.source_coordinate_system, "missing_boundary_after_liftover"))
            continue

        start_chrom, start_start, start_end, start_strand = mapped[start_key]
        end_chrom, end_start, end_end, end_strand = mapped[end_key]
        if start_chrom != end_chrom:
            counts["boundary_chromosomes_disagree"] += 1
            failed.append(failure_row(item.row, args.source_coordinate_system, "boundary_chromosomes_disagree"))
            continue
        if start_strand != end_strand:
            counts["boundary_strands_disagree"] += 1
            failed.append(failure_row(item.row, args.source_coordinate_system, "boundary_strands_disagree"))
            continue

        new_start = min(start_start, end_start) + 1
        new_end = max(start_end, end_end)
        if new_start > new_end:
            counts["invalid_interval_after_liftover"] += 1
            failed.append(failure_row(item.row, args.source_coordinate_system, "invalid_interval_after_liftover"))
            continue

        out = provenance(item.row, item.coord, args.source_coordinate_system, f"lifted_from_{args.source_coordinate_system}_to_{args.target_assembly}_1based")
        out[args.chrom_col] = start_chrom
        out[args.start_col] = str(new_start)
        out[args.end_col] = str(new_end)
        converted.append(out)
        counts["converted_rows"] += 1

    return converted, failed, counts


def convert_liftover(
    header: list[str],
    rows: list[dict[str, str]],
    args: argparse.Namespace,
    source_system: CoordinateSystem,
) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
    if args.chain is None:
        raise ValueError("--chain is required when source assembly differs from target assembly")

    counts: Counter[str] = Counter(total_rows=len(rows))
    valid_row_coords: list[RowCoordinate] = []
    unmapped: list[dict[str, str]] = []

    for row_index, row in enumerate(rows):
        try:
            coord = parse_coordinate(row, args.chrom_col, args.start_col, args.end_col, args.strand_col, source_system)
            boundary_intervals(coord, source_system)
        except ValueError as exc:
            counts["invalid_rows"] += 1
            unmapped.append(failure_row(row, args.source_coordinate_system, str(exc)))
            continue
        valid_row_coords.append(RowCoordinate(row_index, row, coord))

    def convert_in_workdir(workdir: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], Counter[str]]:
        workdir.mkdir(parents=True, exist_ok=True)
        input_bed = workdir / f"circrna_bsj_boundaries_{source_system.assembly}.bed"
        mapped_bed = workdir / f"circrna_bsj_boundaries_{args.target_assembly}.bed"
        unmapped_bed = workdir / "circrna_bsj_boundaries_unmapped.bed"
        write_boundary_bed(valid_row_coords, source_system, input_bed)
        run_liftover(args.liftover_bin, input_bed, args.chain, mapped_bed, unmapped_bed)
        return rebuild_lifted_rows(
            valid_row_coords,
            parse_mapped_bed(mapped_bed),
            parse_unmapped_bed(unmapped_bed),
            args,
        )

    if valid_row_coords:
        if args.workdir:
            converted, failed, liftover_counts = convert_in_workdir(args.workdir)
        else:
            with tempfile.TemporaryDirectory(prefix=f"circrna_{source_system.assembly}_to_{args.target_assembly}_") as tmpdir:
                converted, failed, liftover_counts = convert_in_workdir(Path(tmpdir))
        unmapped.extend(failed)
        counts.update(liftover_counts)
    else:
        converted = []

    counts["unmapped_rows"] = len(unmapped)
    counts["failed_liftover_rows"] = counts["unmapped_rows"] - counts["invalid_rows"]
    return converted, unmapped, counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a unified circRNA coordinate table to target assembly 1-based closed coordinates.",
    )
    parser.add_argument("--input", required=True, type=Path, help="Headered input TSV.")
    parser.add_argument(
        "--source-coordinate-system",
        required=True,
        help="Input coordinate system, formatted like hg19_0based, hg38_1based, mm9_0based, or mm10_1based.",
    )
    parser.add_argument("--target-assembly", required=True, help="Target assembly label for 1-based output, e.g. hg38 or mm10.")
    parser.add_argument("--output", required=True, type=Path, help="Converted target 1-based output TSV.")
    parser.add_argument("--unmapped-output", required=True, type=Path, help="Invalid/unmapped rows TSV.")
    parser.add_argument("--summary-output", required=True, type=Path, help="One-row conversion summary TSV.")
    parser.add_argument("--chrom-col", default="chrom", help="Input chromosome column name. Default: chrom.")
    parser.add_argument("--start-col", default="start", help="Input start coordinate column name. Default: start.")
    parser.add_argument("--end-col", default="end", help="Input end coordinate column name. Default: end.")
    parser.add_argument("--strand-col", default="strand", help="Input strand column name. Use empty string to ignore strand. Default: strand.")
    parser.add_argument("--delimiter", default="\t", help="Input/output delimiter. Default: tab.")
    parser.add_argument("--chain", type=Path, help="liftOver chain file. Required if source assembly differs from target assembly.")
    parser.add_argument("--liftover-bin", default="liftOver", help="liftOver executable. Default: liftOver.")
    parser.add_argument("--workdir", type=Path, help="Directory for temporary liftOver BED files.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.strand_col = args.strand_col or None
    source_system = CoordinateSystem.parse(args.source_coordinate_system)
    header, rows = read_table(args.input, args.delimiter)
    validate_required_columns(header, [args.chrom_col, args.start_col, args.end_col])
    if args.strand_col:
        validate_required_columns(header, [args.strand_col])

    if source_system.assembly == args.target_assembly:
        converted, unmapped, counts = convert_direct(header, rows, args, source_system)
    else:
        converted, unmapped, counts = convert_liftover(header, rows, args, source_system)

    counts["converted_rows"] = len(converted)
    counts["unmapped_rows"] = len(unmapped)
    counts["failed_liftover_rows"] = counts["unmapped_rows"] - counts["invalid_rows"]
    write_table(args.output, header, converted, args.delimiter, PROVENANCE_COLUMNS)
    write_table(args.unmapped_output, header, unmapped, args.delimiter, UNMAPPED_COLUMNS)
    write_summary(args.summary_output, args, source_system, counts)


if __name__ == "__main__":
    main()
