#!/usr/bin/env python3
"""Convert a mouse mm10 1-based circRNA table to hg38 1-based coordinates."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


csv.field_size_limit(sys.maxsize)

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "mm10ToHg38.over.chain.gz"
COORD_COLUMNS = ("chrom", "start", "end", "strand")
PROVENANCE_COLUMNS = (
    "original_chrom", "original_start", "original_end", "original_strand",
    "original_coordinate_system", "conversion_status", "conversion_note",
)
UNMAPPED_COLUMNS = ("failure_reason", "original_strand", "original_coordinate_system")
SUMMARY_COLUMNS = (
    "input_file", "output_file", "unmapped_output_file", "source_coordinate_system",
    "source_assembly", "target_assembly", "total_rows", "converted_rows",
    "unmapped_rows", "reused_old_mapping_rows", "partial_mapping_rows",
    "old_unmapped_rows", "new_unmapped_rows", "invalid_rows",
    "failed_liftover_rows", "chain", "liftover_bin",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--unmapped-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--liftover-bin", default="liftOver")
    parser.add_argument("--workdir", type=Path)
    parser.add_argument(
        "--previous-mapped",
        type=Path,
        help="Optional verified mapped table used for an incremental rebuild.",
    )
    parser.add_argument(
        "--previous-unmapped",
        type=Path,
        help="Optional verified unmapped table used for an incremental rebuild.",
    )
    parser.add_argument(
        "--partial-mapping-overrides",
        type=Path,
        help=(
            "Optional frozen mappings for historical incremental rebuilds. "
            "Rows not found in previous tables or this manifest remain unmapped."
        ),
    )
    parser.add_argument(
        "--liftover-rows-without-history",
        action="store_true",
        help=(
            "Fresh-build mode: apply partial mapping overrides first, then run liftOver "
            "for rows absent from previous tables."
        ),
    )
    return parser.parse_args()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return list(reader.fieldnames), list(reader)


def validate_header(path: Path, header: list[str], required: tuple[str, ...]) -> None:
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")


def source_key(row: dict[str, str], original: bool = False) -> tuple[str, int, int, str]:
    prefix = "original_" if original else ""
    chrom = (row.get(f"{prefix}chrom") or "").strip()
    start = int((row.get(f"{prefix}start") or "").strip())
    end = int((row.get(f"{prefix}end") or "").strip())
    strand = (row.get(f"{prefix}strand") or "NA").strip() or "NA"
    if not chrom or start < 1 or start > end:
        raise ValueError(f"invalid mm10_1based coordinate: {chrom}:{start}-{end}({strand})")
    return chrom, start, end, strand


def coordinate_key(key: tuple[str, int, int, str]) -> tuple[str, int, int]:
    return key[:3]


def mapped_target(row: dict[str, str]) -> tuple[str, int, int]:
    return row["chrom"], int(row["start"]), int(row["end"])


def load_previous_mapped(
    path: Path | None,
) -> tuple[dict[tuple[str, int, int, str], tuple[str, int, int]], set[tuple[str, int, int]]]:
    exact: dict[tuple[str, int, int, str], tuple[str, int, int]] = {}
    coordinates: set[tuple[str, int, int]] = set()
    if path is None:
        return exact, coordinates
    header, rows = read_tsv(path)
    validate_header(path, header, (*COORD_COLUMNS, "original_chrom", "original_start", "original_end", "original_strand"))
    for row in rows:
        key = source_key(row, original=True)
        exact[key] = mapped_target(row)
        coordinates.add(coordinate_key(key))
    return exact, coordinates


def load_previous_unmapped(
    path: Path | None,
) -> tuple[dict[tuple[str, int, int, str], str], set[tuple[str, int, int]]]:
    exact: dict[tuple[str, int, int, str], str] = {}
    coordinates: set[tuple[str, int, int]] = set()
    if path is None:
        return exact, coordinates
    header, rows = read_tsv(path)
    validate_header(path, header, (*COORD_COLUMNS, "failure_reason"))
    for row in rows:
        key = source_key(row)
        exact[key] = row.get("failure_reason") or "previously_unmapped"
        coordinates.add(coordinate_key(key))
    return exact, coordinates


def load_partial_mapping_overrides(
    path: Path | None,
) -> dict[tuple[str, int, int, str], tuple[str, int, int]]:
    overrides: dict[tuple[str, int, int, str], tuple[str, int, int]] = {}
    if path is None:
        return overrides
    header, rows = read_tsv(path)
    validate_header(
        path,
        header,
        ("original_chrom", "original_start", "original_end", "original_strand", "chrom", "start", "end"),
    )
    for row in rows:
        key = source_key(row, original=True)
        if key in overrides:
            raise ValueError(f"{path}: duplicate original coordinate: {key}")
        overrides[key] = mapped_target(row)
    return overrides


def write_boundary_bed(rows: list[tuple[int, tuple[str, int, int, str]]], path: Path) -> None:
    with path.open("w") as handle:
        for row_index, (chrom, start, end, strand) in rows:
            bed_strand = strand if strand in {"+", "-"} else "NA"
            handle.write(f"{chrom}\t{start - 1}\t{start}\t{row_index}|start\t0\t{bed_strand}\n")
            handle.write(f"{chrom}\t{end - 1}\t{end}\t{row_index}|end\t0\t{bed_strand}\n")


def run_liftover(args: argparse.Namespace, pending: list[tuple[int, tuple[str, int, int, str]]], workdir: Path) -> tuple[dict[tuple[int, str], tuple[str, int, int, str]], set[tuple[int, str]]]:
    if shutil.which(args.liftover_bin) is None:
        raise FileNotFoundError(f"liftOver executable not found: {args.liftover_bin}")
    if not args.chain.is_file():
        raise FileNotFoundError(f"liftOver chain not found: {args.chain}")
    workdir.mkdir(parents=True, exist_ok=True)
    source_bed = workdir / "mouse_v4_boundaries.mm10.bed"
    mapped_bed = workdir / "mouse_v4_boundaries.hg38.bed"
    unmapped_bed = workdir / "mouse_v4_boundaries.unmapped.bed"
    write_boundary_bed(pending, source_bed)
    subprocess.run(
        [args.liftover_bin, str(source_bed), str(args.chain), str(mapped_bed), str(unmapped_bed)],
        check=True,
        capture_output=True,
        text=True,
    )
    mapped: dict[tuple[int, str], tuple[str, int, int, str]] = {}
    with mapped_bed.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            row_text, boundary = fields[3].split("|", 1)
            mapped[(int(row_text), boundary)] = (fields[0], int(fields[1]), int(fields[2]), fields[5])
    unmapped: set[tuple[int, str]] = set()
    with unmapped_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and "|" in fields[3]:
                row_text, boundary = fields[3].split("|", 1)
                unmapped.add((int(row_text), boundary))
    return mapped, unmapped


def mapped_row(row: dict[str, str], key: tuple[str, int, int, str], target: tuple[str, int, int], note: str) -> dict[str, str]:
    chrom, start, end, strand = key
    out = dict(row)
    out.update(
        {
            "chrom": target[0], "start": str(target[1]), "end": str(target[2]),
            "original_chrom": chrom, "original_start": str(start),
            "original_end": str(end), "original_strand": strand,
            "original_coordinate_system": "mm10_1based",
            "conversion_status": "lifted_from_mm10_1based_to_hg38_1based",
            "conversion_note": note,
        }
    )
    return out


def unmapped_row(row: dict[str, str], key: tuple[str, int, int, str], reason: str) -> dict[str, str]:
    out = dict(row)
    out.update(
        {
            "failure_reason": reason,
            "original_strand": key[3],
            "original_coordinate_system": "mm10_1based",
        }
    )
    return out


def write_tsv(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, args: argparse.Namespace, counts: Counter[str]) -> None:
    row = {
        "input_file": str(args.input), "output_file": str(args.output),
        "unmapped_output_file": str(args.unmapped_output),
        "source_coordinate_system": "mm10_1based", "source_assembly": "mm10",
        "target_assembly": "hg38", "total_rows": str(counts["total_rows"]),
        "converted_rows": str(counts["converted_rows"]),
        "unmapped_rows": str(counts["unmapped_rows"]),
        "reused_old_mapping_rows": str(counts["reused_old_mapping_rows"]),
        "partial_mapping_rows": str(counts["partial_mapping_rows"]),
        "old_unmapped_rows": str(counts["old_unmapped_rows"]),
        "new_unmapped_rows": str(counts["new_unmapped_rows"]),
        "invalid_rows": str(counts["invalid_rows"]),
        "failed_liftover_rows": str(counts["unmapped_rows"] - counts["invalid_rows"]),
        "chain": str(args.chain), "liftover_bin": args.liftover_bin,
    }
    write_tsv(path, list(SUMMARY_COLUMNS), [row])


def main() -> int:
    args = parse_args()
    for attr in (
        "input", "output", "unmapped_output", "summary_output", "chain",
        "previous_mapped", "previous_unmapped", "partial_mapping_overrides", "workdir",
    ):
        value = getattr(args, attr)
        if value is not None:
            setattr(args, attr, value.expanduser().resolve())
    header, input_rows = read_tsv(args.input)
    validate_header(args.input, header, COORD_COLUMNS)
    previous_mapped, previous_mapped_coordinates = load_previous_mapped(args.previous_mapped)
    previous_unmapped, previous_unmapped_coordinates = load_previous_unmapped(args.previous_unmapped)
    partial_mapping_overrides = load_partial_mapping_overrides(args.partial_mapping_overrides)

    counts: Counter[str] = Counter(total_rows=len(input_rows))
    keys: list[tuple[str, int, int, str] | None] = []
    pending: list[tuple[int, tuple[str, int, int, str]]] = []
    disposition: list[str] = []
    for index, row in enumerate(input_rows):
        try:
            key = source_key(row)
        except (TypeError, ValueError):
            keys.append(None)
            disposition.append("invalid")
            counts["invalid_rows"] += 1
            continue
        keys.append(key)
        coord = coordinate_key(key)
        if key in previous_mapped:
            disposition.append("reuse_mapped")
        elif key in previous_unmapped:
            disposition.append("reuse_unmapped")
        elif key in partial_mapping_overrides:
            disposition.append("override")
        elif coord in previous_mapped_coordinates:
            if args.partial_mapping_overrides:
                disposition.append("new_unmapped")
            else:
                disposition.append("partial")
                pending.append((index, key))
        elif coord in previous_unmapped_coordinates:
            disposition.append("new_unmapped")
        else:
            if args.partial_mapping_overrides and not args.liftover_rows_without_history:
                disposition.append("new_unmapped")
            else:
                disposition.append("partial")
                pending.append((index, key))

    if pending:
        if args.workdir:
            mapped, boundary_unmapped = run_liftover(args, pending, args.workdir)
        else:
            with tempfile.TemporaryDirectory(prefix="mouse_mm10_to_hg38_") as tmpdir:
                mapped, boundary_unmapped = run_liftover(args, pending, Path(tmpdir))
    else:
        mapped, boundary_unmapped = {}, set()

    converted: list[dict[str, str]] = []
    unmapped: list[dict[str, str]] = []
    for index, (row, key, state) in enumerate(zip(input_rows, keys, disposition)):
        if key is None:
            fallback = ((row.get("chrom") or "NA"), 0, 0, (row.get("strand") or "NA"))
            unmapped.append(unmapped_row(row, fallback, "invalid_mm10_1based_coordinate"))
            continue
        if state == "reuse_mapped":
            converted.append(mapped_row(row, key, previous_mapped[key], "reused_previous_verified_mapping"))
            counts["reused_old_mapping_rows"] += 1
            continue
        if state == "reuse_unmapped":
            unmapped.append(unmapped_row(row, key, previous_unmapped[key]))
            counts["old_unmapped_rows"] += 1
            continue
        if state == "override":
            converted.append(
                mapped_row(
                    row,
                    key,
                    partial_mapping_overrides[key],
                    "new_partial_direct_mapping",
                )
            )
            counts["partial_mapping_rows"] += 1
            continue
        if state == "new_unmapped":
            unmapped.append(unmapped_row(row, key, "no_mapping_available_for_new_coordinate"))
            counts["new_unmapped_rows"] += 1
            continue

        start_key = (index, "start")
        end_key = (index, "end")
        if start_key in boundary_unmapped or end_key in boundary_unmapped or start_key not in mapped or end_key not in mapped:
            reason = "no_mapping_available_for_new_coordinate" if args.previous_mapped else "one_or_both_boundaries_unmapped"
            unmapped.append(unmapped_row(row, key, reason))
            counts["new_unmapped_rows"] += 1
            continue
        start_hit = mapped[start_key]
        end_hit = mapped[end_key]
        if start_hit[0] != end_hit[0]:
            unmapped.append(unmapped_row(row, key, "boundary_chromosomes_disagree"))
            counts["new_unmapped_rows"] += 1
            continue
        if args.previous_mapped:
            target = (start_hit[0], start_hit[2], end_hit[2])
            note = "new_partial_direct_mapping"
        else:
            target = (start_hit[0], min(start_hit[1], end_hit[1]) + 1, max(start_hit[2], end_hit[2]))
            note = "NA"
        converted.append(mapped_row(row, key, target, note))
        counts["partial_mapping_rows"] += 1

    counts["converted_rows"] = len(converted)
    counts["unmapped_rows"] = len(unmapped)
    mapped_header = [*header, *(column for column in PROVENANCE_COLUMNS if column not in header)]
    unmapped_header = [*header, *(column for column in UNMAPPED_COLUMNS if column not in header)]
    write_tsv(args.output, mapped_header, converted)
    write_tsv(args.unmapped_output, unmapped_header, unmapped)
    write_summary(args.summary_output, args, counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
