#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


COORD_COLUMNS = ("chrom", "start", "end", "strand")
INPUT_COLUMNS = (*COORD_COLUMNS, "circRNA_id")
MASTER_COLUMNS = (*INPUT_COLUMNS, "source")
PRIMARY_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
CHROM_ORDER = {chrom: idx for idx, chrom in enumerate(PRIMARY_CHROMS)}
PRIMARY_CHROM_SET = set(PRIMARY_CHROMS)

Coord = Tuple[str, int, int, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build master_circRNAs.txt plus 1/0 and circRNA_id presence matrices "
            "from standardized circRNA database files."
        )
    )
    parser.add_argument(
        "--config",
        default="database_sources.tsv",
        help="TSV config with columns: source, path. Default: database_sources.tsv",
    )
    parser.add_argument(
        "--master-out",
        default="master_circRNAs.txt",
        help="Output path for merged long-format table. Default: master_circRNAs.txt",
    )
    parser.add_argument(
        "--presence-out",
        default="circRNA_presence_master.tsv",
        help="Output path for 1/0 presence matrix. Default: circRNA_presence_master.tsv",
    )
    parser.add_argument(
        "--ids-out",
        default="circRNA_presence_master_circ_ids.tsv",
        help=(
            "Output path for circRNA_id presence matrix. "
            "Default: circRNA_presence_master_circ_ids.tsv"
        ),
    )
    parser.add_argument(
        "--compact-out",
        default="circRNA_presence_master_compact.tsv",
        help=(
            "Output path for compact database-name table. "
            "Default: circRNA_presence_master_compact.tsv"
        ),
    )
    parser.add_argument(
        "--compact-ids-out",
        default="circRNA_presence_master_compact_circ_ids.tsv",
        help=(
            "Output path for compact database:circRNA_id table. "
            "Default: circRNA_presence_master_compact_circ_ids.tsv"
        ),
    )
    return parser.parse_args()


def read_config(config_path: Path) -> List[Tuple[str, Path]]:
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    sources: List[Tuple[str, Path]] = []
    seen_sources: Set[str] = set()
    with config_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["source", "path"]:
            raise ValueError(
                f"{config_path} must have exactly these columns: source, path"
            )

        for line_no, row in enumerate(reader, start=2):
            source = (row.get("source") or "").strip()
            raw_path = (row.get("path") or "").strip()
            if not source or not raw_path:
                raise ValueError(f"{config_path}:{line_no} has empty source or path")
            if source in seen_sources:
                raise ValueError(f"{config_path}:{line_no} has duplicate source: {source}")
            seen_sources.add(source)

            input_path = Path(raw_path)
            if not input_path.is_absolute():
                input_path = config_path.parent / input_path
            sources.append((source, input_path))

    if not sources:
        raise ValueError(f"No database sources found in {config_path}")
    return sources


def validate_header(path: Path, fieldnames: Iterable[str] | None) -> None:
    if fieldnames is None:
        raise ValueError(f"{path} is empty")
    missing = [column for column in INPUT_COLUMNS if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def parse_coord(path: Path, line_no: int, row: Dict[str, str]) -> Coord | None:
    values = {column: (row.get(column) or "").strip() for column in INPUT_COLUMNS}
    if not all(values[column] for column in COORD_COLUMNS):
        return None

    try:
        start = int(values["start"])
        end = int(values["end"])
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_no} has non-integer start/end: "
            f"{values['start']!r}, {values['end']!r}"
        ) from exc

    return values["chrom"], start, end, values["strand"]


def unique_join(ids: Set[str]) -> str:
    return ",".join(sorted(ids)) if ids else "0"


def escape_compact_id(circ_id: str) -> str:
    return circ_id.replace(",", "%2C")


def prefixed_id_join(by_source: Dict[str, Set[str]], source_names: List[str]) -> str:
    return ",".join(
        f"{source}:{escape_compact_id(circ_id)}"
        for source in source_names
        for circ_id in sorted(by_source.get(source, set()))
    )


def sort_key(coord: Coord) -> Tuple[int, int, int, str]:
    chrom, start, end, strand = coord
    return CHROM_ORDER[chrom], start, end, strand


def build_outputs(
    sources: List[Tuple[str, Path]],
    master_out: Path,
    presence_out: Path,
    ids_out: Path,
    compact_out: Path,
    compact_ids_out: Path,
) -> None:
    source_names = [source for source, _ in sources]
    coord_to_sources: Dict[Coord, Dict[str, Set[str]]] = {}

    with master_out.open("w", newline="") as master_handle:
        master_writer = csv.writer(master_handle, delimiter="\t", lineterminator="\n")
        master_writer.writerow(MASTER_COLUMNS)

        for source, input_path in sources:
            if not input_path.exists():
                raise FileNotFoundError(f"Input file not found for {source}: {input_path}")

            with input_path.open(newline="") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                validate_header(input_path, reader.fieldnames)

                for line_no, row in enumerate(reader, start=2):
                    coord = parse_coord(input_path, line_no, row)
                    circ_id = (row.get("circRNA_id") or "").strip()
                    if coord is None or not circ_id:
                        continue

                    chrom, start, end, strand = coord
                    master_writer.writerow([chrom, start, end, strand, circ_id, source])

                    if chrom not in PRIMARY_CHROM_SET:
                        continue
                    coord_to_sources.setdefault(coord, {}).setdefault(source, set()).add(circ_id)

    ordered_coords = sorted(coord_to_sources, key=sort_key)
    presence_header = [*COORD_COLUMNS, "#dbs", *source_names]
    ids_header = [*COORD_COLUMNS, "#dbs", *source_names]

    compact_header = [*COORD_COLUMNS, "#dbs", "databases"]
    compact_ids_header = [*COORD_COLUMNS, "#dbs", "circRNA_ids"]

    with presence_out.open("w", newline="") as presence_handle, ids_out.open(
        "w", newline=""
    ) as ids_handle, compact_out.open("w", newline="") as compact_handle, compact_ids_out.open(
        "w", newline=""
    ) as compact_ids_handle:
        presence_writer = csv.writer(
            presence_handle, delimiter="\t", lineterminator="\n"
        )
        ids_writer = csv.writer(ids_handle, delimiter="\t", lineterminator="\n")
        compact_writer = csv.writer(
            compact_handle, delimiter="\t", lineterminator="\n"
        )
        compact_ids_writer = csv.writer(
            compact_ids_handle, delimiter="\t", lineterminator="\n"
        )
        presence_writer.writerow(presence_header)
        ids_writer.writerow(ids_header)
        compact_writer.writerow(compact_header)
        compact_ids_writer.writerow(compact_ids_header)

        for coord in ordered_coords:
            chrom, start, end, strand = coord
            by_source = coord_to_sources[coord]
            base_row = [chrom, start, end, strand]
            present_sources = [source for source in source_names if source in by_source]
            presence_writer.writerow(
                base_row
                + [len(by_source)]
                + ["1" if source in by_source else "0" for source in source_names]
            )
            ids_writer.writerow(
                base_row
                + [len(by_source)]
                + [unique_join(by_source.get(source, set())) for source in source_names]
            )
            compact_writer.writerow(
                base_row + [len(by_source), ",".join(present_sources)]
            )
            compact_ids_writer.writerow(
                base_row + [len(by_source), prefixed_id_join(by_source, source_names)]
            )


def main() -> None:
    args = parse_args()
    config_path = Path(args.config)
    sources = read_config(config_path)
    build_outputs(
        sources=sources,
        master_out=Path(args.master_out),
        presence_out=Path(args.presence_out),
        ids_out=Path(args.ids_out),
        compact_out=Path(args.compact_out),
        compact_ids_out=Path(args.compact_ids_out),
    )


if __name__ == "__main__":
    main()
