#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


csv.field_size_limit(sys.maxsize)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_CONSTRUCTION_DIR = SCRIPT_DIR.parent / "base_construction"
DEFAULT_PRESENCE = BASE_CONSTRUCTION_DIR / "circRNA_presence_master.tsv"
DEFAULT_IDS = BASE_CONSTRUCTION_DIR / "circRNA_presence_master_circ_ids.tsv"
REQUIRED_INPUT_COLUMNS = ("circRNA_id", "chrom", "start", "end")
REQUIRED_BATCH_COLUMNS = ("label", "input", "output", "audit")
OUTPUT_STRAND_COLUMN = "strand_presence_master"
COORD_COLUMNS = ("chrom", "start", "end", "strand")
KEY_COLUMNS = ("chrom", "start", "end")
VALID_STRANDS = {"+", "-"}
EMPTY_VALUES = {"", ".", "*", "-", "NA", "N/A", "None", "none", "null", "NULL"}
BATCH_SUMMARY_COUNT_COLUMNS = (
    "input_rows",
    "filled_strand",
    "resolved_from_presence_master",
    "unresolved_conflicting_base_strands",
    "unresolved_invalid_coordinate",
    "unresolved_no_base_hit",
    "strand_+",
    "strand_-",
    "strand_NA",
)


@dataclass
class StrandSupport:
    strands: set[str] = field(default_factory=set)
    dbs_by_strand: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    ids_by_strand_db: dict[str, dict[str, set[str]]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(set))
    )

    def add(self, strand: str, dbs: Iterable[str], ids_by_db: dict[str, set[str]] | None) -> None:
        self.strands.add(strand)
        for db in dbs:
            self.dbs_by_strand[strand].add(db)
            if ids_by_db and db in ids_by_db:
                self.ids_by_strand_db[strand][db].update(ids_by_db[db])

    def all_dbs(self) -> set[str]:
        dbs: set[str] = set()
        for values in self.dbs_by_strand.values():
            dbs.update(values)
        return dbs

    def all_database_ids(self) -> dict[str, set[str]]:
        ids_by_db: dict[str, set[str]] = defaultdict(set)
        for by_db in self.ids_by_strand_db.values():
            for db, ids in by_db.items():
                ids_by_db[db].update(ids)
        return ids_by_db


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Infer strand_presence_master for a unified circRNA TSV using "
            "circRNA_presence_master.tsv as coordinate-level strand support. "
            "Input coordinates must be hg38 1-based closed."
        )
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        type=Path,
        help="Unified input TSV with circRNA_id/chrom/start/end and optional strand.",
    )
    input_group.add_argument(
        "--batch-manifest",
        type=Path,
        help="Batch manifest TSV with label/input/output/audit columns.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=f"Output TSV with {OUTPUT_STRAND_COLUMN} inferred where possible.",
    )
    parser.add_argument(
        "--presence",
        type=Path,
        default=DEFAULT_PRESENCE,
        help=f"Presence master TSV. Default: {DEFAULT_PRESENCE}",
    )
    parser.add_argument(
        "--ids",
        type=Path,
        default=None,
        help=(
            "Optional circRNA ID presence master TSV. IDs are audit-only and "
            "are reported with database context."
        ),
    )
    parser.add_argument(
        "--audit",
        type=Path,
        help="Per-row audit TSV.",
    )
    parser.add_argument(
        "--batch-summary",
        type=Path,
        help="Batch summary TSV. Required with --batch-manifest.",
    )
    parser.add_argument(
        "--exclude-db",
        action="append",
        default=[],
        help=(
            "Database column to exclude from reference support. May be provided "
            "multiple times, for example --exclude-db CircFunBase."
        ),
    )
    return parser.parse_args()


def clean_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if value in VALID_STRANDS:
        return value
    if value in EMPTY_VALUES:
        return "NA"
    return "NA"


def clean(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    if value in EMPTY_VALUES:
        return "NA"
    return value


def validate_columns(path: Path, fieldnames: list[str] | None, required: Iterable[str]) -> None:
    if fieldnames is None:
        raise ValueError(f"{path} is empty or missing a header")
    missing = [column for column in required if column not in fieldnames]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def parse_key(path: Path, line_no: int, row: dict[str, str]) -> tuple[str, int, int] | None:
    chrom = (row.get("chrom") or "").strip()
    start_text = (row.get("start") or "").strip()
    end_text = (row.get("end") or "").strip()
    if not chrom or not start_text or not end_text:
        return None

    try:
        return chrom, int(start_text), int(end_text)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_no} has non-integer start/end: {start_text!r}, {end_text!r}"
        ) from exc


def database_columns(fieldnames: list[str], exclude_dbs: set[str]) -> list[str]:
    fixed = {*COORD_COLUMNS, "#dbs"}
    dbs = [column for column in fieldnames if column not in fixed]
    unknown = sorted(exclude_dbs.difference(dbs))
    if unknown:
        raise ValueError(f"Unknown --exclude-db value(s): {', '.join(unknown)}")
    return [db for db in dbs if db not in exclude_dbs]


def parse_ids_by_db(row: dict[str, str], dbs: Iterable[str]) -> dict[str, set[str]]:
    ids_by_db: dict[str, set[str]] = {}
    for db in dbs:
        value = (row.get(db) or "").strip()
        if not value or value == "0":
            continue
        ids = {part.strip() for part in value.split(",") if part.strip()}
        if ids:
            ids_by_db[db] = ids
    return ids_by_db


def load_ids_lookup(
    ids_path: Path | None,
    expected_header: list[str],
    exclude_dbs: set[str],
) -> dict[tuple[str, int, int, str], dict[str, set[str]]]:
    if ids_path is None:
        return {}
    if not ids_path.exists():
        raise FileNotFoundError(f"IDs file not found: {ids_path}")

    lookup: dict[tuple[str, int, int, str], dict[str, set[str]]] = {}
    with ids_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(ids_path, reader.fieldnames, COORD_COLUMNS)
        if reader.fieldnames != expected_header:
            raise ValueError(
                f"{ids_path} header does not match the presence header; "
                "the files must come from the same build."
            )
        dbs = database_columns(reader.fieldnames, exclude_dbs)
        for line_no, row in enumerate(reader, start=2):
            key = parse_key(ids_path, line_no, row)
            if key is None:
                continue
            strand = clean_strand(row.get("strand"))
            if strand not in VALID_STRANDS:
                continue
            ids_by_db = parse_ids_by_db(row, dbs)
            if ids_by_db:
                lookup[(*key, strand)] = ids_by_db
    return lookup


def load_presence_support(
    presence_path: Path,
    ids_path: Path | None,
    exclude_dbs: set[str],
) -> dict[tuple[str, int, int], StrandSupport]:
    if not presence_path.exists():
        raise FileNotFoundError(f"Presence file not found: {presence_path}")

    support_by_key: dict[tuple[str, int, int], StrandSupport] = defaultdict(StrandSupport)
    with presence_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(presence_path, reader.fieldnames, (*COORD_COLUMNS, "#dbs"))
        dbs = database_columns(reader.fieldnames, exclude_dbs)
        ids_lookup = load_ids_lookup(ids_path, reader.fieldnames, exclude_dbs)

        for line_no, row in enumerate(reader, start=2):
            key = parse_key(presence_path, line_no, row)
            if key is None:
                continue
            strand = clean_strand(row.get("strand"))
            if strand not in VALID_STRANDS:
                continue

            supporting_dbs = {db for db in dbs if (row.get(db) or "").strip() == "1"}
            if not supporting_dbs:
                continue

            ids_by_db = ids_lookup.get((*key, strand))
            support_by_key[key].add(strand, supporting_dbs, ids_by_db)

    return dict(support_by_key)


def join_values(values: Iterable[str]) -> str:
    values = sorted(set(values))
    return ",".join(values) if values else "NA"


def format_database_ids(ids_by_db: dict[str, set[str]]) -> str:
    if not ids_by_db:
        return "NA"
    parts = []
    for db in sorted(ids_by_db):
        ids = ",".join(sorted(ids_by_db[db]))
        parts.append(f"{db}:{ids}")
    return ";".join(parts)


def resolve_strand(
    support_by_key: dict[tuple[str, int, int], StrandSupport],
    key: tuple[str, int, int] | None,
) -> tuple[str, str, StrandSupport | None]:
    if key is None:
        return "NA", "unresolved_invalid_coordinate", None
    support = support_by_key.get(key)
    if support is None or not support.strands:
        return "NA", "unresolved_no_base_hit", None
    if len(support.strands) == 1:
        return next(iter(support.strands)), "resolved_from_presence_master", support
    return "NA", "unresolved_conflicting_base_strands", support


def add_strands_with_support(
    input_path: Path,
    output_path: Path,
    audit_path: Path,
    support_by_key: dict[tuple[str, int, int], StrandSupport],
) -> Counter[str]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    with input_path.open(newline="") as input_handle, output_path.open(
        "w", newline=""
    ) as output_handle, audit_path.open("w", newline="") as audit_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        validate_columns(input_path, reader.fieldnames, REQUIRED_INPUT_COLUMNS)
        assert reader.fieldnames is not None
        output_fields = list(reader.fieldnames)
        if OUTPUT_STRAND_COLUMN not in output_fields:
            output_fields.append(OUTPUT_STRAND_COLUMN)

        writer = csv.DictWriter(
            output_handle,
            fieldnames=output_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        audit_fields = [
            "circRNA_id",
            "chrom",
            "start",
            "end",
            "original_strand",
            "inferred_strand",
            "status",
            "reference_strands",
            "reference_dbs",
            "reference_database_ids",
        ]
        audit_writer = csv.DictWriter(
            audit_handle,
            fieldnames=audit_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        audit_writer.writeheader()

        for line_no, row in enumerate(reader, start=2):
            counts["input_rows"] += 1
            original_strand = clean_strand(row.get("strand"))

            try:
                key = parse_key(input_path, line_no, row)
            except ValueError:
                key = None

            inferred_strand, status, support = resolve_strand(support_by_key, key)
            if inferred_strand in VALID_STRANDS:
                counts["filled_strand"] += 1

            output_row = {field: row.get(field, "") for field in output_fields}
            output_row[OUTPUT_STRAND_COLUMN] = (
                inferred_strand if inferred_strand in VALID_STRANDS else "NA"
            )
            counts[status] += 1
            counts[f"strand_{output_row[OUTPUT_STRAND_COLUMN]}"] += 1
            writer.writerow(output_row)

            support_strands = support.strands if support else set()
            support_dbs = support.all_dbs() if support else set()
            support_ids = support.all_database_ids() if support else {}
            audit_writer.writerow(
                {
                    "circRNA_id": clean(row.get("circRNA_id")),
                    "chrom": clean(row.get("chrom")),
                    "start": clean(row.get("start")),
                    "end": clean(row.get("end")),
                    "original_strand": original_strand,
                    "inferred_strand": inferred_strand if inferred_strand in VALID_STRANDS else "NA",
                    "status": status,
                    "reference_strands": join_values(support_strands),
                    "reference_dbs": join_values(support_dbs),
                    "reference_database_ids": format_database_ids(support_ids),
                }
            )

    return counts


def add_strands(
    input_path: Path,
    output_path: Path,
    presence_path: Path,
    ids_path: Path | None,
    audit_path: Path,
    exclude_dbs: set[str],
) -> Counter[str]:
    support_by_key = load_presence_support(presence_path, ids_path, exclude_dbs)
    return add_strands_with_support(
        input_path=input_path,
        output_path=output_path,
        audit_path=audit_path,
        support_by_key=support_by_key,
    )


def resolve_manifest_path(path: Path, base_dir: Path) -> Path:
    path = path.expanduser()
    if path.is_absolute():
        return path.resolve()
    return (base_dir / path).resolve()


def read_batch_manifest(manifest_path: Path) -> list[dict[str, Path | str]]:
    manifest_path = manifest_path.expanduser().resolve()
    rows: list[dict[str, Path | str]] = []
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(manifest_path, reader.fieldnames, REQUIRED_BATCH_COLUMNS)
        base_dir = manifest_path.parent
        for line_no, row in enumerate(reader, start=2):
            label = (row.get("label") or "").strip()
            if not label:
                raise ValueError(f"{manifest_path}:{line_no} has empty label")
            for column in ("input", "output", "audit"):
                if not (row.get(column) or "").strip():
                    raise ValueError(f"{manifest_path}:{line_no} has empty {column}")
            rows.append(
                {
                    "label": label,
                    "input": resolve_manifest_path(Path(row["input"]), base_dir),
                    "output": resolve_manifest_path(Path(row["output"]), base_dir),
                    "audit": resolve_manifest_path(Path(row["audit"]), base_dir),
                }
            )
    if not rows:
        raise ValueError(f"{manifest_path} has no batch rows")
    return rows


def write_batch_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["label", "input", "output", "audit", *BATCH_SUMMARY_COUNT_COLUMNS]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def add_strands_batch(
    manifest_path: Path,
    batch_summary_path: Path,
    presence_path: Path,
    ids_path: Path | None,
    exclude_dbs: set[str],
) -> list[dict[str, str]]:
    jobs = read_batch_manifest(manifest_path)
    support_by_key = load_presence_support(presence_path, ids_path, exclude_dbs)
    summary_rows: list[dict[str, str]] = []
    for job in jobs:
        label = str(job["label"])
        input_path = job["input"]
        output_path = job["output"]
        audit_path = job["audit"]
        assert isinstance(input_path, Path)
        assert isinstance(output_path, Path)
        assert isinstance(audit_path, Path)
        counts = add_strands_with_support(
            input_path=input_path,
            output_path=output_path,
            audit_path=audit_path,
            support_by_key=support_by_key,
        )
        summary_rows.append(
            {
                "label": label,
                "input": str(input_path),
                "output": str(output_path),
                "audit": str(audit_path),
                **{key: str(counts[key]) for key in BATCH_SUMMARY_COUNT_COLUMNS},
            }
        )
    write_batch_summary(batch_summary_path, summary_rows)
    return summary_rows


def main() -> None:
    args = parse_args()
    if args.batch_manifest:
        if args.output or args.audit:
            raise ValueError("--output and --audit are only valid with single-file --input mode")
        if not args.batch_summary:
            raise ValueError("--batch-summary is required with --batch-manifest")
        summary_rows = add_strands_batch(
            manifest_path=args.batch_manifest,
            batch_summary_path=args.batch_summary,
            presence_path=args.presence,
            ids_path=args.ids,
            exclude_dbs=set(args.exclude_db),
        )
        sys.stderr.write(f"batch_manifest: {args.batch_manifest}\n")
        sys.stderr.write(f"batch_summary: {args.batch_summary}\n")
        sys.stderr.write(f"presence: {args.presence}\n")
        if args.ids:
            sys.stderr.write(f"ids: {args.ids}\n")
        if args.exclude_db:
            sys.stderr.write(f"excluded_dbs: {','.join(sorted(set(args.exclude_db)))}\n")
        sys.stderr.write(f"batch_jobs: {len(summary_rows)}\n")
        for row in summary_rows:
            sys.stderr.write(
                f"{row['label']}: input_rows={row['input_rows']} "
                f"filled_strand={row['filled_strand']} "
                f"strand_NA={row['strand_NA']}\n"
            )
        return

    if args.output is None or args.audit is None:
        raise ValueError("--output and --audit are required with --input")

    counts = add_strands(
        input_path=args.input,
        output_path=args.output,
        presence_path=args.presence,
        ids_path=args.ids,
        audit_path=args.audit,
        exclude_dbs=set(args.exclude_db),
    )

    sys.stderr.write(f"input: {args.input}\n")
    sys.stderr.write(f"output: {args.output}\n")
    sys.stderr.write(f"presence: {args.presence}\n")
    if args.ids:
        sys.stderr.write(f"ids: {args.ids}\n")
    sys.stderr.write(f"audit: {args.audit}\n")
    if args.exclude_db:
        sys.stderr.write(f"excluded_dbs: {','.join(sorted(set(args.exclude_db)))}\n")
    for key in sorted(counts):
        sys.stderr.write(f"{key}: {counts[key]}\n")


if __name__ == "__main__":
    main()
