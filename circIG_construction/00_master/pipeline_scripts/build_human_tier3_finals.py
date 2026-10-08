#!/usr/bin/env python3
"""Rebuild human tier 3 final tables using locked coordinate systems."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from add_strand_from_presence_master import add_strands_with_support, load_presence_support


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
DEFAULT_MANIFEST = ROOT / "00_master" / "pipeline_configs" / "human_tier3_coordinate_systems.tsv"
DEFAULT_PRESENCE = ROOT / "00_master" / "base_construction" / "circRNA_presence_master.v2.tsv"
DEFAULT_IDS = ROOT / "00_master" / "base_construction" / "circRNA_presence_master_circ_ids.v2.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "hg19ToHg38.over.chain.gz"
LIFTOVER_SCRIPT = SCRIPT_DIR / "liftover_to_hg38_1based.py"
VALID_COORDINATE_SYSTEMS = {"hg19_0based", "hg19_1based", "hg38_0based", "hg38_1based"}
MANIFEST_COLUMNS = ["database", "subset", "clean_input", "coordinate_system", "deduplicate_final"]
FINAL_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
AUDIT_COLUMNS = [
    "subset",
    "circRNA_id",
    "chrom",
    "start",
    "end",
    "gene_symbol",
    "original_strand",
    "strand_presence_master",
    "final_strand",
    "strand_decision_method",
]
VALID_STRANDS = {"+", "-"}


@dataclass(frozen=True)
class Entry:
    database: str
    subset: str
    clean_input: str
    coordinate_system: str
    deduplicate_final: bool

    @property
    def label(self) -> str:
        return self.subset if self.subset != "NA" else self.database


@dataclass(frozen=True)
class Paths:
    converted: Path
    unmapped: Path
    conversion_summary: Path
    presence: Path
    presence_audit: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--input-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--presence", type=Path, default=DEFAULT_PRESENCE)
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--liftover-bin", default="liftOver")
    parser.add_argument("--jobs", type=int, choices=(1, 2), default=1)
    parser.add_argument("--database", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> list[Entry]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != MANIFEST_COLUMNS:
            raise ValueError(f"Manifest columns do not match: {reader.fieldnames}")
        entries = []
        for line_no, row in enumerate(reader, start=2):
            if row["coordinate_system"] not in VALID_COORDINATE_SYSTEMS:
                raise ValueError(f"{path}:{line_no} has an invalid coordinate system: {row['coordinate_system']}")
            if row["deduplicate_final"] not in {"0", "1"}:
                raise ValueError(f"{path}:{line_no} deduplicate_final must be 0 or 1")
            entries.append(
                Entry(
                    database=row["database"],
                    subset=row["subset"],
                    clean_input=row["clean_input"],
                    coordinate_system=row["coordinate_system"],
                    deduplicate_final=row["deduplicate_final"] == "1",
                )
            )
    if not entries:
        raise ValueError(f"Manifest contains no data rows: {path}")
    labels = [(entry.database, entry.subset) for entry in entries]
    if len(labels) != len(set(labels)):
        raise ValueError("Manifest contains duplicate database/subset pairs")
    for database, group in group_entries(entries).items():
        if len({entry.deduplicate_final for entry in group}) != 1:
            raise ValueError(f"Inconsistent deduplicate_final settings for {database}")
    return entries


def group_entries(entries: list[Entry]) -> dict[str, list[Entry]]:
    grouped: dict[str, list[Entry]] = defaultdict(list)
    for entry in entries:
        grouped[entry.database].append(entry)
    return dict(grouped)


def select_entries(entries: list[Entry], requested: list[str]) -> list[Entry]:
    if not requested:
        return entries
    requested_set = set(requested)
    available = {entry.database for entry in entries}
    unknown = sorted(requested_set - available)
    if unknown:
        raise ValueError(f"Requested databases are not in the manifest: {', '.join(unknown)}")
    return [entry for entry in entries if entry.database in requested_set]


def entry_paths(output_root: Path, entry: Entry) -> Paths:
    processed = output_root / entry.database / "processed"
    if entry.subset != "NA":
        processed = processed / "web_annotation_tables"
        stem = f"{entry.subset}.circRNAs"
    else:
        stem = f"{entry.database}.circRNAs"
    return Paths(
        converted=processed / f"{stem}.hg38_1based.tsv",
        unmapped=processed / f"{stem}.hg38_1based.unmapped.tsv",
        conversion_summary=processed / f"{stem}.hg38_1based.summary.tsv",
        presence=processed / f"{stem}.with_strand_presence_master.tsv",
        presence_audit=processed / f"{stem}.with_strand_presence_master.audit.tsv",
    )


def conversion_command(entry: Entry, args: argparse.Namespace) -> tuple[str, ...]:
    paths = entry_paths(args.output_root, entry)
    return (
        sys.executable,
        str(LIFTOVER_SCRIPT),
        "--input",
        str(args.input_root / entry.clean_input),
        "--coordinate-system",
        entry.coordinate_system,
        "--output",
        str(paths.converted),
        "--unmapped-output",
        str(paths.unmapped),
        "--summary-output",
        str(paths.conversion_summary),
        "--chain",
        str(args.chain),
        "--liftover-bin",
        args.liftover_bin,
    )


def validate_environment(entries: list[Entry], args: argparse.Namespace) -> None:
    missing = [str(args.input_root / entry.clean_input) for entry in entries if not (args.input_root / entry.clean_input).is_file()]
    if missing:
        raise FileNotFoundError("Clean input not found:\n" + "\n".join(missing))
    for path in (LIFTOVER_SCRIPT, args.presence, args.ids):
        if not path.is_file():
            raise FileNotFoundError(f"Required file not found: {path}")
    if any(entry.coordinate_system.startswith("hg19") for entry in entries) and not args.dry_run:
        if not args.chain.is_file():
            raise FileNotFoundError(f"Chain file not found: {args.chain}")
        if shutil.which(args.liftover_bin) is None and not Path(args.liftover_bin).is_file():
            raise FileNotFoundError(f"liftOver executable not found: {args.liftover_bin}")


def run_conversion(entry: Entry, args: argparse.Namespace) -> tuple[str, str]:
    command = conversion_command(entry, args)
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        return entry.label, (exc.stderr or exc.stdout or str(exc)).strip()
    return entry.label, ""


def clean_strand(value: str | None) -> str:
    value = (value or "").strip()
    return value if value in VALID_STRANDS else "NA"


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_final(database: str, entries: list[Entry], output_root: Path) -> Path:
    deduplicate = entries[0].deduplicate_final
    final_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()
    counts: Counter[str] = Counter()

    for entry in entries:
        source = entry_paths(output_root, entry).presence
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"{source} has no header")
            missing = [column for column in FINAL_COLUMNS if column not in reader.fieldnames]
            if missing:
                raise ValueError(f"{source} is missing columns: {', '.join(missing)}")
            for row in reader:
                counts["input_rows"] += 1
                original = clean_strand(row.get("strand"))
                presence = clean_strand(row.get("strand_presence_master"))
                if original in VALID_STRANDS:
                    final_strand, method = original, "raw"
                elif presence in VALID_STRANDS:
                    final_strand, method = presence, "presence_master"
                else:
                    final_strand, method = "NA", "unresolved"
                final_row = {column: row.get(column, "NA") or "NA" for column in FINAL_COLUMNS}
                final_row["strand"] = final_strand
                key = tuple(final_row[column] for column in FINAL_COLUMNS)
                if deduplicate and key in seen:
                    counts["duplicate_full_rows_removed"] += 1
                    continue
                seen.add(key)
                final_rows.append(final_row)
                audit_rows.append(
                    {
                        "subset": entry.subset,
                        "circRNA_id": final_row["circRNA_id"],
                        "chrom": final_row["chrom"],
                        "start": final_row["start"],
                        "end": final_row["end"],
                        "gene_symbol": final_row["gene_symbol"],
                        "original_strand": original,
                        "strand_presence_master": presence,
                        "final_strand": final_strand,
                        "strand_decision_method": method,
                    }
                )
                counts[f"method_{method}"] += 1

    processed = output_root / database / "processed"
    final_path = processed / f"{database}.circRNAs.final.tsv"
    write_tsv(final_path, FINAL_COLUMNS, final_rows)
    write_tsv(processed / f"{database}.circRNAs.final.audit.tsv", AUDIT_COLUMNS, audit_rows)
    summary = {
        "database": database,
        "subsets": str(len(entries)),
        "input_rows": str(counts["input_rows"]),
        "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
        "final_rows": str(len(final_rows)),
        "raw_rows": str(counts["method_raw"]),
        "presence_master_rows": str(counts["method_presence_master"]),
        "unresolved_rows": str(counts["method_unresolved"]),
        "final_output": str(final_path),
    }
    write_tsv(processed / f"{database}.circRNAs.final.summary.tsv", list(summary), [summary])
    return final_path


def main() -> int:
    args = parse_args()
    for attr in ("manifest", "input_root", "output_root", "presence", "ids", "chain"):
        setattr(args, attr, getattr(args, attr).expanduser().resolve())
    entries = select_entries(read_manifest(args.manifest), args.database)
    validate_environment(entries, args)

    if args.dry_run:
        for entry in entries:
            print(f"{entry.database}\t{entry.subset}\t{shlex.join(conversion_command(entry, args))}")
        return 0

    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_conversion, entry, args): entry for entry in entries}
        for future in as_completed(futures):
            label, error = future.result()
            print(f"convert\t{label}\t{'failed' if error else 'completed'}", flush=True)
            if error:
                errors.append(f"{label}: {error}")
    if errors:
        raise RuntimeError("Coordinate conversion failed:\n" + "\n".join(errors))

    support = load_presence_support(args.presence, args.ids, set())
    for entry in entries:
        paths = entry_paths(args.output_root, entry)
        counts = add_strands_with_support(paths.converted, paths.presence, paths.presence_audit, support)
        print(f"strand\t{entry.label}\trows={counts['input_rows']}", flush=True)

    for database, database_entries in group_entries(entries).items():
        build_final(database, database_entries, args.output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
