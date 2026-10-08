#!/usr/bin/env python3
"""Rebuild human circRNA final tables from the tier 2 manifest."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
DEFAULT_MANIFEST = ROOT / "00_master" / "pipeline_configs" / "human_tier2_coordinate_systems.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "hg19ToHg38.over.chain.gz"
LIFTOVER_SCRIPT = SCRIPT_DIR / "liftover_to_hg38_1based.py"
WRITE_FINAL_SCRIPT = SCRIPT_DIR / "write_final.py"
VALID_COORDINATE_SYSTEMS = {"hg19_0based", "hg19_1based", "hg38_0based", "hg38_1based"}
REQUIRED_COLUMNS = [
    "database",
    "clean_input",
    "coordinate_system",
    "decision_method",
    "evidence_path",
    "converted_output",
    "unmapped_output",
    "conversion_summary",
    "final_output",
    "notes",
]
@dataclass(frozen=True)
class ManifestEntry:
    database: str
    clean_input: str
    coordinate_system: str
    decision_method: str
    evidence_path: str
    converted_output: str
    unmapped_output: str
    conversion_summary: str
    final_output: str
    notes: str


@dataclass(frozen=True)
class Task:
    entry: ManifestEntry
    commands: tuple[tuple[str, ...], ...]
    rebuilt_final: Path


@dataclass(frozen=True)
class TaskResult:
    database: str
    success: bool
    error: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild hg38_1based and six-column final tables from the human tier 2 manifest.",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=ROOT,
        help="Root directory for relative clean_input paths. Defaults to the repository root.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Staging output root; active outputs referenced by the manifest are not overwritten.",
    )
    parser.add_argument(
        "--database",
        action="append",
        default=[],
        help="Process only the specified database. May be repeated. Defaults to all databases.",
    )
    parser.add_argument("--jobs", type=int, choices=(1, 2), default=1)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--liftover-bin", default="liftOver")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> list[ManifestEntry]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REQUIRED_COLUMNS:
            raise ValueError(
                f"Manifest columns do not match: expected {REQUIRED_COLUMNS}, got {reader.fieldnames}"
            )
        rows = [ManifestEntry(**row) for row in reader]

    if not rows:
        raise ValueError(f"Manifest contains no data rows: {path}")
    databases = [row.database for row in rows]
    duplicates = sorted({name for name in databases if databases.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate databases in manifest: {', '.join(duplicates)}")
    invalid_systems = sorted(
        {row.coordinate_system for row in rows if row.coordinate_system not in VALID_COORDINATE_SYSTEMS}
    )
    if invalid_systems:
        raise ValueError(f"Invalid coordinate_system values in manifest: {', '.join(invalid_systems)}")
    return rows


def select_entries(entries: list[ManifestEntry], requested: list[str]) -> list[ManifestEntry]:
    if not requested:
        return entries
    requested_set = set(requested)
    available = {entry.database for entry in entries}
    unknown = sorted(requested_set - available)
    if unknown:
        raise ValueError(f"Requested databases are not in the manifest: {', '.join(unknown)}")
    return [entry for entry in entries if entry.database in requested_set]


def input_path(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def output_path(root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"Manifest output must be a relative path: {value}")
    return root / path


def make_task(entry: ManifestEntry, args: argparse.Namespace) -> Task:
    clean_input = input_path(args.input_root, entry.clean_input)
    rebuilt_final = output_path(args.output_root, entry.final_output)
    if entry.decision_method == "reference_lookup":
        commands = (
            (
                sys.executable,
                str(WRITE_FINAL_SCRIPT),
                "--input",
                str(clean_input),
                "--output",
                str(rebuilt_final),
            ),
        )
    else:
        for field_name in ("converted_output", "unmapped_output", "conversion_summary"):
            if getattr(entry, field_name) == "NA":
                raise ValueError(f"{entry.database} is missing {field_name}")
        converted = output_path(args.output_root, entry.converted_output)
        unmapped = output_path(args.output_root, entry.unmapped_output)
        summary = output_path(args.output_root, entry.conversion_summary)
        commands = (
            (
                sys.executable,
                str(LIFTOVER_SCRIPT),
                "--input",
                str(clean_input),
                "--coordinate-system",
                entry.coordinate_system,
                "--output",
                str(converted),
                "--unmapped-output",
                str(unmapped),
                "--summary-output",
                str(summary),
                "--chain",
                str(args.chain),
                "--liftover-bin",
                args.liftover_bin,
            ),
            (
                sys.executable,
                str(WRITE_FINAL_SCRIPT),
                "--input",
                str(converted),
                "--output",
                str(rebuilt_final),
            ),
        )
    return Task(entry, commands, rebuilt_final)


def validate_environment(tasks: list[Task], args: argparse.Namespace) -> None:
    for script in (LIFTOVER_SCRIPT, WRITE_FINAL_SCRIPT):
        if not script.is_file():
            raise FileNotFoundError(f"Shared script not found: {script}")
    missing_inputs = sorted(
        {
            command[command.index("--input") + 1]
            for task in tasks
            for command in task.commands[:1]
            if not Path(command[command.index("--input") + 1]).is_file()
        }
    )
    if missing_inputs:
        raise FileNotFoundError("Clean input not found:\n" + "\n".join(missing_inputs))

    requires_liftover = any(
        task.entry.coordinate_system.startswith("hg19")
        and task.entry.decision_method != "reference_lookup"
        for task in tasks
    )
    if requires_liftover and not args.dry_run:
        if not args.chain.is_file():
            raise FileNotFoundError(f"liftOver chain not found: {args.chain}")
        if shutil.which(args.liftover_bin) is None:
            raise FileNotFoundError(f"liftOver executable not found: {args.liftover_bin}")


def run_task(task: Task) -> TaskResult:
    for command in task.commands:
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            details = (exc.stderr or exc.stdout or str(exc)).strip()
            return TaskResult(task.entry.database, False, details)
    return TaskResult(task.entry.database, True)


def main() -> int:
    args = parse_args()
    args.manifest = args.manifest.expanduser().resolve()
    args.input_root = args.input_root.expanduser().resolve()
    args.output_root = args.output_root.expanduser().resolve()
    args.chain = args.chain.expanduser().resolve()

    entries = select_entries(read_manifest(args.manifest), args.database)
    tasks = [make_task(entry, args) for entry in entries]
    validate_environment(tasks, args)

    if args.dry_run:
        for task in tasks:
            for command in task.commands:
                print(f"{task.entry.database}\t{shlex.join(command)}")
        return 0

    results: dict[str, TaskResult] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_task = {executor.submit(run_task, task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover - defensive worker boundary
                result = TaskResult(task.entry.database, False, str(exc))
            results[result.database] = result
            print(f"{result.database}\t{'completed' if result.success else 'failed'}", flush=True)

    return 1 if any(not result.success for result in results.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
