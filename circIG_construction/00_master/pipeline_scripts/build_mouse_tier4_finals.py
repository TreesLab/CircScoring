#!/usr/bin/env python3
"""Rebuild mouse circRNA final tables from the tier 4 manifest."""

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
DEFAULT_MANIFEST = ROOT / "00_master" / "pipeline_configs" / "mouse_tier4_coordinate_systems.tsv"
DEFAULT_PRESENCE = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master.v3.tsv"
DEFAULT_IDS = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master_circ_ids.v3.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "mm9ToMm10.over.chain.gz"
REQUIRED_COLUMNS = [
    "database", "clean_input", "coordinate_strategy", "included_coordinate_classes", "coordinate_check",
    "coordinate_summary", "runner", "legacy_policy", "final_output", "notes",
]
VALID_STRATEGIES = {"precomputed_row_level", "recompute_row_level"}
@dataclass(frozen=True)
class ManifestEntry:
    database: str
    clean_input: str
    coordinate_strategy: str
    included_coordinate_classes: str
    coordinate_check: str
    coordinate_summary: str
    runner: str
    legacy_policy: str
    final_output: str
    notes: str


@dataclass(frozen=True)
class Task:
    entry: ManifestEntry
    command: tuple[str, ...]
    required_inputs: tuple[Path, ...]
    rebuilt_final: Path


@dataclass(frozen=True)
class TaskResult:
    database: str
    success: bool
    error: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild mixed-coordinate final tables with missing strands from the mouse tier 4 manifest."
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--input-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--presence", type=Path, default=DEFAULT_PRESENCE)
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--liftover-bin", default="liftOver")
    parser.add_argument("--database", action="append", default=[])
    parser.add_argument("--jobs", type=int, choices=(1, 2), default=1)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def read_manifest(path: Path) -> list[ManifestEntry]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REQUIRED_COLUMNS:
            raise ValueError(f"Manifest columns do not match: expected {REQUIRED_COLUMNS}, got {reader.fieldnames}")
        entries = [ManifestEntry(**row) for row in reader]
    if not entries:
        raise ValueError(f"Manifest contains no data rows: {path}")
    names = [entry.database for entry in entries]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"Duplicate databases in manifest: {', '.join(duplicates)}")
    invalid = sorted({entry.coordinate_strategy for entry in entries} - VALID_STRATEGIES)
    if invalid:
        raise ValueError(f"Invalid coordinate_strategy values: {', '.join(invalid)}")
    return entries


def select_entries(entries: list[ManifestEntry], requested: list[str]) -> list[ManifestEntry]:
    if not requested:
        return entries
    requested_set = set(requested)
    unknown = sorted(requested_set - {entry.database for entry in entries})
    if unknown:
        raise ValueError(f"Requested databases are not in the manifest: {', '.join(unknown)}")
    return [entry for entry in entries if entry.database in requested_set]


def rooted(root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def staged(root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        raise ValueError(f"Manifest output must be a relative path: {value}")
    return root / path


def make_task(entry: ManifestEntry, args: argparse.Namespace) -> Task:
    clean_input = rooted(args.input_root, entry.clean_input)
    runner = rooted(args.input_root, entry.runner)
    final_output = staged(args.output_root, entry.final_output)
    processed_dir = final_output.parent
    command = [
        sys.executable, str(runner), "--input", str(clean_input),
        "--presence", str(args.presence), "--ids", str(args.ids),
        "--chain", str(args.chain), "--liftover-bin", args.liftover_bin,
        "--processed-dir", str(processed_dir),
    ]
    required = [clean_input, runner]
    if entry.coordinate_strategy == "precomputed_row_level":
        coordinate_check = rooted(args.input_root, entry.coordinate_check)
        command.extend(["--coordinate-check", str(coordinate_check)])
        required.append(coordinate_check)
    else:
        coordinate_dir = args.output_root / entry.database / "mouse" / "coordinate_check"
        command.extend([
            "--base-dataset", str(args.presence),
            "--coordinate-dir", str(coordinate_dir),
        ])
    if entry.legacy_policy == "prefer_mm10_for_ambiguous":
        command.append("--prefer-mm10-for-ambiguous")
    elif entry.legacy_policy not in {"", "NA"}:
        raise ValueError(f"{entry.database}: invalid legacy_policy: {entry.legacy_policy}")
    if entry.included_coordinate_classes not in {"", "NA"}:
        command.extend(
            argument
            for coordinate_class in entry.included_coordinate_classes.split(",")
            for argument in ("--include-coordinate-class", coordinate_class)
        )
    return Task(entry, tuple(command), tuple(required), final_output)


def validate_environment(tasks: list[Task], args: argparse.Namespace) -> None:
    if args.output_root == args.input_root:
        raise ValueError("--output-root must differ from --input-root to avoid overwriting active outputs")
    required = [args.presence, args.ids, args.chain]
    required.extend(path for task in tasks for path in task.required_inputs)
    missing = sorted(str(path) for path in required if not path.is_file())
    if missing:
        raise FileNotFoundError("Required input or script not found:\n" + "\n".join(missing))
    if not args.dry_run and shutil.which(args.liftover_bin) is None:
        raise FileNotFoundError(f"liftOver executable not found: {args.liftover_bin}")


def run_task(task: Task) -> TaskResult:
    try:
        subprocess.run(task.command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or str(exc)).strip()
        return TaskResult(task.entry.database, False, details)
    return TaskResult(task.entry.database, True)


def main() -> int:
    args = parse_args()
    for attr in ("manifest", "input_root", "output_root", "presence", "ids", "chain"):
        setattr(args, attr, getattr(args, attr).expanduser().resolve())
    entries = select_entries(read_manifest(args.manifest), args.database)
    tasks = [make_task(entry, args) for entry in entries]
    validate_environment(tasks, args)

    if args.dry_run:
        for task in tasks:
            print(f"{task.entry.database}\t{shlex.join(task.command)}")
        return 0

    results: dict[str, TaskResult] = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        future_to_task = {executor.submit(run_task, task): task for task in tasks}
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover
                result = TaskResult(task.entry.database, False, str(exc))
            results[result.database] = result
            print(f"{result.database}\t{'completed' if result.success else 'failed'}", flush=True)

    return 1 if any(not result.success for result in results.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
