#!/usr/bin/env python3
"""Rebuild mouse circRNA final tables from the tier 3 manifest."""

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
DEFAULT_MANIFEST = ROOT / "00_master" / "pipeline_configs" / "mouse_tier3_coordinate_systems.tsv"
DEFAULT_PRESENCE = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master.v2.tsv"
DEFAULT_IDS = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master_circ_ids.v2.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "mm9ToMm10.over.chain.gz"
LIFTOVER_SCRIPT = SCRIPT_DIR / "liftover_to_1based.py"
ADD_STRAND_SCRIPT = SCRIPT_DIR / "add_strand_from_presence_master.py"
VALID_STRATEGIES = {"uniform", "per_row"}
VALID_COORDINATE_SYSTEMS = {"mm9_0based", "mm9_1based", "mm10_0based", "mm10_1based", "mixed"}
REQUIRED_COLUMNS = [
    "database", "clean_input", "coordinate_strategy", "coordinate_system",
    "included_coordinate_classes",
    "coordinate_check", "coordinate_summary", "decision_script",
    "converted_output", "unmapped_output", "conversion_summary",
    "presence_output", "presence_audit", "final_output", "final_audit",
    "conflict_output", "final_summary", "notes",
]
@dataclass(frozen=True)
class ManifestEntry:
    database: str
    clean_input: str
    coordinate_strategy: str
    coordinate_system: str
    included_coordinate_classes: str
    coordinate_check: str
    coordinate_summary: str
    decision_script: str
    converted_output: str
    unmapped_output: str
    conversion_summary: str
    presence_output: str
    presence_audit: str
    final_output: str
    final_audit: str
    conflict_output: str
    final_summary: str
    notes: str


@dataclass(frozen=True)
class Task:
    entry: ManifestEntry
    commands: tuple[tuple[str, ...], ...]
    required_inputs: tuple[Path, ...]
    rebuilt_final: Path


@dataclass(frozen=True)
class TaskResult:
    database: str
    success: bool
    error: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild mm10_1based tables, infer strands, and write final tables from the mouse tier 3 manifest."
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
    for entry in entries:
        if entry.coordinate_strategy not in VALID_STRATEGIES:
            raise ValueError(f"{entry.database}: invalid coordinate_strategy: {entry.coordinate_strategy}")
        if entry.coordinate_system not in VALID_COORDINATE_SYSTEMS:
            raise ValueError(f"{entry.database}: invalid coordinate_system: {entry.coordinate_system}")
        if entry.coordinate_strategy == "uniform" and entry.coordinate_system == "mixed":
            raise ValueError(f"{entry.database}: the uniform strategy cannot use mixed")
        if entry.coordinate_strategy == "per_row" and entry.coordinate_system != "mixed":
            raise ValueError(f"{entry.database}: the per_row strategy requires mixed")
    return entries


def select_entries(entries: list[ManifestEntry], requested: list[str]) -> list[ManifestEntry]:
    if not requested:
        return entries
    requested_set = set(requested)
    unknown = sorted(requested_set - {entry.database for entry in entries})
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


def make_uniform_task(entry: ManifestEntry, args: argparse.Namespace) -> Task:
    clean_input = input_path(args.input_root, entry.clean_input)
    decision_script = input_path(args.input_root, entry.decision_script)
    converted = output_path(args.output_root, entry.converted_output)
    unmapped = output_path(args.output_root, entry.unmapped_output)
    conversion_summary = output_path(args.output_root, entry.conversion_summary)
    presence_output = output_path(args.output_root, entry.presence_output)
    presence_audit = output_path(args.output_root, entry.presence_audit)
    final_output = output_path(args.output_root, entry.final_output)
    final_audit = output_path(args.output_root, entry.final_audit)
    conflict_output = output_path(args.output_root, entry.conflict_output)
    final_summary = output_path(args.output_root, entry.final_summary)

    conversion = [
        sys.executable, str(LIFTOVER_SCRIPT), "--input", str(clean_input),
        "--source-coordinate-system", entry.coordinate_system,
        "--target-assembly", "mm10", "--output", str(converted),
        "--unmapped-output", str(unmapped), "--summary-output", str(conversion_summary),
    ]
    if entry.coordinate_system.startswith("mm9_"):
        conversion.extend(["--chain", str(args.chain), "--liftover-bin", args.liftover_bin])
    add_strand = (
        sys.executable, str(ADD_STRAND_SCRIPT), "--input", str(converted),
        "--output", str(presence_output), "--presence", str(args.presence),
        "--ids", str(args.ids), "--audit", str(presence_audit),
    )
    decide = (
        sys.executable, str(decision_script), "--input", str(presence_output),
        "--output", str(final_output), "--audit-output", str(final_audit),
        "--conflict-output", str(conflict_output), "--summary-output", str(final_summary),
    )
    return Task(entry, (tuple(conversion), add_strand, decide), (clean_input, decision_script), final_output)


def make_per_row_task(entry: ManifestEntry, args: argparse.Namespace) -> Task:
    clean_input = input_path(args.input_root, entry.clean_input)
    coordinate_check = input_path(args.input_root, entry.coordinate_check)
    runner = input_path(args.input_root, entry.decision_script)
    final_output = output_path(args.output_root, entry.final_output)
    processed_dir = final_output.parent
    command = (
        sys.executable, str(runner), "--input", str(clean_input),
        "--coordinate-check", str(coordinate_check), "--presence", str(args.presence),
        "--ids", str(args.ids), "--chain", str(args.chain),
        "--liftover-bin", args.liftover_bin, "--processed-dir", str(processed_dir),
    )
    if entry.included_coordinate_classes not in {"", "NA"}:
        command += tuple(
            argument
            for coordinate_class in entry.included_coordinate_classes.split(",")
            for argument in ("--include-coordinate-class", coordinate_class)
        )
    return Task(entry, (command,), (clean_input, coordinate_check, runner), final_output)


def make_task(entry: ManifestEntry, args: argparse.Namespace) -> Task:
    return make_uniform_task(entry, args) if entry.coordinate_strategy == "uniform" else make_per_row_task(entry, args)


def validate_environment(tasks: list[Task], args: argparse.Namespace) -> None:
    if args.output_root == args.input_root:
        raise ValueError("--output-root must differ from --input-root to avoid overwriting active outputs")
    required = [LIFTOVER_SCRIPT, ADD_STRAND_SCRIPT, args.presence, args.ids]
    required.extend(path for task in tasks for path in task.required_inputs)
    missing = sorted(str(path) for path in required if not path.is_file())
    if missing:
        raise FileNotFoundError("Required input or script not found:\n" + "\n".join(missing))
    needs_liftover = any(
        task.entry.coordinate_strategy == "per_row" or task.entry.coordinate_system.startswith("mm9_")
        for task in tasks
    )
    if needs_liftover and not args.dry_run:
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
    for attr in ("manifest", "input_root", "output_root", "presence", "ids", "chain"):
        setattr(args, attr, getattr(args, attr).expanduser().resolve())
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
            except Exception as exc:  # pragma: no cover
                result = TaskResult(task.entry.database, False, str(exc))
            results[result.database] = result
            print(f"{result.database}\t{'completed' if result.success else 'failed'}", flush=True)

    return 1 if any(not result.success for result in results.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
