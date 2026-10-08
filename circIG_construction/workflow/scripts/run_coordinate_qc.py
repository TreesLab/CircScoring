#!/usr/bin/env python3
"""Rerun optional coordinate-system QC from a locked tier manifest."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", choices=("human", "mouse"), required=True)
    parser.add_argument("--tier", type=int, choices=(2, 3, 4), required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--base-dataset", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--liftover-bin", default="liftOver")
    return parser.parse_args()


def safe_label(database: str, subset: str) -> str:
    return database if subset in {"", "NA"} else f"{database}.{subset}"


def main() -> int:
    args = parse_args()
    manifest = args.manifest.expanduser().resolve()
    input_root = args.input_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checker = input_root / "00_master" / "pipeline_scripts" / (
        "check_circrna_coordinate_system.py"
        if args.species == "human"
        else "check_mouse_circrna_coordinate_system.py"
    )
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    completed = 0
    for row in rows:
        if row.get("decision_method") == "reference_lookup":
            continue
        source = Path(row["clean_input"])
        source = source if source.is_absolute() else input_root / source
        label = safe_label(row["database"], row.get("subset", "NA"))
        command = [
            sys.executable, str(checker), "--input", str(source),
            "--base-dataset", str(args.base_dataset),
            "--output", str(output_dir / f"{label}.coordinate_check.tsv"),
            "--summary-tsv-output", str(output_dir / f"{label}.coordinate_check.summary.tsv"),
            "--chain", str(args.chain), "--liftover-bin", args.liftover_bin,
            "--base-id-col", "",
        ]
        if args.tier >= 3:
            command.extend(["--strand-col", "", "--base-strand-col", ""])
        subprocess.run(command, check=True)
        completed += 1
    (output_dir / ".complete").write_text(f"completed\t{completed}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
