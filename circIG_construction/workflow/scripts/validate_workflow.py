#!/usr/bin/env python3
"""Static validation for workflow configuration and repository inputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--allow-missing-clean-inputs",
        action="store_true",
        help="Validate workflow structure only; do not fail for clean inputs that have not been generated.",
    )
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(path: Path, errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f"File not found: {path}")


def source_names(path: Path) -> list[str]:
    return [row["source"] for row in read_tsv(path)]


def validate_tier_growth(paths: list[Path], label: str, errors: list[str]) -> None:
    previous: set[str] = set()
    for tier, path in enumerate(paths, start=1):
        current = set(source_names(path))
        if not previous.issubset(current):
            errors.append(f"{label} v{tier} is missing sources from the previous version: {sorted(previous - current)}")
        previous = current


def main() -> int:
    args = parse_args()
    repo = args.repo.expanduser().resolve()
    config_path = args.config.expanduser().resolve()
    errors: list[str] = []
    try:
        config = json.loads(config_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"Unable to read workflow config: {exc}")

    required_scripts = (
        "00_master/base_construction/build_circRNA_presence.py",
        "00_master/pipeline_scripts/liftover_to_hg38_1based.py",
        "00_master/pipeline_scripts/liftover_to_1based.py",
        "00_master/pipeline_scripts/build_human_tier2_finals.py",
        "00_master/pipeline_scripts/build_human_tier3_finals.py",
        "00_master/pipeline_scripts/build_human_tier4_finals.py",
        "00_master/pipeline_scripts/build_mouse_tier2_finals.py",
        "00_master/pipeline_scripts/build_mouse_tier3_finals.py",
        "00_master/pipeline_scripts/build_mouse_tier4_finals.py",
        "00_master/pipeline_scripts/build_mouse_conservation.py",
    )
    for value in required_scripts:
        require(repo / value, errors)
    for value in config["chains"].values():
        require(repo / value, errors)
    require(repo / config["mouse_partial_mapping_overrides"], errors)

    human_sources = [repo / f"00_master/base_construction/database_sources.v{tier}.tsv" for tier in range(1, 5)]
    mouse_sources = [repo / f"00_master/base_construction/mouse/database_sources.v{tier}.tsv" for tier in range(1, 5)]
    for path in (*human_sources, *mouse_sources):
        require(path, errors)
    if not errors:
        validate_tier_growth(human_sources, "human", errors)
        validate_tier_growth(mouse_sources, "mouse", errors)

    for manifest in (
        "workflow/config/human_tier1.tsv", "workflow/config/mouse_tier1.tsv",
        "00_master/pipeline_configs/human_tier2_coordinate_systems.tsv",
        "00_master/pipeline_configs/human_tier3_coordinate_systems.tsv",
        "00_master/pipeline_configs/human_tier4_coordinate_systems.tsv",
        "00_master/pipeline_configs/mouse_tier2_coordinate_systems.tsv",
        "workflow/config/mouse_tier3_coordinate_systems.tsv",
        "workflow/config/mouse_tier4_coordinate_systems.tsv",
    ):
        path = repo / manifest
        require(path, errors)
        if not path.is_file():
            continue
        rows = read_tsv(path)
        labels = [(row.get("database"), row.get("subset", "NA")) for row in rows]
        if len(labels) != len(set(labels)):
            errors.append(f"Manifest contains duplicate database/subset pairs: {path}")
        for row in rows:
            if row["database"] in {"21_CircTarget", "43_circRNADisease"} and "human_tier2" in manifest:
                continue
            clean = repo / row["clean_input"]
            if not clean.is_file() and not args.allow_missing_clean_inputs:
                errors.append(f"Clean input has not been generated: {clean}")

    if not args.allow_missing_clean_inputs:
        require(repo / "21_CircTarget/clean/21_CircTarget.circRNA_ids.clean.tsv", errors)
        require(repo / "43_circRNADisease/clean/43_circRNADisease.circRNA_ids.clean.tsv", errors)

    if errors:
        print("\n".join(errors))
        return 1
    print("workflow static validation: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
