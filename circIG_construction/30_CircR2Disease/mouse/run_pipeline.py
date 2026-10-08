#!/usr/bin/env python3
"""Run CircR2Disease mouse post-clean pipeline with mixed coordinate handling."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
ROOT = BASE_DIR.parents[1]
PIPELINE_SCRIPTS = ROOT / "00_master" / "pipeline_scripts"
DEFAULT_INPUT = BASE_DIR / "clean" / "30_CircR2Disease.circRNAs.clean.tsv"
DEFAULT_COORDINATE_CHECK = BASE_DIR / "coordinate_check" / "30_CircR2Disease.coordinate_check.tsv"
DEFAULT_PRESENCE = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master.v2.tsv"
DEFAULT_IDS = ROOT / "00_master" / "base_construction" / "mouse" / "circRNA_presence_master_circ_ids.v2.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "mm9ToMm10.over.chain.gz"
DEFAULT_LIFTOVER_BIN = "liftOver"

OUTPUT_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
VALID_STRANDS = {"+", "-"}
CLASS_TO_SOURCE_SYSTEM = {
    "mm10_1based": "mm10_1based",
    "mm10_0based": "mm10_0based",
    "mm9_0based": "mm9_0based",
}
CONVERTIBLE_CLASSES = set(CLASS_TO_SOURCE_SYSTEM)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--coordinate-check", type=Path, default=DEFAULT_COORDINATE_CHECK)
    parser.add_argument("--presence", type=Path, default=DEFAULT_PRESENCE)
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--liftover-bin", default=DEFAULT_LIFTOVER_BIN)
    parser.add_argument("--processed-dir", type=Path, default=BASE_DIR / "processed")
    parser.add_argument(
        "--include-coordinate-class",
        action="append",
        choices=sorted(CONVERTIBLE_CLASSES),
        help=(
            "Only convert the selected row-level coordinate class. May be repeated. "
            "The default converts all supported classes."
        ),
    )
    return parser.parse_args()


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return list(reader.fieldnames), list(reader)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str]:
    return tuple(row.get(column, "") for column in OUTPUT_COLUMNS)


def load_coordinate_check(path: Path) -> dict[tuple[str, str, str, str, str, str], dict[str, str]]:
    header, rows = read_tsv(path)
    missing = [column for column in OUTPUT_COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")
    lookup: dict[tuple[str, str, str, str, str, str], dict[str, str]] = {}
    for row in rows:
        key = normalize_key(row)
        if key in lookup:
            raise ValueError(f"{path}: duplicate coordinate-check key: {key}")
        lookup[key] = row
    return lookup


def split_by_coordinate_class(
    input_path: Path,
    coordinate_check_path: Path,
    subset_dir: Path,
    unmapped_path: Path,
    included_classes: set[str],
) -> dict[str, Path]:
    input_header, input_rows = read_tsv(input_path)
    missing = [column for column in OUTPUT_COLUMNS if column not in input_header]
    if missing:
        raise ValueError(f"{input_path}: missing required columns: {', '.join(missing)}")

    coordinate_check = load_coordinate_check(coordinate_check_path)
    subset_fields = list(input_header)
    for column in ("input_order", "coordinate_check_classification", "coordinate_check_matched_candidate_types"):
        if column not in subset_fields:
            subset_fields.append(column)

    subsets: dict[str, list[dict[str, str]]] = {name: [] for name in included_classes}
    unmapped_rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()

    for index, row in enumerate(input_rows):
        checked = coordinate_check.get(normalize_key(row))
        if checked is None:
            raise ValueError(f"{coordinate_check_path}: missing coordinate-check row for input row {index + 2}")
        classification = checked.get("coordinate_classification", "")
        matched_types = checked.get("matched_candidate_types", "")
        counts[f"class_{classification or 'NA'}"] += 1

        source_class = classification
        note = "from_coordinate_check_classification"
        if classification == "ambiguous" and matched_types == "mm10_1based,mm9_1based":
            source_class = "mm10_1based"
            note = "liftOver_equivalent_ambiguous_treated_as_mm10_1based"

        enriched = dict(row)
        enriched["input_order"] = str(index)
        enriched["coordinate_check_classification"] = classification or "NA"
        enriched["coordinate_check_matched_candidate_types"] = matched_types or "NA"
        if source_class in subsets:
            enriched["conversion_note"] = note
            subsets[source_class].append(enriched)
        else:
            unmapped_rows.append(
                {
                    **enriched,
                    "failure_reason": f"not_convertible_from_coordinate_check:{classification or 'NA'}",
                    "coordinate_check_note": checked.get("notes") or "NA",
                }
            )

    subset_paths: dict[str, Path] = {}
    for source_class, rows in subsets.items():
        path = subset_dir / f"30_CircR2Disease.{source_class}.subset.tsv"
        write_tsv(path, subset_fields, rows)
        subset_paths[source_class] = path

    unmapped_fields = [
        *subset_fields,
        "failure_reason",
        "coordinate_check_note",
    ]
    write_tsv(unmapped_path, unmapped_fields, unmapped_rows)
    return subset_paths


def convert_subsets(args: argparse.Namespace, subset_paths: dict[str, Path], workdir: Path) -> list[Path]:
    converted_paths: list[Path] = []
    for source_class in sorted(subset_paths):
        subset_path = subset_paths[source_class]
        source_system = CLASS_TO_SOURCE_SYSTEM[source_class]
        converted_path = workdir / f"30_CircR2Disease.{source_class}.mm10_1based.tsv"
        unmapped_path = workdir / f"30_CircR2Disease.{source_class}.mm10_1based.unmapped.tsv"
        summary_path = workdir / f"30_CircR2Disease.{source_class}.mm10_1based.summary.tsv"
        cmd = [
            sys.executable,
            str(PIPELINE_SCRIPTS / "liftover_to_1based.py"),
            "--input",
            str(subset_path),
            "--source-coordinate-system",
            source_system,
            "--target-assembly",
            "mm10",
            "--output",
            str(converted_path),
            "--unmapped-output",
            str(unmapped_path),
            "--summary-output",
            str(summary_path),
            "--liftover-bin",
            args.liftover_bin,
            "--workdir",
            str(workdir / f"liftover_{source_class}"),
        ]
        if source_system.startswith("mm9_"):
            cmd.extend(["--chain", str(args.chain)])
        run(cmd)
        converted_paths.append(converted_path)
    return converted_paths


def merge_converted(
    converted_paths: list[Path],
    output_path: Path,
    unmapped_path: Path,
    summary_path: Path,
    initial_unmapped_path: Path,
) -> None:
    converted_rows: list[dict[str, str]] = []
    output_fields: list[str] | None = None
    for path in converted_paths:
        header, rows = read_tsv(path)
        if output_fields is None:
            output_fields = header
        else:
            for column in header:
                if column not in output_fields:
                    output_fields.append(column)
        converted_rows.extend(rows)

    converted_rows.sort(key=lambda row: int(row.get("input_order", "0")))
    if output_fields is None:
        output_fields = OUTPUT_COLUMNS
    write_tsv(output_path, output_fields, converted_rows)

    _, initial_unmapped_rows = read_tsv(initial_unmapped_path)
    all_unmapped_rows = list(initial_unmapped_rows)
    for path in converted_paths:
        liftover_unmapped_path = path.with_name(path.name.replace(".tsv", ".unmapped.tsv"))
        header, rows = read_tsv(liftover_unmapped_path)
        for row in rows:
            all_unmapped_rows.append({**row, "failure_reason": row.get("failure_reason", "liftover_unmapped")})
    unmapped_fields = sorted({key for row in all_unmapped_rows for key in row}) or OUTPUT_COLUMNS
    write_tsv(unmapped_path, unmapped_fields, all_unmapped_rows)

    counts = Counter()
    counts["input_rows"] = len(converted_rows) + len(initial_unmapped_rows)
    counts["converted_rows"] = len(converted_rows)
    counts["unmapped_rows"] = len(all_unmapped_rows)
    counts["initial_unmapped_rows"] = len(initial_unmapped_rows)
    write_tsv(
        summary_path,
        ["metric", "value"],
        [{"metric": key, "value": str(value)} for key, value in sorted(counts.items())],
    )


def clean_strand(value: str | None) -> str:
    if value is None:
        return "NA"
    value = value.strip()
    return value if value in VALID_STRANDS else "NA"


def decide_strand(row: dict[str, str]) -> tuple[str, str]:
    original = clean_strand(row.get("strand"))
    presence_master = clean_strand(row.get("strand_presence_master"))
    if original in VALID_STRANDS:
        return original, "original"
    if presence_master in VALID_STRANDS:
        return presence_master, "presence_master"
    return "NA", "unresolved"


def build_final(input_path: Path, output_path: Path, conflict_path: Path, summary_path: Path) -> None:
    header, rows = read_tsv(input_path)
    missing = [column for column in OUTPUT_COLUMNS if column not in header]
    if missing:
        raise ValueError(f"{input_path}: missing required columns: {', '.join(missing)}")

    counts: Counter[str] = Counter()
    final_rows: list[dict[str, str]] = []
    conflict_rows: list[dict[str, str]] = []
    seen: set[tuple[str, ...]] = set()

    for row in rows:
        counts["input_rows"] += 1
        final_strand, method = decide_strand(row)
        output_row = {column: row.get(column, "NA") or "NA" for column in OUTPUT_COLUMNS}
        output_row["strand"] = final_strand
        key = tuple(output_row[column] for column in OUTPUT_COLUMNS)
        if key in seen:
            counts["duplicate_full_rows_removed"] += 1
            continue
        seen.add(key)
        final_rows.append(output_row)
        counts["final_rows"] += 1
        counts[f"method_{method}"] += 1
        if final_strand in VALID_STRANDS:
            counts["final_stranded_rows"] += 1
        else:
            counts["final_unstranded_rows"] += 1
        if method == "conflict":
            conflict_rows.append({**row, "final_strand": final_strand, "strand_decision_method": method})

    write_tsv(output_path, OUTPUT_COLUMNS, final_rows)
    write_tsv(conflict_path, [*header, "final_strand", "strand_decision_method"], conflict_rows)
    summary = {
        "input_rows": str(counts["input_rows"]),
        "duplicate_full_rows_removed": str(counts["duplicate_full_rows_removed"]),
        "final_rows": str(counts["final_rows"]),
        "final_stranded_rows": str(counts["final_stranded_rows"]),
        "final_unstranded_rows": str(counts["final_unstranded_rows"]),
        "original_rows": str(counts["method_original"]),
        "presence_master_rows": str(counts["method_presence_master"]),
        "conflict_rows": str(counts["method_conflict"]),
        "unresolved_rows": str(counts["method_unresolved"]),
        "final_output": str(output_path),
        "conflict_output": str(conflict_path),
    }
    write_tsv(summary_path, list(summary), [summary])


def main() -> int:
    args = parse_args()
    processed_dir = args.processed_dir.expanduser().resolve()
    workdir = processed_dir / "mixed_coordinate_work"
    subset_dir = workdir / "subsets"
    processed_dir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    initial_unmapped = workdir / "30_CircR2Disease.not_convertible_from_coordinate_check.tsv"
    subset_paths = split_by_coordinate_class(
        input_path=args.input.expanduser().resolve(),
        coordinate_check_path=args.coordinate_check.expanduser().resolve(),
        subset_dir=subset_dir,
        unmapped_path=initial_unmapped,
        included_classes=set(args.include_coordinate_class or CONVERTIBLE_CLASSES),
    )
    converted_paths = convert_subsets(args, subset_paths, workdir)

    converted_output = processed_dir / "30_CircR2Disease.circRNAs.mm10_1based.tsv"
    unmapped_output = processed_dir / "30_CircR2Disease.circRNAs.mm10_1based.unmapped.tsv"
    conversion_summary = processed_dir / "30_CircR2Disease.circRNAs.mm10_1based.summary.tsv"
    merge_converted(converted_paths, converted_output, unmapped_output, conversion_summary, initial_unmapped)

    presence_output = processed_dir / "30_CircR2Disease.circRNAs.with_strand_presence_master.tsv"
    presence_audit = processed_dir / "30_CircR2Disease.circRNAs.with_strand_presence_master.audit.tsv"
    run(
        [
            sys.executable,
            str(PIPELINE_SCRIPTS / "add_strand_from_presence_master.py"),
            "--input",
            str(converted_output),
            "--output",
            str(presence_output),
            "--presence",
            str(args.presence.expanduser().resolve()),
            "--ids",
            str(args.ids.expanduser().resolve()),
            "--audit",
            str(presence_audit),
        ]
    )

    build_final(
        input_path=presence_output,
        output_path=processed_dir / "30_CircR2Disease.circRNAs.final.tsv",
        conflict_path=processed_dir / "30_CircR2Disease.circRNAs.strand_conflicts.tsv",
        summary_path=processed_dir / "30_CircR2Disease.circRNAs.final.summary.tsv",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
