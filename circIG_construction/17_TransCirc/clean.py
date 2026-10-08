#!/usr/bin/env python3
"""Validate TransCirc metadata, select representative isoforms, and clean BSJs."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import TextIO


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE_DIR / "raw" / "transcirc_metadata.tsv"
DEFAULT_PREPROCESSED_OUTPUT = BASE_DIR / "preprocessed" / "17_TransCirc.raw.uniq_bsj.tsv"
DEFAULT_OUTPUT = BASE_DIR / "clean" / "17_TransCirc.circRNAs.clean.tsv"
DEFAULT_SUMMARY_OUTPUT = BASE_DIR / "clean" / "17_TransCirc.circRNAs.clean.summary.tsv"

RAW_COLUMNS = [
    "TransCirc_ID",
    "xref",
    "species",
    "gene",
    "strand",
    "gene_id",
    "chrom",
    "start",
    "end",
    "evidences_num",
    "evidences_score",
    "MS_Qvalue",
    "m6A",
    "RP/RP",
    "TIS",
    "ORF",
    "IRES",
    "Peptide composition",
    "MS",
]
EVIDENCE_COLUMNS = ["m6A", "RP/RP", "TIS", "ORF", "IRES", "Peptide composition", "MS"]
PREPROCESSED_COLUMNS = [*RAW_COLUMNS, "isoform_count", "all_TransCirc_IDs"]
CLEAN_COLUMNS = ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"]
SUMMARY_COLUMNS = [
    "input_file",
    "preprocessed_output",
    "clean_output",
    "raw_rows",
    "unique_transcirc_ids",
    "unique_gene_symbols",
    "unique_bsj_rows",
    "duplicate_bsj_groups",
    "rows_in_duplicate_bsj_groups",
    "rows_removed_by_bsj_selection",
    "representative_score_tie_groups",
    "gene_transcirc_id_mismatches",
    "evidence_formula_errors",
    "ms_qvalue_relation_errors",
    "stranded_rows",
    "unstranded_rows",
]

ID_PATTERN = re.compile(r"^TC-hsa-(.+)_([0-9]{4})$")
BSJKey = tuple[str, str, str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate TransCirc metadata and generate a six-column unique-BSJ clean table."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help=f"Raw TSV. Default: {DEFAULT_INPUT}")
    parser.add_argument(
        "--preprocessed-output",
        type=Path,
        default=DEFAULT_PREPROCESSED_OUTPUT,
        help=f"Representative-isoform audit TSV. Default: {DEFAULT_PREPROCESSED_OUTPUT}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Six-column clean TSV. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Summary TSV. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def normalize_missing(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return text if text else "NA"


def normalize_chrom(value: str, line_number: int) -> str:
    chrom = normalize_missing(value)
    if chrom == "NA":
        raise ValueError(f"Missing chrom at line {line_number}")
    if chrom.lower().startswith("chr"):
        return "chr" + chrom[3:]
    return f"chr{chrom}"


def parse_number(value: str, column: str, line_number: int) -> float:
    text = normalize_missing(value)
    if text == "NA":
        raise ValueError(f"Missing {column} at line {line_number}")
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(f"Invalid {column} at line {line_number}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Non-finite {column} at line {line_number}: {value!r}")
    return number


def parse_integer(value: str, column: str, line_number: int) -> int:
    number = parse_number(value, column, line_number)
    if not number.is_integer():
        raise ValueError(f"Non-integer {column} at line {line_number}: {value!r}")
    return int(number)


def validate_header(fieldnames: list[str] | None) -> None:
    if fieldnames != RAW_COLUMNS:
        raise ValueError(f"Unexpected raw TSV header: {fieldnames!r}; expected {RAW_COLUMNS!r}")


def normalized_raw_row(row: dict[str, str], line_number: int) -> dict[str, str]:
    if None in row:
        raise ValueError(f"Extra TSV field(s) at line {line_number}: {row[None]!r}")
    normalized = {column: normalize_missing(row.get(column)) for column in RAW_COLUMNS}
    normalized["chrom"] = normalize_chrom(normalized["chrom"], line_number)
    return normalized


def validate_row(row: dict[str, str], line_number: int) -> tuple[BSJKey, float, int]:
    circ_id = row["TransCirc_ID"]
    match = ID_PATTERN.fullmatch(circ_id)
    if match is None:
        raise ValueError(f"Invalid TransCirc_ID at line {line_number}: {circ_id!r}")
    if row["species"] != "Homo sapiens":
        raise ValueError(f"Unexpected species at line {line_number}: {row['species']!r}")
    if row["gene"] == "NA":
        raise ValueError(f"Missing gene at line {line_number}")
    if row["gene"] != match.group(1):
        raise ValueError(
            f"Gene/TransCirc_ID mismatch at line {line_number}: "
            f"gene={row['gene']!r}, ID gene={match.group(1)!r}"
        )
    if row["strand"] not in {"+", "-"}:
        raise ValueError(f"Invalid strand at line {line_number}: {row['strand']!r}")

    start = parse_integer(row["start"], "start", line_number)
    end = parse_integer(row["end"], "end", line_number)
    if start < 0 or end < 0:
        raise ValueError(f"Negative coordinate at line {line_number}: {start}-{end}")
    if start > end:
        raise ValueError(f"start > end at line {line_number}: {start} > {end}")
    row["start"] = str(start)
    row["end"] = str(end)

    evidence_values = [parse_number(row[column], column, line_number) for column in EVIDENCE_COLUMNS]
    evidence_num = parse_integer(row["evidences_num"], "evidences_num", line_number)
    evidence_score = parse_number(row["evidences_score"], "evidences_score", line_number)
    calculated_num = sum(value > 0 for value in evidence_values)
    if evidence_num != calculated_num:
        raise ValueError(
            f"evidences_num formula mismatch at line {line_number}: "
            f"reported={evidence_num}, calculated={calculated_num}"
        )
    calculated_score = sum(evidence_values)
    if not math.isclose(evidence_score, calculated_score, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(
            f"evidences_score formula mismatch at line {line_number}: "
            f"reported={evidence_score}, calculated={calculated_score}"
        )

    ms = parse_number(row["MS"], "MS", line_number)
    if ms not in {0, 1}:
        raise ValueError(f"Invalid MS indicator at line {line_number}: {row['MS']!r}")
    qvalue_present = row["MS_Qvalue"] != "NA"
    if (ms == 1) != qvalue_present:
        raise ValueError(
            f"MS/MS_Qvalue relation mismatch at line {line_number}: "
            f"MS={row['MS']!r}, MS_Qvalue={row['MS_Qvalue']!r}"
        )
    if qvalue_present:
        parse_number(row["MS_Qvalue"], "MS_Qvalue", line_number)

    return (row["chrom"], row["start"], row["end"], row["strand"]), evidence_score, evidence_num


def temporary_output(target: Path) -> tuple[TextIO, Path]:
    target.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        prefix=f".{target.name}.",
        suffix=".tmp",
        dir=target.parent,
        delete=False,
    )
    return handle, Path(handle.name)


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError("Output exists; use --force to overwrite: " + ", ".join(existing))


def inspect_raw(input_path: Path) -> tuple[dict[BSJKey, list[object]], dict[BSJKey, list[str]], Counter[str]]:
    counts: Counter[str] = Counter()
    seen_ids: set[str] = set()
    genes: set[str] = set()
    # Each key maps to isoform count, best score, evidence count, source order,
    # representative ID, and the number of tied top scores.
    states: dict[BSJKey, list[object]] = {}
    duplicate_ids: dict[BSJKey, list[str]] = {}

    with input_path.open(newline="", encoding="utf-8", errors="strict") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(reader.fieldnames)
        for source_order, raw_row in enumerate(reader, start=1):
            line_number = source_order + 1
            counts["raw_rows"] += 1
            row = normalized_raw_row(raw_row, line_number)
            key, score, evidence_num = validate_row(row, line_number)

            circ_id = row["TransCirc_ID"]
            if circ_id in seen_ids:
                raise ValueError(f"Duplicate TransCirc_ID at line {line_number}: {circ_id!r}")
            seen_ids.add(circ_id)
            genes.add(row["gene"])

            state = states.get(key)
            if state is None:
                states[key] = [1, score, evidence_num, source_order, circ_id, 1]
                continue

            state[0] = int(state[0]) + 1
            if state[0] == 2:
                duplicate_ids[key] = [str(state[4]), circ_id]
            else:
                duplicate_ids[key].append(circ_id)

            if score > float(state[1]):
                state[5] = 1
            elif score == float(state[1]):
                state[5] = int(state[5]) + 1

            candidate_rank = (score, evidence_num, -source_order)
            representative_rank = (float(state[1]), int(state[2]), -int(state[3]))
            if candidate_rank > representative_rank:
                state[1:5] = [score, evidence_num, source_order, circ_id]

    counts["unique_transcirc_ids"] = len(seen_ids)
    counts["unique_gene_symbols"] = len(genes)
    counts["unique_bsj_rows"] = len(states)
    counts["duplicate_bsj_groups"] = len(duplicate_ids)
    counts["rows_in_duplicate_bsj_groups"] = sum(int(states[key][0]) for key in duplicate_ids)
    counts["rows_removed_by_bsj_selection"] = counts["raw_rows"] - len(states)
    counts["representative_score_tie_groups"] = sum(int(states[key][5]) > 1 for key in duplicate_ids)
    counts["gene_transcirc_id_mismatches"] = 0
    counts["evidence_formula_errors"] = 0
    counts["ms_qvalue_relation_errors"] = 0
    counts["stranded_rows"] = len(states)
    counts["unstranded_rows"] = 0
    return states, duplicate_ids, counts


def write_outputs(
    input_path: Path,
    preprocessed_path: Path,
    clean_path: Path,
    summary_path: Path,
    states: dict[BSJKey, list[object]],
    duplicate_ids: dict[BSJKey, list[str]],
    counts: Counter[str],
) -> None:
    temporary_paths: list[Path] = []
    try:
        preprocessed_handle, preprocessed_tmp = temporary_output(preprocessed_path)
        clean_handle, clean_tmp = temporary_output(clean_path)
        temporary_paths.extend([preprocessed_tmp, clean_tmp])
        written_keys: set[BSJKey] = set()

        with (
            input_path.open(newline="", encoding="utf-8", errors="strict") as input_handle,
            preprocessed_handle,
            clean_handle,
        ):
            reader = csv.DictReader(input_handle, delimiter="\t")
            validate_header(reader.fieldnames)
            preprocessed_writer = csv.DictWriter(
                preprocessed_handle,
                fieldnames=PREPROCESSED_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            clean_writer = csv.DictWriter(
                clean_handle,
                fieldnames=CLEAN_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            preprocessed_writer.writeheader()
            clean_writer.writeheader()

            for source_order, raw_row in enumerate(reader, start=1):
                line_number = source_order + 1
                row = normalized_raw_row(raw_row, line_number)
                key, _score, _evidence_num = validate_row(row, line_number)
                state = states[key]
                if row["TransCirc_ID"] != state[4]:
                    continue
                if key in written_keys:
                    raise ValueError(f"Representative BSJ written more than once: {key!r}")
                written_keys.add(key)

                preprocessed_row = {
                    **row,
                    "isoform_count": str(state[0]),
                    "all_TransCirc_IDs": ",".join(duplicate_ids.get(key, [row["TransCirc_ID"]])),
                }
                preprocessed_writer.writerow(preprocessed_row)
                clean_writer.writerow(
                    {
                        "circRNA_id": row["TransCirc_ID"],
                        "chrom": row["chrom"],
                        "start": row["start"],
                        "end": row["end"],
                        "strand": row["strand"],
                        "gene_symbol": row["gene"],
                    }
                )

        if len(written_keys) != len(states):
            raise ValueError(
                f"Representative output count mismatch: wrote {len(written_keys)}, expected {len(states)}"
            )

        summary_handle, summary_tmp = temporary_output(summary_path)
        temporary_paths.append(summary_tmp)
        with summary_handle:
            writer = csv.DictWriter(
                summary_handle,
                fieldnames=SUMMARY_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            summary = {column: str(counts[column]) for column in SUMMARY_COLUMNS}
            summary.update(
                {
                    "input_file": str(input_path),
                    "preprocessed_output": str(preprocessed_path),
                    "clean_output": str(clean_path),
                }
            )
            writer.writerow(summary)

        for temporary_path, target_path in [
            (preprocessed_tmp, preprocessed_path),
            (clean_tmp, clean_path),
            (summary_tmp, summary_path),
        ]:
            temporary_path.chmod(0o644)
            temporary_path.replace(target_path)
    except Exception:
        for temporary_path in temporary_paths:
            temporary_path.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    preprocessed_path = args.preprocessed_output.expanduser().resolve()
    clean_path = args.output.expanduser().resolve()
    summary_path = args.summary_output.expanduser().resolve()
    output_paths = [preprocessed_path, clean_path, summary_path]

    if not input_path.is_file():
        raise FileNotFoundError(f"Raw TSV does not exist: {input_path}")
    if len(set(output_paths)) != len(output_paths):
        raise ValueError("Output paths must be distinct")
    refuse_overwrite(output_paths, args.force)

    states, duplicate_ids, counts = inspect_raw(input_path)
    write_outputs(
        input_path=input_path,
        preprocessed_path=preprocessed_path,
        clean_path=clean_path,
        summary_path=summary_path,
        states=states,
        duplicate_ids=duplicate_ids,
        counts=counts,
    )

    for key in [
        "raw_rows",
        "unique_transcirc_ids",
        "unique_bsj_rows",
        "duplicate_bsj_groups",
        "rows_removed_by_bsj_selection",
        "representative_score_tie_groups",
    ]:
        print(f"{key}: {counts[key]}")
    print(f"preprocessed_output: {preprocessed_path}")
    print(f"clean_output: {clean_path}")
    print(f"summary_output: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
