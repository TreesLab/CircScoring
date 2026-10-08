import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import (
    CHECKAA_COLUMNS,
    DETAILED_FEATURE_COLUMNS,
    PREDICTION_COLUMNS,
    PREDICTCS_FEATURE_COLUMNS,
    RAW_COLUMNS,
    atomic_text_writer,
    finish_atomic_write,
)
from scripts.score_rank_schema import RANK_COLUMNS


RESULT_KEY_COLUMNS = ["chr", "pos1", "pos2", "strand", "circRNA_id"]


def open_table(path, expected_header, label):
    handle = open(path, newline="", encoding="utf-8")
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader, None)
    if header != expected_header:
        handle.close()
        raise ValueError(
            f"Unexpected {label} header in {path}: expected "
            f"{expected_header!r}, got {header!r}"
        )
    return handle, reader


def next_row(reader, expected_columns, label, row_number):
    row = next(reader, None)
    if row is None:
        raise ValueError(f"Missing {label} row {row_number}")
    if len(row) != len(expected_columns):
        raise ValueError(
            f"Invalid {label} row {row_number}: expected "
            f"{len(expected_columns)} columns, got {len(row)}"
        )
    return row


def ensure_exhausted(reader, label):
    extra = next(reader, None)
    if extra is not None:
        raise ValueError(f"Unexpected extra {label} row: {extra!r}")


def merge_analysis_results(detailed, output_handle):
    feature_columns = (
        DETAILED_FEATURE_COLUMNS if detailed else PREDICTCS_FEATURE_COLUMNS
    )
    score_columns = ["circRNA_id", *PREDICTION_COLUMNS]
    feature_handle, feature_reader = open_table(
        snakemake.input.features, feature_columns, "feature"
    )
    score_handle, score_reader = open_table(
        snakemake.input.scores, score_columns, "predictCS score"
    )
    rank_columns = ["circRNA_id", *RANK_COLUMNS]
    rank_handle, rank_reader = open_table(
        snakemake.input.ranks, rank_columns, "predictCS score rank"
    )

    try:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                *RESULT_KEY_COLUMNS,
                *feature_columns[1:],
                *PREDICTION_COLUMNS,
                *RANK_COLUMNS,
            ]
        )
        with open(
            snakemake.input.circRNAs, newline="", encoding="utf-8"
        ) as circ_handle:
            circ_reader = csv.reader(circ_handle, delimiter="\t")
            for row_number, circ_row in enumerate(circ_reader, start=1):
                if len(circ_row) != len(RAW_COLUMNS):
                    raise ValueError(
                        f"Invalid active circRNA row {row_number}: expected "
                        f"{len(RAW_COLUMNS)} columns, got {len(circ_row)}"
                    )
                feature_row = next_row(
                    feature_reader, feature_columns, "feature", row_number
                )
                score_row = next_row(
                    score_reader, score_columns, "predictCS score", row_number
                )
                rank_row = next_row(
                    rank_reader, rank_columns, "predictCS score rank", row_number
                )
                circ_id = circ_row[4]
                if (
                    feature_row[0] != circ_id
                    or score_row[0] != circ_id
                    or rank_row[0] != circ_id
                ):
                    raise ValueError(
                        f"Result ID mismatch at row {row_number}: expected "
                        f"{circ_id!r}, feature={feature_row[0]!r}, "
                        f"predictCS={score_row[0]!r}, rank={rank_row[0]!r}"
                    )
                writer.writerow(
                    [*circ_row, *feature_row[1:], *score_row[1:], *rank_row[1:]]
                )
        ensure_exhausted(feature_reader, "feature")
        ensure_exhausted(score_reader, "predictCS score")
        ensure_exhausted(rank_reader, "predictCS score rank")
    finally:
        feature_handle.close()
        score_handle.close()
        rank_handle.close()


def append_checkaa_results(detailed, output_handle):
    feature_columns = (
        DETAILED_FEATURE_COLUMNS if detailed else PREDICTCS_FEATURE_COLUMNS
    )
    base_columns = [
        *RESULT_KEY_COLUMNS,
        *feature_columns[1:],
        *PREDICTION_COLUMNS,
        *RANK_COLUMNS,
    ]
    base_handle, base_reader = open_table(
        snakemake.input.base, base_columns, "base result"
    )
    checkaa_handle, checkaa_reader = open_table(
        snakemake.input.checkAA, CHECKAA_COLUMNS, "checkAA"
    )

    try:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([*base_columns, *CHECKAA_COLUMNS[1:]])
        for row_number, base_row in enumerate(base_reader, start=1):
            if len(base_row) != len(base_columns):
                raise ValueError(
                    f"Invalid base result row {row_number}: expected "
                    f"{len(base_columns)} columns, got {len(base_row)}"
                )
            checkaa_row = next_row(
                checkaa_reader, CHECKAA_COLUMNS, "checkAA", row_number
            )
            if checkaa_row[0] != base_row[4]:
                raise ValueError(
                    f"checkAA ID mismatch at row {row_number}: expected "
                    f"{base_row[4]!r}, got {checkaa_row[0]!r}"
                )
            writer.writerow([*base_row, *checkaa_row[1:]])
        ensure_exhausted(checkaa_reader, "checkAA")
    finally:
        base_handle.close()
        checkaa_handle.close()


mode = str(snakemake.params.mode)
detailed = bool(snakemake.params.detailed)
output_handle, temporary_path, output_path = atomic_text_writer(
    snakemake.output[0], newline=""
)

try:
    if mode == "analysis":
        merge_analysis_results(detailed, output_handle)
    elif mode == "checkAA":
        append_checkaa_results(detailed, output_handle)
    else:
        raise ValueError(f"Unsupported result merge mode: {mode!r}")
    finish_atomic_write(output_handle, temporary_path, output_path)
except Exception:
    output_handle.close()
    temporary_path.unlink(missing_ok=True)
    raise
