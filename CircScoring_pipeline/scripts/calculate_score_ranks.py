import csv
import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np


REFERENCE_VERSION = 3
MISSING_VALUES = {"", ".", "NA", "NaN", "nan"}


def iter_chunks(reader, chunk_size):
    chunk = []
    for row in reader:
        chunk.append(row)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def calculate_chunk(rows, score_columns, sorted_scores, row_count, first_line):
    ranks = [["."] * len(score_columns) for _ in rows]
    for column_index, column in enumerate(score_columns):
        positions = []
        values = []
        for row_index, row in enumerate(rows):
            raw_value = row[column_index + 1].strip()
            if raw_value in MISSING_VALUES:
                continue
            try:
                value = float(raw_value)
            except ValueError as error:
                raise ValueError(
                    f"Invalid {column} value at score line "
                    f"{first_line + row_index}: {raw_value!r}"
                ) from error
            if not math.isfinite(value):
                raise ValueError(
                    f"Non-finite {column} value at score line "
                    f"{first_line + row_index}: {raw_value!r}"
                )
            positions.append(row_index)
            values.append(value)

        if not values:
            continue
        query = np.asarray(values, dtype=np.float64)
        left = np.searchsorted(sorted_scores[column_index], query, side="left")
        right = np.searchsorted(sorted_scores[column_index], query, side="right")
        equal_counts = right - left
        rank_positions = left.astype(np.float64)
        matched = equal_counts > 0
        rank_positions[matched] += 0.5 * (equal_counts[matched] + 1)
        percentiles = (rank_positions / row_count) * 100.0
        for row_index, percentile in zip(positions, percentiles):
            ranks[row_index][column_index] = f"{percentile:.10f}"
    return ranks


score_columns = list(snakemake.params.score_columns)
rank_columns = list(snakemake.params.rank_columns)
array_paths = [Path(path) for path in snakemake.input.arrays]
chunk_size = int(snakemake.params.chunk_size)

if len(score_columns) != len(rank_columns) or len(score_columns) != len(array_paths):
    raise ValueError("Score, rank, and reference-array counts must match")
if chunk_size <= 0:
    raise ValueError("Score-rank chunk size must be positive")

with open(snakemake.input.metadata, encoding="utf-8") as handle:
    metadata = json.load(handle)
if metadata.get("version") != REFERENCE_VERSION:
    raise ValueError(f"Unsupported score-rank reference metadata: {metadata!r}")
if metadata.get("score_columns") != score_columns:
    raise ValueError(
        "Score-rank reference columns do not match the requested score columns"
    )
row_count = metadata.get("row_count")
if not isinstance(row_count, int) or row_count <= 0:
    raise ValueError(f"Invalid score-rank reference row count: {row_count!r}")

sorted_scores = [np.load(path, mmap_mode="r", allow_pickle=False) for path in array_paths]
for column, values in zip(score_columns, sorted_scores):
    if values.dtype != np.float64 or values.shape != (row_count,):
        raise ValueError(
            f"Invalid sorted reference for {column}: "
            f"dtype={values.dtype}, shape={values.shape}"
        )

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
temporary_path = None
total_rows = 0
missing_counts = [0] * len(score_columns)
try:
    with open(snakemake.input.scores, newline="", encoding="utf-8") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        expected_header = ["circRNA_id", *score_columns]
        header = next(reader, None)
        if header != expected_header:
            raise ValueError(
                f"Unexpected predictCS score header: expected "
                f"{expected_header!r}, got {header!r}"
            )

        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=output_path.parent,
            prefix=f".{output_path.name}.", suffix=".tmp", delete=False,
        ) as output_handle:
            temporary_path = Path(output_handle.name)
            writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["circRNA_id", *rank_columns])
            for rows in iter_chunks(reader, chunk_size):
                first_line = total_rows + 2
                for offset, row in enumerate(rows):
                    if len(row) != len(expected_header):
                        raise ValueError(
                            f"Invalid predictCS score line {first_line + offset}: "
                            f"expected {len(expected_header)} columns, got {len(row)}"
                        )
                    if not row[0]:
                        raise ValueError(
                            f"Empty circRNA_id at score line {first_line + offset}"
                        )
                ranks = calculate_chunk(
                    rows, score_columns, sorted_scores, row_count, first_line
                )
                for row, rank_values in zip(rows, ranks):
                    writer.writerow([row[0], *rank_values])
                    for index, value in enumerate(rank_values):
                        if value == ".":
                            missing_counts[index] += 1
                total_rows += len(rows)
            output_handle.flush()
            os.fsync(output_handle.fileno())
    os.replace(temporary_path, output_path)
    temporary_path = None
finally:
    if temporary_path is not None:
        temporary_path.unlink(missing_ok=True)

log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w", encoding="utf-8") as log:
    log.write(f"reference_rows={row_count}\n")
    log.write(f"score_rows={total_rows}\n")
    for column, missing_count in zip(score_columns, missing_counts):
        log.write(f"missing_{column}={missing_count}\n")
