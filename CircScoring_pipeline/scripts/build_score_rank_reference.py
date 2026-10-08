import csv
import json
import logging
import os
import tempfile
from pathlib import Path

import numpy as np


REFERENCE_VERSION = 3


def atomic_save_array(path, values):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.",
            suffix=".tmp", delete=False,
        ) as output:
            temporary_path = Path(output.name)
            np.save(output, values, allow_pickle=False)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def atomic_write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as output:
            temporary_path = Path(output.name)
            json.dump(payload, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def read_header(path):
    with open(path, newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if header is None:
        raise ValueError(f"Empty score-rank reference table: {path}")
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate columns in score-rank reference table: {header!r}")
    return header


log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=log_path,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger("build_score_rank_reference")

input_path = Path(snakemake.input.table)
score_columns = list(snakemake.params.score_columns)
array_paths = [Path(path) for path in snakemake.output.arrays]
expected_rows = int(snakemake.params.expected_rows)

if len(score_columns) != len(array_paths):
    raise ValueError(
        "The number of score columns must match the number of output arrays"
    )

logger.info("Reading score-rank reference table: %s", input_path)
header = read_header(input_path)
missing_columns = [column for column in score_columns if column not in header]
if missing_columns:
    raise ValueError(
        f"Missing score columns in {input_path}: {missing_columns!r}"
    )
use_columns = tuple(header.index(column) for column in score_columns)

values = np.loadtxt(
    input_path,
    delimiter="\t",
    dtype=np.float64,
    skiprows=1,
    usecols=use_columns,
    ndmin=2,
)
if values.shape != (expected_rows, len(score_columns)):
    raise ValueError(
        f"Unexpected score-rank reference shape: expected "
        f"{(expected_rows, len(score_columns))!r}, got {values.shape!r}"
    )
if not np.isfinite(values).all():
    raise ValueError("Score-rank reference contains non-finite score values")

logger.info("Sorting %d rows for %d score columns", *values.shape)
values.sort(axis=0)
for column_index, (column, path) in enumerate(zip(score_columns, array_paths)):
    atomic_save_array(path, values[:, column_index])
    logger.info("Wrote sorted index for %s: %s", column, path)

metadata = {
    "version": REFERENCE_VERSION,
    "row_count": expected_rows,
    "dtype": "float64",
    "score_columns": score_columns,
    "rank_formula": (
        "if count_equal > 0: ((count_less + 0.5 * (count_equal + 1)) / N) "
        "* 100; else: (count_less / N) * 100"
    ),
}
atomic_write_json(snakemake.output.metadata, metadata)
logger.info("Score-rank reference index completed")
