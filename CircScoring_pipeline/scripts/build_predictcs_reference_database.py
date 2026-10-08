#!/usr/bin/env python3

import csv
import os
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import parse_coordinate


REFERENCE_SCHEMA_VERSION = 1
DATABASE_COLUMNS = (
    "CIRCpedia", "ExoRBase", "CircNet", "CIRIonco",
    "CircASbase", "TSCD", "MiOncoCirc",
)


def require_columns(fieldnames, columns, label):
    if fieldnames is None:
        raise ValueError(f"{label} has no header")
    missing = [column for column in columns if column not in fieldnames]
    if missing:
        raise ValueError(f"{label} is missing: " + ", ".join(missing))


def parse_binary(value, column, label):
    value = value.strip()
    if value not in {"0", "1"}:
        raise ValueError(f"Invalid binary {column} in {label}: {value!r}")
    return int(value)


def load_human_support(connection, path):
    connection.execute(
        "CREATE TABLE human_support (chr TEXT, pos1 INTEGER, pos2 INTEGER, "
        "strand TEXT, seven_db INTEGER, fl_circas INTEGER, "
        "PRIMARY KEY(chr,pos1,pos2,strand)) WITHOUT ROWID"
    )
    source_rows = 0
    missing_strand_rows = 0
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ("chrom", "start", "end", "strand", *DATABASE_COLUMNS, "FL-circAS")
        require_columns(reader.fieldnames, required, "human presence table")
        batch = []
        for line_number, row in enumerate(reader, start=2):
            source_rows += 1
            if row["strand"].strip().upper() == "NA":
                missing_strand_rows += 1
                continue
            label = f"human presence line {line_number}"
            key = parse_coordinate(
                row["chrom"], row["start"], row["end"], row["strand"], label
            )
            seven_db = sum(
                parse_binary(row[column], column, label) for column in DATABASE_COLUMNS
            )
            fl_circas = parse_binary(row["FL-circAS"], "FL-circAS", label)
            batch.append((*key, seven_db, fl_circas))
            if len(batch) >= 100_000:
                connection.executemany(
                    "INSERT INTO human_support VALUES (?, ?, ?, ?, ?, ?)", batch
                )
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO human_support VALUES (?, ?, ?, ?, ?, ?)", batch
            )
    stored_rows = connection.execute("SELECT COUNT(*) FROM human_support").fetchone()[0]
    return source_rows, stored_rows, missing_strand_rows


def load_coordinate_support(connection, table, path, columns, label):
    connection.execute(
        f"CREATE TABLE {table} (chr TEXT, pos1 INTEGER, pos2 INTEGER, "
        "strand TEXT, PRIMARY KEY(chr,pos1,pos2,strand)) WITHOUT ROWID"
    )
    source_rows = 0
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, columns, label)
        batch = []
        for line_number, row in enumerate(reader, start=2):
            key = parse_coordinate(
                row[columns[0]], row[columns[1]], row[columns[2]], row[columns[3]],
                f"{label} line {line_number}",
            )
            batch.append(key)
            source_rows += 1
            if len(batch) >= 100_000:
                connection.executemany(
                    f"INSERT OR IGNORE INTO {table} VALUES (?, ?, ?, ?)", batch
                )
                batch.clear()
        if batch:
            connection.executemany(
                f"INSERT OR IGNORE INTO {table} VALUES (?, ?, ?, ?)", batch
            )
    stored_rows = connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    return source_rows, stored_rows


output_path = Path(snakemake.output.database)
output_path.parent.mkdir(parents=True, exist_ok=True)
fd, temporary_name = tempfile.mkstemp(
    prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
)
os.close(fd)
temporary_path = Path(temporary_name)
started = time.perf_counter()
timings = []
connection = None

try:
    connection = sqlite3.connect(temporary_path)
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = FILE")
    connection.execute(
        "CREATE TABLE reference_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL) "
        "WITHOUT ROWID"
    )

    phase_started = time.perf_counter()
    human_source_rows, human_rows, missing_human_strand_rows = load_human_support(
        connection, snakemake.input.human_presence
    )
    timings.append(("human_support_seconds", time.perf_counter() - phase_started))

    phase_started = time.perf_counter()
    mouse_source_rows, mouse_rows = load_coordinate_support(
        connection,
        "mouse_support",
        snakemake.input.mouse,
        ("chrom", "start", "end", "strand"),
        "lifted mouse table",
    )
    timings.append(("mouse_support_seconds", time.perf_counter() - phase_started))

    phase_started = time.perf_counter()
    circfl_source_rows, circfl_rows = load_coordinate_support(
        connection,
        "circfl_support",
        snakemake.input.circfl,
        ("chr", "pos1", "pos2", "strand"),
        "circFL table",
    )
    timings.append(("circfl_support_seconds", time.perf_counter() - phase_started))

    metadata = {
        "schema_version": REFERENCE_SCHEMA_VERSION,
        "human_source_rows": human_source_rows,
        "human_support_rows": human_rows,
        "human_missing_strand_rows": missing_human_strand_rows,
        "mouse_source_rows": mouse_source_rows,
        "mouse_support_rows": mouse_rows,
        "circfl_source_rows": circfl_source_rows,
        "circfl_support_rows": circfl_rows,
    }
    connection.executemany(
        "INSERT INTO reference_metadata VALUES (?, ?)",
        ((key, str(value)) for key, value in metadata.items()),
    )
    connection.commit()
    connection.execute("ANALYZE")
    connection.commit()
    connection.close()
    connection = None

    with temporary_path.open("rb") as handle:
        os.fsync(handle.fileno())
    os.replace(temporary_path, output_path)
except Exception:
    if connection is not None:
        connection.close()
    temporary_path.unlink(missing_ok=True)
    raise

timings.append(("total_seconds", time.perf_counter() - started))
with open(snakemake.log[0], "w", encoding="utf-8") as log:
    for key, value in metadata.items():
        log.write(f"{key}={value}\n")
    for key, value in timings:
        log.write(f"{key}={value:.3f}\n")
