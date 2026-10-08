import csv
import json
import math
import sqlite3
import sys
from itertools import zip_longest
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import (
    DETAILED_FEATURE_COLUMNS,
    FEATURE_COLUMNS,
    PREDICTCS_FEATURE_COLUMNS,
    RAW_COLUMNS,
    REQUEST_COLUMNS,
    atomic_text_writer,
    finish_atomic_write,
    validate_header,
)


REFERENCE_SCHEMA_VERSION = 1
ANNOTATION_FEATURES = (
    "donor_site_at_the_annotated_boundary",
    "acceptor_site_at_the_annotated_boundary",
    "donor_acceptor_sites_at_the_same_transcript_isoform",
)
AS_FEATURES = ("has_AS_event(donor)", "has_AS_event(acceptor)")
MISSING_VALUES = {"", ".", "na", "n/a", "nan"}


def parse_binary(value, column, circ_id):
    value = value.strip()
    if value not in {"0", "1"}:
        raise ValueError(f"Invalid binary {column} for {circ_id}: {value!r}")
    return int(value)


def normalize_numeric(value, column, circ_id):
    text = str(value).strip()
    if text.lower() in MISSING_VALUES:
        return ""
    try:
        number = float(text)
    except ValueError:
        raise ValueError(
            f"Nonnumeric {column} for {circ_id}: {text!r}"
        ) from None
    return text if math.isfinite(number) else ""


def load_active(connection, active_path, request_path):
    connection.execute(
        "CREATE TEMP TABLE active (row_num INTEGER PRIMARY KEY, chr TEXT, "
        "pos1 INTEGER, pos2 INTEGER, strand TEXT, circ_id TEXT, "
        "event_id TEXT UNIQUE)"
    )
    with open(active_path, newline="", encoding="utf-8") as active_handle, open(
        request_path, newline="", encoding="utf-8"
    ) as request_handle:
        active_rows = csv.reader(active_handle, delimiter="\t")
        requests = csv.DictReader(request_handle, delimiter="\t")
        validate_header(requests, REQUEST_COLUMNS, "feature requests")
        batch = []
        for row_num, pair in enumerate(
            zip_longest(active_rows, requests), start=1
        ):
            active, request = pair
            if active is None or request is None:
                raise ValueError("Input and feature request row counts differ")
            if len(active) != len(RAW_COLUMNS):
                raise ValueError(f"Invalid active circRNA row {row_num}")
            if active[:4] != [request[column] for column in RAW_COLUMNS[:4]]:
                raise ValueError(
                    f"Input/request coordinate mismatch at row {row_num}"
                )
            batch.append((row_num - 1, *active, request["circ_id"]))
            if len(batch) >= 100_000:
                connection.executemany(
                    "INSERT INTO active VALUES (?, ?, ?, ?, ?, ?, ?)", batch
                )
                batch.clear()
        if batch:
            connection.executemany(
                "INSERT INTO active VALUES (?, ?, ?, ?, ?, ?, ?)", batch
            )


def load_features(connection, group, path, key_column):
    columns = FEATURE_COLUMNS[group]
    table = f"feature_{group}"
    connection.execute(
        f"CREATE TEMP TABLE {table} (event_id TEXT PRIMARY KEY, payload TEXT)"
    )
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Feature table has no header: {path}")
        missing = [
            column for column in (key_column, *columns)
            if column not in reader.fieldnames
        ]
        if missing:
            raise ValueError(f"Missing {group} columns: {', '.join(missing)}")
        connection.executemany(
            f"INSERT INTO {table} VALUES (?, ?)",
            (
                (
                    row[key_column],
                    json.dumps(
                        {column: row[column] for column in columns},
                        ensure_ascii=True,
                        separators=(",", ":"),
                    ),
                )
                for row in reader
            ),
        )


reference_path = Path(snakemake.input.reference).resolve()
connection = sqlite3.connect(
    f"{reference_path.as_uri()}?mode=ro&immutable=1", uri=True
)
connection.execute("PRAGMA temp_store = FILE")
version = connection.execute(
    "SELECT value FROM reference_metadata WHERE key='schema_version'"
).fetchone()
if version is None or version[0] != str(REFERENCE_SCHEMA_VERSION):
    connection.close()
    raise ValueError(f"Incompatible predictCS reference database: {version}")

try:
    load_active(connection, snakemake.input.circRNAs, snakemake.input.requests)
    load_features(connection, "annotation", snakemake.input.annotation, "circ_id")
    load_features(
        connection,
        "alternative_splicing",
        snakemake.input.alternative_splicing,
        "circ_id",
    )
    load_features(connection, "conservation", snakemake.input.conservation, "event_id")
    load_features(connection, "splicing", snakemake.input.splicing, "event_id")

    query = """
        SELECT q.circ_id,
               COALESCE(h.seven_db, 0), COALESCE(h.fl_circas, 0),
               CASE WHEN m.chr IS NULL THEN 0 ELSE 1 END,
               CASE WHEN f.chr IS NULL THEN 0 ELSE 1 END,
               a.payload, s.payload, c.payload, x.payload
        FROM active AS q
        LEFT JOIN human_support AS h USING (chr,pos1,pos2,strand)
        LEFT JOIN mouse_support AS m USING (chr,pos1,pos2,strand)
        LEFT JOIN circfl_support AS f USING (chr,pos1,pos2,strand)
        LEFT JOIN feature_annotation AS a USING (event_id)
        LEFT JOIN feature_alternative_splicing AS s USING (event_id)
        LEFT JOIN feature_conservation AS c USING (event_id)
        LEFT JOIN feature_splicing AS x USING (event_id)
        ORDER BY q.row_num
    """

    feature_handle, feature_tmp, feature_output = atomic_text_writer(
        snakemake.output.features, newline=""
    )
    detailed_handle, detailed_tmp, detailed_output = atomic_text_writer(
        snakemake.output.detailed, newline=""
    )
    try:
        feature_writer = csv.writer(
            feature_handle, delimiter="\t", lineterminator="\n"
        )
        detailed_writer = csv.writer(
            detailed_handle, delimiter="\t", lineterminator="\n"
        )
        feature_writer.writerow(PREDICTCS_FEATURE_COLUMNS)
        detailed_writer.writerow(DETAILED_FEATURE_COLUMNS)
        for row in connection.execute(query):
            circ_id, seven_db, fl_circas, mouse, circfl, *payloads = row
            if any(payload is None for payload in payloads):
                groups = (
                    "annotation",
                    "alternative_splicing",
                    "conservation",
                    "splicing",
                )
                missing = [
                    group for group, payload in zip(groups, payloads)
                    if payload is None
                ]
                raise ValueError(
                    f"Missing feature result for {circ_id}: {', '.join(missing)}"
                )
            annotation, alternative_splicing, conservation, splicing = map(
                json.loads, payloads
            )
            feature_writer.writerow([
                circ_id,
                *(parse_binary(annotation[c], c, circ_id) for c in ANNOTATION_FEATURES),
                *(parse_binary(alternative_splicing[c], c, circ_id) for c in AS_FEATURES),
                *(
                    normalize_numeric(conservation[c], c, circ_id)
                    for c in FEATURE_COLUMNS["conservation"]
                ),
                *(
                    normalize_numeric(splicing[c], c, circ_id)
                    for c in FEATURE_COLUMNS["splicing"]
                ),
                seven_db,
                int(bool(fl_circas or circfl)),
                mouse,
            ])
            detailed_writer.writerow([
                circ_id,
                *(annotation[c] for c in FEATURE_COLUMNS["annotation"]),
                *(
                    alternative_splicing[c]
                    for c in FEATURE_COLUMNS["alternative_splicing"]
                ),
                *(
                    normalize_numeric(conservation[c], c, circ_id)
                    for c in FEATURE_COLUMNS["conservation"]
                ),
                *(
                    normalize_numeric(splicing[c], c, circ_id)
                    for c in FEATURE_COLUMNS["splicing"]
                ),
                seven_db,
                int(bool(fl_circas or circfl)),
                mouse,
            ])
        finish_atomic_write(feature_handle, feature_tmp, feature_output)
        finish_atomic_write(detailed_handle, detailed_tmp, detailed_output)
    except Exception:
        if not feature_handle.closed:
            feature_handle.close()
        if not detailed_handle.closed:
            detailed_handle.close()
        feature_tmp.unlink(missing_ok=True)
        detailed_tmp.unlink(missing_ok=True)
        raise
finally:
    connection.close()
