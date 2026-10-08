import csv
import os
import sqlite3
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import circ_id_for_key, parse_coordinate


HEADER_ALIASES = (
    {"chr", "chrom", "chromosome", "chrm"},
    {"start", "pos1"},
    {"end", "pos2"},
    {"strand"},
)


def is_header(row):
    if len(row) != len(HEADER_ALIASES):
        return False
    normalized = [value.strip().lower() for value in row]
    return all(value in aliases for value, aliases in zip(normalized, HEADER_ALIASES))


output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
output_handle = tempfile.NamedTemporaryFile(
    mode="w", newline="", encoding="utf-8", dir=output_path.parent,
    prefix=f".{output_path.name}.", suffix=".tmp", delete=False,
)
temporary_output = Path(output_handle.name)
fd, database_name = tempfile.mkstemp(
    prefix=".input-circrnas.", suffix=".sqlite", dir=output_path.parent
)
os.close(fd)
database_path = Path(database_name)

try:
    connection = sqlite3.connect(database_path)
    connection.execute(
        "CREATE TABLE seen (chr TEXT, pos1 INTEGER, pos2 INTEGER, strand TEXT, "
        "PRIMARY KEY (chr, pos1, pos2, strand)) WITHOUT ROWID"
    )
    writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
    row_count = 0
    with open(snakemake.input[0], newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 4:
                raise ValueError(
                    f"Input circRNA line {line_number} must have four columns, "
                    f"got {len(row)}"
                )
            if line_number == 1 and is_header(row):
                continue
            key = parse_coordinate(*row, f"input circRNA line {line_number}")
            try:
                connection.execute("INSERT INTO seen VALUES (?, ?, ?, ?)", key)
            except sqlite3.IntegrityError:
                raise ValueError(
                    f"Duplicated input circRNA coordinate at line "
                    f"{line_number}: {key}"
                ) from None
            writer.writerow([*key, circ_id_for_key(key)])
            row_count += 1
    if row_count == 0:
        raise ValueError("Input circRNA table contains no records")
    connection.close()
    output_handle.flush()
    os.fsync(output_handle.fileno())
    output_handle.close()
    os.replace(temporary_output, output_path)
finally:
    if not output_handle.closed:
        output_handle.close()
    temporary_output.unlink(missing_ok=True)
    database_path.unlink(missing_ok=True)
