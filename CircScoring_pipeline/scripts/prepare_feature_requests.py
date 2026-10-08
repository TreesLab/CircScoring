import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import (
    RAW_COLUMNS,
    REQUEST_COLUMNS,
    atomic_text_writer,
    event_id_for_key,
    finish_atomic_write,
    get_donor_acceptor,
    parse_coordinate,
)


handle, temporary_path, output_path = atomic_text_writer(
    snakemake.output[0], newline=""
)
try:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(REQUEST_COLUMNS)
    with open(snakemake.input.circRNAs, newline="", encoding="utf-8") as source:
        for line_number, row in enumerate(csv.reader(source, delimiter="\t"), start=1):
            if len(row) != len(RAW_COLUMNS):
                raise ValueError(
                    f"Expected five input columns at line {line_number}, got {len(row)}"
                )
            key = parse_coordinate(*row[:4], f"circRNA input line {line_number}")
            donor, acceptor = get_donor_acceptor(key[1], key[2], key[3])
            writer.writerow([*key, event_id_for_key(key), donor, acceptor])
    finish_atomic_write(handle, temporary_path, output_path)
except Exception:
    if not handle.closed:
        handle.close()
    temporary_path.unlink(missing_ok=True)
    raise
