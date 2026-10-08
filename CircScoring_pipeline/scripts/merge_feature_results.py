import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import (
    FEATURE_COLUMNS,
    REQUEST_COLUMNS,
    atomic_text_writer,
    finish_atomic_write,
    validate_header,
)


group = str(snakemake.params.group)
columns = FEATURE_COLUMNS[group]
result_key = "event_id"
expected_header = [result_key, *columns]
results = {}

for label, path in (
    ("calculated", snakemake.input.calculated),
    ("unavailable", snakemake.input.unavailable),
):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_header(reader, expected_header, f"query {group} {label} results")
        for row in reader:
            event_id = row[result_key]
            if event_id in results:
                raise ValueError(f"Duplicated {group} result: {event_id}")
            results[event_id] = [row[column] for column in columns]

handle, temporary_path, output_path = atomic_text_writer(snakemake.output[0], newline="")
try:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(expected_header)
    request_count = 0
    with open(snakemake.input.requests, newline="", encoding="utf-8") as requests:
        reader = csv.DictReader(requests, delimiter="\t")
        validate_header(reader, REQUEST_COLUMNS, "query feature requests")
        for row in reader:
            event_id = row["circ_id"]
            if event_id not in results:
                raise ValueError(f"Missing {group} result: {event_id}")
            writer.writerow([event_id, *results.pop(event_id)])
            request_count += 1
    if results:
        raise ValueError(f"Unexpected {group} result: {next(iter(results))}")
    finish_atomic_write(handle, temporary_path, output_path)
except Exception:
    handle.close()
    temporary_path.unlink(missing_ok=True)
    raise
