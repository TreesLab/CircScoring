import csv
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import CHECKAA_COLUMNS


RESULT_HEADER = CHECKAA_COLUMNS


def expected_ids(path):
    with open(path, newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(row) != 5:
                raise ValueError(
                    f"Active circRNA line {line_number} must have five columns"
                )
            yield row[4]


expected = iter(expected_ids(snakemake.input.circRNAs))
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
fd, temporary_name = tempfile.mkstemp(
    prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
)
os.close(fd)
temporary_path = Path(temporary_name)
row_count = 0

try:
    with temporary_path.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(RESULT_HEADER)
        for result_path in snakemake.input.results:
            with open(result_path, newline="", encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, None)
                if header != RESULT_HEADER:
                    raise ValueError(
                        f"Unexpected checkAA result header in {result_path}: {header!r}"
                    )
                for line_number, row in enumerate(reader, start=2):
                    if len(row) != len(RESULT_HEADER):
                        raise ValueError(
                            f"Invalid checkAA result at {result_path}:{line_number}"
                        )
                    expected_id = next(expected, None)
                    if expected_id is None or row[0] != expected_id:
                        raise ValueError(
                            f"checkAA result order mismatch: expected {expected_id!r}, "
                            f"got {row[0]!r}"
                        )
                    values = row[1:]
                    if any(value not in {"0", "1"} for value in values):
                        raise ValueError(
                            f"Non-binary checkAA result at {result_path}:{line_number}"
                        )
                    if int(values[2]) != int(values[0] == "1" or values[1] == "1"):
                        raise ValueError(
                            f"Invalid ambiguity OR at {result_path}:{line_number}"
                        )
                    writer.writerow(row)
                    row_count += 1
        extra_expected = next(expected, None)
        if extra_expected is not None:
            raise ValueError(
                f"Missing checkAA results beginning with {extra_expected}"
            )
    os.replace(temporary_path, output_path)
except Exception:
    temporary_path.unlink(missing_ok=True)
    raise
