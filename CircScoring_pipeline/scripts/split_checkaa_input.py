import csv
import json
import os
import shutil
import tempfile
from pathlib import Path


chunk_size = int(snakemake.params.chunk_size)
if chunk_size < 1:
    raise ValueError("checkAA chunk_size must be positive")

output_directory = Path(snakemake.output.chunks)
output_directory.parent.mkdir(parents=True, exist_ok=True)
temporary_directory = Path(tempfile.mkdtemp(
    prefix=f".{output_directory.name}.", dir=output_directory.parent
))

chunk_handle = None
chunk_writer = None
chunk_count = 0
row_count = 0
chunk_rows = 0

try:
    with open(snakemake.input.circRNAs, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 5:
                raise ValueError(
                    f"Active circRNA line {line_number} must have five columns"
                )
            if chunk_handle is None or chunk_rows == chunk_size:
                if chunk_handle is not None:
                    chunk_handle.close()
                chunk_count += 1
                chunk_rows = 0
                chunk_handle = open(
                    temporary_directory / f"part_{chunk_count:06d}.tsv",
                    "w", newline="", encoding="utf-8",
                )
                chunk_writer = csv.writer(
                    chunk_handle, delimiter="\t", lineterminator="\n"
                )
            chunk_writer.writerow(row[:4])
            chunk_rows += 1
            row_count += 1

    if chunk_handle is not None:
        chunk_handle.close()
        chunk_handle = None
    if row_count == 0:
        raise ValueError("Active circRNA table contains no records")

    manifest = {
        "chunk_size": chunk_size,
        "chunks": chunk_count,
        "rows": row_count,
    }
    (temporary_directory / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    if output_directory.exists():
        shutil.rmtree(output_directory)
    os.replace(temporary_directory, output_directory)
except Exception:
    if chunk_handle is not None:
        chunk_handle.close()
    shutil.rmtree(temporary_directory, ignore_errors=True)
    raise
