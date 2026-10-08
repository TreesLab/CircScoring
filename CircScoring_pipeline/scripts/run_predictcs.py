import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import PREDICTION_COLUMNS


input_path = Path(snakemake.input[0])
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

with input_path.open(newline="", encoding="utf-8") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader, None)
    has_records = next(reader, None) is not None

if header is None or not header or header[0] != "circRNA_id":
    raise ValueError(f"Invalid predictCS input header: {input_path}")

if not has_records:
    with output_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["circRNA_id", *PREDICTION_COLUMNS])
    Path(snakemake.log[0]).write_text(
        "No predictCS input rows; predictor was not executed.\n",
        encoding="utf-8",
    )
else:
    with open(snakemake.log[0], "w", encoding="utf-8") as log:
        subprocess.run(
            [
                sys.executable,
                str(snakemake.params.predictor),
                str(input_path),
                str(output_path),
                "--chunk-size",
                str(snakemake.params.chunk_size),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
