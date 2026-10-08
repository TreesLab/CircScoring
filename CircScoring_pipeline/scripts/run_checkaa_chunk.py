import csv
import logging
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path


release_root = Path(snakemake.scriptdir).parent
checkaa_root = release_root / "scripts" / "checkaa"
sys.path.insert(0, str(checkaa_root))

from checkAA_reads import (
    RESULT_HEADER,
    AmbAlnChecker,
    read_fasta_records,
    write_fasta_records,
)


def read_chunk(path):
    mappings = OrderedDict()
    with open(path, newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if len(row) != 4:
                raise ValueError(
                    f"checkAA chunk line {line_number} must have four columns"
                )
            chrom, pos1_text, pos2_text, strand = row
            try:
                pos1 = int(pos1_text)
                pos2 = int(pos2_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid checkAA coordinates at line {line_number}: {row!r}"
                ) from exc
            if pos1 < 1 or pos1 > pos2 or strand not in {"+", "-"}:
                raise ValueError(
                    f"Invalid checkAA coordinate at line {line_number}: {row!r}"
                )
            standard_id = f"{chrom}:{pos1}|{pos2}({strand})"
            internal_id = (
                f"{chrom}:{pos2}|{pos1}(+)" if strand == "+" else standard_id
            )
            if internal_id in mappings:
                raise ValueError(f"Duplicate checkAA ID in chunk: {internal_id}")
            mappings[internal_id] = standard_id
    if not mappings:
        raise ValueError(f"Empty checkAA chunk: {path}")
    return mappings


def write_result(path, records, mappings, computed):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent
    )
    os.close(fd)
    temporary_path = Path(temporary_name)
    try:
        with temporary_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(RESULT_HEADER)
            for record in records:
                values = computed.get(record.id)
                if values is None:
                    raise ValueError(f"Missing checkAA result for {record.id}")
                writer.writerow([mappings[record.id], *values])
        os.replace(temporary_path, output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


output_path = Path(snakemake.output.result)
output_path.parent.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=log_path,
    filemode="w",
    format="{asctime} - {message}",
    level=logging.INFO,
    style="{",
    force=True,
)

mappings = read_chunk(snakemake.input.table)
flank_length = int(snakemake.params.flank_length)
if flank_length < 1:
    raise ValueError("checkAA flank_length must be positive")

with tempfile.TemporaryDirectory(
    prefix="checkAA.", dir=output_path.parent
) as temporary_directory:
    temporary_directory = Path(temporary_directory)
    reads_path = temporary_directory / "pseudo_bsj.fa"
    with log_path.open("a", encoding="utf-8") as log_handle:
        subprocess.run(
            [
                sys.executable,
                str(checkaa_root / "get_flanking_seq.py"),
                str(snakemake.input.genome),
                str(snakemake.input.table),
                str(reads_path),
                str(flank_length),
            ],
            stdout=log_handle,
            stderr=log_handle,
            check=True,
        )

    records = read_fasta_records(reads_path)
    record_ids = [record.id for record in records]
    if record_ids != list(mappings):
        expected = set(mappings)
        observed = set(record_ids)
        raise ValueError(
            "Pseudo-BSJ FASTA IDs do not match the input chunk: "
            f"missing={sorted(expected - observed)[:5]}, "
            f"extra={sorted(observed - expected)[:5]}"
        )
    expected_length = 2 * flank_length
    invalid_lengths = [
        (record.id, len(record.sequence))
        for record in records
        if len(record.sequence) != expected_length
    ]
    if invalid_lengths:
        raise ValueError(
            f"Pseudo-BSJ sequences must be {expected_length} bp: "
            f"{invalid_lengths[:5]}"
        )

    input_fasta = temporary_directory / "checkaa_input.fa"
    write_fasta_records(records, input_fasta)
    checker = AmbAlnChecker(
        ref_genome=str(snakemake.input.genome),
        ref_others=str(snakemake.input.other_references),
        work_dir=str(temporary_directory),
        num_proc=int(snakemake.threads),
        aligner="pblat",
        pblat_bin=str(snakemake.params.pblat),
    )
    computed = checker.check(input_fasta).result
    write_result(output_path, records, mappings, computed)
