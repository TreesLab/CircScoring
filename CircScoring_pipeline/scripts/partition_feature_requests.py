import csv
import faulthandler
import logging
import sys
import time
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
if group not in {"conservation", "splicing"}:
    raise ValueError(f"Unsupported query partition group: {group!r}")
result_key = "event_id"
result_columns = [result_key, *FEATURE_COLUMNS[group]]
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)
logger = logging.getLogger(f"partition_query_feature_requests.{group}")
logger.setLevel(logging.INFO)
logger.propagate = False
handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
handler.setFormatter(logging.Formatter(
    "%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
))
logger.handlers[:] = [handler]
faulthandler.enable(file=handler.stream, all_threads=True)


def bigwig_chromosomes():
    import pyBigWig

    references = []
    for path in (snakemake.input.phyloP, snakemake.input.phastCons):
        logger.info("Opening bigWig reference: %s", path)
        handle = pyBigWig.open(path)
        if handle is None or not handle.isBigWig():
            raise ValueError(f"Invalid bigWig file: {path}")
        chromosomes = handle.chroms()
        logger.info("Loaded %d chromosome entries from %s", len(chromosomes), path)
        references.append((handle, chromosomes))
    for handle, _ in references:
        handle.close()
    return [chromosomes for _, chromosomes in references]


def fasta_lengths(path):
    index_path = Path(f"{path}.fai")
    if index_path.exists():
        lengths = {}
        with index_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    raise ValueError(f"Invalid FASTA index line {line_number}: {index_path}")
                lengths[fields[0]] = int(fields[1])
        return lengths

    lengths = {}
    name = None
    length = 0
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                if not name or name in lengths:
                    raise ValueError(f"Invalid FASTA header at line {line_number}: {path}")
                length = 0
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence before header at line {line_number}: {path}")
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def is_conservation_calculable(row, references):
    chrom = row["chr"]
    pos1 = int(row["pos1"])
    pos2 = int(row["pos2"])
    windows = ((pos1 - 11, pos1 + 9), (pos2 - 10, pos2 + 10))
    return all(
        chrom in chromosomes
        and all(start >= 0 and end <= chromosomes[chrom] for start, end in windows)
        for chromosomes in references
    )


def is_splicing_calculable(row, lengths):
    chrom = row["chr"]
    if chrom not in lengths:
        return False
    donor = int(row["donor"])
    acceptor = int(row["acceptor"])
    if row["strand"] == "+":
        windows = ((donor - 3, donor + 6), (acceptor - 21, acceptor + 2))
    else:
        windows = ((donor - 7, donor + 2), (acceptor - 3, acceptor + 20))
    return all(start >= 0 and end <= lengths[chrom] for start, end in windows)


def main():
    started_at = time.perf_counter()
    logger.info("Starting %s request partition", group)
    logger.info("Request input: %s", snakemake.input.requests)
    reference = (
        bigwig_chromosomes()
        if group == "conservation"
        else fasta_lengths(snakemake.input.genome)
    )
    is_calculable = (
        is_conservation_calculable
        if group == "conservation"
        else is_splicing_calculable
    )
    valid_handle, valid_temporary, valid_output = atomic_text_writer(
        snakemake.output.calculable, newline=""
    )
    missing_handle, missing_temporary, missing_output = atomic_text_writer(
        snakemake.output.unavailable, newline=""
    )

    try:
        valid_writer = csv.DictWriter(
            valid_handle,
            fieldnames=REQUEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        missing_writer = csv.writer(
            missing_handle, delimiter="\t", lineterminator="\n"
        )
        valid_writer.writeheader()
        missing_writer.writerow(result_columns)
        valid_count = 0
        missing_count = 0
        with open(snakemake.input.requests, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            validate_header(reader, REQUEST_COLUMNS, "query feature requests")
            for row_count, row in enumerate(reader, start=1):
                if is_calculable(row, reference):
                    valid_writer.writerow(row)
                    valid_count += 1
                else:
                    missing_writer.writerow(
                        [row["circ_id"], *([""] * len(FEATURE_COLUMNS[group]))]
                    )
                    missing_count += 1
                if row_count % 100_000 == 0:
                    logger.info(
                        "Processed %d requests: calculable=%d unavailable=%d",
                        row_count,
                        valid_count,
                        missing_count,
                    )
        finish_atomic_write(valid_handle, valid_temporary, valid_output)
        finish_atomic_write(missing_handle, missing_temporary, missing_output)
    except BaseException:
        valid_handle.close()
        missing_handle.close()
        valid_temporary.unlink(missing_ok=True)
        missing_temporary.unlink(missing_ok=True)
        raise
    logger.info(
        "Partitioned %s query requests in %.2f seconds: calculable=%d unavailable=%d",
        group, time.perf_counter() - started_at, valid_count, missing_count,
    )


try:
    main()
except BaseException:
    logger.exception("Failed to partition %s query requests", group)
    handler.flush()
    raise

handler.flush()
