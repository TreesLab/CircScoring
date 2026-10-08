import atexit
import csv
import logging
import multiprocessing
import os
import shutil
import statistics
import tempfile
import time
from functools import lru_cache

import pyBigWig


IO_BUFFER_SIZE = 8 * 1024 * 1024
WRITE_BATCH_SIZE = 10_000
REQUIRED_COLUMNS = ("chr", "pos1", "pos2", "strand", "circ_id")
OUTPUT_COLUMNS = (
    "event_id",
    "phyloP(acceptor_in)",
    "phyloP(acceptor_out)",
    "phyloP(donor_in)",
    "phyloP(donor_out)",
    "phastCons(acceptor_in)",
    "phastCons(acceptor_out)",
    "phastCons(donor_in)",
    "phastCons(donor_out)",
)
POS1_WINDOW = "pos1"
POS2_WINDOW = "pos2"


_WORKER_PHYLOP = None
_WORKER_PHASTCONS = None
_WORKER_GET_ENDPOINT_SCORES = None
_WORKER_COLUMN_INDEXES = None


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def to_ucsc_chromosome_name(chromosome):
    chromosome = chromosome.strip()
    if chromosome == "MT":
        return "chrM"
    if chromosome.startswith("chr"):
        return chromosome
    return f"chr{chromosome}"


def close_worker_bigwigs():
    global _WORKER_PHYLOP
    global _WORKER_PHASTCONS

    if _WORKER_PHYLOP is not None:
        _WORKER_PHYLOP.close()
        _WORKER_PHYLOP = None
    if _WORKER_PHASTCONS is not None:
        _WORKER_PHASTCONS.close()
        _WORKER_PHASTCONS = None


def validate_bigwig(path, label):
    bigwig = None
    try:
        bigwig = pyBigWig.open(path)
        if bigwig is None or not bigwig.isBigWig():
            raise ValueError(f"{label} is not a valid bigWig file: {path}")
    except RuntimeError as exc:
        raise ValueError(f"Cannot open {label} bigWig file {path}: {exc}") \
            from exc
    finally:
        if bigwig is not None:
            bigwig.close()


def initialize_worker(
        phylop_path,
        phastcons_path,
        cache_size,
        column_indexes):
    global _WORKER_PHYLOP
    global _WORKER_PHASTCONS
    global _WORKER_GET_ENDPOINT_SCORES
    global _WORKER_COLUMN_INDEXES

    _WORKER_PHYLOP = pyBigWig.open(phylop_path)
    _WORKER_PHASTCONS = pyBigWig.open(phastcons_path)
    _WORKER_GET_ENDPOINT_SCORES = lru_cache(maxsize=cache_size)(
        get_endpoint_scores_uncached
    )
    _WORKER_COLUMN_INDEXES = column_indexes
    atexit.register(close_worker_bigwigs)


def get_endpoint_scores_uncached(chr_, position, window_type):
    if window_type == POS1_WINDOW:
        start = position - 11
        end = position + 9
    elif window_type == POS2_WINDOW:
        start = position - 10
        end = position + 10
    else:
        raise ValueError(f"Unsupported endpoint window type: {window_type!r}")

    if start < 0:
        raise ValueError(
            f"Conservation window starts before chromosome position 1: "
            f"({chr_}, {position}, {window_type})"
        )

    phylop_values = _WORKER_PHYLOP.values(chr_, start, end)
    phastcons_values = _WORKER_PHASTCONS.values(chr_, start, end)
    if len(phylop_values) != 20 or len(phastcons_values) != 20:
        raise ValueError(
            f"Expected 20 conservation values for "
            f"({chr_}, {position}, {window_type})"
        )

    return (
        statistics.fmean(phylop_values[:10]),
        statistics.fmean(phylop_values[10:]),
        statistics.fmean(phastcons_values[:10]),
        statistics.fmean(phastcons_values[10:]),
    )


def get_column_indexes(header):
    missing_columns = [
        column for column in REQUIRED_COLUMNS if column not in header
    ]
    if missing_columns:
        raise ValueError(
            "Missing required input column(s): "
            + ", ".join(missing_columns)
        )
    return {column: header.index(column) for column in REQUIRED_COLUMNS}


def get_required_values(row, row_label):
    missing_columns = [
        column
        for column, index in _WORKER_COLUMN_INDEXES.items()
        if index >= len(row) or row[index] == ""
    ]
    if missing_columns:
        raise ValueError(
            f"Missing required column(s) at {row_label}: "
            + ", ".join(missing_columns)
        )
    return {
        column: row[index]
        for column, index in _WORKER_COLUMN_INDEXES.items()
    }


def parse_position(value, column, row_label):
    try:
        position = int(value)
    except ValueError:
        raise ValueError(
            f"Invalid {column} coordinate at {row_label}: {value!r}"
        ) from None
    if position < 1:
        raise ValueError(
            f"Invalid {column} coordinate at {row_label}: {position}"
        )
    return position


def calculate_circrna_scores(values, row_label):
    chr_ = to_ucsc_chromosome_name(values["chr"])
    pos1 = parse_position(values["pos1"], "pos1", row_label)
    pos2 = parse_position(values["pos2"], "pos2", row_label)
    pos1, pos2 = sorted((pos1, pos2))
    strand = values["strand"]
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid strand at {row_label}: {strand!r}")

    try:
        pos1_scores = _WORKER_GET_ENDPOINT_SCORES(
            chr_,
            pos1,
            POS1_WINDOW,
        )
        pos2_scores = _WORKER_GET_ENDPOINT_SCORES(
            chr_,
            pos2,
            POS2_WINDOW,
        )
    except (RuntimeError, ValueError) as exc:
        raise ValueError(
            f"Cannot read conservation scores at {row_label}: {exc}"
        ) from exc

    phylop_pos1_left, phylop_pos1_right, \
        phastcons_pos1_left, phastcons_pos1_right = pos1_scores
    phylop_pos2_left, phylop_pos2_right, \
        phastcons_pos2_left, phastcons_pos2_right = pos2_scores

    if strand == "+":
        phylop_scores = (
            phylop_pos1_right,
            phylop_pos1_left,
            phylop_pos2_left,
            phylop_pos2_right,
        )
        phastcons_scores = (
            phastcons_pos1_right,
            phastcons_pos1_left,
            phastcons_pos2_left,
            phastcons_pos2_right,
        )
    else:
        phylop_scores = (
            phylop_pos2_left,
            phylop_pos2_right,
            phylop_pos1_right,
            phylop_pos1_left,
        )
        phastcons_scores = (
            phastcons_pos2_left,
            phastcons_pos2_right,
            phastcons_pos1_right,
            phastcons_pos1_left,
        )

    return (values["circ_id"], *phylop_scores, *phastcons_scores)


def read_header_and_data_start(path):
    with open(path, "rb") as f_in:
        header_line = f_in.readline()
        if not header_line:
            raise ValueError("Input circRNA file is empty")
        data_start = f_in.tell()

    try:
        header_text = header_line.decode("utf-8")
    except UnicodeDecodeError:
        raise ValueError("Input circRNA header is not valid UTF-8") from None

    header = next(csv.reader([header_text], delimiter="\t"))
    return header, data_start


def plan_byte_range_chunks(path, data_start, requested_chunks):
    file_size = os.path.getsize(path)
    if data_start >= file_size:
        return []

    boundaries = [data_start]
    data_size = file_size - data_start
    with open(path, "rb") as f_in:
        for chunk_index in range(1, requested_chunks):
            target = data_start + data_size * chunk_index // requested_chunks
            f_in.seek(target)
            f_in.readline()
            boundary = f_in.tell()
            if boundary < file_size and boundary > boundaries[-1]:
                boundaries.append(boundary)
    boundaries.append(file_size)

    return [
        (chunk_id, start, end)
        for chunk_id, (start, end) in enumerate(
            zip(boundaries, boundaries[1:])
        )
        if start < end
    ]


def process_chunk(task):
    chunk_id, start, end, input_path, part_path = task
    started = time.perf_counter()
    processed_rows = 0
    output_batch = []
    cache_before = _WORKER_GET_ENDPOINT_SCORES.cache_info()

    with open(input_path, "rb") as f_in, open(
            part_path,
            "w",
            newline="",
            encoding="utf-8",
            buffering=IO_BUFFER_SIZE) as f_out:
        f_in.seek(start)
        writer = csv.writer(
            f_out,
            delimiter="\t",
            lineterminator="\n",
        )
        record_num = 0

        while f_in.tell() < end:
            byte_offset = f_in.tell()
            raw_line = f_in.readline()
            if not raw_line:
                break
            record_num += 1
            try:
                line = raw_line.decode("utf-8").rstrip("\r\n")
            except UnicodeDecodeError:
                raise ValueError(
                    f"Input row at byte offset {byte_offset} is not valid UTF-8"
                ) from None
            if not line:
                continue

            row_label = (
                f"chunk {chunk_id}, record {record_num}, "
                f"byte offset {byte_offset}"
            )
            row = line.split("\t")
            values = get_required_values(row, row_label)
            output_batch.append(calculate_circrna_scores(values, row_label))
            processed_rows += 1

            if len(output_batch) >= WRITE_BATCH_SIZE:
                writer.writerows(output_batch)
                output_batch.clear()

        if output_batch:
            writer.writerows(output_batch)

    cache_after = _WORKER_GET_ENDPOINT_SCORES.cache_info()
    return {
        "chunk_id": chunk_id,
        "part_path": part_path,
        "processed_rows": processed_rows,
        "elapsed": time.perf_counter() - started,
        "pid": os.getpid(),
        "cache_hits": cache_after.hits - cache_before.hits,
        "cache_misses": cache_after.misses - cache_before.misses,
    }


def concatenate_parts(output_path, part_paths):
    with open(output_path, "wb", buffering=IO_BUFFER_SIZE) as f_out:
        f_out.write(("\t".join(OUTPUT_COLUMNS) + "\n").encode("utf-8"))
        for part_path in part_paths:
            with open(part_path, "rb", buffering=IO_BUFFER_SIZE) as f_part:
                shutil.copyfileobj(
                    f_part,
                    f_out,
                    length=IO_BUFFER_SIZE,
                )


def main():
    total_started = time.perf_counter()
    part_directory = None
    logging.info("Starting conservation-score calculation")

    try:
        header, data_start = read_header_and_data_start(
            snakemake.input.circRNAs
        )
        column_indexes = get_column_indexes(header)
        worker_count = max(1, int(snakemake.threads))
        chunks_per_worker = int(snakemake.params.chunks_per_worker)
        cache_size = int(snakemake.params.cache_size)
        if chunks_per_worker < 1:
            raise ValueError("chunks_per_worker must be at least 1")
        if cache_size < 0:
            raise ValueError("conservation_cache_size cannot be negative")
        validate_bigwig(snakemake.input.phyloP, "phyloP")
        validate_bigwig(snakemake.input.phastCons, "phastCons")
        chunks = plan_byte_range_chunks(
            snakemake.input.circRNAs,
            data_start,
            worker_count * chunks_per_worker,
        )
        logging.info(
            "Planned %d byte-range chunks for %d workers with cache "
            "size %d per worker",
            len(chunks),
            worker_count,
            cache_size,
        )

        output_directory = os.path.dirname(snakemake.output[0]) or "."
        os.makedirs(output_directory, exist_ok=True)
        part_directory = tempfile.mkdtemp(
            prefix=".conservation.parts.",
            dir=output_directory,
        )
        tasks = [
            (
                chunk_id,
                start,
                end,
                snakemake.input.circRNAs,
                os.path.join(part_directory, f"part.{chunk_id:05d}.tsv"),
            )
            for chunk_id, start, end in chunks
        ]

        processing_started = time.perf_counter()
        results = []
        if tasks:
            context = multiprocessing.get_context("fork")
            with context.Pool(
                    processes=worker_count,
                    initializer=initialize_worker,
                    initargs=(
                        snakemake.input.phyloP,
                        snakemake.input.phastCons,
                        cache_size,
                        column_indexes,
                    )) as pool:
                for result in pool.imap_unordered(
                        process_chunk,
                        tasks,
                        chunksize=1):
                    results.append(result)
                    logging.info(
                        "Worker %d finished chunk %d: %d circRNAs in "
                        "%.2f seconds (%.2f rows/second), cache hits=%d "
                        "misses=%d",
                        result["pid"],
                        result["chunk_id"],
                        result["processed_rows"],
                        result["elapsed"],
                        (
                            result["processed_rows"] / result["elapsed"]
                            if result["elapsed"]
                            else 0
                        ),
                        result["cache_hits"],
                        result["cache_misses"],
                    )

        results.sort(key=lambda result: result["chunk_id"])
        staged_output = os.path.join(
            part_directory,
            "circRNAs.conservation_scores.tsv",
        )
        concatenate_started = time.perf_counter()
        concatenate_parts(
            staged_output,
            [result["part_path"] for result in results],
        )
        os.replace(staged_output, snakemake.output[0])

        processed_rows = sum(result["processed_rows"] for result in results)
        cache_hits = sum(result["cache_hits"] for result in results)
        cache_misses = sum(result["cache_misses"] for result in results)
        processing_elapsed = time.perf_counter() - processing_started
        logging.info(
            "Concatenated %d parts in %.2f seconds",
            len(results),
            time.perf_counter() - concatenate_started,
        )
        logging.info(
            "Endpoint cache totals: %d hits, %d misses (%.2f%% hit "
            "rate); bigWig calls=%d",
            cache_hits,
            cache_misses,
            (
                100 * cache_hits / (cache_hits + cache_misses)
                if cache_hits + cache_misses
                else 0
            ),
            cache_misses * 2,
        )
        logging.info(
            "Finished processing %d circRNAs in %.2f seconds "
            "(%.2f rows/second)",
            processed_rows,
            processing_elapsed,
            (
                processed_rows / processing_elapsed
                if processing_elapsed
                else 0
            ),
        )
        logging.info(
            "Conservation-score calculation completed in %.2f seconds",
            time.perf_counter() - total_started,
        )
    except Exception:
        logging.exception(
            "Conservation-score calculation failed after %.2f seconds",
            time.perf_counter() - total_started,
        )
        raise
    finally:
        if part_directory is not None:
            shutil.rmtree(part_directory, ignore_errors=True)


main()
