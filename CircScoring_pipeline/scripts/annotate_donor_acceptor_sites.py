import csv
import gc
import logging
import multiprocessing
import os
import shutil
import sqlite3
import tempfile
import time
from contextlib import ExitStack
from functools import lru_cache
from pathlib import Path

from circmimi.annotation import Annotation
from circmimi.models import (
    AcceptorSite,
    Chromosome,
    DonorSite,
    Exon,
    Gene,
    Strand,
    Transcript,
    TranscriptExon,
)


logging.basicConfig(
    filename=snakemake.log[0],
    encoding='utf-8',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)


IO_BUFFER_SIZE = 8 * 1024 * 1024
WRITE_BATCH_SIZE = 10_000
EMPTY_TRANSCRIPT_IDS = frozenset()
EMPTY_BOUNDARY_RECORD = (EMPTY_TRANSCRIPT_IDS, '', '')
EMIT_DETAILS = bool(getattr(snakemake.params, 'emit_details', True))


SUMMARY_COLUMNS = [
    'chr',
    'pos1',
    'pos2',
    'strand',
    'circ_id',
    'donor_site_at_the_annotated_boundary',
    'acceptor_site_at_the_annotated_boundary',
    'both_donor_acceptor_at_annotated_boundary',
    'donor_acceptor_sites_at_the_same_transcript_isoform',
    'the_number_of_transcripts_of_donor_acceptor_sites_at_the_same_transcript_isoform',
]


DETAIL_COLUMNS = [
    'chr',
    'pos1',
    'pos2',
    'strand',
    'circ_id',
    'donor_site',
    'acceptor_site',
    'donor_boundary_gene_ids',
    'donor_boundary_transcript_ids',
    'acceptor_boundary_gene_ids',
    'acceptor_boundary_transcript_ids',
    'shared_boundary_gene_ids',
    'shared_boundary_gene_symbols',
    'shared_boundary_transcript_ids',
]


REQUIRED_COLUMNS = (
    'chr',
    'pos1',
    'pos2',
    'strand',
    'circ_id',
    'donor',
    'acceptor',
)


_WORKER_DONOR_BOUNDARY_INDEX = None
_WORKER_ACCEPTOR_BOUNDARY_INDEX = None
_WORKER_TRANSCRIPT_GENE_INDEX = None
_WORKER_COLUMN_INDEXES = None
_WORKER_FORMAT_SHARED_RECORD = None


def build_boundary_transcript_index(
        annotation_db,
        SiteType,
        site_id_column,
        transcript_gene_index,
        emit_details):
    boundary_index = {}

    query_started = time.perf_counter()
    site_rows = annotation_db.session.query(
        Chromosome.name,
        Strand.name,
        SiteType.junc_site,
    ).select_from(
        SiteType
    ).join(
        Chromosome, SiteType.chr_id == Chromosome.id
    ).join(
        Strand, SiteType.strand_id == Strand.id
    ).all()
    logging.info(
        "Loaded %d %s site rows in %.2f seconds",
        len(site_rows),
        SiteType.__name__,
        time.perf_counter() - query_started,
    )

    for chr_, strand, site in site_rows:
        boundary_index[(chr_, strand, site)] = set()
    del site_rows

    query_started = time.perf_counter()
    transcript_rows = annotation_db.session.query(
        Chromosome.name,
        Strand.name,
        SiteType.junc_site,
        Transcript.transcript_id,
    ).select_from(
        SiteType
    ).join(
        Exon, site_id_column == SiteType.id
    ).join(
        TranscriptExon, TranscriptExon.eid == Exon.id
    ).join(
        Transcript, Transcript.id == TranscriptExon.tid
    ).join(
        Chromosome, SiteType.chr_id == Chromosome.id
    ).join(
        Strand, SiteType.strand_id == Strand.id
    ).all()
    logging.info(
        "Loaded %d %s transcript rows in %.2f seconds",
        len(transcript_rows),
        SiteType.__name__,
        time.perf_counter() - query_started,
    )

    for chr_, strand, site, transcript_id in transcript_rows:
        boundary_index.setdefault((chr_, strand, site), set()).add(
            transcript_id
        )
    del transcript_rows

    format_started = time.perf_counter()
    for key, transcript_ids in boundary_index.items():
        transcript_ids = frozenset(transcript_ids)
        if emit_details:
            gene_ids, transcript_ids_text = format_gene_and_transcript_ids(
                transcript_ids,
                transcript_gene_index,
            )
            boundary_index[key] = (
                transcript_ids,
                gene_ids,
                transcript_ids_text,
            )
        else:
            boundary_index[key] = transcript_ids
    logging.info(
        "Finalized %d immutable %s boundary records in %.2f seconds "
        "(details=%s)",
        len(boundary_index),
        SiteType.__name__,
        time.perf_counter() - format_started,
        emit_details,
    )

    return boundary_index


def build_transcript_gene_index(annotation_db):
    rows = annotation_db.session.query(
        Transcript.transcript_id,
        Gene.gene_id,
        Gene.gene_symbol,
    ).join(
        Gene, Transcript.gid == Gene.id
    ).all()

    return {
        transcript_id: (gene_id, gene_symbol)
        for transcript_id, gene_id, gene_symbol in rows
    }


def format_gene_and_transcript_ids(
        transcript_ids,
        transcript_gene_index):
    if not transcript_ids:
        return '', ''

    transcripts_by_gene = {}
    for transcript_id in sorted(transcript_ids):
        try:
            gene_id = transcript_gene_index[transcript_id][0]
        except KeyError:
            raise ValueError(
                f"Transcript ID missing gene mapping: {transcript_id}"
            ) from None
        transcripts_by_gene.setdefault(gene_id, []).append(transcript_id)

    gene_ids = sorted(transcripts_by_gene)
    return (
        ';'.join(gene_ids),
        ';'.join(
            ','.join(transcripts_by_gene[gene_id])
            for gene_id in gene_ids
        ),
    )


def format_shared_boundary_record_uncached(transcript_ids):
    if not transcript_ids:
        return '', '', ''

    transcripts_by_gene = {}
    gene_symbols_by_id = {}
    for transcript_id in sorted(transcript_ids):
        try:
            gene_id, gene_symbol = \
                _WORKER_TRANSCRIPT_GENE_INDEX[transcript_id]
        except KeyError:
            raise ValueError(
                f"Transcript ID missing gene mapping: {transcript_id}"
            ) from None
        transcripts_by_gene.setdefault(gene_id, []).append(transcript_id)
        gene_symbols_by_id[gene_id] = gene_symbol or ''

    gene_ids = sorted(transcripts_by_gene)
    return (
        ';'.join(gene_ids),
        ';'.join(gene_symbols_by_id[gene_id] for gene_id in gene_ids),
        ';'.join(
            ','.join(transcripts_by_gene[gene_id])
            for gene_id in gene_ids
        ),
    )


def initialize_worker(cache_size):
    global _WORKER_FORMAT_SHARED_RECORD
    _WORKER_FORMAT_SHARED_RECORD = lru_cache(maxsize=cache_size)(
        format_shared_boundary_record_uncached
    )


def get_required_column_indexes(header):
    column_indexes = {}
    for column in REQUIRED_COLUMNS:
        try:
            column_indexes[column] = header.index(column)
        except ValueError:
            raise ValueError(
                f"Missing required input column: {column}"
            ) from None
    return column_indexes


def get_required_values(row, column_indexes, row_label):
    missing_columns = [
        column
        for column, index in column_indexes.items()
        if index >= len(row) or row[index] == ''
    ]
    if missing_columns:
        raise ValueError(
            f"Missing required column(s) at {row_label}: "
            + ", ".join(missing_columns)
        )
    return tuple(row[index] for index in column_indexes.values())


def format_circRNA_annotation(values, donor_record, acceptor_record):
    (
        chr_,
        pos1,
        pos2,
        strand,
        circ_id,
        donor,
        acceptor,
    ) = values
    donor_site = int(donor)
    acceptor_site = int(acceptor)

    donor_at_boundary = int(donor_record is not None)
    acceptor_at_boundary = int(acceptor_record is not None)
    if EMIT_DETAILS:
        if donor_record is None:
            donor_record = EMPTY_BOUNDARY_RECORD
        if acceptor_record is None:
            acceptor_record = EMPTY_BOUNDARY_RECORD
        (
            donor_transcript_ids,
            donor_gene_ids,
            donor_transcript_ids_text,
        ) = donor_record
        (
            acceptor_transcript_ids,
            acceptor_gene_ids,
            acceptor_transcript_ids_text,
        ) = acceptor_record
    else:
        donor_transcript_ids = donor_record or EMPTY_TRANSCRIPT_IDS
        acceptor_transcript_ids = acceptor_record or EMPTY_TRANSCRIPT_IDS

    if donor_transcript_ids and acceptor_transcript_ids:
        common_transcript_ids = (
            donor_transcript_ids & acceptor_transcript_ids
        )
    else:
        common_transcript_ids = EMPTY_TRANSCRIPT_IDS

    both_at_boundary = int(donor_at_boundary and acceptor_at_boundary)
    same_transcript_isoform_count = len(common_transcript_ids)
    same_transcript_isoform = int(same_transcript_isoform_count > 0)

    summary_row = [
        chr_,
        pos1,
        pos2,
        strand,
        circ_id,
        donor_at_boundary,
        acceptor_at_boundary,
        both_at_boundary,
        same_transcript_isoform,
        same_transcript_isoform_count,
    ]
    details_row = None
    if EMIT_DETAILS:
        (
            shared_gene_ids,
            shared_gene_symbols,
            shared_transcript_ids_text,
        ) = _WORKER_FORMAT_SHARED_RECORD(common_transcript_ids)
        details_row = [
            chr_,
            pos1,
            pos2,
            strand,
            circ_id,
            donor_site,
            acceptor_site,
            donor_gene_ids,
            donor_transcript_ids_text,
            acceptor_gene_ids,
            acceptor_transcript_ids_text,
            shared_gene_ids,
            shared_gene_symbols,
            shared_transcript_ids_text,
        ]
    return summary_row, details_row


def check_circRNA(values):
    chr_, _, _, strand, _, donor, acceptor = values
    donor_record = _WORKER_DONOR_BOUNDARY_INDEX.get(
        (chr_, strand, int(donor))
    )
    acceptor_record = _WORKER_ACCEPTOR_BOUNDARY_INDEX.get(
        (chr_, strand, int(acceptor))
    )
    return format_circRNA_annotation(values, donor_record, acceptor_record)


def read_header_and_data_start(path):
    with open(path, 'rb') as f_in:
        header_line = f_in.readline()
        if not header_line:
            raise ValueError("Input circRNA file is empty")
        data_start = f_in.tell()

    try:
        header_text = header_line.decode('utf-8')
    except UnicodeDecodeError:
        raise ValueError("Input circRNA header is not valid UTF-8") from None
    try:
        header = next(csv.reader([header_text], delimiter='\t'))
    except StopIteration:
        raise ValueError("Input circRNA header is empty") from None

    return header, data_start


def count_records_up_to(path, data_start, limit):
    if limit <= 0:
        return None
    count = 0
    with open(path, 'rb') as f_in:
        f_in.seek(data_start)
        for raw_line in f_in:
            if not raw_line.rstrip(b'\r\n'):
                continue
            count += 1
            if count > limit:
                return None
    return count


def open_annotation_sqlite(path):
    database_path = Path(path).resolve()
    connection = sqlite3.connect(
        f'{database_path.as_uri()}?mode=ro&immutable=1',
        uri=True,
    )
    connection.execute('PRAGMA query_only = ON')
    return connection


def load_targeted_boundary(
        connection,
        chromosome_ids,
        strand_ids,
        site_type,
        chromosome,
        strand,
        site):
    chromosome_id = chromosome_ids.get(chromosome)
    strand_id = strand_ids.get(strand)
    if chromosome_id is None or strand_id is None:
        return None

    if site_type == 'donor':
        site_table = 'donor_site'
        site_index = 'donor_site_index'
        exon_column = 'donor_id'
        exon_index = 'ix_exon_donor_id'
    elif site_type == 'acceptor':
        site_table = 'acceptor_site'
        site_index = 'acceptor_site_index'
        exon_column = 'acceptor_id'
        exon_index = 'ix_exon_acceptor_id'
    else:
        raise ValueError(f'Unsupported site type: {site_type}')

    query = f'''
        SELECT s.id, te.tid
        FROM {site_table} AS s INDEXED BY {site_index}
        LEFT JOIN exon AS e INDEXED BY {exon_index}
          ON e.{exon_column} = s.id
        LEFT JOIN transcript_exon AS te INDEXED BY ix_transcript_exon_eid
          ON te.eid = e.id
        WHERE s.chr_id = ? AND s.junc_site = ? AND s.strand_id = ?
    '''
    rows = connection.execute(
        query,
        (chromosome_id, int(site), strand_id),
    )
    boundary_found = False
    transcript_ids = set()
    for _site_id, transcript_id in rows:
        boundary_found = True
        if transcript_id is not None:
            transcript_ids.add(transcript_id)
    return frozenset(transcript_ids) if boundary_found else None


def process_targeted_queries(input_path, output_path, column_indexes):
    started = time.perf_counter()
    connection = open_annotation_sqlite(snakemake.input.anno_db)
    chromosome_ids = dict(connection.execute('SELECT name, id FROM chromosome'))
    strand_ids = dict(connection.execute('SELECT name, id FROM strand'))
    boundary_cache = {}
    processed_rows = 0

    output_directory = os.path.dirname(output_path) or '.'
    os.makedirs(output_directory, exist_ok=True)
    output_handle = tempfile.NamedTemporaryFile(
        mode='w',
        newline='',
        encoding='utf-8',
        dir=output_directory,
        prefix=f'.{os.path.basename(output_path)}.',
        suffix='.tmp',
        delete=False,
    )
    temporary_output = output_handle.name
    try:
        writer = csv.writer(output_handle, delimiter='\t', lineterminator='\n')
        writer.writerow(SUMMARY_COLUMNS)
        with open(input_path, newline='', encoding='utf-8') as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            next(reader)
            for line_number, row in enumerate(reader, start=2):
                if not row:
                    continue
                values = get_required_values(
                    row,
                    column_indexes,
                    f'input line {line_number}',
                )
                chromosome, _, _, strand, _, donor, acceptor = values
                donor_key = ('donor', chromosome, strand, int(donor))
                acceptor_key = ('acceptor', chromosome, strand, int(acceptor))
                if donor_key not in boundary_cache:
                    boundary_cache[donor_key] = load_targeted_boundary(
                        connection,
                        chromosome_ids,
                        strand_ids,
                        *donor_key,
                    )
                if acceptor_key not in boundary_cache:
                    boundary_cache[acceptor_key] = load_targeted_boundary(
                        connection,
                        chromosome_ids,
                        strand_ids,
                        *acceptor_key,
                    )
                summary_row, _ = format_circRNA_annotation(
                    values,
                    boundary_cache[donor_key],
                    boundary_cache[acceptor_key],
                )
                writer.writerow(summary_row)
                processed_rows += 1
        output_handle.flush()
        os.fsync(output_handle.fileno())
        output_handle.close()
        os.replace(temporary_output, output_path)
    except Exception:
        if not output_handle.closed:
            output_handle.close()
        try:
            os.unlink(temporary_output)
        except FileNotFoundError:
            pass
        raise
    finally:
        connection.close()

    elapsed = time.perf_counter() - started
    logging.info(
        'Finished targeted annotation for %d circRNAs with %d unique '
        'boundaries in %.2f seconds (%.2f rows/second)',
        processed_rows,
        len(boundary_cache),
        elapsed,
        processed_rows / elapsed if elapsed else 0,
    )
    return processed_rows


def plan_byte_range_chunks(path, data_start, requested_chunks):
    file_size = os.path.getsize(path)
    if data_start >= file_size:
        return []

    boundaries = [data_start]
    data_size = file_size - data_start
    with open(path, 'rb') as f_in:
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
    (
        chunk_id,
        start,
        end,
        input_path,
        summary_part,
        details_part,
    ) = task
    started = time.perf_counter()
    processed_rows = 0
    summary_batch = []
    details_batch = []
    cache_before = _WORKER_FORMAT_SHARED_RECORD.cache_info()

    with ExitStack() as stack:
        f_in = stack.enter_context(open(input_path, 'rb'))
        f_summary = stack.enter_context(open(
            summary_part,
            'w',
            newline='',
            encoding='utf-8',
            buffering=IO_BUFFER_SIZE,
        ))
        f_details = None
        if EMIT_DETAILS:
            f_details = stack.enter_context(open(
                details_part,
                'w',
                newline='',
                encoding='utf-8',
                buffering=IO_BUFFER_SIZE,
            ))
        f_in.seek(start)
        summary_writer = csv.writer(
            f_summary,
            delimiter='\t',
            lineterminator='\n',
        )
        details_writer = None
        if f_details is not None:
            details_writer = csv.writer(
                f_details,
                delimiter='\t',
                lineterminator='\n',
            )

        record_num = 0
        while f_in.tell() < end:
            byte_offset = f_in.tell()
            raw_line = f_in.readline()
            if not raw_line:
                break
            record_num += 1
            try:
                line = raw_line.decode('utf-8').rstrip('\r\n')
            except UnicodeDecodeError:
                raise ValueError(
                    f"Input row at byte offset {byte_offset} is not "
                    "valid UTF-8"
                ) from None
            if not line:
                continue
            row = line.split('\t')

            row_label = (
                f"chunk {chunk_id}, record {record_num}, "
                f"byte offset {byte_offset}"
            )
            values = get_required_values(
                row,
                _WORKER_COLUMN_INDEXES,
                row_label,
            )
            try:
                summary_row, details_row = check_circRNA(values)
            except ValueError as exc:
                raise ValueError(f"Invalid circRNA at {row_label}: {exc}") \
                    from exc
            summary_batch.append(summary_row)
            if details_writer is not None:
                details_batch.append(details_row)
            processed_rows += 1

            if len(summary_batch) >= WRITE_BATCH_SIZE:
                summary_writer.writerows(summary_batch)
                if details_writer is not None:
                    details_writer.writerows(details_batch)
                summary_batch.clear()
                details_batch.clear()

        if summary_batch:
            summary_writer.writerows(summary_batch)
            if details_writer is not None:
                details_writer.writerows(details_batch)

    cache_after = _WORKER_FORMAT_SHARED_RECORD.cache_info()
    return {
        'chunk_id': chunk_id,
        'summary_part': summary_part,
        'details_part': details_part,
        'processed_rows': processed_rows,
        'elapsed': time.perf_counter() - started,
        'pid': os.getpid(),
        'cache_hits': cache_after.hits - cache_before.hits,
        'cache_misses': cache_after.misses - cache_before.misses,
    }


def concatenate_parts(output_path, columns, part_paths):
    with open(output_path, 'wb', buffering=IO_BUFFER_SIZE) as f_out:
        f_out.write(('\t'.join(columns) + '\n').encode('utf-8'))
        for part_path in part_paths:
            with open(
                    part_path,
                    'rb',
                    buffering=IO_BUFFER_SIZE) as f_part:
                shutil.copyfileobj(
                    f_part,
                    f_out,
                    length=IO_BUFFER_SIZE,
                )


def main():
    global _WORKER_DONOR_BOUNDARY_INDEX
    global _WORKER_ACCEPTOR_BOUNDARY_INDEX
    global _WORKER_TRANSCRIPT_GENE_INDEX
    global _WORKER_COLUMN_INDEXES

    total_started = time.perf_counter()
    logging.info(
        "Starting donor/acceptor boundary annotation (details=%s)",
        EMIT_DETAILS,
    )
    part_directory = None
    pool_started = False

    try:
        header, data_start = read_header_and_data_start(
            snakemake.input.circRNAs
        )
        column_indexes = get_required_column_indexes(header)
        targeted_threshold = int(
            getattr(snakemake.params, 'targeted_threshold', 0)
        )
        targeted_row_count = None
        if not EMIT_DETAILS:
            targeted_row_count = count_records_up_to(
                snakemake.input.circRNAs,
                data_start,
                targeted_threshold,
            )
        if targeted_row_count is not None:
            logging.info(
                "Using targeted SQLite annotation for %d circRNAs "
                "(threshold=%d)",
                targeted_row_count,
                targeted_threshold,
            )
            process_targeted_queries(
                snakemake.input.circRNAs,
                snakemake.output.summary,
                column_indexes,
            )
            logging.info(
                "Boundary annotation completed in %.2f seconds",
                time.perf_counter() - total_started,
            )
            return

        logging.info(
            "Using full in-memory annotation index "
            "(targeted threshold=%d)",
            targeted_threshold,
        )
        phase_started = time.perf_counter()
        annotation_db = Annotation(snakemake.input.anno_db)
        logging.info(
            "Opened annotation DB in %.2f seconds",
            time.perf_counter() - phase_started,
        )

        transcript_gene_index = {}
        if EMIT_DETAILS:
            phase_started = time.perf_counter()
            transcript_gene_index = build_transcript_gene_index(annotation_db)
            logging.info(
                "Built transcript-to-gene index with %d transcripts "
                "in %.2f seconds",
                len(transcript_gene_index),
                time.perf_counter() - phase_started,
            )
        else:
            logging.info("Skipped transcript-to-gene index in summary-only mode")

        phase_started = time.perf_counter()
        donor_boundary_index = build_boundary_transcript_index(
            annotation_db,
            DonorSite,
            Exon.donor_id,
            transcript_gene_index,
            EMIT_DETAILS,
        )
        logging.info(
            "Built donor boundary index with %d sites in %.2f seconds",
            len(donor_boundary_index),
            time.perf_counter() - phase_started,
        )

        phase_started = time.perf_counter()
        acceptor_boundary_index = build_boundary_transcript_index(
            annotation_db,
            AcceptorSite,
            Exon.acceptor_id,
            transcript_gene_index,
            EMIT_DETAILS,
        )
        logging.info(
            "Built acceptor boundary index with %d sites in %.2f seconds",
            len(acceptor_boundary_index),
            time.perf_counter() - phase_started,
        )
        annotation_db.session.close()
        del annotation_db

        worker_count = max(1, int(snakemake.threads))
        chunks_per_worker = int(snakemake.params.chunks_per_worker)
        cache_size = int(snakemake.params.shared_cache_size)
        chunks = plan_byte_range_chunks(
            snakemake.input.circRNAs,
            data_start,
            worker_count * chunks_per_worker,
        )
        logging.info(
            "Planned %d byte-range chunks for %d workers",
            len(chunks),
            worker_count,
        )

        output_directory = os.path.dirname(snakemake.output.summary) or '.'
        os.makedirs(output_directory, exist_ok=True)
        part_directory = tempfile.mkdtemp(
            prefix='.annotation_boundaries.parts.',
            dir=output_directory,
        )
        tasks = []
        for chunk_id, start, end in chunks:
            tasks.append((
                chunk_id,
                start,
                end,
                snakemake.input.circRNAs,
                os.path.join(
                    part_directory,
                    f'summary.{chunk_id:05d}.tsv',
                ),
                os.path.join(
                    part_directory,
                    f'details.{chunk_id:05d}.tsv',
                ),
            ))

        _WORKER_DONOR_BOUNDARY_INDEX = donor_boundary_index
        _WORKER_ACCEPTOR_BOUNDARY_INDEX = acceptor_boundary_index
        _WORKER_TRANSCRIPT_GENE_INDEX = transcript_gene_index
        _WORKER_COLUMN_INDEXES = column_indexes

        processing_started = time.perf_counter()
        results = []
        if tasks:
            gc.collect()
            gc.freeze()
            pool_started = True
            context = multiprocessing.get_context('fork')
            with context.Pool(
                    processes=worker_count,
                    initializer=initialize_worker,
                    initargs=(cache_size,)) as pool:
                for result in pool.imap_unordered(
                        process_chunk,
                        tasks,
                        chunksize=1):
                    results.append(result)
                    logging.info(
                        "Worker %d finished chunk %d: %d circRNAs "
                        "in %.2f seconds (%.2f rows/second), cache "
                        "hits=%d misses=%d",
                        result['pid'],
                        result['chunk_id'],
                        result['processed_rows'],
                        result['elapsed'],
                        (
                            result['processed_rows'] / result['elapsed']
                            if result['elapsed']
                            else 0
                        ),
                        result['cache_hits'],
                        result['cache_misses'],
                    )
            gc.unfreeze()
            pool_started = False

        results.sort(key=lambda result: result['chunk_id'])
        concatenate_started = time.perf_counter()
        staged_summary = os.path.join(
            part_directory,
            'annotation_boundaries.final.tsv',
        )
        concatenate_parts(
            staged_summary,
            SUMMARY_COLUMNS,
            [result['summary_part'] for result in results],
        )
        os.replace(staged_summary, snakemake.output.summary)
        if EMIT_DETAILS:
            staged_details = os.path.join(
                part_directory,
                'annotation_boundary_transcripts.final.tsv',
            )
            concatenate_parts(
                staged_details,
                DETAIL_COLUMNS,
                [result['details_part'] for result in results],
            )
            os.replace(staged_details, snakemake.output.details)

        processed_rows = sum(
            result['processed_rows']
            for result in results
        )
        cache_hits = sum(result['cache_hits'] for result in results)
        cache_misses = sum(result['cache_misses'] for result in results)
        processing_elapsed = time.perf_counter() - processing_started
        logging.info(
            "Concatenated %d annotation chunks in %.2f seconds "
            "(details=%s)",
            len(results),
            time.perf_counter() - concatenate_started,
            EMIT_DETAILS,
        )
        logging.info(
            "Shared boundary cache totals: %d hits, %d misses "
            "(%.2f%% hit rate)",
            cache_hits,
            cache_misses,
            (
                100 * cache_hits / (cache_hits + cache_misses)
                if cache_hits + cache_misses
                else 0
            ),
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
            "Boundary annotation completed in %.2f seconds",
            time.perf_counter() - total_started,
        )
    except Exception:
        if pool_started:
            gc.unfreeze()
        logging.exception(
            "Boundary annotation failed after %.2f seconds",
            time.perf_counter() - total_started,
        )
        raise
    finally:
        if part_directory is not None:
            shutil.rmtree(part_directory, ignore_errors=True)


main()
