import csv
import logging
import os
import tempfile
import time
from bisect import bisect_left, bisect_right
from itertools import groupby

from circmimi.annotation import Annotation
from circmimi.models import (
    Chromosome,
    Exon,
    Gene,
    Strand,
    Transcript,
    TranscriptExon,
)


IO_BUFFER_SIZE = 8 * 1024 * 1024
QUERY_BATCH_SIZE = 100_000
PROGRESS_INTERVAL = 10_000
START_BOUNDARY = 1
END_BOUNDARY = 2
NO_AS_EVENT = ""
EXON_TO_EXON = "exon_to_exon"
EXON_TO_INTRON = "exon_to_intron"
AS_EVENT_TYPE_PRIORITY = {
    NO_AS_EVENT: 0,
    EXON_TO_EXON: 1,
    EXON_TO_INTRON: 2,
}

OUTPUT_COLUMNS = (
    "chr",
    "pos",
    "strand",
    "is_annotated_splicing_site",
    "has_AS_event",
    "AS_event_type",
)


logging.basicConfig(
    filename=snakemake.log[0],
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(set(intervals)):
        if merged and start <= merged[-1][1] + 1:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def add_interval_coverage(differences, positions, start, end):
    left = bisect_left(positions, start)
    right = bisect_right(positions, end)
    if left < right:
        differences[left] += 1
        differences[right] -= 1


def interval_contains(intervals, starts, position):
    interval_index = bisect_right(starts, position) - 1
    return (
        interval_index >= 0
        and intervals[interval_index][1] >= position
    )


def group_exons_by_transcript(rows):
    transcripts = []
    for _, transcript_rows in groupby(rows, key=lambda row: row[3]):
        exons = sorted({
            (int(row[4]), int(row[5])) for row in transcript_rows
        })
        if not exons:
            continue
        transcripts.append((
            min(start for start, _ in exons),
            max(end for _, end in exons),
            exons,
            merge_intervals(exons),
        ))
    return transcripts


def annotate_gene_boundaries(transcripts):
    boundary_sides = {}
    all_exons = []

    for transcript_start, transcript_end, exons, merged_exons in transcripts:
        all_exons.extend(merged_exons)

        for exon_start, exon_end in exons:
            if exon_start != transcript_start:
                boundary_sides[exon_start] = (
                    boundary_sides.get(exon_start, 0) | START_BOUNDARY
                )
            if exon_end != transcript_end:
                boundary_sides[exon_end] = (
                    boundary_sides.get(exon_end, 0) | END_BOUNDARY
                )

    if not boundary_sides:
        return {}

    positions = sorted(boundary_sides)
    active_transcript_differences = [0] * (len(positions) + 1)
    exonic_transcript_differences = [0] * (len(positions) + 1)

    # Case 1 以差分計數找出「位於 transcript span、但不在 exon」的座標。
    for transcript_start, transcript_end, _, merged_exons in transcripts:
        add_interval_coverage(
            active_transcript_differences,
            positions,
            transcript_start,
            transcript_end,
        )
        for exon_start, exon_end in merged_exons:
            add_interval_coverage(
                exonic_transcript_differences,
                positions,
                exon_start,
                exon_end,
            )

    # Case 2 查詢 internal boundary 外側相鄰鹼基是否仍被同 gene exon 涵蓋。
    merged_gene_exons = merge_intervals(all_exons)
    merged_gene_exon_starts = [start for start, _ in merged_gene_exons]
    active_transcripts = 0
    exonic_transcripts = 0
    results = {}

    for index, position in enumerate(positions):
        active_transcripts += active_transcript_differences[index]
        exonic_transcripts += exonic_transcript_differences[index]
        has_case_1 = active_transcripts > exonic_transcripts

        sides = boundary_sides[position]
        has_case_2 = (
            bool(sides & START_BOUNDARY)
            and interval_contains(
                merged_gene_exons,
                merged_gene_exon_starts,
                position - 1,
            )
        ) or (
            bool(sides & END_BOUNDARY)
            and interval_contains(
                merged_gene_exons,
                merged_gene_exon_starts,
                position + 1,
            )
        )
        if has_case_1:
            event_type = EXON_TO_INTRON
        elif has_case_2:
            event_type = EXON_TO_EXON
        else:
            event_type = NO_AS_EVENT
        results[position] = event_type

    return results


def iter_annotation_rows(annotation_db):
    return annotation_db.session.query(
        Chromosome.name,
        Strand.name,
        Gene.id,
        Transcript.id,
        Exon.start,
        Exon.end,
    ).select_from(
        Transcript
    ).join(
        Gene, Transcript.gid == Gene.id
    ).join(
        TranscriptExon, Transcript.id == TranscriptExon.tid
    ).join(
        Exon, TranscriptExon.eid == Exon.id
    ).join(
        Chromosome, Exon.chr_id == Chromosome.id
    ).join(
        Strand, Exon.strand_id == Strand.id
    ).order_by(
        Chromosome.name,
        Strand.name,
        Gene.id,
        Transcript.id,
        Exon.start,
        Exon.end,
    ).yield_per(
        QUERY_BATCH_SIZE
    ).execution_options(
        stream_results=True
    )


def build_site_index(annotation_db):
    site_index = {}
    gene_count = 0
    started = time.perf_counter()
    rows = iter_annotation_rows(annotation_db)

    for group_key, gene_rows in groupby(
            rows, key=lambda row: (row[0], row[1], row[2])):
        chr_, strand, _ = group_key
        transcripts = group_exons_by_transcript(gene_rows)
        gene_results = annotate_gene_boundaries(transcripts)
        region_index = site_index.setdefault((chr_, strand), {})
        for position, event_type in gene_results.items():
            previous_type = region_index.get(position)
            if (
                    previous_type is None
                    or AS_EVENT_TYPE_PRIORITY[event_type]
                    > AS_EVENT_TYPE_PRIORITY[previous_type]):
                region_index[position] = event_type

        gene_count += 1
        if gene_count % PROGRESS_INTERVAL == 0:
            site_count = sum(len(region) for region in site_index.values())
            logging.info(
                "Processed %d gene regions and collected %d AS sites in "
                "%.2f seconds",
                gene_count,
                site_count,
                time.perf_counter() - started,
            )

    site_count = sum(len(region) for region in site_index.values())
    logging.info(
        "Built AS index for %d gene regions with %d unique sites in "
        "%.2f seconds",
        gene_count,
        site_count,
        time.perf_counter() - started,
    )
    return site_index, gene_count, site_count


def write_site_index(output_path, site_index):
    output_directory = os.path.dirname(output_path) or "."
    os.makedirs(output_directory, exist_ok=True)
    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
                mode="w",
                newline="",
                encoding="utf-8",
                buffering=IO_BUFFER_SIZE,
                prefix=".alternative-splicing-index.",
                suffix=".tsv",
                dir=output_directory,
                delete=False) as f_out:
            temporary_path = f_out.name
            writer = csv.writer(
                f_out,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writerow(OUTPUT_COLUMNS)
            for chr_, strand in sorted(site_index):
                for position in sorted(site_index[(chr_, strand)]):
                    event_type = site_index[(chr_, strand)][position]
                    writer.writerow((
                        chr_,
                        position,
                        strand,
                        1,
                        int(bool(event_type)),
                        event_type,
                    ))
        os.replace(temporary_path, output_path)
    except Exception:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise


def main():
    total_started = time.perf_counter()
    logging.info("Starting alternative-splicing site index construction")
    annotation_db = Annotation(snakemake.input.anno_db)

    try:
        site_index, gene_count, site_count = build_site_index(annotation_db)
    finally:
        annotation_db.session.close()

    write_started = time.perf_counter()
    write_site_index(snakemake.output.site_index, site_index)
    logging.info(
        "Wrote %d AS sites from %d gene regions in %.2f seconds",
        site_count,
        gene_count,
        time.perf_counter() - write_started,
    )
    logging.info(
        "Alternative-splicing site index completed in %.2f seconds",
        time.perf_counter() - total_started,
    )


main()
