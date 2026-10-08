import csv
import logging
import os
import tempfile
import time


IO_BUFFER_SIZE = 8 * 1024 * 1024
PROGRESS_INTERVAL = 500_000
REQUIRED_CIRCRNA_COLUMNS = ("chr", "donor", "acceptor", "strand")
INDEX_COLUMNS = (
    "chr",
    "pos",
    "strand",
    "is_annotated_splicing_site",
    "has_AS_event",
    "AS_event_type",
)
FEATURE_COLUMNS = (
    "is_annotated_splicing_site(donor)",
    "has_AS_event(donor)",
    "is_annotated_splicing_site(acceptor)",
    "has_AS_event(acceptor)",
    "AS_event_type(donor)",
    "AS_event_type(acceptor)",
)
NO_AS_EVENT = ""
EXON_TO_EXON = "exon_to_exon"
EXON_TO_INTRON = "exon_to_intron"
VALID_AS_EVENT_TYPES = {
    NO_AS_EVENT,
    EXON_TO_EXON,
    EXON_TO_INTRON,
}
NO_ANNOTATION = (0, 0, NO_AS_EVENT)
ANNOTATED_WITHOUT_AS = (1, 0, NO_AS_EVENT)


logging.basicConfig(
    filename=snakemake.log[0],
    encoding="utf-8",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def get_column_indexes(header, required_columns, source_label):
    missing_columns = [
        column for column in required_columns if column not in header
    ]
    if missing_columns:
        raise ValueError(
            f"Missing required column(s) in {source_label}: "
            + ", ".join(missing_columns)
        )
    return {column: header.index(column) for column in required_columns}


def parse_binary_value(value, column, line_num):
    if value not in {"0", "1"}:
        raise ValueError(
            f"Invalid {column} value at AS index line {line_num}: "
            f"{value!r}"
        )
    return int(value)


def load_site_index(index_path):
    started = time.perf_counter()
    site_index = {}

    with open(
            index_path,
            newline="",
            encoding="utf-8",
            buffering=IO_BUFFER_SIZE) as f_index:
        reader = csv.reader(f_index, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError("Alternative-splicing site index is empty") from None

        indexes = get_column_indexes(header, INDEX_COLUMNS, "AS site index")
        for line_num, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"Expected {len(header)} columns at AS index line "
                    f"{line_num}, got {len(row)}"
                )
            chr_ = row[indexes["chr"]]
            strand = row[indexes["strand"]]
            try:
                position = int(row[indexes["pos"]])
            except ValueError:
                raise ValueError(
                    f"Invalid position at AS index line {line_num}: "
                    f"{row[indexes['pos']]!r}"
                ) from None

            is_annotated = parse_binary_value(
                row[indexes["is_annotated_splicing_site"]],
                "is_annotated_splicing_site",
                line_num,
            )
            has_as_event = parse_binary_value(
                row[indexes["has_AS_event"]],
                "has_AS_event",
                line_num,
            )
            event_type = row[indexes["AS_event_type"]]
            if event_type not in VALID_AS_EVENT_TYPES:
                raise ValueError(
                    f"Invalid AS_event_type at AS index line {line_num}: "
                    f"{event_type!r}"
                )
            if not is_annotated:
                raise ValueError(
                    "AS site index must contain annotated sites only; "
                    f"found 0 at line {line_num}"
                )
            if has_as_event != int(bool(event_type)):
                raise ValueError(
                    "Inconsistent has_AS_event and AS_event_type at AS "
                    f"index line {line_num}: {has_as_event!r}, "
                    f"{event_type!r}"
                )
            region_index = site_index.setdefault((chr_, strand), {})
            if position in region_index:
                raise ValueError(
                    "Duplicate coordinate in AS site index at line "
                    f"{line_num}: ({chr_}, {position}, {strand})"
                )
            region_index[position] = event_type

    site_count = sum(len(region) for region in site_index.values())
    logging.info(
        "Loaded %d AS sites across %d chromosome/strand regions in "
        "%.2f seconds",
        site_count,
        len(site_index),
        time.perf_counter() - started,
    )
    return site_index


def parse_position(value, column, line_num):
    try:
        position = int(value)
    except ValueError:
        raise ValueError(
            f"Invalid {column} position at circRNA line {line_num}: "
            f"{value!r}"
        ) from None
    if position < 1:
        raise ValueError(
            f"Invalid {column} position at circRNA line {line_num}: "
            f"{position}"
        )
    return position


def get_site_annotation(region_index, position):
    if region_index is None:
        return NO_ANNOTATION
    event_type = region_index.get(position)
    if event_type is None:
        return NO_ANNOTATION
    if not event_type:
        return ANNOTATED_WITHOUT_AS
    return (1, 1, event_type)


def annotate_circrnas(circrna_path, output_path, site_index):
    output_directory = os.path.dirname(output_path) or "."
    os.makedirs(output_directory, exist_ok=True)
    temporary_path = None
    started = time.perf_counter()
    row_count = 0
    donor_annotated_count = 0
    acceptor_annotated_count = 0

    try:
        with open(
                circrna_path,
                newline="",
                encoding="utf-8",
                buffering=IO_BUFFER_SIZE) as f_in, tempfile.NamedTemporaryFile(
                    mode="w",
                    newline="",
                    encoding="utf-8",
                    buffering=IO_BUFFER_SIZE,
                    prefix=".alternative-splicing.",
                    suffix=".tsv",
                    dir=output_directory,
                    delete=False) as f_out:
            temporary_path = f_out.name
            reader = csv.reader(f_in, delimiter="\t")
            writer = csv.writer(
                f_out,
                delimiter="\t",
                lineterminator="\n",
            )

            try:
                header = next(reader)
            except StopIteration:
                raise ValueError("circRNA input table is empty") from None

            indexes = get_column_indexes(
                header,
                REQUIRED_CIRCRNA_COLUMNS,
                "circRNA table",
            )
            duplicated_features = [
                column for column in FEATURE_COLUMNS if column in header
            ]
            if duplicated_features:
                raise ValueError(
                    "circRNA table already contains AS feature column(s): "
                    + ", ".join(duplicated_features)
                )
            writer.writerow([*header, *FEATURE_COLUMNS])

            for line_num, row in enumerate(reader, start=2):
                if len(row) != len(header):
                    raise ValueError(
                        f"Expected {len(header)} columns at circRNA line "
                        f"{line_num}, got {len(row)}"
                    )
                chr_ = row[indexes["chr"]]
                strand = row[indexes["strand"]]
                if strand not in {"+", "-"}:
                    raise ValueError(
                        f"Invalid strand at circRNA line {line_num}: "
                        f"{strand!r}"
                    )
                donor = parse_position(
                    row[indexes["donor"]],
                    "donor",
                    line_num,
                )
                acceptor = parse_position(
                    row[indexes["acceptor"]],
                    "acceptor",
                    line_num,
                )
                region_index = site_index.get((chr_, strand))
                donor_annotation = get_site_annotation(region_index, donor)
                acceptor_annotation = get_site_annotation(
                    region_index,
                    acceptor,
                )

                writer.writerow((
                    *row,
                    donor_annotation[0],
                    donor_annotation[1],
                    acceptor_annotation[0],
                    acceptor_annotation[1],
                    donor_annotation[2],
                    acceptor_annotation[2],
                ))
                row_count += 1
                donor_annotated_count += donor_annotation[0]
                acceptor_annotated_count += acceptor_annotation[0]

                if row_count % PROGRESS_INTERVAL == 0:
                    elapsed = time.perf_counter() - started
                    logging.info(
                        "Processed %d circRNAs in %.2f seconds "
                        "(%.2f rows/second)",
                        row_count,
                        elapsed,
                        row_count / elapsed,
                    )

        os.replace(temporary_path, output_path)
    except Exception:
        if temporary_path and os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise

    elapsed = time.perf_counter() - started
    logging.info(
        "Annotated %d circRNAs in %.2f seconds (%.2f rows/second); "
        "annotated donor sites: %d; annotated acceptor sites: %d",
        row_count,
        elapsed,
        row_count / elapsed if elapsed else 0,
        donor_annotated_count,
        acceptor_annotated_count,
    )


def main():
    total_started = time.perf_counter()
    logging.info("Starting circRNA alternative-splicing annotation")
    site_index = load_site_index(snakemake.input.site_index)
    annotate_circrnas(
        snakemake.input.circRNAs,
        snakemake.output[0],
        site_index,
    )
    logging.info(
        "circRNA alternative-splicing annotation completed in %.2f seconds",
        time.perf_counter() - total_started,
    )


main()
