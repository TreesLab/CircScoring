#!/usr/bin/env python

import argparse
import csv
import logging
import subprocess
import tempfile
import time
from pathlib import Path


REQUIRED_COLUMNS = ("chr", "BSJ_ID", "strand")
READ_COUNT_PREFIX = "#read_count("
VALID_STRANDS = {"+", "-"}


def add_reason(group, reason):
    if reason not in group["reasons"]:
        group["reasons"].append(reason)


def parse_bsj_id(value):
    fields = value.split("|")
    if len(fields) != 3 or not fields[0]:
        raise ValueError("BSJ_ID must use chrom|pos1|pos2")

    chrom, pos1_text, pos2_text = fields
    try:
        pos1 = int(pos1_text)
        pos2 = int(pos2_text)
    except ValueError as exc:
        raise ValueError("BSJ_ID coordinates must be integers") from exc

    if pos1 < 1 or pos2 < 1:
        raise ValueError("BSJ_ID coordinates must be positive")
    if pos1 > pos2:
        raise ValueError("BSJ_ID must satisfy pos1 <= pos2")

    return chrom, pos1, pos2


def new_group(bsj_id, strand, sample_count):
    group = {
        "source_bsj_id": bsj_id,
        "source_chr": "",
        "source_pos1": "",
        "source_pos2": "",
        "source_strand": strand,
        "counts": [0] * sample_count,
        "seen_numeric": [False] * sample_count,
        "reasons": [],
    }

    if strand not in VALID_STRANDS:
        add_reason(group, "strand must be + or -")

    try:
        chrom, pos1, pos2 = parse_bsj_id(bsj_id)
        group["source_chr"] = chrom
        group["source_pos1"] = pos1
        group["source_pos2"] = pos2
    except ValueError as exc:
        add_reason(group, str(exc))

    return group


def add_read_counts(group, row, read_count_columns):
    for index, column in enumerate(read_count_columns):
        value = row[column].strip()
        if value == "NA":
            continue

        try:
            count = int(value)
        except ValueError:
            add_reason(group, f"{column} must be a non-negative integer or NA")
            continue

        if count < 0:
            add_reason(group, f"{column} must be a non-negative integer or NA")
            continue

        group["counts"][index] += count
        group["seen_numeric"][index] = True


def read_and_group_table(input_path):
    groups = []
    group_indexes = {}

    with open(input_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("circFL table is empty and has no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("circFL table contains duplicate column names")

        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            raise ValueError(
                "circFL table is missing required column(s): " + ", ".join(missing)
            )

        read_count_columns = [
            column
            for column in reader.fieldnames
            if column.startswith(READ_COUNT_PREFIX)
        ]
        if not read_count_columns:
            raise ValueError(
                f"circFL table has no columns beginning with {READ_COUNT_PREFIX}"
            )

        row_count = 0
        for row_count, row in enumerate(reader, start=1):
            bsj_id = row["BSJ_ID"].strip()
            strand = row["strand"].strip()
            key = (bsj_id, strand)
            group_index = group_indexes.get(key)
            if group_index is None:
                group_index = len(groups)
                group_indexes[key] = group_index
                groups.append(new_group(bsj_id, strand, len(read_count_columns)))

            group = groups[group_index]
            row_chrom = row["chr"].strip()
            if group["source_chr"] and row_chrom != group["source_chr"]:
                add_reason(group, "chr column differs from BSJ_ID chromosome")
            add_read_counts(group, row, read_count_columns)

    logging.info("Input isoform rows: %s", row_count)
    logging.info("Source BSJ-strand groups: %s", len(groups))
    return groups, read_count_columns


def write_source_bed(groups, bed_path):
    valid_group_count = 0
    with open(bed_path, "w", encoding="utf-8") as handle:
        for group_index, group in enumerate(groups):
            if group["reasons"]:
                continue

            chrom = group["source_chr"]
            pos1 = group["source_pos1"]
            pos2 = group["source_pos2"]
            strand = group["source_strand"]
            handle.write(
                f"{chrom}\t{pos1 - 1}\t{pos1}\tg{group_index}|left\t0\t{strand}\n"
            )
            handle.write(
                f"{chrom}\t{pos2 - 1}\t{pos2}\tg{group_index}|right\t0\t{strand}\n"
            )
            valid_group_count += 1

    return valid_group_count


def run_liftover(executable, input_bed, chain_path, mapped_bed, unmapped_bed):
    subprocess.run(
        [
            executable,
            str(input_bed),
            str(chain_path),
            str(mapped_bed),
            str(unmapped_bed),
        ],
        check=True,
    )


def parse_mapped_bed(mapped_bed, group_count):
    mappings = {}
    with open(mapped_bed, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(
                    f"mapped BED line {line_number} has fewer than 6 columns"
                )

            name_fields = fields[3].split("|")
            if (
                len(name_fields) != 2
                or not name_fields[0].startswith("g")
                or name_fields[1] not in {"left", "right"}
            ):
                raise ValueError(
                    f"mapped BED line {line_number} has an invalid record name"
                )

            try:
                group_index = int(name_fields[0][1:])
                start0 = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"mapped BED line {line_number} contains invalid coordinates"
                ) from exc

            if group_index < 0 or group_index >= group_count:
                raise ValueError(
                    f"mapped BED line {line_number} refers to an unknown group"
                )

            key = (group_index, name_fields[1])
            mappings.setdefault(key, []).append(
                {
                    "chrom": fields[0],
                    "start0": start0,
                    "end": end,
                    "strand": fields[5],
                }
            )

    return mappings


def classify_mapping(group_index, mappings):
    left_mappings = mappings.get((group_index, "left"), [])
    right_mappings = mappings.get((group_index, "right"), [])

    if len(left_mappings) > 1 or len(right_mappings) > 1:
        return None, "ambiguous_mapping", "a BSJ endpoint has multiple mappings"
    if not left_mappings or not right_mappings:
        missing = []
        if not left_mappings:
            missing.append("left")
        if not right_mappings:
            missing.append("right")
        return None, "unmapped", f"missing mapped endpoint: {','.join(missing)}"

    left = left_mappings[0]
    right = right_mappings[0]
    if left["end"] - left["start0"] != 1 or right["end"] - right["start0"] != 1:
        return None, "inconsistent_mapping", "mapped endpoint is not 1 bp"
    if left["chrom"] != right["chrom"]:
        return None, "inconsistent_mapping", "mapped endpoints are on different chromosomes"
    if (
        left["strand"] not in VALID_STRANDS
        or right["strand"] not in VALID_STRANDS
        or left["strand"] != right["strand"]
    ):
        return None, "inconsistent_mapping", "mapped endpoint strands conflict"

    pos1, pos2 = sorted((left["end"], right["end"]))

    target = {
        "chr": left["chrom"],
        "pos1": pos1,
        "pos2": pos2,
        "strand": left["strand"],
    }
    target["BSJ_ID"] = f"{target['chr']}|{pos1}|{pos2}"
    return target, "", ""


def merge_counts(target_group, source_group):
    for index, seen in enumerate(source_group["seen_numeric"]):
        if seen:
            target_group["counts"][index] += source_group["counts"][index]
            target_group["seen_numeric"][index] = True


def format_counts(group):
    return [
        str(value) if group["seen_numeric"][index] else "NA"
        for index, value in enumerate(group["counts"])
    ]


def evaluate_groups(groups, mappings):
    target_groups = []
    target_indexes = {}
    rejected = []

    for group_index, group in enumerate(groups):
        if group["reasons"]:
            rejected.append((group, "invalid_input", "; ".join(group["reasons"])))
            continue

        target, status, reason = classify_mapping(group_index, mappings)
        if target is None:
            rejected.append((group, status, reason))
            continue

        target_key = (target["BSJ_ID"], target["strand"])
        target_index = target_indexes.get(target_key)
        if target_index is None:
            target_index = len(target_groups)
            target_indexes[target_key] = target_index
            target_group = {
                **target,
                "counts": [0] * len(group["counts"]),
                "seen_numeric": [False] * len(group["counts"]),
            }
            target_groups.append(target_group)

        merge_counts(target_groups[target_index], group)

    return target_groups, rejected


def write_outputs(
    output_path,
    rejected_path,
    read_count_columns,
    target_groups,
    rejected_groups,
):
    output_path = Path(output_path)
    rejected_path = Path(rejected_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rejected_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["BSJ_ID", "chr", "pos1", "pos2", "strand", *read_count_columns])
        for group in target_groups:
            writer.writerow(
                [
                    group["BSJ_ID"],
                    group["chr"],
                    group["pos1"],
                    group["pos2"],
                    group["strand"],
                    *format_counts(group),
                ]
            )

    rejected_header = [
        "source_BSJ_ID",
        "source_chr",
        "source_pos1",
        "source_pos2",
        "source_strand",
        *read_count_columns,
        "liftover_status",
        "liftover_reason",
    ]
    with open(rejected_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(rejected_header)
        for group, status, reason in rejected_groups:
            writer.writerow(
                [
                    group["source_bsj_id"],
                    group["source_chr"],
                    group["source_pos1"],
                    group["source_pos2"],
                    group["source_strand"],
                    *format_counts(group),
                    status,
                    reason,
                ]
            )


def lift_bsj_table(input_path, chain_path, output_path, rejected_path, liftover):
    started = time.perf_counter()
    groups, read_count_columns = read_and_group_table(input_path)

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_path = Path(temporary_directory)
        input_bed = temporary_path / "circFL_bsj.hg19.bed"
        mapped_bed = temporary_path / "circFL_bsj.hg38.bed"
        unmapped_bed = temporary_path / "circFL_bsj.unmapped.bed"

        valid_group_count = write_source_bed(groups, input_bed)
        logging.info("BSJ groups sent to liftOver: %s", valid_group_count)
        if valid_group_count:
            liftover_started = time.perf_counter()
            run_liftover(
                liftover,
                input_bed,
                chain_path,
                mapped_bed,
                unmapped_bed,
            )
            logging.info(
                "liftOver execution: %.2f seconds",
                time.perf_counter() - liftover_started,
            )
            mappings = parse_mapped_bed(mapped_bed, len(groups))
        else:
            mappings = {}

    target_groups, rejected_groups = evaluate_groups(groups, mappings)
    write_outputs(
        output_path,
        rejected_path,
        read_count_columns,
        target_groups,
        rejected_groups,
    )
    logging.info("Lifted target BSJ-strand groups: %s", len(target_groups))
    logging.info("Rejected source BSJ-strand groups: %s", len(rejected_groups))
    logging.info("Total execution: %.2f seconds", time.perf_counter() - started)


def setup_logging(log_path=None):
    handlers = [logging.StreamHandler()]
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        handlers = [logging.FileHandler(log_path)]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=handlers,
        force=True,
    )


def create_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Group circFL isoforms by BSJ and strand, lift BSJ endpoints, "
            "and sum sample read counts."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rejected", required=True)
    parser.add_argument("--liftover", default="liftOver")
    parser.add_argument("--log")
    return parser


def main():
    if "snakemake" in globals():
        setup_logging(str(snakemake.log[0]) if snakemake.log else None)
        lift_bsj_table(
            input_path=snakemake.input.table,
            chain_path=snakemake.input.chain,
            output_path=snakemake.output.table,
            rejected_path=snakemake.output.unmapped,
            liftover=snakemake.params.liftover,
        )
        return

    args = create_parser().parse_args()
    setup_logging(args.log)
    lift_bsj_table(
        input_path=args.input,
        chain_path=args.chain,
        output_path=args.output,
        rejected_path=args.rejected,
        liftover=args.liftover,
    )


if __name__ == "__main__":
    main()
