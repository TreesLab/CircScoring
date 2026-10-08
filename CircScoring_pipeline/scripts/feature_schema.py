import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path


RAW_COLUMNS = ["chr", "pos1", "pos2", "strand", "circ_id"]
REQUEST_COLUMNS = [*RAW_COLUMNS, "donor", "acceptor"]

FEATURE_COLUMNS = {
    "annotation": [
        "donor_site_at_the_annotated_boundary",
        "acceptor_site_at_the_annotated_boundary",
        "both_donor_acceptor_at_annotated_boundary",
        "donor_acceptor_sites_at_the_same_transcript_isoform",
        "the_number_of_transcripts_of_donor_acceptor_sites_at_the_same_transcript_isoform",
    ],
    "alternative_splicing": [
        "is_annotated_splicing_site(donor)",
        "has_AS_event(donor)",
        "is_annotated_splicing_site(acceptor)",
        "has_AS_event(acceptor)",
        "AS_event_type(donor)",
        "AS_event_type(acceptor)",
    ],
    "conservation": [
        "phyloP(acceptor_in)",
        "phyloP(acceptor_out)",
        "phyloP(donor_in)",
        "phyloP(donor_out)",
        "phastCons(acceptor_in)",
        "phastCons(acceptor_out)",
        "phastCons(donor_in)",
        "phastCons(donor_out)",
    ],
    "splicing": ["MAXENT(donor)", "MAXENT(acceptor)"],
}

PREDICTION_COLUMNS = [
    "XGBoost_CS-R_score",
    "XGBoost_CS-C_score",
    "ElasticNet_CS-R_score",
    "ElasticNet_CS-C_score",
]

REFERENCE_SUPPORT_COLUMNS = [
    "7 DB",
    "FL-circAS or circFL_seq",
    "mouse_conserved",
]

PREDICTCS_FEATURE_COLUMNS = [
    "circRNA_id",
    "donor_site_at_the_annotated_boundary",
    "acceptor_site_at_the_annotated_boundary",
    "donor_acceptor_sites_at_the_same_transcript_isoform",
    "has_AS_event(donor)",
    "has_AS_event(acceptor)",
    *FEATURE_COLUMNS["conservation"],
    *FEATURE_COLUMNS["splicing"],
    *REFERENCE_SUPPORT_COLUMNS,
]

DETAILED_FEATURE_COLUMNS = [
    "circRNA_id",
    *(column for group in FEATURE_COLUMNS.values() for column in group),
    *REFERENCE_SUPPORT_COLUMNS,
]

CHECKAA_COLUMNS = [
    "circRNA_id",
    "with an alternative co-linear explanation",
    "with multiple_hits",
    "alignment ambiguity (with an alternative co-linear explanation or multiple hits)",
]


def normalize_chromosome(value):
    chromosome = str(value).strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    if chromosome.upper() == "MT":
        chromosome = "M"
    if not chromosome:
        raise ValueError("Empty chromosome")
    return f"chr{chromosome}"


def parse_coordinate(chrom, pos1, pos2, strand, label):
    chrom = normalize_chromosome(chrom)
    strand = str(strand).strip()
    if strand not in {"+", "-"}:
        raise ValueError(f"Unsupported strand in {label}: {strand!r}")
    try:
        pos1 = int(pos1)
        pos2 = int(pos2)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid coordinates in {label}: {pos1!r}, {pos2!r}"
        ) from exc
    if pos1 < 1 or pos2 < 1 or pos1 > pos2:
        raise ValueError(f"Invalid coordinate range in {label}: {pos1}, {pos2}")
    return chrom, pos1, pos2, strand


def get_donor_acceptor(pos1, pos2, strand):
    if strand == "+":
        return pos2, pos1
    if strand == "-":
        return pos1, pos2
    raise ValueError(f"Unsupported strand: {strand!r}")


def circ_id_for_key(key):
    chrom, pos1, pos2, strand = key
    return f"{chrom}:{pos1}|{pos2}({strand})"


def event_id_for_key(key):
    encoded = json.dumps(key, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_text_writer(path, newline=""):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline=newline,
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    return handle, Path(handle.name), path


def finish_atomic_write(handle, temporary_path, output_path):
    handle.flush()
    os.fsync(handle.fileno())
    handle.close()
    os.replace(temporary_path, output_path)


def validate_header(reader, expected, label):
    if reader.fieldnames != expected:
        raise ValueError(
            f"Unexpected columns in {label}: expected {expected!r}, "
            f"got {reader.fieldnames!r}"
        )
