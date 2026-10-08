#!/usr/bin/env python3
"""Known coordinate exceptions for riboCIRC raw tables."""

from __future__ import annotations


# The literature table generally uses hg38_0based, but this row's circBase ID
# is supported by multiple databases as having an originally 1-based start.
# Convert it here to the 0-based start used consistently by the clean table.
RAW_COORDINATE_EXCEPTIONS = {
    (
        "literature_reported",
        "hsa_circ_0002301",
        "hg38",
        "chr14",
        31127785,
        31133675,
        "-",
    ): 31127784,
}


def normalize_hg38_0based_start(
    *,
    source_name: str,
    circbase_id: str | None,
    assembly: str | None,
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> tuple[int, bool]:
    """Return an hg38_0based start and whether a known exception was applied."""

    key = (
        source_name,
        (circbase_id or "").strip(),
        (assembly or "").strip(),
        chrom,
        start,
        end,
        strand,
    )
    corrected_start = RAW_COORDINATE_EXCEPTIONS.get(key)
    if corrected_start is None:
        return start, False
    return corrected_start, True
