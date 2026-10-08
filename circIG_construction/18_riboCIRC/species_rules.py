#!/usr/bin/env python3
"""Species classification rules for riboCIRC raw tables."""

from __future__ import annotations

from collections.abc import Mapping


MISSING_VALUES = {"", "-", "NA", "N/A", "NULL", "NONE", "NAN"}
HUMAN_ID_FIELDS = ("riboCIRC_ID", "circBase_ID", "circRNA_name")


def _is_human_id(value: str | None) -> bool:
    text = (value or "").strip()
    if text.upper() in MISSING_VALUES:
        return False
    lowered = text.lower()
    return lowered.startswith(("hsa_", "hsa-"))


def classify_literature_human(
    row: Mapping[str, str | None],
) -> tuple[bool, str | None]:
    """Determine whether a literature row provides human circRNA coordinates."""

    species = (row.get("Species") or "").strip()
    if species == "Homo sapiens" or species.startswith("Homo sapiens/"):
        return True, "species"

    assembly = (row.get("Genome_assembly") or "").strip()
    has_human_id = any(
        _is_human_id(row.get(field)) for field in HUMAN_ID_FIELDS
    )
    if assembly == "hg38" and has_human_id:
        return True, "hg38_hsa_id"

    return False, None
