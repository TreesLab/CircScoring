#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "resolve_circrna_ids.py"
SPEC = importlib.util.spec_from_file_location("resolve_circrna_ids", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ResolveCircRNAIdsTests(unittest.TestCase):
    def test_unique_shared_unmatched_and_ambiguous_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            ids = root / "ids.tsv"
            ref_a = root / "ref_a.tsv"
            ref_b = root / "ref_b.tsv"
            output = root / "resolved.tsv"
            audit = root / "audit.tsv"
            summary = root / "summary.tsv"

            write_tsv(
                ids,
                ["circRNA_id"],
                [
                    {"circRNA_id": "id_unique"},
                    {"circRNA_id": "id_shared"},
                    {"circRNA_id": "id_ambiguous"},
                    {"circRNA_id": "id_missing"},
                    {"circRNA_id": "id_unique"},
                    {"circRNA_id": "NA"},
                ],
            )
            write_tsv(
                ref_a,
                MODULE.COLUMNS,
                [
                    {"circRNA_id": "id_unique", "chrom": "chr2", "start": "20", "end": "30", "strand": "+", "gene_symbol": "G2"},
                    {"circRNA_id": "id_shared", "chrom": "chr1", "start": "10", "end": "15", "strand": "-", "gene_symbol": "G1"},
                    {"circRNA_id": "id_ambiguous", "chrom": "chr3", "start": "40", "end": "50", "strand": "+", "gene_symbol": "G3"},
                ],
            )
            write_tsv(
                ref_b,
                MODULE.COLUMNS,
                [
                    {"circRNA_id": "id_shared", "chrom": "chr1", "start": "10", "end": "15", "strand": "-", "gene_symbol": "G1"},
                    {"circRNA_id": "id_ambiguous", "chrom": "chr4", "start": "60", "end": "70", "strand": "-", "gene_symbol": "G4"},
                ],
            )

            counts = MODULE.resolve_ids(
                ids,
                [MODULE.Reference("ref_a", ref_a), MODULE.Reference("ref_b", ref_b)],
                output,
                audit,
                summary,
                False,
            )

            self.assertEqual([row["circRNA_id"] for row in read_tsv(output)], ["id_shared", "id_unique"])
            audit_by_id = {row["circRNA_id"]: row for row in read_tsv(audit)}
            self.assertEqual(audit_by_id["id_missing"]["issue"], "unmatched")
            self.assertEqual(audit_by_id["id_ambiguous"]["issue"], "ambiguous_reference")
            self.assertEqual(audit_by_id["id_ambiguous"]["candidate_count"], "2")
            self.assertEqual(counts["duplicate_input_ids"], 1)
            self.assertEqual(counts["invalid_input_rows"], 1)
            self.assertEqual(counts["resolved_ids"], 2)

    def test_duplicate_reference_names_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            path = Path(temporary_dir) / "reference.tsv"
            write_tsv(path, MODULE.COLUMNS, [])
            with self.assertRaisesRegex(ValueError, "duplicate reference name"):
                MODULE.load_references([MODULE.Reference("same", path), MODULE.Reference("same", path)])


if __name__ == "__main__":
    unittest.main()
