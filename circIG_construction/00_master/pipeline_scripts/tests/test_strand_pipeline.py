#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_DIR / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {filename}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


presence_module = load_module(
    "add_strand_from_presence_master",
    "add_strand_from_presence_master.py",
)
audit_module = load_module(
    "write_final_strand_audit",
    "write_final_strand_audit.py",
)


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class PresenceMasterTests(unittest.TestCase):
    def test_resolution_conflict_missing_and_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            presence = root / "presence.tsv"
            ids = root / "ids.tsv"
            input_path = root / "input.tsv"
            output = root / "output.tsv"
            audit = root / "audit.tsv"

            fields = ["chrom", "start", "end", "strand", "#dbs", "DB1", "DB2"]
            write_tsv(
                presence,
                fields,
                [
                    {"chrom": "chr1", "start": "10", "end": "20", "strand": "+", "#dbs": "1", "DB1": "1", "DB2": "0"},
                    {"chrom": "chr1", "start": "30", "end": "40", "strand": "+", "#dbs": "1", "DB1": "1", "DB2": "0"},
                    {"chrom": "chr1", "start": "30", "end": "40", "strand": "-", "#dbs": "1", "DB1": "0", "DB2": "1"},
                    {"chrom": "chr1", "start": "50", "end": "60", "strand": "-", "#dbs": "1", "DB1": "0", "DB2": "1"},
                ],
            )
            write_tsv(
                ids,
                fields,
                [
                    {"chrom": "chr1", "start": "10", "end": "20", "strand": "+", "#dbs": "1", "DB1": "id1", "DB2": "0"},
                    {"chrom": "chr1", "start": "30", "end": "40", "strand": "+", "#dbs": "1", "DB1": "id2", "DB2": "0"},
                    {"chrom": "chr1", "start": "30", "end": "40", "strand": "-", "#dbs": "1", "DB1": "0", "DB2": "id3"},
                    {"chrom": "chr1", "start": "50", "end": "60", "strand": "-", "#dbs": "1", "DB1": "0", "DB2": "id4"},
                ],
            )
            write_tsv(
                input_path,
                ["circRNA_id", "chrom", "start", "end", "strand", "gene_symbol"],
                [
                    {"circRNA_id": "a", "chrom": "chr1", "start": "10", "end": "20", "strand": "NA", "gene_symbol": "A"},
                    {"circRNA_id": "b", "chrom": "chr1", "start": "30", "end": "40", "strand": "+", "gene_symbol": "B"},
                    {"circRNA_id": "c", "chrom": "chr1", "start": "70", "end": "80", "strand": "NA", "gene_symbol": "C"},
                ],
            )

            counts = presence_module.add_strands(
                input_path,
                output,
                presence,
                ids,
                audit,
                set(),
            )
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["strand_presence_master"] for row in rows], ["+", "NA", "NA"])
            self.assertEqual(counts["resolved_from_presence_master"], 1)
            self.assertEqual(counts["unresolved_conflicting_base_strands"], 1)
            self.assertEqual(counts["unresolved_no_base_hit"], 1)

            support = presence_module.load_presence_support(presence, ids, {"DB2"})
            inferred, status, _ = presence_module.resolve_strand(support, ("chr1", 50, 60))
            self.assertEqual((inferred, status), ("NA", "unresolved_no_base_hit"))


class FinalDecisionTests(unittest.TestCase):
    def test_gene_annotation_is_audit_only(self) -> None:
        self.assertEqual(
            audit_module.decide_strand(
                {
                    "strand": "NA",
                    "strand_gene_annotation": "+",
                    "strand_presence_master": "-",
                },
                "raw",
            ),
            ("-", "presence_master"),
        )
        self.assertEqual(
            audit_module.decide_strand(
                {
                    "strand": "NA",
                    "strand_gene_annotation": "+",
                    "strand_presence_master": "NA",
                },
                "raw",
            ),
            ("NA", "unresolved"),
        )
        self.assertEqual(
            audit_module.decide_strand(
                {
                    "strand": "+",
                    "strand_gene_annotation": "-",
                    "strand_presence_master": "-",
                },
                "raw",
            ),
            ("+", "raw"),
        )


if __name__ == "__main__":
    unittest.main()
