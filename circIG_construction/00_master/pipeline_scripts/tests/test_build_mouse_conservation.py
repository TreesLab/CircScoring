from __future__ import annotations

import csv
import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "build_mouse_conservation.py"
SPEC = importlib.util.spec_from_file_location("build_mouse_conservation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_tsv(path: Path, columns: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


class MouseConservationTests(unittest.TestCase):
    def test_preserves_na_strand_and_collapses_mouse_ids(self) -> None:
        with self.subTest("temporary tables"):
            import tempfile

            with tempfile.TemporaryDirectory() as tmpdir:
                root = Path(tmpdir)
                human = root / "human.tsv"
                mouse = root / "mouse.tsv"
                output = root / "output.tsv"
                summary = root / "summary.tsv"
                write_tsv(
                    human,
                    ["chrom", "start", "end", "strand", "#dbs"],
                    [["chr1", "10", "20", "NA", "1"], ["chr2", "30", "40", "+", "2"]],
                )
                write_tsv(
                    mouse,
                    [
                        "chrom", "start", "end", "strand", "original_chrom", "original_start",
                        "original_end", "original_strand",
                    ],
                    [
                        ["chr1", "10", "20", "nan", "chr3", "100", "110", "NA"],
                        ["chr1", "10", "20", "NA", "chr4", "200", "210", "-"],
                        ["chr1", "10", "20", "NA", "chr3", "100", "110", "NA"],
                    ],
                )

                counts = MODULE.build(human, mouse, output, summary)

                with output.open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(rows[0]["circ_id"], "chr1:10|20(NA)")
                self.assertEqual(
                    rows[0]["mouse_circ_id"],
                    "chr3:100|110(NA),chr4:200|210(-)",
                )
                self.assertEqual(rows[1]["mouse_circ_id"], "")
                self.assertEqual(counts["mouse_grouped_circ_ids"], 1)
                self.assertEqual(counts["conserved_human_rows"], 1)
