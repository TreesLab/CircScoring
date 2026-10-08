#!/usr/bin/env python3
"""Clean CircTarget raw interactions and extract human circRNA IDs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW = BASE_DIR / "raw" / "all.txt"
DEFAULT_OUT_DIR = BASE_DIR / "clean"

ID_OUTPUT = "21_CircTarget.circRNA_ids.clean.tsv"
INTERACTION_OUTPUT = "21_CircTarget.interactions.human.clean.tsv"
SUMMARY_OUTPUT = "21_CircTarget.clean.summary.tsv"

INTERACTION_COLUMNS = [
    "circRNA_id",
    "target_ensembl_id",
    "target_gene_name",
    "target_gene_type",
    "chimeric_read_count",
    "p_value",
    "cell_line_tissue",
    "species",
    "interaction_type",
    "detected_method",
]

RAW_COLUMNS = {
    "circRNA ID",
    "Ensembl ID",
    "Gene name",
    "Gene type",
    "Chimeric read count",
    "P-value",
    "Cell line/Tissue",
    "Species",
    "Interaction type",
    "Detected method",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW, help=f"CircTarget raw CSV. Default: {DEFAULT_RAW}")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help=f"Clean output directory. Default: {DEFAULT_OUT_DIR}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def resolve(path: Path) -> Path:
    return path.expanduser().resolve()


def normalize(value: str | None) -> str:
    if value is None:
        return "NA"
    text = value.strip()
    return text if text else "NA"


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError("Output exists; use --force to overwrite: " + ", ".join(existing))


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def clean(raw_path: Path, output_dir: Path, force: bool) -> Counter[str]:
    id_output = output_dir / ID_OUTPUT
    interaction_output = output_dir / INTERACTION_OUTPUT
    summary_output = output_dir / SUMMARY_OUTPUT
    refuse_overwrite([id_output, interaction_output, summary_output], force)

    counts: Counter[str] = Counter()
    circ_ids: set[str] = set()
    interaction_rows: list[dict[str, str]] = []

    with raw_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = RAW_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{raw_path} missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            counts["raw_rows"] += 1
            species = normalize(row.get("Species"))
            counts[f"raw_species_rows:{species}"] += 1
            if species != "Human":
                continue
            counts["human_interaction_rows"] += 1

            circ_id = normalize(row.get("circRNA ID"))
            if circ_id == "NA":
                counts["human_rows_without_circRNA_id"] += 1
                continue

            circ_ids.add(circ_id)
            interaction_rows.append(
                {
                    "circRNA_id": circ_id,
                    "target_ensembl_id": normalize(row.get("Ensembl ID")),
                    "target_gene_name": normalize(row.get("Gene name")),
                    "target_gene_type": normalize(row.get("Gene type")),
                    "chimeric_read_count": normalize(row.get("Chimeric read count")),
                    "p_value": normalize(row.get("P-value")),
                    "cell_line_tissue": normalize(row.get("Cell line/Tissue")),
                    "species": species,
                    "interaction_type": normalize(row.get("Interaction type")),
                    "detected_method": normalize(row.get("Detected method")),
                }
            )

    counts["clean_human_interaction_rows"] = len(interaction_rows)
    counts["clean_unique_circRNA_ids"] = len(circ_ids)
    write_tsv(id_output, ["circRNA_id"], [{"circRNA_id": circ_id} for circ_id in sorted(circ_ids)])
    write_tsv(interaction_output, INTERACTION_COLUMNS, interaction_rows)
    write_tsv(
        summary_output,
        ["metric", "value"],
        [{"metric": key, "value": str(value)} for key, value in counts.items()],
    )
    return counts


def main() -> int:
    args = parse_args()
    counts = clean(resolve(args.raw), resolve(args.out_dir), args.force)
    for key, value in counts.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
