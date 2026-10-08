#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
DEFAULT_BASE_DATASET = ROOT / "02_circAtlas" / "human_browse_table_standardized.tsv"
DEFAULT_CHAIN = ROOT / "00_master" / "pipeline_resources" / "hg19ToHg38.over.chain.gz"
MATCHED_TYPES = ("hg38_1based", "hg38_0based", "hg19_1based", "hg19_0based")
APPENDED_COLUMNS = [
    "coordinate_classification",
    "matched_candidate_types",
    "best_candidate_type",
    "matched_base_ids",
    "matched_base_coordinates",
    "matched_base_database_counts",
    "matched_base_databases",
    "liftover_status",
    "notes",
]
UNIQUE_CLASSIFICATIONS = set(MATCHED_TYPES)


@dataclass(frozen=True)
class Coordinate:
    chrom: str
    start: int
    end: int
    strand: str = "."

    def key(self, use_strand: bool) -> tuple[str, int, int] | tuple[str, int, int, str]:
        if use_strand:
            return (self.chrom, self.start, self.end, self.strand)
        return (self.chrom, self.start, self.end)

    def label(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}:{self.strand}"


@dataclass(frozen=True)
class Candidate:
    name: str
    coord: Coordinate | None
    status: str = "ok"
    note: str = ""


@dataclass(frozen=True)
class BaseMatch:
    base_id: str
    coord: Coordinate
    databases: tuple[str, ...] = ()


def normalize_chrom(value: str) -> str:
    chrom = value.strip()
    if not chrom:
        return chrom
    if chrom.lower().startswith("chr"):
        return "chr" + chrom[3:]
    return f"chr{chrom}"


def normalize_strand(value: str | None) -> str:
    if value is None:
        return "."
    strand = value.strip()
    return strand if strand in {"+", "-", "."} else "."


def parse_int(value: str, column: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Column {column!r} contains a non-integer coordinate: {value!r}") from exc


def read_table(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Input table has no header: {path}")
        return list(reader.fieldnames), rows


def load_base_index(
    path: Path,
    delimiter: str,
    chrom_col: str,
    start_col: str,
    end_col: str,
    strand_col: str | None,
    id_col: str | None,
) -> tuple[dict[tuple, list[BaseMatch]], bool, bool]:
    header, rows = read_table(path, delimiter)
    use_strand = strand_col is not None
    database_columns = header[header.index("#dbs") + 1 :] if "#dbs" in header else []
    has_database_annotations = bool(database_columns) or "databases" in header
    index: dict[tuple, list[BaseMatch]] = defaultdict(list)

    for row_number, row in enumerate(rows, start=2):
        try:
            coord = Coordinate(
                chrom=normalize_chrom(row[chrom_col]),
                start=parse_int(row[start_col], start_col),
                end=parse_int(row[end_col], end_col),
                strand=normalize_strand(row.get(strand_col) if strand_col else None),
            )
        except KeyError as exc:
            raise KeyError(f"Missing required base dataset column {exc.args[0]!r}") from exc
        except ValueError as exc:
            raise ValueError(f"{path}:{row_number}: {exc}") from exc

        if coord.start > coord.end:
            raise ValueError(f"{path}:{row_number}: start is greater than end: {coord.label()}")

        base_id = row.get(id_col, "") if id_col else ""
        if not base_id:
            base_id = f"base_row_{row_number}"
        if "databases" in header:
            databases = tuple(sorted({item.strip() for item in row.get("databases", "").split(",") if item.strip() and item.strip() != "NA"}))
        else:
            databases = tuple(sorted(column for column in database_columns if row.get(column, "").strip() not in {"", "0", "NA", "nan"}))
        match = BaseMatch(base_id, coord, databases)
        if match not in index[coord.key(use_strand)]:
            index[coord.key(use_strand)].append(match)

    return dict(index), use_strand, has_database_annotations


def make_hg38_candidates(coord: Coordinate) -> list[Candidate]:
    candidates = [
        Candidate("hg38_1based", coord),
        Candidate(
            "hg38_0based",
            Coordinate(coord.chrom, coord.start + 1, coord.end, coord.strand),
        ),
    ]
    return [candidate for candidate in candidates if candidate.coord and candidate.coord.start <= candidate.coord.end]


def boundary_bed_intervals(coord: Coordinate, coordinate_system: str) -> tuple[tuple[int, int], tuple[int, int]]:
    if coordinate_system == "1based":
        start0 = coord.start - 1
        end0 = coord.end
    elif coordinate_system == "0based":
        start0 = coord.start
        end0 = coord.end
    else:
        raise ValueError(f"Unsupported coordinate system: {coordinate_system}")

    if start0 < 0 or start0 >= end0:
        raise ValueError(f"Invalid {coordinate_system} interval for liftOver: {coord.label()}")

    return (start0, start0 + 1), (end0 - 1, end0)


def write_liftover_input(rows: list[tuple[str, Coordinate]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        for row_key, coord in rows:
            for system in ("1based", "0based"):
                try:
                    start_interval, end_interval = boundary_bed_intervals(coord, system)
                except ValueError:
                    continue
                for boundary_name, interval in (("start", start_interval), ("end", end_interval)):
                    handle.write(
                        "\t".join(
                            [
                                coord.chrom,
                                str(interval[0]),
                                str(interval[1]),
                                f"{row_key}|hg19_{system}|{boundary_name}",
                                "0",
                                coord.strand,
                            ]
                        )
                        + "\n"
                    )


def parse_liftover_output(path: Path) -> dict[str, tuple[str, int, int, str]]:
    mapped: dict[str, tuple[str, int, int, str]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start, end, name, _score, strand = fields[:6]
            mapped[name] = (normalize_chrom(chrom), int(start), int(end), normalize_strand(strand))
    return mapped


def run_liftover(
    row_coords: list[tuple[str, Coordinate]],
    liftover_bin: str,
    chain: Path,
    workdir: Path,
) -> dict[str, Candidate]:
    candidates: dict[str, Candidate] = {}
    if not row_coords:
        return candidates

    if shutil.which(liftover_bin) is None:
        note = f"liftOver executable not found: {liftover_bin}"
        for row_key, _coord in row_coords:
            for name in ("hg19_1based", "hg19_0based"):
                candidates[f"{row_key}|{name}"] = Candidate(name, None, "not_evaluated", note)
        return candidates

    if not chain.exists():
        note = f"chain file not found: {chain}"
        for row_key, _coord in row_coords:
            for name in ("hg19_1based", "hg19_0based"):
                candidates[f"{row_key}|{name}"] = Candidate(name, None, "not_evaluated", note)
        return candidates

    workdir.mkdir(parents=True, exist_ok=True)
    input_bed = workdir / "circrna_coordinate_check_hg19_candidates.bed"
    mapped_bed = workdir / "circrna_coordinate_check_hg19_candidates.hg38.bed"
    unmapped_bed = workdir / "circrna_coordinate_check_hg19_candidates.unmapped.bed"
    write_liftover_input(row_coords, input_bed)

    try:
        subprocess.run(
            [liftover_bin, str(input_bed), str(chain), str(mapped_bed), str(unmapped_bed)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        note = "liftOver failed"
        if exc.stderr:
            note = f"{note}: {exc.stderr.strip()}"
        for row_key, _coord in row_coords:
            for name in ("hg19_1based", "hg19_0based"):
                candidates[f"{row_key}|{name}"] = Candidate(name, None, "failed", note)
        return candidates

    mapped = parse_liftover_output(mapped_bed)
    for row_key, _coord in row_coords:
        for system in ("1based", "0based"):
            name = f"hg19_{system}"
            start_key = f"{row_key}|{name}|start"
            end_key = f"{row_key}|{name}|end"
            candidate_key = f"{row_key}|{name}"

            if start_key not in mapped or end_key not in mapped:
                candidates[candidate_key] = Candidate(name, None, "unmapped", "one or both BSJ boundaries did not liftOver")
                continue

            start_chrom, start_start, start_end, start_strand = mapped[start_key]
            end_chrom, end_start, end_end, end_strand = mapped[end_key]
            if start_chrom != end_chrom:
                candidates[candidate_key] = Candidate(name, None, "failed", "lifted BSJ boundaries mapped to different chromosomes")
                continue
            if start_strand != end_strand:
                candidates[candidate_key] = Candidate(name, None, "failed", "lifted BSJ boundaries mapped to different strands")
                continue

            lifted_start = min(start_start, end_start) + 1
            lifted_end = max(start_end, end_end)
            if lifted_start > lifted_end:
                candidates[candidate_key] = Candidate(name, None, "failed", "invalid interval after liftOver rebuild")
                continue

            candidates[candidate_key] = Candidate(
                name,
                Coordinate(start_chrom, lifted_start, lifted_end, start_strand),
            )

    return candidates


def classify_matches(matches: dict[str, list[BaseMatch]]) -> tuple[str, str, str, str, str]:
    matched_types = [name for name in MATCHED_TYPES if matches.get(name)]
    if not matched_types:
        return "unmatched", "", "", "", ""
    if len(matched_types) == 1:
        classification = matched_types[0]
        best = matched_types[0]
    else:
        classification = "ambiguous"
        best = ""

    matched_ids: list[str] = []
    matched_coords: list[str] = []
    for name in matched_types:
        for match in matches[name]:
            matched_ids.append(f"{name}:{match.base_id}")
            matched_coords.append(f"{name}:{match.coord.label()}")

    return classification, ",".join(matched_types), best, ";".join(matched_ids), ";".join(matched_coords)


def format_database_support(matches: dict[str, list[BaseMatch]]) -> tuple[str, str]:
    counts: list[str] = []
    databases: list[str] = []
    for name in MATCHED_TYPES:
        matched_databases = sorted({database for match in matches.get(name, []) for database in match.databases})
        if matched_databases:
            counts.append(f"{name}:{len(matched_databases)}")
            databases.append(f"{name}:{','.join(matched_databases)}")
    return ";".join(counts) or "NA", ";".join(databases) or "NA"


def append_results(
    header: list[str],
    rows: list[dict[str, str]],
    input_coords: list[tuple[str, Coordinate]],
    base_index: dict[tuple, list[BaseMatch]],
    use_strand: bool,
    hg19_candidates: dict[str, Candidate],
    has_database_annotations: bool,
) -> list[dict[str, str]]:
    coord_by_key = dict(input_coords)
    output_rows: list[dict[str, str]] = []

    for row_index, row in enumerate(rows):
        row_key = str(row_index)
        coord = coord_by_key.get(row_key)
        matches: dict[str, list[BaseMatch]] = {}
        notes: list[str] = []
        liftover_statuses: list[str] = []

        if coord is None:
            notes.append("invalid input coordinate")
        else:
            for candidate in make_hg38_candidates(coord):
                if candidate.coord is None:
                    continue
                labels = base_index.get(candidate.coord.key(use_strand), [])
                if labels:
                    matches[candidate.name] = labels

            for name in ("hg19_1based", "hg19_0based"):
                candidate = hg19_candidates.get(f"{row_key}|{name}", Candidate(name, None, "not_evaluated", "missing liftOver candidate"))
                liftover_statuses.append(f"{name}:{candidate.status}")
                if candidate.note:
                    notes.append(f"{name}:{candidate.note}")
                if candidate.coord is None:
                    continue
                labels = base_index.get(candidate.coord.key(use_strand), [])
                if labels:
                    matches[candidate.name] = labels

        classification, matched_types, best, matched_ids, matched_coords = classify_matches(matches)
        database_counts, matched_databases = format_database_support(matches)
        if not has_database_annotations:
            database_counts = "NA"
            matched_databases = "NA"
        output_row = dict(row)
        output_row.update(
            {
                "coordinate_classification": classification,
                "matched_candidate_types": matched_types,
                "best_candidate_type": best,
                "matched_base_ids": matched_ids,
                "matched_base_coordinates": matched_coords,
                "matched_base_database_counts": database_counts,
                "matched_base_databases": matched_databases,
                "liftover_status": ";".join(liftover_statuses),
                "notes": ";".join(dict.fromkeys(notes)),
            }
        )
        output_rows.append(output_row)

    return output_rows


def write_table(path: Path, header: list[str], rows: list[dict[str, str]], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    output_header = list(header) + [column for column in APPENDED_COLUMNS if column not in header]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_header, delimiter=delimiter, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_matched_coordinates(value: str) -> dict[str, set[str]]:
    coords_by_type: dict[str, set[str]] = defaultdict(set)
    if not value:
        return {}

    for item in value.split(";"):
        if not item:
            continue
        candidate_type, separator, coord = item.partition(":")
        if not separator or not coord:
            continue
        coords_by_type[candidate_type].add(coord)

    return dict(coords_by_type)


def effective_summary_classification(row: dict[str, str]) -> tuple[str, str]:
    classification = row.get("coordinate_classification", "")
    if classification in UNIQUE_CLASSIFICATIONS:
        return classification, ""

    matched_types = row.get("matched_candidate_types", "")
    if matched_types == "hg38_1based,hg19_1based":
        coords_by_type = split_matched_coordinates(row.get("matched_base_coordinates", ""))
        hg38_coords = coords_by_type.get("hg38_1based", set())
        hg19_coords = coords_by_type.get("hg19_1based", set())
        if hg38_coords and hg38_coords == hg19_coords:
            return "hg38_1based", "liftover_equivalent_ambiguous"

    return "", ""


def percentage(count: int, total: int) -> str:
    if total == 0:
        return "0.00%"
    return f"{count / total * 100:.2f}%"


def write_count_section(handle, title: str, counts: Counter[str], total: int) -> None:
    handle.write(f"\n{title}\n")
    if not counts:
        handle.write("  (none)\n")
        return
    for key, count in counts.most_common():
        label = key if key else "(blank)"
        handle.write(f"  {label}\t{count}\t{percentage(count, total)}\n")


def dominant_classification(effective_counts: Counter[str]) -> tuple[str, int, int]:
    matched_total = sum(effective_counts.values())
    if not effective_counts:
        return "unknown", 0, matched_total
    dominant, count = effective_counts.most_common(1)[0]
    return dominant, count, matched_total


def decide_overall_judgment(effective_counts: Counter[str], total_rows: int, judgment_mode: str) -> str:
    if judgment_mode == "strict":
        unique_effective = set(effective_counts)
        if len(unique_effective) == 1:
            return next(iter(unique_effective))
        if len(unique_effective) > 1:
            return "mixed"
        return "no_unique_coordinate_system"

    dominant, dominant_count, matched_total = dominant_classification(effective_counts)
    if matched_total == 0:
        return "no_unique_coordinate_system"
    if len(effective_counts) == 1:
        return dominant
    if dominant_count / matched_total >= 0.995:
        return dominant
    if dominant_count / matched_total >= 0.95 and dominant_count / max(total_rows, 1) >= 0.05:
        return dominant
    return "mixed"


def summarize_output_rows(
    rows: list[dict[str, str]],
    input_path: Path,
    base_dataset: Path,
    output_path: Path,
    summary_path: Path | None = None,
    judgment_mode: str = "loose",
) -> dict[str, str]:
    total_rows = len(rows)
    raw_classification_counts: Counter[str] = Counter()
    effective_classification_counts: Counter[str] = Counter()
    matched_candidate_counts: Counter[str] = Counter()
    liftover_status_counts: Counter[str] = Counter()
    unresolved_ambiguous_count = 0
    liftover_equivalent_ambiguous_count = 0

    for row in rows:
        raw_classification = row.get("coordinate_classification", "")
        raw_classification_counts[raw_classification] += 1
        matched_candidate_counts[row.get("matched_candidate_types", "") or "(none)"] += 1
        liftover_status_counts[row.get("liftover_status", "") or "(blank)"] += 1

        effective_classification, effective_note = effective_summary_classification(row)
        if effective_classification:
            effective_classification_counts[effective_classification] += 1
        if effective_note == "liftover_equivalent_ambiguous":
            liftover_equivalent_ambiguous_count += 1
        elif raw_classification == "ambiguous":
            unresolved_ambiguous_count += 1

    overall_judgment = decide_overall_judgment(effective_classification_counts, total_rows, judgment_mode)
    dominant, dominant_count, matched_total = dominant_classification(effective_classification_counts)
    dominant_matched_percent = percentage(dominant_count, matched_total)

    summary = {
        "input_file": str(input_path),
        "base_dataset": str(base_dataset),
        "row_level_output": str(output_path),
        "total_rows": str(total_rows),
        "judgment_mode": judgment_mode,
        "overall_judgment": overall_judgment,
        "dominant_effective_classification": dominant,
        "dominant_effective_count": str(dominant_count),
        "matched_total": str(matched_total),
        "dominant_matched_percent": dominant_matched_percent,
        "unique_classified_rows_for_summary": str(sum(effective_classification_counts.values())),
        "hg38_1based_count": str(effective_classification_counts.get("hg38_1based", 0)),
        "hg38_0based_count": str(effective_classification_counts.get("hg38_0based", 0)),
        "hg19_1based_count": str(effective_classification_counts.get("hg19_1based", 0)),
        "hg19_0based_count": str(effective_classification_counts.get("hg19_0based", 0)),
        "ambiguous_count": str(raw_classification_counts.get("ambiguous", 0)),
        "unresolved_ambiguous_count": str(unresolved_ambiguous_count),
        "liftover_equivalent_ambiguous_count": str(liftover_equivalent_ambiguous_count),
        "unmatched_count": str(raw_classification_counts.get("unmatched", 0)),
        "invalid_input_count": str(raw_classification_counts.get("invalid_input", 0)),
    }

    if summary_path:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with summary_path.open("w", newline="") as handle:
            handle.write("circRNA coordinate system check summary\n")
            handle.write(f"Input table: {input_path}\n")
            handle.write(f"Base dataset: {base_dataset}\n")
            handle.write(f"Row-level output: {output_path}\n")
            handle.write(f"Total rows: {total_rows}\n")
            handle.write(f"Judgment mode: {judgment_mode}\n")
            handle.write(f"Overall judgment: {overall_judgment}\n")
            handle.write(f"Dominant effective classification: {dominant}\n")
            handle.write(f"Dominant effective count: {dominant_count}\n")
            handle.write(f"Matched total: {matched_total}\n")
            handle.write(f"Dominant matched percent: {dominant_matched_percent}\n")
            handle.write(f"Unique-classified rows for summary: {sum(effective_classification_counts.values())}\n")
            handle.write(f"LiftOver-equivalent ambiguous rows: {liftover_equivalent_ambiguous_count}\n")
            handle.write(f"Unresolved ambiguous rows: {unresolved_ambiguous_count}\n")
            handle.write(f"Unmatched rows: {raw_classification_counts.get('unmatched', 0)}\n")
            handle.write(f"Invalid input rows: {raw_classification_counts.get('invalid_input', 0)}\n")

            write_count_section(handle, "Raw coordinate_classification counts", raw_classification_counts, total_rows)
            write_count_section(handle, "Effective classification counts for summary", effective_classification_counts, total_rows)
            write_count_section(handle, "Matched candidate type counts", matched_candidate_counts, total_rows)
            write_count_section(handle, "LiftOver status counts", liftover_status_counts, total_rows)

            handle.write("\nRule notes\n")
            if judgment_mode == "strict":
                handle.write("  Overall judgment uses strict unique effective classification.\n")
            else:
                handle.write("  Overall judgment uses loose dominant effective classification.\n")
                handle.write("  Loose mode accepts a dominant class if it is >=99.5% of matched rows, or >=95% of matched rows and >=5% of total rows.\n")
            handle.write("  Unmatched, invalid_input, and unresolved ambiguous rows are excluded from the overall coordinate-system decision.\n")
            handle.write("  matched_candidate_types=hg38_1based,hg19_1based with identical matched coordinates is treated as hg38_1based for summary only.\n")

    return summary


def write_summary_tsv(path: Path, summary: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "input_file",
        "base_dataset",
        "row_level_output",
        "total_rows",
        "judgment_mode",
        "overall_judgment",
        "dominant_effective_classification",
        "dominant_effective_count",
        "matched_total",
        "dominant_matched_percent",
        "unique_classified_rows_for_summary",
        "hg38_1based_count",
        "hg38_0based_count",
        "hg19_1based_count",
        "hg19_0based_count",
        "ambiguous_count",
        "unresolved_ambiguous_count",
        "liftover_equivalent_ambiguous_count",
        "unmatched_count",
        "invalid_input_count",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(summary)


def parse_input_coordinates(
    rows: list[dict[str, str]],
    chrom_col: str,
    start_col: str,
    end_col: str,
    strand_col: str | None,
) -> tuple[list[tuple[str, Coordinate]], dict[int, str]]:
    coords: list[tuple[str, Coordinate]] = []
    errors: dict[int, str] = {}
    for row_index, row in enumerate(rows):
        try:
            coord = Coordinate(
                chrom=normalize_chrom(row[chrom_col]),
                start=parse_int(row[start_col], start_col),
                end=parse_int(row[end_col], end_col),
                strand=normalize_strand(row.get(strand_col) if strand_col else None),
            )
            if coord.start > coord.end:
                raise ValueError(f"start is greater than end: {coord.label()}")
            coords.append((str(row_index), coord))
        except (KeyError, ValueError) as exc:
            errors[row_index] = str(exc)
    return coords, errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Classify circRNA coordinates as hg38/hg19 and 1-based/0-based by matching to a known hg38 1-based base dataset.",
    )
    parser.add_argument("--input", required=True, type=Path, help="Headered input circRNA table to check.")
    parser.add_argument(
        "--base-dataset",
        type=Path,
        default=DEFAULT_BASE_DATASET,
        help=f"Known hg38 1-based closed circRNA table. Default: {DEFAULT_BASE_DATASET}",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output table with appended classification columns.")
    parser.add_argument("--summary-output", type=Path, help="Optional TXT summary output path.")
    parser.add_argument("--summary-tsv-output", type=Path, help="Optional one-row TSV summary output path for pipeline use.")
    parser.add_argument("--chrom-col", default="chrom", help="Input chromosome column name. Default: chrom.")
    parser.add_argument("--start-col", default="start", help="Input start coordinate column name. Default: start.")
    parser.add_argument("--end-col", default="end", help="Input end coordinate column name. Default: end.")
    parser.add_argument(
        "--strand-col",
        default="strand",
        help="Input strand column name. Use empty string to ignore strand. Default: strand.",
    )
    parser.add_argument("--id-col", default="circRNA_id", help="Input circRNA ID column name. Preserved only; matching does not require it. Default: circRNA_id.")
    parser.add_argument("--base-chrom-col", default="chrom", help="Base dataset chromosome column name.")
    parser.add_argument("--base-start-col", default="start", help="Base dataset 1-based start column name.")
    parser.add_argument("--base-end-col", default="end", help="Base dataset 1-based end column name.")
    parser.add_argument("--base-strand-col", default="strand", help="Base dataset strand column name. Use empty string to ignore strand.")
    parser.add_argument("--base-id-col", default="circRNA_id", help="Base dataset ID column name. Use empty string to synthesize row IDs.")
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN, help="hg19ToHg38 liftOver chain file.")
    parser.add_argument("--liftover-bin", default="liftOver", help="liftOver executable.")
    parser.add_argument("--workdir", type=Path, help="Directory for temporary liftOver BED files.")
    parser.add_argument("--delimiter", default="\t", help="Input/output delimiter. Default: tab.")
    parser.add_argument("--base-delimiter", default="\t", help="Base dataset delimiter. Default: tab.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Use strict overall judgment: report mixed whenever more than one effective coordinate class is present. Default uses loose dominant-class judgment.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    strand_col = args.strand_col or None
    base_strand_col = args.base_strand_col or None
    base_id_col = args.base_id_col or None

    base_index, use_strand, has_database_annotations = load_base_index(
        args.base_dataset,
        args.base_delimiter,
        args.base_chrom_col,
        args.base_start_col,
        args.base_end_col,
        base_strand_col,
        base_id_col,
    )
    header, rows = read_table(args.input, args.delimiter)
    input_coords, input_errors = parse_input_coordinates(rows, args.chrom_col, args.start_col, args.end_col, strand_col)

    if args.workdir:
        hg19_candidates = run_liftover(input_coords, args.liftover_bin, args.chain, args.workdir)
    else:
        with tempfile.TemporaryDirectory(prefix="circrna_coordinate_check_") as tmpdir:
            hg19_candidates = run_liftover(input_coords, args.liftover_bin, args.chain, Path(tmpdir))

    output_rows = append_results(header, rows, input_coords, base_index, use_strand, hg19_candidates, has_database_annotations)
    for row_index, error in input_errors.items():
        output_rows[row_index].update(
            {
                "coordinate_classification": "invalid_input",
                "matched_candidate_types": "",
                "best_candidate_type": "",
                "matched_base_ids": "",
                "matched_base_coordinates": "",
                "matched_base_database_counts": "",
                "matched_base_databases": "",
                "liftover_status": "",
                "notes": error,
            }
        )

    write_table(args.output, header, output_rows, args.delimiter)
    summary: dict[str, str] | None = None
    if args.summary_output:
        summary = summarize_output_rows(
            output_rows,
            args.input,
            args.base_dataset,
            args.output,
            args.summary_output,
            "strict" if args.strict else "loose",
        )
    if args.summary_tsv_output:
        if summary is None:
            summary = summarize_output_rows(
                output_rows,
                args.input,
                args.base_dataset,
                args.output,
                judgment_mode="strict" if args.strict else "loose",
            )
        write_summary_tsv(args.summary_tsv_output, summary)


if __name__ == "__main__":
    main()
