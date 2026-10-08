#!/usr/bin/env python3
"""Download the deepBase human circRNA browser coordinate table."""

from __future__ import annotations

import argparse
import csv
import html
import json
import plistlib
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path
from typing import Any


BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / "raw"
DEFAULT_ENDPOINT = "https://rna.sysu.edu.cn/deepbase3/subpages/server_processing_circRNA_browse_deepBase.php"
DEFAULT_BROWSE_URL = "https://rna.sysu.edu.cn/deepbase3/subpages/browse_circRNA.php?SClade=mammal&SOrganism=hg19"
DEFAULT_OUTPUT = RAW_DIR / "deepbase_human_circRNA_browser_api.tsv"
DEFAULT_SUMMARY_OUTPUT = RAW_DIR / "deepbase_human_circRNA_browser_api.download_summary.tsv"
DEFAULT_WEBARCHIVE = RAW_DIR / "DeepBase.webarchive"
EXPECTED_COUNT = 918
EXPECTED_HEADER = [
    "Name",
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Source Gene ID",
    "Source Gene Symbol",
]


class TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self._in_table = False
        self._in_cell = False
        self._current_table: list[list[str]] = []
        self._current_row: list[str] = []
        self._cell_parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table":
            self._in_table = True
            self._current_table = []
        elif self._in_table and tag == "tr":
            self._current_row = []
        elif self._in_table and tag in {"td", "th"}:
            self._in_cell = True
            self._cell_parts = []

    def handle_endtag(self, tag: str) -> None:
        if self._in_table and tag in {"td", "th"} and self._in_cell:
            self._current_row.append(normalize_text("".join(self._cell_parts)))
            self._in_cell = False
        elif self._in_table and tag == "tr":
            if any(cell for cell in self._current_row):
                self._current_table.append(self._current_row)
        elif self._in_table and tag == "table":
            self.tables.append(self._current_table)
            self._in_table = False

    def handle_data(self, data: str) -> None:
        if self._in_cell:
            self._cell_parts.append(data)


def normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", html.unescape(value)).strip()


def normalize_missing(value: Any) -> str:
    if value is None:
        return "NA"
    text = normalize_text(str(value))
    return text if text else "NA"


def request_json(url: str, timeout: int) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "data-downloader/1.0",
            "Accept": "application/json,text/javascript,*/*;q=0.8",
            "Referer": DEFAULT_BROWSE_URL,
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        body = response.read().decode("utf-8", errors="replace")
    try:
        parsed = json.loads(body)
    except json.JSONDecodeError as exc:
        raise ValueError(f"API did not return JSON: {body[:200]!r}") from exc
    if not isinstance(parsed, dict):
        raise ValueError("API JSON root is not an object")
    return parsed


def api_url(endpoint: str, start: int, length: int) -> str:
    params = {
        "db": "hg19",
        "draw": "1",
        "start": str(start),
        "length": str(length),
    }
    return f"{endpoint}?{urllib.parse.urlencode(params)}"


def parse_api_record(record: Any) -> dict[str, str]:
    if isinstance(record, dict):
        values = [record.get(column, record.get(str(index), "")) for index, column in enumerate(EXPECTED_HEADER)]
    elif isinstance(record, list):
        values = record[: len(EXPECTED_HEADER)]
    else:
        raise ValueError(f"unsupported API record type: {type(record).__name__}")
    if len(values) < len(EXPECTED_HEADER):
        raise ValueError(f"API record has too few columns: {values!r}")
    return {column: normalize_missing(value) for column, value in zip(EXPECTED_HEADER, values)}


def fetch_api_rows(endpoint: str, expected_count: int, timeout: int) -> tuple[list[dict[str, str]], dict[str, str]]:
    payload = request_json(api_url(endpoint, 0, expected_count), timeout)
    raw_rows = payload.get("data", payload.get("aaData", []))
    if not isinstance(raw_rows, list):
        raise ValueError("API JSON does not contain a list in data/aaData")
    rows = [parse_api_record(record) for record in raw_rows[:expected_count]]
    metadata = {
        "source": "api",
        "endpoint": endpoint,
        "requested_start": "0",
        "requested_length": str(expected_count),
        "records_total": str(payload.get("recordsTotal", payload.get("iTotalRecords", "NA"))),
        "records_filtered": str(payload.get("recordsFiltered", payload.get("iTotalDisplayRecords", "NA"))),
        "downloaded_rows": str(len(rows)),
        "fallback_reason": "NA",
    }
    return rows, metadata


def read_webarchive(path: Path) -> tuple[str, str]:
    with path.open("rb") as handle:
        archive = plistlib.load(handle)
    main_resource = archive.get("WebMainResource")
    if not isinstance(main_resource, dict):
        raise ValueError(f"{path} does not contain a WebMainResource")
    data = main_resource.get("WebResourceData")
    if not isinstance(data, bytes):
        raise ValueError(f"{path} WebMainResource does not contain byte data")
    encoding = main_resource.get("WebResourceTextEncodingName") or "utf-8"
    source_url = str(main_resource.get("WebResourceURL") or path)
    return data.decode(encoding, errors="replace"), source_url


def find_circrna_table(markup: str) -> list[list[str]]:
    parser = TableParser()
    parser.feed(markup)
    for table in parser.tables:
        for index, row in enumerate(table):
            if row[: len(EXPECTED_HEADER)] == EXPECTED_HEADER:
                return table[index:]
    raise ValueError("could not find the deepBase circRNA browser table")


def rows_from_table(table: list[list[str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for raw in table[1:]:
        if len(raw) < len(EXPECTED_HEADER):
            continue
        row = {column: normalize_missing(value) for column, value in zip(EXPECTED_HEADER, raw[: len(EXPECTED_HEADER)])}
        if row["Name"] != "NA" and row["Chromosome"] != "NA" and row["Start"] != "NA" and row["End"] != "NA":
            rows.append(row)
    return rows


def fetch_webarchive_rows(path: Path, fallback_reason: str, expected_count: int) -> tuple[list[dict[str, str]], dict[str, str]]:
    markup, source_url = read_webarchive(path)
    rows = rows_from_table(find_circrna_table(markup))[:expected_count]
    metadata = {
        "source": "webarchive",
        "endpoint": source_url,
        "requested_start": "0",
        "requested_length": str(expected_count),
        "records_total": "NA",
        "records_filtered": "NA",
        "downloaded_rows": str(len(rows)),
        "fallback_reason": fallback_reason,
    }
    return rows, metadata


def validate_rows(rows: list[dict[str, str]], expected_count: int) -> None:
    if len(rows) != expected_count:
        raise ValueError(f"expected {expected_count} rows, found {len(rows)}")
    seen: set[str] = set()
    for row in rows:
        circ_id = row["Name"]
        if circ_id in seen:
            raise ValueError(f"duplicate circRNA ID: {circ_id}")
        seen.add(circ_id)
        if not re.fullmatch(r"chr[\w.]+", row["Chromosome"]):
            raise ValueError(f"invalid chromosome for {circ_id}: {row['Chromosome']}")
        if not row["Start"].isdigit() or not row["End"].isdigit():
            raise ValueError(f"invalid coordinates for {circ_id}: {row['Start']}-{row['End']}")
        if int(row["Start"]) > int(row["End"]):
            raise ValueError(f"start > end for {circ_id}: {row['Start']}>{row['End']}")
        if row["Strand"] not in {"+", "-"}:
            raise ValueError(f"invalid strand for {circ_id}: {row['Strand']}")


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPECTED_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, metadata: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in [
            "source",
            "endpoint",
            "requested_start",
            "requested_length",
            "records_total",
            "records_filtered",
            "downloaded_rows",
            "fallback_reason",
        ]:
            writer.writerow([key, metadata.get(key, "NA")])


def refuse_overwrite(paths: list[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(f"Output exists; use --force to overwrite: {', '.join(existing)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT, help=f"DataTables API endpoint. Default: {DEFAULT_ENDPOINT}")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"Raw output TSV. Default: {DEFAULT_OUTPUT}")
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help=f"Download summary TSV. Default: {DEFAULT_SUMMARY_OUTPUT}",
    )
    parser.add_argument("--expected-count", type=int, default=EXPECTED_COUNT, help=f"Expected row count. Default: {EXPECTED_COUNT}")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds. Default: 60.")
    parser.add_argument(
        "--webarchive-fallback",
        type=Path,
        default=DEFAULT_WEBARCHIVE,
        help=f"Fallback Safari webarchive. Default: {DEFAULT_WEBARCHIVE}",
    )
    parser.add_argument("--no-fallback", action="store_true", help="Fail instead of using webarchive fallback when API fetch fails.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing outputs.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    refuse_overwrite([args.output, args.summary_output], args.force)

    fallback_reason = "NA"
    try:
        rows, metadata = fetch_api_rows(args.endpoint, args.expected_count, args.timeout)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError) as exc:
        fallback_reason = f"API fetch failed: {exc}"
        if args.no_fallback:
            raise
        if not args.webarchive_fallback.exists():
            raise FileNotFoundError(f"{fallback_reason}; fallback file not found: {args.webarchive_fallback}") from exc
        rows, metadata = fetch_webarchive_rows(args.webarchive_fallback, fallback_reason, args.expected_count)

    validate_rows(rows, args.expected_count)
    write_rows(args.output, rows)
    write_summary(args.summary_output, metadata)
    print(f"wrote {len(rows)} deepBase circRNA rows to {args.output}", file=sys.stderr)
    if metadata["source"] == "webarchive":
        print(f"used webarchive fallback: {fallback_reason}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
