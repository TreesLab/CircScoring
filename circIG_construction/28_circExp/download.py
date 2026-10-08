#!/usr/bin/env python3
"""Download circExp raw data and website annotation tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import shutil
import sys
import tempfile
import time
import zipfile
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin
from urllib.request import Request, urlopen


BASE_URL = "https://bioinfo-minzhao.org/soft/circexp/"
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

OFFICIAL_RAW_FILES = {
    "dataset.tsv": "dataset.tsv",
    "annot.zip": "annot.zip",
    "diff_exp.zip": "diff_exp.zip",
    "exp_matrix.zip": "exp_matrix.zip",
}

WEB_OUTPUT_COLUMNS = [
    "source_file",
    "geo_dataset",
    "original_id",
    "circbase_id",
    "genomic_location",
    "strand",
    "feature",
    "parental_gene",
]

SUMMARY_COLUMNS = [
    "geo_dataset",
    "html_url",
    "ajax_endpoint",
    "records_total",
    "records_filtered",
    "downloaded_rows",
    "status",
    "message",
]

MANIFEST_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "validation_status",
]


class TextExtractor(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)

    def text(self) -> str:
        return "".join(self.parts)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty files.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds. Default: 120")
    parser.add_argument("--retries", type=int, default=3, help="Network retry count. Default: 3")
    parser.add_argument(
        "--pause",
        type=float,
        default=1.0,
        help="Seconds between API pages and GSE datasets; retry backoff is automatic. Default: 1.0",
    )
    parser.add_argument("--page-size", type=int, default=10000, help="DataTables page size. Default: 10000")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive GSE downloads fail. Default: 5",
    )
    parser.add_argument("--only", nargs="*", help="Optional GSE IDs to download from the website annotation tables.")
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads without writing files.")
    args = parser.parse_args()
    if args.pause < 0:
        parser.error("--pause must be >= 0")
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    if args.max_consecutive_errors < 1:
        parser.error("--max-consecutive-errors must be >= 1")
    return args


def retry_after_seconds(error: HTTPError) -> float | None:
    value = error.headers.get("Retry-After") if error.headers else None
    if not value:
        return None
    try:
        return max(float(value), 0.0)
    except ValueError:
        try:
            retry_at = parsedate_to_datetime(value)
        except (TypeError, ValueError, OverflowError):
            return None
        if retry_at.tzinfo is None:
            retry_at = retry_at.replace(tzinfo=timezone.utc)
        return max((retry_at - datetime.now(timezone.utc)).total_seconds(), 0.0)


def retry_delay(error: BaseException, attempt: int) -> float | None:
    if isinstance(error, HTTPError):
        if error.code not in RETRYABLE_HTTP_STATUS:
            return None
        server_delay = retry_after_seconds(error)
    elif isinstance(error, (URLError, TimeoutError, OSError)):
        server_delay = None
    else:
        return None
    return max(min(float(2**attempt), 60.0), server_delay or 0.0)


def request_bytes(url: str, *, timeout: int, retries: int) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last_error = exc
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                break
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RuntimeError(f"failed to fetch {url}: {last_error}")


def request_text(url: str, *, timeout: int, retries: int) -> str:
    return request_bytes(url, timeout=timeout, retries=retries).decode("utf-8", errors="replace")


def request_json(url: str, *, timeout: int, retries: int) -> dict:
    text = request_text(url, timeout=timeout, retries=retries)
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"invalid JSON from {url}: {exc}") from exc


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, output: Path, *, force: bool, timeout: int, retries: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False) as handle:
        tmp_path = Path(handle.name)

    try:
        for attempt in range(1, retries + 1):
            request = Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urlopen(request, timeout=timeout) as response:
                    with tmp_path.open("wb") as handle:
                        shutil.copyfileobj(response, handle)
                if tmp_path.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
                tmp_path.replace(output)
                return "downloaded"
            except (HTTPError, URLError, TimeoutError, OSError) as exc:
                tmp_path.unlink(missing_ok=True)
                delay = retry_delay(exc, attempt)
                if delay is None or attempt == retries:
                    raise
                print(
                    f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                    file=sys.stderr,
                )
                time.sleep(delay)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    raise RuntimeError("unreachable retry state")


def validate_official_raw(path: Path) -> str:
    if path.name == "dataset.tsv":
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
        if not header or header[0] != "GSE dataset":
            raise ValueError(f"unexpected dataset.tsv header: {header!r}")
        return "dataset_tsv_ok"

    with zipfile.ZipFile(path) as archive:
        members = archive.namelist()
        if not members:
            raise ValueError(f"zip archive is empty: {path}")
        bad_member = archive.testzip()
        if bad_member is not None:
            raise ValueError(f"zip archive failed validation at member: {bad_member}")
    return "zip_ok"


def download_official_raw(args: argparse.Namespace, raw_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    today = date.today().isoformat()
    for filename, url_path in OFFICIAL_RAW_FILES.items():
        url = urljoin(BASE_URL, url_path)
        path = raw_dir / filename
        status = download_file(url, path, force=args.force, timeout=args.timeout, retries=args.retries)
        validation_status = validate_official_raw(path)
        rows.append(
            {
                "download_date": today,
                "file": filename,
                "source_url": url,
                "status": status,
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256sum(path),
                "validation_status": validation_status,
            }
        )
        print(f"[{status}] {filename}", flush=True)
    return rows


def clean_cell(value: object) -> str:
    text = "" if value is None else str(value)
    if "<" in text and ">" in text:
        parser = TextExtractor()
        parser.feed(text)
        text = parser.text()
    return html.unescape(text).strip()


def read_gse_ids(dataset_path: Path) -> list[str]:
    with dataset_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row["GSE dataset"].strip() for row in reader if row.get("GSE dataset", "").strip()]


def parse_ajax_endpoint(html_text: str, gse: str) -> str:
    match = re.search(r'"ajax"\s*:\s*"([^"]+)"', html_text)
    if not match:
        raise RuntimeError(f"no DataTables ajax endpoint found for {gse}")
    return html.unescape(match.group(1))


def parse_table_headers(html_text: str) -> list[str]:
    match = re.search(r'<table[^>]+id="example"[^>]*>.*?<thead>\s*<tr>(.*?)</tr>', html_text, re.S)
    if not match:
        return []
    return [clean_cell(cell) for cell in re.findall(r"<th[^>]*>(.*?)</th>", match.group(1), re.S)]


def page_url(endpoint: str, draw: int, start: int, length: int) -> str:
    query = urlencode({"draw": draw, "start": start, "length": length})
    return f"{urljoin(BASE_URL, endpoint)}?{query}"


def normalize_web_row(gse: str, endpoint: str, row: list[object]) -> dict[str, str]:
    values = [clean_cell(value) for value in row]
    values.extend([""] * (6 - len(values)))
    return {
        "source_file": endpoint,
        "geo_dataset": gse,
        "original_id": values[0],
        "circbase_id": values[1],
        "genomic_location": values[2],
        "strand": values[3],
        "feature": values[4],
        "parental_gene": values[5],
    }


def read_web_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{column: row.get(column, "") for column in WEB_OUTPUT_COLUMNS} for row in reader]


def write_web_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WEB_OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def download_web_dataset(
    gse: str,
    *,
    page_size: int,
    timeout: int,
    retries: int,
    pause: float,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    html_url = urljoin(BASE_URL, f"{gse}.annot.html")
    html_text = request_text(html_url, timeout=timeout, retries=retries)
    endpoint = parse_ajax_endpoint(html_text, gse)
    headers = parse_table_headers(html_text)
    expected_headers = ["Original ID", "circBase ID", "Genomic location", "Strand", "Feature", "Parental gene"]
    if headers and headers != expected_headers:
        raise RuntimeError(f"unexpected table headers for {gse}: {headers}")

    first = request_json(page_url(endpoint, 1, 0, page_size), timeout=timeout, retries=retries)
    total = int(first.get("recordsTotal", 0))
    filtered = int(first.get("recordsFiltered", total))
    rows = [normalize_web_row(gse, endpoint, row) for row in first.get("data", [])]

    draw = 2
    for start in range(page_size, total, page_size):
        data = request_json(page_url(endpoint, draw, start, page_size), timeout=timeout, retries=retries)
        rows.extend(normalize_web_row(gse, endpoint, row) for row in data.get("data", []))
        draw += 1
        if pause:
            time.sleep(pause)

    status = {
        "geo_dataset": gse,
        "html_url": html_url,
        "ajax_endpoint": endpoint,
        "records_total": total,
        "records_filtered": filtered,
        "downloaded_rows": len(rows),
        "status": "ok" if len(rows) == total else "incomplete",
        "message": "" if len(rows) == total else f"downloaded {len(rows)} of {total}",
    }
    return rows, status


def write_web_summary(raw_dir: Path, summary: list[dict[str, object]], failures: list[dict[str, str]]) -> None:
    with (raw_dir / "web_annotation_download_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)

    with (raw_dir / "web_annotation_download_failures.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["geo_dataset", "stage", "message"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(failures)


def download_web_annotations(args: argparse.Namespace, raw_dir: Path) -> tuple[int, int]:
    dataset_path = raw_dir / "dataset.tsv"
    table_dir = raw_dir / "web_annotation_tables"
    table_dir.mkdir(parents=True, exist_ok=True)

    gse_ids = read_gse_ids(dataset_path)
    if args.only:
        wanted = set(args.only)
        gse_ids = [gse for gse in gse_ids if gse in wanted]
        missing = wanted.difference(gse_ids)
        if missing:
            print(f"warning: requested GSE IDs not in dataset.tsv: {','.join(sorted(missing))}", file=sys.stderr)

    merged_rows: list[dict[str, str]] = []
    summary: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    consecutive_errors = 0

    for index, gse in enumerate(gse_ids, start=1):
        output = table_dir / f"{gse}.web_annotation.tsv"
        print(f"[{index}/{len(gse_ids)}] {gse}", flush=True)
        request_attempted = False
        try:
            if output.exists() and output.stat().st_size > 0 and not args.force:
                rows = read_web_rows(output)
                summary.append(
                    {
                        "geo_dataset": gse,
                        "html_url": urljoin(BASE_URL, f"{gse}.annot.html"),
                        "ajax_endpoint": "NA",
                        "records_total": len(rows),
                        "records_filtered": len(rows),
                        "downloaded_rows": len(rows),
                        "status": "skipped_existing_file",
                        "message": "",
                    }
                )
            else:
                request_attempted = True
                rows, status = download_web_dataset(
                    gse,
                    page_size=args.page_size,
                    timeout=args.timeout,
                    retries=args.retries,
                    pause=args.pause,
                )
                write_web_rows(output, rows)
                summary.append(status)
                if status["status"] != "ok":
                    failures.append({"geo_dataset": gse, "stage": "download", "message": str(status["message"])})
                    consecutive_errors += 1
                else:
                    consecutive_errors = 0
            merged_rows.extend(rows)
        except Exception as exc:
            failures.append({"geo_dataset": gse, "stage": "download", "message": str(exc)})
            print(f"  failed: {exc}", file=sys.stderr, flush=True)
            consecutive_errors += 1
        if consecutive_errors >= args.max_consecutive_errors:
            print(
                f"[stop] {consecutive_errors} consecutive GSE downloads failed",
                file=sys.stderr,
            )
            break
        if request_attempted and index < len(gse_ids) and args.pause > 0:
            time.sleep(args.pause)

    write_web_rows(raw_dir / "web_circRNA_annotation.tsv", merged_rows)
    write_web_summary(raw_dir, summary, failures)
    return len(merged_rows), len(failures)


def write_manifest(raw_dir: Path, rows: list[dict[str, str]]) -> None:
    with (raw_dir / "download_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def dry_run(args: argparse.Namespace) -> int:
    raw_dir = args.raw_dir.expanduser()
    print("official_raw")
    for filename, url_path in OFFICIAL_RAW_FILES.items():
        print(f"{urljoin(BASE_URL, url_path)}\t{raw_dir / filename}")

    print("web_annotation_tables")
    dataset_path = raw_dir / "dataset.tsv"
    if dataset_path.exists():
        gse_ids = read_gse_ids(dataset_path)
        if args.only:
            wanted = set(args.only)
            gse_ids = [gse for gse in gse_ids if gse in wanted]
        for gse in gse_ids:
            print(f"{urljoin(BASE_URL, f'{gse}.annot.html')}\t{raw_dir / 'web_annotation_tables' / f'{gse}.web_annotation.tsv'}")
    else:
        print(f"dataset list will be read from {raw_dir / 'dataset.tsv'} after official raw download")
    return 0


def main() -> int:
    args = parse_args()
    args.raw_dir = args.raw_dir.expanduser().resolve()

    if args.dry_run:
        return dry_run(args)

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = download_official_raw(args, args.raw_dir)
    write_manifest(args.raw_dir, manifest_rows)

    merged_rows, failure_count = download_web_annotations(args, args.raw_dir)
    print(f"web_annotation_rows\t{merged_rows}")
    print(f"web_annotation_failures\t{failure_count}")
    return 1 if failure_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
