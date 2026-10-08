#!/usr/bin/env python3
"""Download CIRIonco human circRNA raw API data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
import tempfile
import time
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import urlencode
import urllib.error
import urllib.request


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
API_URL = "https://ngdc.cncb.ac.cn/cirionco/api/circrnas/"
USER_AGENT = "data-downloader/1.0"
PAGE_SIZE = 1000
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

RAW_COLUMNS = [
    "CIRC_ID",
    "regulation_status",
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Gene_Name",
    "CIRC_type",
    "Tissue_number",
    "Tissue_cancer_number",
    "Tissue_BSJ_number",
    "Tissue_cancer_BSJ_number",
    "Tissue_names",
    "cancer_tissue_names",
    "start_exon",
    "end_exon",
    "start_exon_boundary",
    "end_exon_boundary",
    "exon_composition",
    "system",
    "tissue",
    "disease",
    "system_regulation",
    "tissue_regulation",
    "disease_regulation",
]

MANIFEST_COLUMNS = [
    "download_date",
    "source_url",
    "page",
    "page_size",
    "output_file",
    "status",
    "result_count",
    "size_bytes",
    "sha256",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CIRIonco human circRNA API pages.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--page-size", type=int, default=PAGE_SIZE, help=f"API page size. Default: {PAGE_SIZE}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty API page files.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds. Default: 120")
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Delay between API page requests. Use at least 1 second for the public site. Default: 1.0",
    )
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive API pages fail. Default: 5",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned download URLs without downloading.")
    args = parser.parse_args()
    if args.sleep < 0:
        parser.error("--sleep must be >= 0")
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    if args.max_consecutive_errors < 1:
        parser.error("--max-consecutive-errors must be >= 1")
    return args


def retry_after_seconds(error: urllib.error.HTTPError) -> float | None:
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
    if isinstance(error, urllib.error.HTTPError):
        if error.code not in RETRYABLE_HTTP_STATUS:
            return None
        server_delay = retry_after_seconds(error)
    elif isinstance(error, (urllib.error.URLError, TimeoutError, OSError)):
        server_delay = None
    else:
        return None
    exponential_delay = min(float(2**attempt), 60.0)
    return max(exponential_delay, server_delay or 0.0)


def api_page_url(page: int, page_size: int) -> str:
    return f"{API_URL}?{urlencode({'page': page, 'page_size': page_size})}"


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_json(url: str, output: Path, *, force: bool, timeout: int, retries: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False) as handle:
        tmp_path = Path(handle.name)

    try:
        for attempt in range(1, retries + 1):
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    with tmp_path.open("wb") as handle:
                        shutil.copyfileobj(response, handle)
                if tmp_path.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
                json.loads(tmp_path.read_text())
                tmp_path.replace(output)
                return "downloaded"
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
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


def load_page(path: Path) -> dict:
    with path.open() as handle:
        obj = json.load(handle)
    if not isinstance(obj, dict):
        raise ValueError(f"API page is not a JSON object: {path}")
    if "results" not in obj or not isinstance(obj["results"], list):
        raise ValueError(f"API page has no results list: {path}")
    return obj


def normalize_value(value) -> str:
    if value is None:
        return ""
    return str(value)


def write_merged_outputs(raw_dir: Path, page_count: int) -> tuple[int, Path, Path]:
    tsv_path = raw_dir / "human_circrnas_api_merged.tsv"
    jsonl_path = raw_dir / "human_circrnas_api_merged.jsonl"
    rows_written = 0

    with (
        tsv_path.open("w", newline="") as tsv_handle,
        jsonl_path.open("w") as jsonl_handle,
    ):
        writer = csv.DictWriter(tsv_handle, fieldnames=RAW_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()

        for page in range(1, page_count + 1):
            page_path = raw_dir / "api_pages" / f"circrnas_page_{page}.json"
            obj = load_page(page_path)
            for record in obj["results"]:
                if not isinstance(record, dict):
                    raise ValueError(f"non-object record in {page_path}")
                writer.writerow({column: normalize_value(record.get(column)) for column in RAW_COLUMNS})
                jsonl_handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
                rows_written += 1

    return rows_written, tsv_path, jsonl_path


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    pages_dir = raw_dir / "api_pages"
    raw_dir.mkdir(parents=True, exist_ok=True)
    pages_dir.mkdir(parents=True, exist_ok=True)

    page_size = args.page_size
    first_url = api_page_url(1, page_size)
    first_page = pages_dir / "circrnas_page_1.json"

    if args.dry_run:
        print(first_url)
        return 0

    today = date.today().isoformat()
    manifest_rows: list[dict[str, str]] = []
    failures: list[str] = []

    try:
        status = fetch_json(
            first_url,
            first_page,
            force=args.force,
            timeout=args.timeout,
            retries=args.retries,
        )
        first_obj = load_page(first_page)
        total_count = int(first_obj["count"])
        page_count = math.ceil(total_count / page_size)
        print(f"[{status}] circrnas_page_1.json count={total_count} pages={page_count}")
    except (OSError, urllib.error.URLError, ValueError, RuntimeError, KeyError) as exc:
        print(f"[failed] circrnas_page_1.json: {exc}", file=sys.stderr)
        return 1

    first_result_count = len(first_obj["results"])
    manifest_rows.append(
        {
            "download_date": today,
            "source_url": first_url,
            "page": "1",
            "page_size": str(page_size),
            "output_file": str(first_page.relative_to(BASE_DIR)),
            "status": status,
            "result_count": str(first_result_count),
            "size_bytes": str(first_page.stat().st_size),
            "sha256": sha256sum(first_page),
        }
    )
    if status == "downloaded" and page_count > 1 and args.sleep > 0:
        time.sleep(args.sleep)

    consecutive_errors = 0
    for page in range(2, page_count + 1):
        url = api_page_url(page, page_size)
        page_path = pages_dir / f"circrnas_page_{page}.json"
        request_attempted = False
        try:
            status = fetch_json(
                url,
                page_path,
                force=args.force,
                timeout=args.timeout,
                retries=args.retries,
            )
            request_attempted = status == "downloaded"
            result_count = len(load_page(page_path)["results"])
            print(f"[{status}] circrnas_page_{page}.json rows={result_count}")
            consecutive_errors = 0
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            result_count = "NA"
            failures.append(f"page {page}: {exc}")
            consecutive_errors += 1
            request_attempted = True
            print(f"[failed] circrnas_page_{page}.json: {exc}", file=sys.stderr)

        manifest_rows.append(
            {
                "download_date": today,
                "source_url": url,
                "page": str(page),
                "page_size": str(page_size),
                "output_file": str(page_path.relative_to(BASE_DIR)),
                "status": status,
                "result_count": str(result_count),
                "size_bytes": str(page_path.stat().st_size if page_path.exists() else 0),
                "sha256": sha256sum(page_path) if page_path.exists() and page_path.stat().st_size > 0 else "NA",
            }
        )
        if consecutive_errors >= args.max_consecutive_errors:
            print(
                f"[stop] {consecutive_errors} consecutive API pages failed",
                file=sys.stderr,
            )
            break
        if request_attempted and page < page_count and args.sleep > 0:
            time.sleep(args.sleep)

    write_manifest(raw_dir / "download_manifest.tsv", manifest_rows)
    if failures:
        print(f"failed_pages\t{len(failures)}", file=sys.stderr)
        return 1

    rows_written, tsv_path, jsonl_path = write_merged_outputs(raw_dir, page_count)
    if rows_written != total_count:
        raise RuntimeError(f"merged row count {rows_written} does not match API count {total_count}")

    print(f"merged_rows\t{rows_written}")
    print(f"merged_tsv\t{tsv_path}")
    print(f"merged_jsonl\t{jsonl_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
