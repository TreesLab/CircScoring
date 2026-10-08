#!/usr/bin/env python3
"""Download CircMine raw circRNA records.

This script intentionally keeps the raw API shape. It downloads the CircMine
project list, fetches one JSON response per HSACM dataset, then writes merged
raw TSV/JSONL files with the API fields only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


BASE_URL = "http://www.biomedical-web.com"
PROJECT_LIST_URL = BASE_URL + "/circmine_v2/tool/home/proj.view.txt"
DATASET_API_URL = BASE_URL + "/circmine_v2/circrnalyer_v2/circRNAs-of-circMine"
RAW_FIELDS = ["circmine", "circrna", "circbank", "hg19", "hg38", "strand", "hostGene"]
USER_AGENT = "data-downloader/1.0"
DATASET_ID_RE = re.compile(r"^HSACM\d{6}$")
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CircMine raw circRNA data.")
    parser.add_argument("--raw-dir", type=Path, default=Path(__file__).resolve().parent / "raw")
    parser.add_argument("--project-list-url", default=PROJECT_LIST_URL)
    parser.add_argument("--dataset-api-url", default=DATASET_API_URL)
    parser.add_argument("--force", action="store_true", help="Redownload existing files.")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Seconds between API requests. Use at least 1 second for the public site. Default: 1.0",
    )
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive dataset requests fail. Default: 5",
    )
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads only.")
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


def request_bytes(
    url: str,
    *,
    data: bytes | None = None,
    timeout: int,
    retries: int,
) -> bytes:
    headers = {"User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "text/plain;charset=UTF-8"
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data else "GET")
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                raise RuntimeError(f"failed to fetch {url}: {exc}") from exc
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable retry state")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_project_list(args: argparse.Namespace, project_path: Path) -> None:
    if project_path.exists() and project_path.stat().st_size > 0 and not args.force:
        return
    content = request_bytes(args.project_list_url, timeout=args.timeout, retries=args.retries)
    project_path.write_bytes(content)


def read_dataset_ids(project_path: Path) -> tuple[list[str], dict[str, int]]:
    ids: list[str] = []
    counts: dict[str, int] = {}
    with project_path.open(newline="") as handle:
        for line in handle:
            dataset_id = line.rstrip("\n\r").split(";", 1)[0].strip()
            if not DATASET_ID_RE.match(dataset_id):
                continue
            counts[dataset_id] = counts.get(dataset_id, 0) + 1
            if dataset_id not in ids:
                ids.append(dataset_id)
    return ids, counts


def validate_dataset_json(path: Path) -> list[dict[str, object]]:
    with path.open() as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict) or "results" not in payload:
        raise ValueError(f"{path} does not contain a results object")
    results = payload["results"]
    if not isinstance(results, list):
        raise ValueError(f"{path} results is not a list")
    for row in results:
        if not isinstance(row, dict):
            raise ValueError(f"{path} contains a non-object result")
    return results


def download_dataset_jsons(args: argparse.Namespace, dataset_ids: list[str], pages_dir: Path) -> None:
    pages_dir.mkdir(parents=True, exist_ok=True)
    total = len(dataset_ids)
    failures: list[str] = []
    consecutive_errors = 0
    for index, dataset_id in enumerate(dataset_ids, start=1):
        path = pages_dir / f"{dataset_id}.json"
        if path.exists() and path.stat().st_size > 0 and not args.force:
            validate_dataset_json(path)
            continue
        if args.dry_run:
            print(f"would download {dataset_id} ({index}/{total})")
            continue
        try:
            content = request_bytes(
                args.dataset_api_url,
                data=dataset_id.encode(),
                timeout=args.timeout,
                retries=args.retries,
            )
            path.write_bytes(content)
            validate_dataset_json(path)
            consecutive_errors = 0
        except (OSError, RuntimeError, ValueError) as exc:
            path.unlink(missing_ok=True)
            failures.append(f"{dataset_id}: {exc}")
            consecutive_errors += 1
            print(f"[failed] {dataset_id}: {exc}", file=sys.stderr)
            if consecutive_errors >= args.max_consecutive_errors:
                print(
                    f"[stop] {consecutive_errors} consecutive dataset requests failed",
                    file=sys.stderr,
                )
                break
        if args.sleep and index < total:
            time.sleep(args.sleep)

    if failures:
        raise RuntimeError(
            f"{len(failures)} dataset request(s) failed; first failure: {failures[0]}"
        )


def norm(value: object) -> str:
    if value is None:
        return "NA"
    text = str(value)
    return text if text else "NA"


def merge_raw_outputs(
    raw_dir: Path,
    dataset_ids: list[str],
    source_counts: dict[str, int],
) -> dict[str, object]:
    pages_dir = raw_dir / "circRNAs_of_circMine_pages"
    merged_tsv = raw_dir / "human_circRNAs_of_circMine.raw_merged.tsv"
    merged_jsonl = raw_dir / "human_circRNAs_of_circMine.raw_merged.jsonl"
    counts_tsv = raw_dir / "human_circRNAs_of_circMine.dataset_counts.tsv"
    manifest_tsv = raw_dir / "download_manifest.tsv"

    total_records = 0
    empty_datasets = []
    with merged_tsv.open("w", newline="") as tsv_handle, merged_jsonl.open("w") as jsonl_handle:
        writer = csv.DictWriter(tsv_handle, fieldnames=RAW_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        with counts_tsv.open("w", newline="") as counts_handle:
            counts_writer = csv.writer(counts_handle, delimiter="\t", lineterminator="\n")
            counts_writer.writerow(["circmine", "records", "source_occurrences"])
            for dataset_id in dataset_ids:
                path = pages_dir / f"{dataset_id}.json"
                results = validate_dataset_json(path)
                counts_writer.writerow([dataset_id, len(results), source_counts.get(dataset_id, 0)])
                if not results:
                    empty_datasets.append(dataset_id)
                for row in results:
                    writer.writerow({field: norm(row.get(field)) for field in RAW_FIELDS})
                    jsonl_handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
                    total_records += 1

    downloaded_at = datetime.now(timezone.utc).date().isoformat()
    manifest_rows = [
        ("project_list_url", PROJECT_LIST_URL),
        ("dataset_api_url", DATASET_API_URL),
        ("downloaded_at_utc_date", downloaded_at),
        ("project_list_rows_with_dataset_id", sum(source_counts.values())),
        ("unique_dataset_ids", len(dataset_ids)),
        ("empty_dataset_ids", ",".join(empty_datasets) if empty_datasets else "NA"),
        ("raw_records", total_records),
        ("merged_tsv_sha256", sha256(merged_tsv)),
        ("merged_jsonl_sha256", sha256(merged_jsonl)),
    ]
    with manifest_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        writer.writerows(manifest_rows)

    return {
        "unique_dataset_ids": len(dataset_ids),
        "raw_records": total_records,
        "empty_dataset_ids": ",".join(empty_datasets) if empty_datasets else "NA",
        "merged_tsv": str(merged_tsv),
        "merged_jsonl": str(merged_jsonl),
    }


def main() -> int:
    args = parse_args()
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    project_path = args.raw_dir / "circmine_v2_proj.view.txt"

    if args.dry_run:
        print(f"project list: {args.project_list_url}")
        print(f"dataset API: {args.dataset_api_url}")
        return 0

    download_project_list(args, project_path)
    dataset_ids, source_counts = read_dataset_ids(project_path)
    if not dataset_ids:
        raise SystemExit("no HSACM dataset IDs found in project list")
    download_dataset_jsons(args, dataset_ids, args.raw_dir / "circRNAs_of_circMine_pages")
    summary = merge_raw_outputs(args.raw_dir, dataset_ids, source_counts)
    for key, value in summary.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
