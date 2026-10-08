#!/usr/bin/env python3
"""Download TCCIA raw circRNA files from Zenodo."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sys
import tempfile
from datetime import date
from pathlib import Path
import urllib.error
import urllib.request


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
ZENODO_RECORD_ID = "10091408"
DEFAULT_RECORD_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"
USER_AGENT = "data-downloader/1.0"

RAW_FILES = {
    "TCCIA_Ensemble_circRNAs.tsv.gz": ["ID", "Gene", "Mean", "Median", "Max", "SD", "N", "Cohort"],
    "TCCIA_4_methods_circRNAs.tsv.gz": ["cohort", "chr", "start", "end", "strand", "tool"],
}
RECORD_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "md5",
    "row_count",
    "validation_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download TCCIA raw circRNA files from Zenodo record 10091408.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--record-url", default=DEFAULT_RECORD_URL, help=f"Zenodo record API URL. Default: {DEFAULT_RECORD_URL}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty files.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds per file. Default: 180")
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads without downloading files.")
    return parser.parse_args()


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, output: Path, *, force: bool, timeout: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False) as handle:
        tmp_path = Path(handle.name)

    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with tmp_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        if tmp_path.stat().st_size == 0:
            raise RuntimeError(f"downloaded file is empty: {url}")
        tmp_path.replace(output)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    return "downloaded"


def load_or_download_record(raw_dir: Path, record_url: str, *, force: bool, timeout: int) -> tuple[dict, str]:
    record_path = raw_dir / f"zenodo_record_{ZENODO_RECORD_ID}.json"
    status = download_file(record_url, record_path, force=force, timeout=timeout)
    with record_path.open() as handle:
        return json.load(handle), status


def zenodo_file_urls(record: dict) -> dict[str, str]:
    urls: dict[str, str] = {}
    for item in record.get("files", []):
        key = item.get("key")
        if key in RAW_FILES:
            urls[key] = item.get("links", {}).get("self", "")
    missing = sorted(set(RAW_FILES) - set(urls))
    if missing:
        raise ValueError(f"Zenodo metadata is missing expected file(s): {', '.join(missing)}")
    return urls


def validate_gzip_table(path: Path, expected_header: list[str]) -> int:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty gzip table: {path}") from exc
        if header != expected_header:
            raise ValueError(f"unexpected header in {path}: {header!r}")
        return sum(1 for _row in reader)


def write_record(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(args.record_url)
        for filename in RAW_FILES:
            print(filename)
        return 0

    failures: list[str] = []
    rows: list[dict[str, str]] = []
    today = date.today().isoformat()

    try:
        record, record_status = load_or_download_record(raw_dir, args.record_url, force=args.force, timeout=args.timeout)
        file_urls = zenodo_file_urls(record)
        record_path = raw_dir / f"zenodo_record_{ZENODO_RECORD_ID}.json"
        rows.append(
            {
                "download_date": today,
                "file": record_path.name,
                "source_url": args.record_url,
                "status": record_status,
                "size_bytes": str(record_path.stat().st_size),
                "md5": md5sum(record_path),
                "row_count": "NA",
                "validation_status": "json_ok",
            }
        )
    except (OSError, urllib.error.URLError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"[failed] Zenodo metadata: {exc}", file=sys.stderr)
        return 1

    for filename, expected_header in RAW_FILES.items():
        path = raw_dir / filename
        url = file_urls[filename]
        try:
            status = download_file(url, path, force=args.force, timeout=args.timeout)
            row_count = validate_gzip_table(path, expected_header)
            validation_status = "gzip_tsv_ok"
            print(f"[{status}] {filename} rows={row_count}")
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            row_count = "NA"
            validation_status = "failed"
            failures.append(f"{filename}: {exc}")
            print(f"[failed] {filename}: {exc}", file=sys.stderr)

        rows.append(
            {
                "download_date": today,
                "file": filename,
                "source_url": url,
                "status": status,
                "size_bytes": str(path.stat().st_size if path.exists() else 0),
                "md5": md5sum(path) if path.exists() and path.stat().st_size > 0 else "NA",
                "row_count": str(row_count),
                "validation_status": validation_status,
            }
        )

    write_record(raw_dir / "download_record.tsv", rows)

    if failures:
        print(f"failed_downloads\t{len(failures)}", file=sys.stderr)
        return 1

    print(f"downloaded_or_validated_files\t{len(RAW_FILES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
