#!/usr/bin/env python3
"""Download CircRNADb raw circRNA dataset."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_URL = "http://reprod.njmu.edu.cn/circrnadb/doc/circRNA_dataset.zip"
ZIP_NAME = "circRNA_dataset.zip"
EXTRACT_DIR_NAME = "circRNA_dataset"
REQUIRED_RAW_FILES = [
    "circRNA_dataset.txt",
    "circRNA_fasta",
    "circRNA_IRES",
    "circRNA_ORF",
]
REQUIRED_COORDINATE_FIELDS = 15
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    return max(min(float(2**attempt), 60.0), server_delay or 0.0)


def download_file(url: str, output: Path, force: bool, timeout: int, retries: int) -> str:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.", suffix=".tmp", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        for attempt in range(1, retries + 1):
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response, tmp_path.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
                if tmp_path.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
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


def extract_zip(zip_path: Path, raw_dir: Path, force: bool) -> str:
    extract_dir = raw_dir / EXTRACT_DIR_NAME
    required_paths = [extract_dir / name for name in REQUIRED_RAW_FILES]
    if extract_dir.exists() and all(path.exists() for path in required_paths) and not force:
        return "skipped_existing_extract"

    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.namelist():
            name = Path(member).name
            if not name:
                continue
            with archive.open(member) as source, (extract_dir / name).open("wb") as target:
                shutil.copyfileobj(source, target)

    missing = [str(path) for path in required_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Extracted zip is missing required file(s): {', '.join(missing)}")
    return "extracted"


def count_coordinate_rows(path: Path) -> tuple[int, int]:
    rows = 0
    invalid_field_count = 0
    with path.open(errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            rows += 1
            if len(row) != REQUIRED_COORDINATE_FIELDS:
                invalid_field_count += 1
    return rows, invalid_field_count


def write_record(path: Path, values: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        for key, value in values.items():
            handle.write(f"{key}: {value}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CircRNADb raw dataset zip.")
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Raw dataset zip URL. Default: {DEFAULT_URL}")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Redownload and re-extract raw data.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds. Default: 180")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    args = parser.parse_args()
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    return args


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    zip_path = raw_dir / ZIP_NAME
    extract_dir = raw_dir / EXTRACT_DIR_NAME

    download_status = download_file(args.url, zip_path, args.force, args.timeout, args.retries)
    extract_status = extract_zip(zip_path, raw_dir, args.force)
    coordinate_path = extract_dir / "circRNA_dataset.txt"
    coordinate_rows, invalid_field_count = count_coordinate_rows(coordinate_path)

    write_record(
        raw_dir / "download_record.txt",
        {
            "download_date": datetime.now().isoformat(timespec="seconds"),
            "source_url": args.url,
            "download_status": download_status,
            "extract_status": extract_status,
            "zip_path": str(zip_path),
            "zip_size_bytes": str(zip_path.stat().st_size),
            "zip_md5": md5sum(zip_path),
            "extract_dir": str(extract_dir),
            "extracted_files": ",".join(REQUIRED_RAW_FILES),
            "coordinate_table": str(coordinate_path),
            "coordinate_rows": str(coordinate_rows),
            "coordinate_rows_with_invalid_field_count": str(invalid_field_count),
            "required_for_clean": "circRNA_dataset/circRNA_dataset.txt",
        },
    )
    print(f"[{download_status}] {zip_path}")
    print(f"[{extract_status}] {extract_dir}")
    print(f"[write] {raw_dir / 'download_record.txt'}")
    print(f"coordinate rows: {coordinate_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
