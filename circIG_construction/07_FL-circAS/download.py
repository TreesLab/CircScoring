#!/usr/bin/env python3
"""Download FL-circAS human circRNA raw coordinate data."""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "raw"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

DOWNLOADS = [
    (
        "BSJ_human.csv",
        "https://cosbi.ee.ncku.edu.tw/st52021_static/download/BSJ_human.csv",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download the FL-circAS Homo sapiens BSJ-level raw circRNA table."
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Output directory for raw files. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even when non-empty destination files already exist.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Download timeout in seconds per file. Default: 120",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Maximum attempts per request. Default: 3",
    )
    args = parser.parse_args()
    if args.retries < 1:
        parser.error("--retries must be >= 1")
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
    return max(min(float(2**attempt), 60.0), server_delay or 0.0)


def download_file(url: str, output_path: Path, timeout: int, retries: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        delete=False,
    ) as tmp_handle:
        tmp_path = Path(tmp_handle.name)

    try:
        for attempt in range(1, retries + 1):
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "data-downloader/1.0"},
            )
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    with tmp_path.open("wb") as tmp_handle:
                        shutil.copyfileobj(response, tmp_handle)

                if tmp_path.stat().st_size == 0:
                    raise RuntimeError("downloaded file is empty")

                tmp_path.replace(output_path)
                return
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
                tmp_path.unlink(missing_ok=True)
                delay = retry_delay(exc, attempt)
                if delay is None or attempt == retries:
                    raise
                print(
                    f"[retry] attempt={attempt + 1}/{retries} "
                    f"wait={delay:.1f}s url={url}",
                    file=sys.stderr,
                )
                time.sleep(delay)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    failed = []
    skipped = 0
    downloaded = 0

    for filename, url in DOWNLOADS:
        output_path = raw_dir / filename
        if output_path.exists() and output_path.stat().st_size > 0 and not args.force:
            print(f"[skip] {filename} exists ({output_path.stat().st_size} bytes)")
            skipped += 1
            continue

        print(f"[download] {filename}")
        try:
            download_file(url, output_path, args.timeout, args.retries)
        except (urllib.error.URLError, TimeoutError, RuntimeError, OSError) as exc:
            print(f"[failed] {filename}: {exc}", file=sys.stderr)
            failed.append((filename, url, str(exc)))
            continue

        print(f"[done] {filename} ({output_path.stat().st_size} bytes)")
        downloaded += 1

    print(
        f"Summary: {downloaded} downloaded, {skipped} skipped, "
        f"{len(failed)} failed"
    )

    if failed:
        print("Failed downloads:", file=sys.stderr)
        for filename, url, error in failed:
            print(f"- {filename}\t{url}\t{error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
