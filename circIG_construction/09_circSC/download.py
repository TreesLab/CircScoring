#!/usr/bin/env python3
"""Download circSC raw data files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import time
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


DOWNLOAD_PAGE = "https://ngdc.cncb.ac.cn/circatlas/circSC/download.html"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}
FILES = [
    (
        "circSC_circRNA_list.zip",
        "https://ngdc.cncb.ac.cn/circatlas/circSC/data/download/circSC_circRNA_list.zip",
    ),
    (
        "sample_detail.zip",
        "https://ngdc.cncb.ac.cn/circatlas/circSC/data/download/sample_detail.zip",
    ),
]
MANIFEST_COLUMNS = [
    "download_date",
    "download_page",
    "source_url",
    "output_file",
    "status",
    "bytes",
    "sha256",
]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=root / "raw")
    parser.add_argument("--force", action="store_true", help="Re-download existing files.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned downloads without writing files.")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between file requests. Default: 1.0")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive file downloads fail. Default: 5",
    )
    args = parser.parse_args()
    if args.sleep < 0:
        parser.error("--sleep must be >= 0")
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    if args.max_consecutive_errors < 1:
        parser.error("--max-consecutive-errors must be >= 1")
    return args


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def download(url: str, output: Path, force: bool, timeout: int, retries: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "existing"

    tmp = output.with_suffix(output.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    try:
        for attempt in range(1, retries + 1):
            request = Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urlopen(request, timeout=timeout) as response, tmp.open("wb") as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                if not tmp.exists() or tmp.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
                tmp.replace(output)
                return "downloaded"
            except (HTTPError, URLError, TimeoutError, OSError) as exc:
                tmp.unlink(missing_ok=True)
                delay = retry_delay(exc, attempt)
                if delay is None or attempt == retries:
                    raise RuntimeError(f"download failed for {url}: {exc}") from exc
                print(
                    f"[retry] attempt={attempt + 1}/{retries} "
                    f"wait={delay:.1f}s url={url}",
                    file=sys.stderr,
                )
                time.sleep(delay)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise

    raise RuntimeError("unreachable retry state")


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()

    if args.dry_run:
        print(f"download_page\t{DOWNLOAD_PAGE}")
        for filename, url in FILES:
            print(f"source_url\t{url}")
            print(f"output\t{raw_dir / filename}")
        return 0

    raw_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    failures: list[str] = []
    consecutive_errors = 0
    request_attempted = False
    for filename, url in FILES:
        output = raw_dir / filename
        needs_download = not (
            output.exists() and output.stat().st_size > 0 and not args.force
        )
        if needs_download and request_attempted and args.sleep > 0:
            time.sleep(args.sleep)
        if needs_download:
            request_attempted = True
        try:
            status = download(url, output, args.force, args.timeout, args.retries)
            size = str(output.stat().st_size)
            digest = sha256(output)
            consecutive_errors = 0
            print(f"{filename}\t{status}\t{size}")
        except (RuntimeError, OSError) as exc:
            status = f"failed: {exc}"
            size = "0"
            digest = "NA"
            failures.append(f"{url}: {exc}")
            consecutive_errors += 1
            print(f"[failed] {filename}: {exc}", file=sys.stderr)
        row = {
            "download_date": date.today().isoformat(),
            "download_page": DOWNLOAD_PAGE,
            "source_url": url,
            "output_file": str(output),
            "status": status,
            "bytes": size,
            "sha256": digest,
        }
        rows.append(row)
        if consecutive_errors >= args.max_consecutive_errors:
            print(
                f"[stop] {consecutive_errors} consecutive downloads failed",
                file=sys.stderr,
            )
            break

    write_manifest(raw_dir / "download_manifest.tsv", rows)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
