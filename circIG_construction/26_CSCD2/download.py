#!/usr/bin/env python3
"""Download CSCD2 human hg38 circRNA raw files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
import sys
import tempfile
import time
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
import urllib.error
import urllib.request


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_BASE_URL = "http://geneyun.net/CSCD2/public/static/download/"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

GROUPS = ["cancer", "normal", "common"]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
EXPECTED_HEADER = [
    "circRNA",
    "circbase_id",
    "type",
    "sample_source",
    "position",
    "host_gene",
    "reads_counts",
    "algorithm",
    "chain",
    "subcellular_location",
    "gene_type",
    "average_read_counts",
    "average_log2SRPTM",
    "log2SRPTM",
    "alternative_splicing",
    "ratios",
    "works_cited",
    "mioncocirc",
    "number_of_algorithm",
    "sequence",
]
RECORD_COLUMNS = [
    "download_date",
    "group",
    "chrom",
    "url",
    "path",
    "status",
    "size_bytes",
    "md5",
    "row_count",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CSCD2 raw hg38 circRNA annotation files.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL, help=f"Download base URL. Default: {DEFAULT_BASE_URL}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty files.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds per file. Default: 180")
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Seconds between file requests. Default: 1.0",
    )
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive file downloads fail. Default: 5",
    )
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads without downloading files.")
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
    return max(min(float(2**attempt), 60.0), server_delay or 0.0)


def normalize_base_url(base_url: str) -> str:
    return base_url.rstrip("/") + "/"


def filename_for(group: str, chrom: str) -> str:
    return f"hg38_{group}_circrna_circ_{chrom}.txt.gz"


def iter_targets(base_url: str, raw_dir: Path) -> list[dict[str, str | Path]]:
    base_url = normalize_base_url(base_url)
    targets: list[dict[str, str | Path]] = []
    for group in GROUPS:
        for chrom in CHROMS:
            filename = filename_for(group, chrom)
            targets.append(
                {
                    "group": group,
                    "chrom": chrom,
                    "url": base_url + filename,
                    "path": raw_dir / group / filename,
                }
            )
    return targets


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gzip_table(path: Path) -> int:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty gzip table: {path}") from exc
        if header != EXPECTED_HEADER:
            raise ValueError(f"unexpected header in {path}: {header!r}")
        return sum(1 for _row in reader)


def download_file(url: str, output: Path, *, force: bool, timeout: int, retries: int) -> str:
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


def write_record(record_path: Path, records: list[dict[str, str]]) -> None:
    record_path.parent.mkdir(parents=True, exist_ok=True)
    with record_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    targets = iter_targets(args.base_url, raw_dir)

    if args.dry_run:
        for target in targets:
            print(f"{target['group']}\t{target['chrom']}\t{target['url']}\t{target['path']}")
        return 0

    records: list[dict[str, str]] = []
    failures: list[str] = []
    consecutive_errors = 0
    today = date.today().isoformat()

    for index, target in enumerate(targets, start=1):
        group = str(target["group"])
        chrom = str(target["chrom"])
        url = str(target["url"])
        path = Path(target["path"])
        request_attempted = args.force or not (path.exists() and path.stat().st_size > 0)
        try:
            status = download_file(
                url,
                path,
                force=args.force,
                timeout=args.timeout,
                retries=args.retries,
            )
            row_count = validate_gzip_table(path)
            size_bytes = path.stat().st_size
            digest = md5sum(path)
            print(f"[{status}] {group}/{path.name} rows={row_count}")
            consecutive_errors = 0
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            size_bytes = path.stat().st_size if path.exists() else 0
            digest = md5sum(path) if path.exists() and path.stat().st_size > 0 else "NA"
            row_count = "NA"
            failures.append(f"{url}: {exc}")
            consecutive_errors += 1
            print(f"[failed] {url}: {exc}", file=sys.stderr)

        records.append(
            {
                "download_date": today,
                "group": group,
                "chrom": chrom,
                "url": url,
                "path": str(path),
                "status": status,
                "size_bytes": str(size_bytes),
                "md5": digest,
                "row_count": str(row_count),
            }
        )
        if consecutive_errors >= args.max_consecutive_errors:
            print(
                f"[stop] {consecutive_errors} consecutive file downloads failed",
                file=sys.stderr,
            )
            break
        if request_attempted and index < len(targets) and args.sleep > 0:
            time.sleep(args.sleep)

    write_record(raw_dir / "download_record.tsv", records)

    if failures:
        print(f"failed_downloads\t{len(failures)}", file=sys.stderr)
        return 1

    print(f"downloaded_or_validated_files\t{len(targets)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
