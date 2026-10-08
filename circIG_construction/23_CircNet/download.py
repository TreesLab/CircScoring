#!/usr/bin/env python3
"""Download CircNet human cancer circRNA CSV files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
BASE_URL = "https://awi.cuhk.edu.cn/~CircNet/file"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

CIRCRNA_FILES = [
    "ACC.csv",
    "BLCA.csv",
    "BRCA.csv",
    "CHOL.csv",
    "COAD.csv",
    "COLO.csv",
    "DLBC.csv",
    "ESCA.csv",
    "GBM.csv",
    "HNSC.csv",
    "KDNY.csv",
    "LAML.csv",
    "LGG.csv",
    "LIHC.csv",
    "LUAD.csv",
    "LUNG.csv",
    "LUSC.csv",
    "LYMP.csv",
    "MBL.csv",
    "MESO.csv",
    "MISC.csv",
    "MM.csv",
    "MPN.csv",
    "NRBL.csv",
    "OV.csv",
    "PAAD.csv",
    "PRAD.csv",
    "READ.csv",
    "RHABDO.csv",
    "SARC.csv",
    "SECR.csv",
    "SKCM.csv",
    "STAD.csv",
    "TGCT.csv",
    "THCA.csv",
    "THYM.csv",
    "UCEC.csv",
]

SUPPLEMENTARY_FILES = [
    "mir_update_0914.csv",
    "sample_0913.csv",
]

EXPECTED_HEADER = ["", "CircID", "Strand", "CircType", "HostGene", "Algorithm", "Sequence"]
SUMMARY_COLUMNS = [
    "file",
    "url",
    "status",
    "size_bytes",
    "md5",
    "row_count",
    "header_status",
]


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_csv_rows(path: Path, expect_circrna_header: bool) -> tuple[int, str]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            return 0, "empty"

        header_status = "not_checked"
        if expect_circrna_header:
            header_status = "ok" if header == EXPECTED_HEADER else f"unexpected:{header!r}"

        count = sum(1 for _ in reader)
    return count, header_status


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
    if output.exists() and output.stat().st_size > 0 and not force:
        return "exists"

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    try:
        for attempt in range(1, retries + 1):
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response, tmp_output.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
                if tmp_output.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
                tmp_output.replace(output)
                return "downloaded"
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
                tmp_output.unlink(missing_ok=True)
                delay = retry_delay(exc, attempt)
                if delay is None or attempt == retries:
                    raise
                print(
                    f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                    file=sys.stderr,
                )
                time.sleep(delay)
    except Exception:
        tmp_output.unlink(missing_ok=True)
        raise
    raise RuntimeError("unreachable retry state")


def write_record(path: Path, rows: list[dict[str, str]], base_url: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        handle.write("# CircNet download record\n")
        handle.write(f"# downloaded_at_utc\t{datetime.now(timezone.utc).isoformat()}\n")
        handle.write(f"# base_url\t{base_url}\n")
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CircNet raw circRNA CSV files.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--base-url", default=BASE_URL, help=f"CircNet file base URL. Default: {BASE_URL}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing downloaded files.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds. Default: 180")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between file requests. Default: 1.0")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive file downloads fail. Default: 5",
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Also download non-coordinate supplementary CSV files.",
    )
    args = parser.parse_args()
    if args.sleep < 0:
        parser.error("--sleep must be >= 0")
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    if args.max_consecutive_errors < 1:
        parser.error("--max-consecutive-errors must be >= 1")
    return args


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    filenames = list(CIRCRNA_FILES)
    if args.include_supplementary:
        filenames = SUPPLEMENTARY_FILES + filenames

    rows: list[dict[str, str]] = []
    failures: list[str] = []
    consecutive_errors = 0
    for index, filename in enumerate(filenames, start=1):
        url = f"{args.base_url.rstrip('/')}/{filename}"
        output = raw_dir / filename
        request_attempted = args.force or not (output.exists() and output.stat().st_size > 0)
        try:
            status = download_file(url, output, args.force, args.timeout, args.retries)
            expect_circrna_header = filename in CIRCRNA_FILES
            row_count, header_status = count_csv_rows(output, expect_circrna_header)
            size_bytes = str(output.stat().st_size)
            digest = md5sum(output)
            consecutive_errors = 0
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            row_count = "NA"
            header_status = "failed"
            size_bytes = "0"
            digest = "NA"
            failures.append(f"{url}: {exc}")
            consecutive_errors += 1
        rows.append(
            {
                "file": filename,
                "url": url,
                "status": status,
                "size_bytes": size_bytes,
                "md5": digest,
                "row_count": str(row_count),
                "header_status": header_status,
            }
        )
        print(f"{status}: {output}", file=sys.stderr)
        if consecutive_errors >= args.max_consecutive_errors:
            print(f"[stop] {consecutive_errors} consecutive file downloads failed", file=sys.stderr)
            break
        if request_attempted and index < len(filenames) and args.sleep > 0:
            time.sleep(args.sleep)

    write_record(raw_dir / "download_record.tsv", rows, args.base_url)
    print(f"wrote: {raw_dir / 'download_record.tsv'}", file=sys.stderr)
    if failures:
        print(f"failed_downloads\t{len(failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
