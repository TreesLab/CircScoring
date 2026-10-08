#!/usr/bin/env python3
"""Download ExoRBase raw circRNA files for the unified pipeline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import tempfile
import time
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import quote, urljoin
import urllib.error
import urllib.request


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_BASE_URL = "http://exorbase.org/"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

ANNOTATION_FILENAME = "circRNAs_anno.csv"
ANNOTATION_PATH = "data/anno/circRNAs_anno.csv"
EXPRESSION_DIRNAME = "raw_expression"
EXPRESSION_SOURCE_FOLDER = "circRNA"
EXPRESSION_SUFFIX = "_circRNAs.txt"

EXPECTED_COLUMNS = {
    "circID",
    "Genomic position",
    "Strand",
    "Gene symbol",
}

EXPRESSION_COHORTS = [
    "EVPs in Benign",
    "EVPs in Bile",
    "EVPs in BRCA",
    "EVPs in CRC",
    "EVPs in CSF",
    "EVPs in GBM",
    "EVPs in GC",
    "EVPs in HCC",
    "EVPs in Healthy",
    "EVPs in KIRC",
    "EVPs in ML",
    "EVPs in MEL",
    "EVPs in NSCLC",
    "EVPs in OV",
    "EVPs in PDAC",
    "EVPs in SCLC",
    "EVPs in Urine",
    "CHB",
    "CRC",
    "CSF",
    "ESCA",
    "HCC",
    "Healthy",
    "Liver_Cirrhosis",
    "NSCLC",
    "MGUS",
    "MM",
    "GC",
    "Urine",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download ExoRBase circRNA annotation raw data. By default this "
            "downloads only circRNAs_anno.csv, which is the raw coordinate "
            "annotation used by clean.py."
        )
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help=f"ExoRBase site base URL. Default: {DEFAULT_BASE_URL}",
    )
    parser.add_argument(
        "--include-expression",
        action="store_true",
        help="Also download the 29 circRNA expression raw matrices.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download files even when non-empty destination files already exist.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=180,
        help="Network timeout in seconds per file. Default: 180",
    )
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between file requests. Default: 1.0")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive expression downloads fail. Default: 5",
    )
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


def annotation_url(base_url: str) -> str:
    return urljoin(normalize_base_url(base_url), ANNOTATION_PATH)


def expression_source_name(cohort: str) -> str:
    return cohort.replace(" ", "_") + EXPRESSION_SUFFIX


def expression_url(base_url: str, cohort: str) -> str:
    folder = quote(EXPRESSION_SOURCE_FOLDER)
    filename = quote(expression_source_name(cohort))
    return urljoin(normalize_base_url(base_url), f"data/Downloads/{folder}/{filename}")


def download_file(url: str, output: Path, *, force: bool, timeout: int, retries: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"[skip] {output}")
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
        delete=False,
    ) as tmp_handle:
        tmp_path = Path(tmp_handle.name)

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
                print(f"[downloaded] {output}")
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


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_annotation_rows(path: Path) -> tuple[int, list[str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty annotation file: {path}") from exc
        missing = sorted(EXPECTED_COLUMNS - set(header))
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")
        rows = sum(1 for _row in reader)
    return rows, header


def expression_output_path(raw_dir: Path, cohort: str) -> Path:
    return raw_dir / EXPRESSION_DIRNAME / expression_source_name(cohort)


def download_expression_files(
    raw_dir: Path,
    base_url: str,
    *,
    force: bool,
    timeout: int,
    retries: int,
    sleep: float,
    max_consecutive_errors: int,
) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    failures: list[str] = []
    consecutive_errors = 0
    for index, cohort in enumerate(EXPRESSION_COHORTS, start=1):
        output = expression_output_path(raw_dir, cohort)
        url = expression_url(base_url, cohort)
        request_attempted = force or not (output.exists() and output.stat().st_size > 0)
        try:
            status = download_file(url, output, force=force, timeout=timeout, retries=retries)
            size_bytes = str(output.stat().st_size)
            digest = md5sum(output)
            consecutive_errors = 0
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            size_bytes = "0"
            digest = "NA"
            failures.append(f"{url}: {exc}")
            consecutive_errors += 1
            print(f"[failed] {url}: {exc}", file=sys.stderr)
        rows.append(
            {
                "file": str(output.relative_to(raw_dir)),
                "source_url": url,
                "status": status,
                "bytes": size_bytes,
                "md5": digest,
            }
        )
        if consecutive_errors >= max_consecutive_errors:
            print(f"[stop] {consecutive_errors} consecutive expression downloads failed", file=sys.stderr)
            break
        if request_attempted and index < len(EXPRESSION_COHORTS) and sleep > 0:
            time.sleep(sleep)
    return rows, failures


def existing_expression_count(raw_dir: Path) -> int:
    expression_dir = raw_dir / EXPRESSION_DIRNAME
    return len(sorted(expression_dir.glob("*.txt"))) if expression_dir.exists() else 0


def write_record(
    raw_dir: Path,
    *,
    base_url: str,
    annotation_status: str,
    annotation_rows: int,
    annotation_header: list[str],
    expression_rows: list[dict[str, str]],
    include_expression: bool,
) -> None:
    annotation = raw_dir / ANNOTATION_FILENAME
    record = raw_dir / "download_record.txt"
    lines = [
        "ExoRBase raw data download record",
        f"checked_date: {date.today().isoformat()}",
        f"base_url: {normalize_base_url(base_url)}",
        f"annotation_source_url: {annotation_url(base_url)}",
        f"annotation_file: {annotation}",
        f"annotation_status: {annotation_status}",
        f"annotation_size_bytes: {annotation.stat().st_size}",
        f"annotation_md5: {md5sum(annotation)}",
        f"annotation_rows_excluding_header: {annotation_rows}",
        f"annotation_columns: {len(annotation_header)}",
        f"include_expression: {str(include_expression).lower()}",
        f"expression_matrix_files_present: {existing_expression_count(raw_dir)}",
        "required_for_clean: circRNAs_anno.csv",
        "notes: download.py downloads raw files only; no cleaning, coordinate conversion, or expression merging is performed.",
        "",
        "annotation_columns_list:",
        *[f"- {column}" for column in annotation_header],
    ]
    if expression_rows:
        lines.extend(["", "expression_downloads:", "file\tsource_url\tstatus\tbytes\tmd5"])
        lines.extend(
            "\t".join(
                [
                    row["file"],
                    row["source_url"],
                    row["status"],
                    row["bytes"],
                    row["md5"],
                ]
            )
            for row in expression_rows
        )
    record.write_text("\n".join(lines) + "\n")
    print(f"[write] {record}")


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)
    base_url = normalize_base_url(args.base_url)

    annotation = raw_dir / ANNOTATION_FILENAME
    expression_rows: list[dict[str, str]] = []
    expression_failures: list[str] = []

    try:
        annotation_status = download_file(
            annotation_url(base_url),
            annotation,
            force=args.force,
            timeout=args.timeout,
            retries=args.retries,
        )
        if args.include_expression:
            expression_rows, expression_failures = download_expression_files(
                raw_dir,
                base_url,
                force=args.force,
                timeout=args.timeout,
                retries=args.retries,
                sleep=args.sleep,
                max_consecutive_errors=args.max_consecutive_errors,
            )
        annotation_rows, annotation_header = count_annotation_rows(annotation)
        write_record(
            raw_dir,
            base_url=base_url,
            annotation_status=annotation_status,
            annotation_rows=annotation_rows,
            annotation_header=annotation_header,
            expression_rows=expression_rows,
            include_expression=args.include_expression,
        )
    except (urllib.error.URLError, TimeoutError, RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"validated annotation rows: {annotation_rows}")
    print(f"expression matrix files present: {existing_expression_count(raw_dir)}")
    if expression_failures:
        print(f"failed_expression_downloads\t{len(expression_failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
