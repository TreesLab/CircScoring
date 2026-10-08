#!/usr/bin/env python3
"""Download TSCD human hg38 circRNA raw data."""

from __future__ import annotations

import argparse
import shutil
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

DOWNLOADS = [
    (
        "hg38_adult_TS_circRNAs.tar.gz",
        "https://www.geneyun.net/TSCD/download/hg38_adult_TS_circRNAs.tar.gz",
    ),
    (
        "hg38_fetal_TS_circRNAs.tar.gz",
        "https://www.geneyun.net/TSCD/download/hg38_fetal_TS_circRNAs.tar.gz",
    ),
]


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


def download_file(
    url: str,
    output: Path,
    force: bool,
    timeout: int,
    retries: int,
) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"exists, skip: {output}")
        return "existing"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
        delete=False,
    ) as tmp_handle:
        tmp_output = Path(tmp_handle.name)

    try:
        for attempt in range(1, retries + 1):
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            context = None
            if url.startswith("https://"):
                import ssl

                context = ssl._create_unverified_context()

            try:
                with urllib.request.urlopen(
                    request,
                    timeout=timeout,
                    context=context,
                ) as response, tmp_output.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
                if tmp_output.stat().st_size == 0:
                    raise RuntimeError(f"downloaded file is empty: {url}")
                tmp_output.replace(output)
                print(f"downloaded: {output}")
                return "downloaded"
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
                tmp_output.unlink(missing_ok=True)
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
        tmp_output.unlink(missing_ok=True)
        raise

    raise RuntimeError("unreachable retry state")


def extract_tarball(path: Path, raw_dir: Path, force: bool) -> None:
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            output = raw_dir / Path(member.name).name
            if output.exists() and output.stat().st_size > 0 and not force:
                print(f"exists, skip: {output}")
                continue
            source = archive.extractfile(member)
            if source is None:
                continue
            with tempfile.NamedTemporaryFile(
                prefix=f".{output.name}.",
                suffix=".tmp",
                dir=output.parent,
                delete=False,
            ) as tmp_handle:
                tmp_output = Path(tmp_handle.name)
            try:
                with source, tmp_output.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
                if tmp_output.stat().st_size == 0:
                    raise RuntimeError(f"extracted file is empty: {member.name}")
                tmp_output.replace(output)
            except Exception:
                tmp_output.unlink(missing_ok=True)
                raise
            print(f"extracted: {output}")


def write_download_links(raw_dir: Path) -> None:
    lines = [
        "TSCD human hg38 circRNA download links",
        "Homepage: https://www.geneyun.net/TSCD/",
        f"Checked date: {date.today().isoformat()}",
        "",
    ]
    for filename, url in DOWNLOADS:
        lines.append(f"{filename}\t{url}")
    lines.extend(
        [
            "",
            "Notes:",
            "- This pipeline uses only human hg38 adult/fetal tissue-specific circRNA tables.",
            "- Download uses an unverified TLS context because the TSCD certificate has previously been expired.",
            "- Download stage only saves and extracts raw files; no cleaning or coordinate conversion is done here.",
        ]
    )
    path = raw_dir / "download_links.txt"
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote: {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download TSCD human hg38 adult/fetal circRNA raw tables.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing downloaded and extracted files.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds. Default: 180")
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


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    failures: list[str] = []
    consecutive_errors = 0
    request_attempted = False

    for filename, url in DOWNLOADS:
        output = raw_dir / filename
        needs_download = not (
            output.exists() and output.stat().st_size > 0 and not args.force
        )
        if needs_download and request_attempted and args.sleep > 0:
            time.sleep(args.sleep)
        if needs_download:
            request_attempted = True
        try:
            download_file(url, output, args.force, args.timeout, args.retries)
            extract_tarball(output, raw_dir, args.force)
            consecutive_errors = 0
        except (urllib.error.URLError, TimeoutError, RuntimeError, OSError, tarfile.TarError) as exc:
            failures.append(f"{url}: {exc}")
            consecutive_errors += 1
            print(f"[failed] {filename}: {exc}", file=sys.stderr)
            if consecutive_errors >= args.max_consecutive_errors:
                print(
                    f"[stop] {consecutive_errors} consecutive downloads failed",
                    file=sys.stderr,
                )
                break

    write_download_links(raw_dir)
    if failures:
        print("Failed downloads:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
