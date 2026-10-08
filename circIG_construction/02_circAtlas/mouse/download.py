#!/usr/bin/env python3
"""Download circAtlas 3.0 mouse circRNA raw data."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "raw"

API_URL = "https://ngdc.cncb.ac.cn/circatlas/browse_species1.php"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}
SPECIES = "mouse"

DOWNLOADS = [
    {
        "name": "mouse_bed_v3.0.zip",
        "url": "https://ngdc.cncb.ac.cn/circatlas/analysis/bed/mouse_bed_v3.0.zip",
        "extract": "mouse_bed_v3.0.txt",
    },
    {
        "name": "mouse_sequence_v3.0.zip",
        "url": "https://ngdc.cncb.ac.cn/circatlas/analysis/bed/mouse_sequence_v3.0.zip",
        "extract": "mouse_sequence_v3.0",
    },
]

API_COLUMNS = [
    ("species", "Species"),
    ("name", "circAtlas ID"),
    ("uid", "Uniform ID"),
    ("pos", "Position"),
    ("strand", "Strand"),
    ("ctype", "circRNA type"),
    ("score", "Multiple Conservation Score(MCS)"),
    ("nspe", "MCS species"),
    ("ntis", "MCS tissue"),
    ("nsam", "MCS sample"),
    ("len", "Length (nt)"),
    ("tlen", "Type length"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download circAtlas 3.0 mouse BED, sequence, and browse table raw data."
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
        "--page-size",
        type=int,
        default=10000,
        help="Rows requested per circAtlas browse API page. Default: 10000",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Download/API timeout in seconds per request. Default: 120",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Seconds to sleep between browse API requests. Default: 1.0",
    )
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    args = parser.parse_args()
    if args.sleep < 0:
        parser.error("--sleep must be >= 0")
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
            request = urllib.request.Request(url, headers={"User-Agent": "data-downloader/1.0"})
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
                print(f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}", file=sys.stderr)
                time.sleep(delay)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def extract_member(zip_path: Path, member_name: str, output_path: Path, force: bool) -> str:
    if output_path.exists() and output_path.stat().st_size > 0 and not force:
        return "skipped"

    with zipfile.ZipFile(zip_path) as archive:
        members = {Path(info.filename).name: info for info in archive.infolist()}
        if member_name not in members:
            available = ", ".join(sorted(members))
            raise RuntimeError(
                f"{member_name} not found in {zip_path.name}; available files: {available}"
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
            delete=False,
        ) as tmp_handle:
            tmp_path = Path(tmp_handle.name)

        try:
            with archive.open(members[member_name]) as source:
                with tmp_path.open("wb") as target:
                    shutil.copyfileobj(source, target)

            if tmp_path.stat().st_size == 0:
                raise RuntimeError(f"extracted file is empty: {member_name}")

            tmp_path.replace(output_path)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

    return "extracted"


def api_payload(start: int, length: int, draw: int) -> bytes:
    payload: list[tuple[str, str]] = [
        ("species", SPECIES),
        ("draw", str(draw)),
        ("start", str(start)),
        ("length", str(length)),
        ("order[0][column]", "0"),
        ("order[0][dir]", "asc"),
    ]
    for idx, (key, _label) in enumerate(API_COLUMNS):
        payload.append((f"columns[{idx}][data]", key))
    return urllib.parse.urlencode(payload).encode()


def fetch_api_page(start: int, length: int, draw: int, timeout: int, retries: int) -> dict[str, Any]:
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            API_URL,
            data=api_payload(start, length, draw),
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "User-Agent": "data-downloader/1.0",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode())
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                raise
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={API_URL}",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable retry state")


def fetch_browse_table(output_path: Path, page_size: int, timeout: int, retries: int, sleep: float) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    first = fetch_api_page(start=0, length=1, draw=1, timeout=timeout, retries=retries)
    total = int(first["recordsTotal"])
    print(f"[api] total {SPECIES} records reported by circAtlas: {total}", file=sys.stderr)

    with tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        delete=False,
    ) as tmp_handle:
        tmp_path = Path(tmp_handle.name)
        writer = csv.writer(tmp_handle, delimiter="\t", lineterminator="\n")
        writer.writerow([label for _key, label in API_COLUMNS])

        written = 0
        seen_names: set[str] = set()
        page_index = 0
        try:
            if total > 0 and sleep > 0:
                time.sleep(sleep)
            while written < total:
                start = page_index * page_size
                page = fetch_api_page(
                    start=start,
                    length=page_size,
                    draw=page_index + 2,
                    timeout=timeout,
                    retries=retries,
                )
                rows = page["data"]
                if not rows:
                    raise RuntimeError(
                        f"No rows returned at start={start}; aborting to avoid silent truncation."
                    )

                for row in rows:
                    circ_id = row.get("name", "")
                    if circ_id in seen_names:
                        continue
                    seen_names.add(circ_id)
                    writer.writerow(
                        [
                            row.get(key) if row.get(key) is not None else "NA"
                            for key, _label in API_COLUMNS
                        ]
                    )
                    written += 1

                print(
                    f"[api] fetched start={start} rows={len(rows)} unique_written={written}",
                    file=sys.stderr,
                )
                page_index += 1
                if written < total and sleep > 0:
                    time.sleep(sleep)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise

    if written != total:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"Expected {total} rows, wrote {written}.")

    tmp_path.replace(output_path)
    return written


def count_tsv_records(path: Path) -> int | str:
    if not path.exists() or path.stat().st_size == 0:
        return "NA"
    with path.open() as handle:
        lines = sum(1 for _line in handle)
    return max(lines - 1, 0)


def write_download_links(raw_dir: Path) -> None:
    output_path = raw_dir / "download_links.txt"
    with output_path.open("w") as handle:
        for item in DOWNLOADS:
            handle.write(f"{item['name']}\t{item['url']}\n")
        handle.write(f"mouse_browse_table.tsv\t{API_URL}\tspecies={SPECIES}\n")


def write_summary(raw_dir: Path, rows: list[dict[str, str]]) -> None:
    output_path = raw_dir / "download_summary.tsv"
    fieldnames = ["file", "source", "status", "records", "bytes", "url"]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    failed: list[tuple[str, str]] = []
    summary_rows: list[dict[str, str]] = []

    for item in DOWNLOADS:
        zip_path = raw_dir / item["name"]
        extracted_path = raw_dir / item["extract"]
        try:
            if zip_path.exists() and zip_path.stat().st_size > 0 and not args.force:
                download_status = "skipped"
                print(f"[skip] {zip_path.name} exists ({zip_path.stat().st_size} bytes)")
            else:
                print(f"[download] {zip_path.name}")
                download_file(item["url"], zip_path, args.timeout, args.retries)
                download_status = "downloaded"
                print(f"[done] {zip_path.name} ({zip_path.stat().st_size} bytes)")

            extract_status = extract_member(
                zip_path=zip_path,
                member_name=item["extract"],
                output_path=extracted_path,
                force=args.force,
            )
            print(f"[{extract_status}] {extracted_path.name} ({extracted_path.stat().st_size} bytes)")

            status = download_status if extract_status == "skipped" else f"{download_status}+{extract_status}"
            summary_rows.append(
                {
                    "file": extracted_path.name,
                    "source": "official_download_zip",
                    "status": status,
                    "records": str(count_tsv_records(extracted_path))
                    if extracted_path.suffix == ".txt"
                    else "NA",
                    "bytes": str(extracted_path.stat().st_size),
                    "url": item["url"],
                }
            )
        except (urllib.error.URLError, TimeoutError, RuntimeError, OSError, zipfile.BadZipFile) as exc:
            print(f"[failed] {item['name']}: {exc}", file=sys.stderr)
            failed.append((item["name"], str(exc)))

    browse_path = raw_dir / "mouse_browse_table.tsv"
    try:
        if browse_path.exists() and browse_path.stat().st_size > 0 and not args.force:
            print(f"[skip] {browse_path.name} exists ({browse_path.stat().st_size} bytes)")
            browse_status = "skipped"
            browse_records = count_tsv_records(browse_path)
        else:
            print(f"[api] fetch {browse_path.name}")
            browse_records = fetch_browse_table(
                output_path=browse_path,
                page_size=args.page_size,
                timeout=args.timeout,
                retries=args.retries,
                sleep=args.sleep,
            )
            browse_status = "downloaded"
            print(f"[done] {browse_path.name} ({browse_records} records)")

        summary_rows.append(
            {
                "file": browse_path.name,
                "source": "browse_api",
                "status": browse_status,
                "records": str(browse_records),
                "bytes": str(browse_path.stat().st_size),
                "url": API_URL,
            }
        )
    except (urllib.error.URLError, TimeoutError, RuntimeError, OSError, json.JSONDecodeError) as exc:
        print(f"[failed] {browse_path.name}: {exc}", file=sys.stderr)
        failed.append((browse_path.name, str(exc)))

    write_download_links(raw_dir)
    write_summary(raw_dir, summary_rows)
    print(f"[summary] wrote {raw_dir / 'download_links.txt'}")
    print(f"[summary] wrote {raw_dir / 'download_summary.tsv'}")

    if failed:
        print("Failed downloads:", file=sys.stderr)
        for filename, error in failed:
            print(f"- {filename}\t{error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
