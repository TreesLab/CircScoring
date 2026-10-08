#!/usr/bin/env python3
"""Download CircFunBase mouse circRNA raw data."""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "raw"

BROWSE_URL = "https://bis.zju.edu.cn/CircFunBase/searchresult.php"
NAME_LIST_URL = "https://bis.zju.edu.cn/CircFunBase/assets/Data/json/Mus_musculus_circ.name.json"
API_URL = "https://bis.zju.edu.cn/CircFunBase/Api/Public/circfunapi/"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}

BROWSE_COLUMNS = [
    "browse_id",
    "count",
    "circRNA",
    "gene_symbol",
    "gene_link",
    "species",
    "function",
    "detail_path",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download CircFunBase Mus musculus browse records and rebuild "
            "single-entry CircRNA.getinfo API responses as raw JSONL."
        )
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Output directory for raw files. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing raw files.")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip API records already present in mouse_api_results.jsonl.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Limit API calls for testing. Browse/name-list files are still downloaded.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=1.0,
        help="Seconds between API calls. Use at least 1 second for the public site. Default: 1.0",
    )
    parser.add_argument("--timeout", type=int, default=60, help="Request timeout in seconds. Default: 60")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    parser.add_argument(
        "--max-consecutive-errors",
        type=int,
        default=5,
        help="Stop after this many consecutive API records fail. Default: 5",
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
    exponential_delay = min(float(2**attempt), 60.0)
    return max(exponential_delay, server_delay or 0.0)


def request(
    url: str,
    *,
    data: bytes | None = None,
    method: str | None = None,
    timeout: int,
    retries: int,
) -> bytes:
    headers = {"User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                raise
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable retry state")


def write_atomic_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        tmp_path.write_bytes(data)
        if tmp_path.stat().st_size == 0:
            raise RuntimeError(f"downloaded file is empty: {path.name}")
        tmp_path.replace(path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def download_name_list(raw_dir: Path, timeout: int, retries: int, force: bool) -> None:
    output = raw_dir / "Mus_musculus_circ.name.json"
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"[skip] {output.name}")
        return
    print(f"[download] {output.name}")
    write_atomic_bytes(output, request(NAME_LIST_URL, timeout=timeout, retries=retries))


def fetch_browse_json(raw_dir: Path, timeout: int, retries: int, force: bool) -> list[dict[str, Any]]:
    output = raw_dir / "mouse_browse_results.json"
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"[skip] {output.name}")
        return json.loads(output.read_text())

    payload = urllib.parse.urlencode({"species": "Mus musculus", "type": "browse"}).encode()
    print(f"[download] {output.name}")
    data = request(BROWSE_URL, data=payload, method="POST", timeout=timeout, retries=retries)
    records = json.loads(data.decode())
    if not isinstance(records, list):
        raise RuntimeError("browse endpoint did not return a JSON list")
    write_atomic_bytes(output, json.dumps(records, ensure_ascii=False, indent=2).encode())
    return records


def normalize(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    return text if text else "NA"


def html_text(value: Any) -> str:
    text = normalize(value)
    if text == "NA":
        return text
    text = re.sub(r"<[^>]*>", "", text)
    text = html.unescape(text).strip()
    return text if text else "NA"


def html_href(value: Any) -> str:
    text = normalize(value)
    if text == "NA":
        return text
    match = re.search(r"""href\s*=\s*["']([^"']+)["']""", text, flags=re.IGNORECASE)
    if not match:
        return "NA"
    href = html.unescape(match.group(1)).strip()
    return href if href else "NA"


def browse_row(record: dict[str, Any]) -> dict[str, str]:
    return {
        "browse_id": normalize(record.get("ID")),
        "count": normalize(record.get("count")),
        "circRNA": normalize(record.get("circRNA")),
        "gene_symbol": html_text(record.get("Gene")),
        "gene_link": normalize(record.get("Gene_link")),
        "species": normalize(record.get("Species")),
        "function": normalize(record.get("Function")),
        "detail_path": html_href(record.get("Detail")),
    }


def write_browse_table(raw_dir: Path, records: list[dict[str, Any]], force: bool) -> Path:
    output = raw_dir / "mouse_browse_table.tsv"
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"[skip] {output.name}")
        return output

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
        delete=False,
    ) as tmp:
        tmp_path = Path(tmp.name)
        writer = csv.DictWriter(tmp, fieldnames=BROWSE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(browse_row(record))

    try:
        tmp_path.replace(output)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise
    print(f"[write] {output.name}")
    return output


def load_existing_api_names(path: Path) -> set[str]:
    if not path.exists():
        return set()
    names: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            name = normalize(record.get("circRNA"))
            if name != "NA":
                names.add(name)
    return names


def api_url(name: str) -> str:
    query = urllib.parse.urlencode({"service": "CircRNA.getinfo", "circRNAname": name})
    return f"{API_URL}?{query}"


def fetch_api_response(name: str, timeout: int, retries: int) -> dict[str, Any]:
    return json.loads(request(api_url(name), timeout=timeout, retries=retries).decode())


def fetch_api_records(
    raw_dir: Path,
    names: list[str],
    *,
    timeout: int,
    retries: int,
    sleep: float,
    max_consecutive_errors: int,
    force: bool,
    resume: bool,
    limit: int | None,
) -> None:
    output = raw_dir / "mouse_api_results.jsonl"
    errors = raw_dir / "mouse_api_errors.tsv"

    if output.exists() and output.stat().st_size > 0 and not force and not resume:
        print(f"[skip] {output.name}")
        return

    completed = set() if force else load_existing_api_names(output) if resume else set()
    mode = "a" if resume and output.exists() and not force else "w"
    selected_names = [name for name in names if name not in completed]
    if limit is not None:
        selected_names = selected_names[:limit]

    raw_dir.mkdir(parents=True, exist_ok=True)
    with output.open(mode) as out_handle, errors.open("w", newline="") as err_handle:
        err_writer = csv.writer(err_handle, delimiter="\t", lineterminator="\n")
        err_writer.writerow(["circRNA", "error"])
        consecutive_errors = 0
        stopped_early = False
        for index, name in enumerate(selected_names, start=1):
            try:
                response = fetch_api_response(name, timeout, retries)
                out_handle.write(json.dumps({"circRNA": name, "response": response}, ensure_ascii=False) + "\n")
                out_handle.flush()
                consecutive_errors = 0
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as exc:
                err_writer.writerow([name, str(exc)])
                err_handle.flush()
                consecutive_errors += 1
                if consecutive_errors >= max_consecutive_errors:
                    print(
                        f"[stop] {consecutive_errors} consecutive API records failed",
                        file=sys.stderr,
                    )
                    stopped_early = True
                    break
            if sleep > 0 and index < len(selected_names):
                time.sleep(sleep)
            if index % 100 == 0:
                print(f"[api] {index}/{len(selected_names)} queried", file=sys.stderr)

    print(f"[write] {output.name}")
    print(f"[write] {errors.name}")
    if stopped_early:
        raise RuntimeError(
            f"stopped after {max_consecutive_errors} consecutive API record failures"
        )


def write_download_links(raw_dir: Path) -> None:
    output = raw_dir / "download_links.txt"
    text = (
        "CircFunBase Mus musculus raw endpoint URLs\n"
        f"browse POST: {BROWSE_URL} ; species=Mus musculus&type=browse\n"
        f"name list: {NAME_LIST_URL}\n"
        f"single-entry API: {API_URL}?service=CircRNA.getinfo&circRNAname=<circRNA>\n"
    )
    output.write_text(text)
    print(f"[write] {output.name}")


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    try:
        download_name_list(raw_dir, args.timeout, args.retries, args.force)
        browse_records = fetch_browse_json(raw_dir, args.timeout, args.retries, args.force)
        browse_table = write_browse_table(raw_dir, browse_records, args.force)
        write_download_links(raw_dir)
        names = []
        with browse_table.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                name = normalize(row.get("circRNA"))
                if name != "NA" and name not in names:
                    names.append(name)
        fetch_api_records(
            raw_dir,
            names,
            timeout=args.timeout,
            retries=args.retries,
            sleep=args.sleep,
            max_consecutive_errors=args.max_consecutive_errors,
            force=args.force,
            resume=args.resume,
            limit=args.limit,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
