#!/usr/bin/env python3
"""Download BloodCircR raw files for the unified pipeline."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
RAW_FILENAME = "CircRNA_Atlas"
SOURCE_URL = "https://mega.nz/file/zKwF2ZBb#pKYs6luAWg3mNjJaEunCoUarXwbSiQWRGlRVVxkE1N8"
USER_AGENT = "data-downloader/1.0"
RETRYABLE_HTTP_STATUS = {429, 500, 502, 503, 504}
EXPECTED_COLUMNS = {
    "BloodCircR_ID",
    "chr",
    "circRNA_start",
    "circRNA_end",
    "strand",
    "host_gene",
}


def b64url_decode(value: str) -> bytes:
    value = value.replace("-", "+").replace("_", "/")
    value += "=" * ((4 - len(value) % 4) % 4)
    return base64.b64decode(value)


def parse_mega_url(url: str) -> tuple[str, str]:
    if "/file/" not in url or "#" not in url:
        raise ValueError("Expected a public MEGA file URL like https://mega.nz/file/<id>#<key>")
    file_id = url.split("/file/", 1)[1].split("#", 1)[0]
    file_key = url.split("#", 1)[1]
    return file_id, file_key


def derive_keys(file_key: str) -> tuple[bytes, bytes]:
    raw = b64url_decode(file_key)
    if len(raw) != 32:
        raise ValueError(f"Expected 32 decoded key bytes, got {len(raw)}")
    words = struct.unpack(">8I", raw)
    aes_words = (
        words[0] ^ words[4],
        words[1] ^ words[5],
        words[2] ^ words[6],
        words[3] ^ words[7],
    )
    aes_key = struct.pack(">4I", *aes_words)
    iv = struct.pack(">4I", words[4], words[5], 0, 0)
    return aes_key, iv


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
    elif isinstance(error, (urllib.error.URLError, TimeoutError, OSError, RuntimeError)):
        server_delay = None
    else:
        return None
    return max(min(float(2**attempt), 60.0), server_delay or 0.0)


def mega_api_request(file_id: str, timeout: int, retries: int) -> dict[str, object]:
    body = json.dumps([{"a": "g", "g": 1, "p": file_id}]).encode("utf-8")
    url = "https://g.api.mega.co.nz/cs?id=1"
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            data=body,
            headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            break
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                raise
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={url}",
                file=sys.stderr,
            )
            time.sleep(delay)
    else:
        raise RuntimeError("unreachable retry state")
    item = payload[0]
    if isinstance(item, int):
        raise RuntimeError(f"MEGA API returned error code {item}")
    if "g" not in item:
        raise RuntimeError(f"MEGA API response did not include a download URL: {item}")
    return item


def download_encrypted(
    download_url: str,
    encrypted_path: Path,
    expected_size: int | None,
    timeout: int,
    retries: int,
) -> int:
    encrypted_path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        bytes_written = 0
        request = urllib.request.Request(download_url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response, encrypted_path.open("wb") as out:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
                    bytes_written += len(chunk)
            if expected_size is not None and bytes_written != expected_size:
                raise RuntimeError(f"Downloaded {bytes_written} encrypted bytes, expected {expected_size}")
            return bytes_written
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError, RuntimeError) as exc:
            encrypted_path.unlink(missing_ok=True)
            delay = retry_delay(exc, attempt)
            if delay is None or attempt == retries:
                raise
            print(
                f"[retry] attempt={attempt + 1}/{retries} wait={delay:.1f}s url={download_url}",
                file=sys.stderr,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable retry state")


def decrypt_with_openssl(encrypted_path: Path, output_path: Path, aes_key: bytes, iv: bytes) -> None:
    openssl = shutil.which("openssl")
    if openssl is None:
        raise FileNotFoundError("openssl executable not found; required to decrypt the MEGA download")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(output_path.name + ".part")
    try:
        subprocess.run(
            [
                openssl,
                "enc",
                "-d",
                "-aes-128-ctr",
                "-nopad",
                "-K",
                aes_key.hex(),
                "-iv",
                iv.hex(),
                "-in",
                str(encrypted_path),
                "-out",
                str(tmp_path),
            ],
            check=True,
        )
        os.replace(tmp_path, output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def download_mega_file(url: str, output_path: Path, timeout: int, retries: int, sleep: float) -> dict[str, object]:
    file_id, file_key = parse_mega_url(url)
    aes_key, iv = derive_keys(file_key)
    metadata = mega_api_request(file_id, timeout, retries)
    expected_size = metadata.get("s")
    if expected_size is not None:
        expected_size = int(expected_size)
    encrypted_path = output_path.with_name(output_path.name + ".enc.part")
    try:
        if sleep > 0:
            time.sleep(sleep)
        bytes_written = download_encrypted(
            str(metadata["g"]),
            encrypted_path,
            expected_size,
            timeout,
            retries,
        )
        decrypt_with_openssl(encrypted_path, output_path, aes_key, iv)
    finally:
        encrypted_path.unlink(missing_ok=True)
    return {"file_id": file_id, "expected_size": expected_size, "downloaded_encrypted_bytes": bytes_written}


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_atlas(path: Path) -> tuple[int, list[str]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"empty raw atlas: {path}") from exc
        missing = sorted(EXPECTED_COLUMNS - set(header))
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")
        rows = sum(1 for _row in reader)
    return rows, header


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download BloodCircR raw atlas for the unified pipeline.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--url", default=SOURCE_URL, help="Official public MEGA URL.")
    parser.add_argument("--force", action="store_true", help="Re-download and overwrite the raw atlas if it already exists.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds. Default: 120")
    parser.add_argument("--sleep", type=float, default=1.0, help="Seconds between metadata and file requests. Default: 1.0")
    parser.add_argument("--retries", type=int, default=3, help="Maximum attempts per request. Default: 3")
    args = parser.parse_args()
    if args.sleep < 0:
        parser.error("--sleep must be >= 0")
    if args.retries < 1:
        parser.error("--retries must be >= 1")
    return args


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    atlas = raw_dir / RAW_FILENAME
    record = raw_dir / "download_record.txt"

    raw_dir.mkdir(parents=True, exist_ok=True)
    download_status = "skipped_existing_file"
    download_metadata: dict[str, object] = {}
    if args.force or not atlas.exists() or atlas.stat().st_size == 0:
        print(f"downloading {args.url}", flush=True)
        download_metadata = download_mega_file(
            args.url,
            atlas,
            args.timeout,
            args.retries,
            args.sleep,
        )
        download_status = "downloaded"
    else:
        print(f"raw atlas already exists; skipping download: {atlas}", flush=True)

    rows, header = inspect_atlas(atlas)
    digest = md5sum(atlas)
    size = atlas.stat().st_size

    record.write_text(
        "\n".join(
            [
                "BloodCircR raw data validation",
                f"checked_date: {date.today().isoformat()}",
                f"source_url: {args.url}",
                f"download_status: {download_status}",
                f"download_metadata: {json.dumps(download_metadata, sort_keys=True)}",
                f"raw_file: {atlas}",
                f"size_bytes: {size}",
                f"rows_excluding_header: {rows}",
                f"md5: {digest}",
                "required_for_clean: yes",
                "notes: This script downloads the official raw atlas when needed, then validates it. It does not clean or transform data.",
                "columns:",
                *[f"- {column}" for column in header],
                "",
            ]
        )
    )

    print(f"validated {atlas}")
    print(f"rows_excluding_header={rows}")
    print(f"md5={digest}")
    print(f"wrote {record}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
