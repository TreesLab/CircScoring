#!/usr/bin/env python3
"""Download MiOncoCirc raw circRNA release data."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import HTTPCookieProcessor, Request, build_opener


SOURCE_PAGE = "https://mioncocirc.github.io/download/"
DRIVE_FILE_ID = "1jVOQqAoB64LWyVEnBza-OeYKt65r3X8i"
DRIVE_OPEN_URL = f"https://drive.google.com/open?id={DRIVE_FILE_ID}"
DRIVE_DOWNLOAD_URL = "https://drive.google.com/uc"
USER_AGENT = "data-downloader/1.0"
OUTPUT_NAME = "v0.1.release.txt"
MANIFEST_COLUMNS = [
    "download_date",
    "source_page",
    "source_url",
    "drive_file_id",
    "output_file",
    "status",
    "bytes",
]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, default=root / "raw")
    parser.add_argument("--force", action="store_true", help="Re-download existing output file.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned download without writing files.")
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def request(opener, url: str, timeout: int):
    return opener.open(
        Request(url, headers={"User-Agent": USER_AGENT}),
        timeout=timeout,
    )


def find_confirm_token(text: str) -> str | None:
    patterns = [
        r"confirm=([0-9A-Za-z_]+)&",
        r'"confirm"\s*,\s*"([0-9A-Za-z_]+)"',
        r"confirm=([0-9A-Za-z_]+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    return None


def is_html(path: Path) -> bool:
    with path.open("rb") as handle:
        prefix = handle.read(512).lower()
    return b"<html" in prefix or b"<!doctype html" in prefix


def download_from_drive(output: Path, force: bool, timeout: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "existing"
    if output.exists() or output.is_symlink():
        output.unlink()

    opener = build_opener(HTTPCookieProcessor())
    params = urlencode({"export": "download", "id": DRIVE_FILE_ID})
    url = f"{DRIVE_DOWNLOAD_URL}?{params}"
    tmp = output.with_suffix(output.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    try:
        response = request(opener, url, timeout)
        content_type = response.headers.get("Content-Type", "")
        data = response.read()
        if "text/html" in content_type.lower():
            html = data.decode("utf-8", errors="replace")
            token = find_confirm_token(html)
            if token:
                params = urlencode({"export": "download", "confirm": token, "id": DRIVE_FILE_ID})
                response = request(opener, f"{DRIVE_DOWNLOAD_URL}?{params}", timeout)
                data = response.read()
        tmp.write_bytes(data)
    except (HTTPError, URLError, TimeoutError) as exc:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"download failed: {exc}") from exc

    if not tmp.exists() or tmp.stat().st_size == 0:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError("downloaded file is empty")
    if is_html(tmp):
        tmp.unlink()
        raise RuntimeError(
            "download returned an HTML page instead of the raw data; "
            "Google Drive may require a browser confirmation or access may be restricted"
        )

    tmp.replace(output)
    return "downloaded"


def write_manifest(path: Path, row: dict[str, str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = raw_dir / OUTPUT_NAME

    if args.dry_run:
        print(f"source_page\t{SOURCE_PAGE}")
        print(f"source_url\t{DRIVE_OPEN_URL}")
        print(f"drive_file_id\t{DRIVE_FILE_ID}")
        print(f"output\t{output}")
        return 0

    raw_dir.mkdir(parents=True, exist_ok=True)
    status = download_from_drive(output, args.force, args.timeout)
    row = {
        "download_date": date.today().isoformat(),
        "source_page": SOURCE_PAGE,
        "source_url": DRIVE_OPEN_URL,
        "drive_file_id": DRIVE_FILE_ID,
        "output_file": str(output),
        "status": status,
        "bytes": str(output.stat().st_size),
    }
    write_manifest(raw_dir / "download_manifest.tsv", row)
    print(f"status\t{status}")
    print(f"output\t{output}")
    print(f"bytes\t{output.stat().st_size}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
