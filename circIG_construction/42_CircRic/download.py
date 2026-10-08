#!/usr/bin/env python3
"""Download CircRic raw circRNA expression data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from datetime import date
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


DOWNLOAD_PAGE = "https://hanlaboratory.com/cRic/download"
SOURCE_URL = "https://hanlaboratory.com/static_UT/download/circRNA_expression.csv.zip"
OUTPUT_NAME = "circRNA_expression.csv.zip"
USER_AGENT = "data-downloader/1.0"
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
    parser.add_argument("--force", action="store_true", help="Re-download existing file.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned download without writing files.")
    parser.add_argument("--timeout", type=int, default=120)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(output: Path, force: bool, timeout: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "existing"
    if output.exists() or output.is_symlink():
        output.unlink()

    tmp = output.with_suffix(output.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    request = Request(SOURCE_URL, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=timeout) as response, tmp.open("wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    except (HTTPError, URLError, TimeoutError) as exc:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError(f"download failed: {exc}") from exc

    if not tmp.exists() or tmp.stat().st_size == 0:
        if tmp.exists():
            tmp.unlink()
        raise RuntimeError("downloaded file is empty")

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
        print(f"download_page\t{DOWNLOAD_PAGE}")
        print(f"source_url\t{SOURCE_URL}")
        print(f"output\t{output}")
        return 0

    raw_dir.mkdir(parents=True, exist_ok=True)
    status = download(output, args.force, args.timeout)
    row = {
        "download_date": date.today().isoformat(),
        "download_page": DOWNLOAD_PAGE,
        "source_url": SOURCE_URL,
        "output_file": str(output),
        "status": status,
        "bytes": str(output.stat().st_size),
        "sha256": sha256(output),
    }
    write_manifest(raw_dir / "download_manifest.tsv", row)
    print(f"status\t{status}")
    print(f"output\t{output}")
    print(f"bytes\t{output.stat().st_size}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
