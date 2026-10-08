#!/usr/bin/env python3
"""Download CircBank human circRNA raw data."""

from __future__ import annotations

import argparse
import shutil
import urllib.request
from datetime import date
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"

HUMAN_BASIC_INFORMATION_URL = "http://update-img.circbank.cn/upload/b6f32a7651f4481bb81bd08f344581e6.gz"
HUMAN_BASIC_INFORMATION_FILENAME = "human_basic_information.tsv.gz"
HUMAN_MATURE_SEQUENCE_URL = "http://update-img.circbank.cn/upload/376b76e8625a44eaad4eabcbbfa24c18.gz"


def download_file(url: str, output: Path, force: bool, timeout: int) -> None:
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"exists, skip: {output}")
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    with urllib.request.urlopen(request, timeout=timeout) as response, tmp_output.open("wb") as handle:
        shutil.copyfileobj(response, handle)

    tmp_output.replace(output)
    print(f"downloaded: {output}")


def write_download_links(raw_dir: Path) -> None:
    path = raw_dir / "download_links.txt"
    text = (
        "CircBank human circRNA download links\n"
        "Source page: http://www.circbank.cn/#/download\n"
        "Source API: http://www.circbank.cn/download/selectForMem\n"
        f"Checked date: {date.today().isoformat()}\n"
        "\n"
        "[Human circRNA data]\n"
        "1. human basic information\n"
        "Format: tsv.gz\n"
        "Version: v1.0\n"
        "Upload date: 2024-08-16\n"
        f"Direct URL: {HUMAN_BASIC_INFORMATION_URL}\n"
        "\n"
        "[Optional human sequence file]\n"
        "2. human mature sequence\n"
        "Format: tsv.gz\n"
        "Version: v1.0\n"
        "Upload date: 2024-08-16\n"
        f"Direct URL: {HUMAN_MATURE_SEQUENCE_URL}\n"
        "\n"
        "[Notes]\n"
        "- The integration pipeline uses only human basic information.\n"
        "- The mature sequence file is not downloaded by this script.\n"
    )
    path.write_text(text)
    print(f"wrote: {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CircBank human basic information raw table.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing downloaded files.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds. Default: 120")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    download_file(
        HUMAN_BASIC_INFORMATION_URL,
        raw_dir / HUMAN_BASIC_INFORMATION_FILENAME,
        args.force,
        args.timeout,
    )
    write_download_links(raw_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
