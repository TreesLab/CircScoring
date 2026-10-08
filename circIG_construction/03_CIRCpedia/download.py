#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
ANNOTATION_URL = "https://bits.fudan.edu.cn/circpediav3/static/download_cache/annotation/Homo-sapiens.txt.zip"
USER_AGENT = "data-downloader/1.0"


def download_file(url: str, output: Path, force: bool, timeout: int) -> None:
    if output.exists() and not force:
        print(f"exists, skip: {output}")
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    with urllib.request.urlopen(request, timeout=timeout) as response, tmp_output.open("wb") as handle:
        shutil.copyfileobj(response, handle)

    tmp_output.replace(output)
    print(f"downloaded: {output}")


def extract_annotation(zip_path: Path, raw_dir: Path, force: bool) -> Path:
    output = raw_dir / "Homo-sapiens.txt"
    if output.exists() and not force:
        print(f"exists, skip: {output}")
        return output

    with zipfile.ZipFile(zip_path) as archive:
        members = [name for name in archive.namelist() if Path(name).name == "Homo-sapiens.txt"]
        if not members:
            raise ValueError(f"Homo-sapiens.txt not found in archive: {zip_path}")
        member = members[0]
        tmp_output = output.with_suffix(output.suffix + ".tmp")
        with archive.open(member) as source, tmp_output.open("wb") as target:
            shutil.copyfileobj(source, target)
        tmp_output.replace(output)

    print(f"extracted: {output}")
    return output


def write_download_links(raw_dir: Path) -> None:
    path = raw_dir / "download_links.txt"
    text = (
        "CIRCpedia v3 Homo sapiens circRNA annotation\n"
        f"{ANNOTATION_URL}\n"
    )
    path.write_text(text)
    print(f"wrote: {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CIRCpedia v3 Homo sapiens raw circRNA annotation.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Overwrite existing downloaded/extracted files.")
    parser.add_argument("--timeout", type=int, default=120, help="Network timeout in seconds. Default: 120")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir
    raw_dir.mkdir(parents=True, exist_ok=True)

    zip_path = raw_dir / "Homo-sapiens.txt.zip"
    download_file(ANNOTATION_URL, zip_path, args.force, args.timeout)
    extract_annotation(zip_path, raw_dir, args.force)
    write_download_links(raw_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
