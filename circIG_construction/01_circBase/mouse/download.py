#!/usr/bin/env python3
"""Download circBase mouse circRNA raw coordinate tables."""

from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = SCRIPT_DIR / "raw"
BASE_URL = "https://www.circbase.org/download"

FILENAMES = [
    "mmu_mm9_circRNA.txt",
    "mmu_mm9_Rybak2015.txt",
    "mmu_mm9_Memczak2013.txt",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download circBase Mus musculus mm9 raw TXT coordinate files."
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
        "--timeout",
        type=int,
        default=60,
        help="Download timeout in seconds per file. Default: 60",
    )
    return parser.parse_args()


def download_file(url: str, output_path: Path, timeout: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        delete=False,
    ) as tmp_handle:
        tmp_path = Path(tmp_handle.name)

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "data-downloader/1.0"},
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with tmp_path.open("wb") as tmp_handle:
                shutil.copyfileobj(response, tmp_handle)

        if tmp_path.stat().st_size == 0:
            raise RuntimeError("downloaded file is empty")

        tmp_path.replace(output_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def write_download_links(raw_dir: Path) -> None:
    links_path = raw_dir / "download_links.txt"
    with links_path.open("w") as handle:
        for filename in FILENAMES:
            handle.write(f"{BASE_URL}/{filename}\n")


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    failed = 0
    skipped = 0
    downloaded = 0

    for filename in FILENAMES:
        url = f"{BASE_URL}/{filename}"
        output_path = raw_dir / filename

        if output_path.exists() and output_path.stat().st_size > 0 and not args.force:
            size = output_path.stat().st_size
            print(f"[skip] {filename} exists ({size} bytes)")
            rows.append((filename, url, "skipped_exists", str(size), ""))
            skipped += 1
            continue

        print(f"[download] {filename}")
        try:
            download_file(url, output_path, args.timeout)
        except (urllib.error.URLError, TimeoutError, RuntimeError, OSError) as exc:
            print(f"[failed] {filename}: {exc}", file=sys.stderr)
            rows.append((filename, url, "failed", "0", str(exc).replace("\t", " ")))
            failed += 1
            continue

        size = output_path.stat().st_size
        print(f"[done] {filename} ({size} bytes)")
        rows.append((filename, url, "downloaded", str(size), ""))
        downloaded += 1

    write_download_links(raw_dir)

    manifest_path = raw_dir / "download_manifest.tsv"
    with manifest_path.open("w") as handle:
        handle.write("filename\turl\tstatus\tsize_bytes\terror\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")

    print(
        f"Summary: {downloaded} downloaded, {skipped} skipped, "
        f"{failed} failed"
    )
    print(f"Manifest: {manifest_path}")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
