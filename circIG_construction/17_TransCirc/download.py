#!/usr/bin/env python3
"""Download the official TransCirc metadata table."""

from __future__ import annotations

import argparse
import shutil
import tempfile
import urllib.request
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
DEFAULT_URL = "https://www.biosino.org/Transcirc/userfiles/download/transcirc_metadata.tsv"
OUTPUT_FILENAME = "transcirc_metadata.tsv"
USER_AGENT = "data-downloader/1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download the official TransCirc metadata TSV.")
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Download URL. Default: {DEFAULT_URL}")
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument("--timeout", type=int, default=300, help="Network timeout in seconds. Default: 300")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing non-empty raw file.")
    return parser.parse_args()


def download_file(url: str, output: Path, *, timeout: int, force: bool) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_path = Path(handle.name)

    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with temporary_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)

        if temporary_path.stat().st_size == 0:
            raise RuntimeError(f"Downloaded file is empty: {url}")

        temporary_path.replace(output)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    return "downloaded"


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = raw_dir / OUTPUT_FILENAME
    status = download_file(args.url, output, timeout=args.timeout, force=args.force)

    print(f"{status}: {output}")
    print(f"size_bytes: {output.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
