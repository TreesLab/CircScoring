#!/usr/bin/env python3
"""Download Circ2Traits raw data files from GitHub."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import tempfile
from datetime import date
from pathlib import Path
import urllib.error
import urllib.request


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"
BASE_URL = "https://raw.githubusercontent.com/shaoli86/circ2Traits/master"

DOWNLOADS = [
    {
        "file": "circRNA_info_all.txt",
        "url": f"{BASE_URL}/circRNA_info_all.txt",
    },
    {
        "file": "all_disease_set.txt",
        "url": f"{BASE_URL}/all_disease_set.txt",
    },
    {
        "file": "gwascatalog.txt",
        "url": f"{BASE_URL}/gwascatalog.txt",
    },
]

MANIFEST_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "line_count",
    "validation_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download Circ2Traits raw text files from GitHub.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty files.")
    parser.add_argument("--timeout", type=int, default=180, help="Network timeout in seconds per file. Default: 180")
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads without downloading files.")
    return parser.parse_args()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, output: Path, *, force: bool, timeout: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent, delete=False) as handle:
        tmp_path = Path(handle.name)

    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with tmp_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        if tmp_path.stat().st_size == 0:
            raise RuntimeError(f"downloaded file is empty: {url}")
        tmp_path.replace(output)
        output.chmod(0o644)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    return "downloaded"


def validate_text_file(path: Path) -> tuple[int, str]:
    with path.open("rb") as handle:
        prefix = handle.read(512)
    if not prefix:
        raise ValueError(f"empty file: {path}")
    lowered = prefix.lower()
    if b"<html" in lowered or b"<!doctype html" in lowered:
        raise ValueError(f"downloaded file looks like HTML, not raw text: {path}")

    with path.open(errors="replace") as handle:
        line_count = sum(1 for _line in handle)
    if line_count == 0:
        raise ValueError(f"file has no readable lines: {path}")
    return line_count, "text_ok"


def write_manifest(raw_dir: Path, rows: list[dict[str, str]]) -> None:
    path = raw_dir / "download_manifest.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()

    if args.dry_run:
        for item in DOWNLOADS:
            print(f"{item['url']}\t{raw_dir / item['file']}")
        return 0

    rows: list[dict[str, str]] = []
    failures: list[str] = []
    today = date.today().isoformat()

    for item in DOWNLOADS:
        path = raw_dir / item["file"]
        try:
            status = download_file(item["url"], path, force=args.force, timeout=args.timeout)
            line_count, validation_status = validate_text_file(path)
            print(f"[{status}] {item['file']} lines={line_count}")
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            line_count = "NA"
            validation_status = "failed"
            failures.append(f"{item['file']}: {exc}")
            print(f"[failed] {item['file']}: {exc}", file=sys.stderr)

        rows.append(
            {
                "download_date": today,
                "file": item["file"],
                "source_url": item["url"],
                "status": status,
                "size_bytes": str(path.stat().st_size if path.exists() else 0),
                "sha256": sha256sum(path) if path.exists() and path.stat().st_size > 0 else "NA",
                "line_count": str(line_count),
                "validation_status": validation_status,
            }
        )

    write_manifest(raw_dir, rows)
    print(f"manifest: {raw_dir / 'download_manifest.tsv'}")

    if failures:
        print(f"failed_downloads\t{len(failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
