#!/usr/bin/env python3
"""Download CircAge human circRNA raw main files."""

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


csv.field_size_limit(sys.maxsize)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"

RAW_FILES = [
    {
        "tissue": "Aortic",
        "filename": "human_Aortic_main.txt",
        "url": "https://file.kiz.ac.cn/circage/circAge/human/human_Aortic_main.txt",
    },
    {
        "tissue": "Lung",
        "filename": "human_Lung_main.txt",
        "url": "https://file.kiz.ac.cn/circage/circAge/human/human_Lung_main.txt",
    },
    {
        "tissue": "Skin",
        "filename": "human_Skin_main.txt",
        "url": "https://file.kiz.ac.cn/circage/circAge/human/human_Skin_main.txt",
    },
    {
        "tissue": "Umbilical",
        "filename": "human_Umbilical_main.txt",
        "url": "https://file.kiz.ac.cn/circage/circAge/human/human_Umbilical_main.txt",
    },
    {
        "tissue": "white_blood_cell",
        "filename": "human_white_blood_cell_main.txt",
        "url": "https://file.kiz.ac.cn/circage/circAge/human/human_white_blood_cell_main.txt",
    },
]

RECORD_COLUMNS = [
    "download_date",
    "tissue",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "row_count",
    "validation_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download CircAge human circRNA raw main files.")
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
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    return "downloaded"


def validate_raw_table(path: Path) -> int:
    row_count = 0
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 11:
                raise ValueError(f"{path} line {line_number} has {len(row)} fields, expected 11")
            row_count += 1
    if row_count == 0:
        raise ValueError(f"empty raw table: {path}")
    return row_count


def write_record(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    raw_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        for item in RAW_FILES:
            print(f"{item['tissue']}\t{item['url']}\t{raw_dir / item['filename']}")
        return 0

    rows: list[dict[str, str]] = []
    failures: list[str] = []
    today = date.today().isoformat()

    for item in RAW_FILES:
        path = raw_dir / item["filename"]
        try:
            status = download_file(item["url"], path, force=args.force, timeout=args.timeout)
            row_count = validate_raw_table(path)
            validation_status = "tsv_11_columns_ok"
            print(f"[{status}] {item['filename']} rows={row_count}")
        except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
            status = f"failed: {exc}"
            row_count = "NA"
            validation_status = "failed"
            failures.append(f"{item['filename']}: {exc}")
            print(f"[failed] {item['filename']}: {exc}", file=sys.stderr)

        rows.append(
            {
                "download_date": today,
                "tissue": item["tissue"],
                "file": item["filename"],
                "source_url": item["url"],
                "status": status,
                "size_bytes": str(path.stat().st_size if path.exists() else 0),
                "sha256": sha256sum(path) if path.exists() and path.stat().st_size > 0 else "NA",
                "row_count": str(row_count),
                "validation_status": validation_status,
            }
        )

    write_record(raw_dir / "download_record.tsv", rows)

    if failures:
        print(f"failed_downloads\t{len(failures)}", file=sys.stderr)
        return 1

    print(f"downloaded_or_validated_files\t{len(RAW_FILES)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
