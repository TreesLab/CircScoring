#!/usr/bin/env python3
"""Download m6A2Circ human circRNA raw data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import tempfile
from datetime import date
from pathlib import Path
import urllib.request


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"

RAW_FILE = {
    "file": "huamn_circRNAs_Information.zip",
    "url": "http://m6a2circ.canceromics.org/api/download/huamn_circRNAs_Information.zip",
    "zip_member": "huamn_circRNAs_Information.txt",
}

MANIFEST_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "zip_member",
    "zip_validation_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download m6A2Circ human circRNA raw zip.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty file.")
    parser.add_argument("--timeout", type=int, default=300, help="Network timeout in seconds. Default: 300")
    parser.add_argument("--dry-run", action="store_true", help="List planned download without downloading.")
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


def validate_zip(path: Path, member: str) -> str:
    list_result = subprocess.run(
        ["unzip", "-Z", "-1", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    names = list_result.stdout.splitlines()
    if member not in names:
        raise ValueError(f"{path} does not contain required member {member!r}; found: {names}")

    subprocess.run(
        ["unzip", "-t", str(path), member],
        check=True,
        capture_output=True,
        text=True,
    )
    return "ok"


def write_manifest(raw_dir: Path, row: dict[str, str]) -> None:
    path = raw_dir / "download_manifest.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = raw_dir / RAW_FILE["file"]

    if args.dry_run:
        print(f"{RAW_FILE['url']}\t{output}")
        return 0

    status = download_file(RAW_FILE["url"], output, force=args.force, timeout=args.timeout)
    zip_status = validate_zip(output, RAW_FILE["zip_member"])
    row = {
        "download_date": date.today().isoformat(),
        "file": RAW_FILE["file"],
        "source_url": RAW_FILE["url"],
        "status": status,
        "size_bytes": str(output.stat().st_size),
        "sha256": sha256sum(output),
        "zip_member": RAW_FILE["zip_member"],
        "zip_validation_status": zip_status,
    }
    write_manifest(raw_dir, row)
    print(f"{status}: {output}")
    print(f"zip_validation_status: {zip_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
