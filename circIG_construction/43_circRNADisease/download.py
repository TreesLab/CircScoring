#!/usr/bin/env python3
"""Download the circRNADisease v3 circRNA detail table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile
from datetime import date
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_RAW_DIR = BASE_DIR / "raw"
USER_AGENT = "data-downloader/1.0"
FILENAME = "circRNADisease_V3_circrna_details.xlsx"
DOWNLOAD_URL = (
    "https://cgga.org.cn/circRNADisease/download/"
    "circRNADisease_V3_circrna_details.xlsx"
)
REQUIRED_XLSX_MEMBERS = {
    "[Content_Types].xml",
    "xl/workbook.xml",
    "xl/worksheets/sheet1.xml",
}
MANIFEST_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "validation_status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download the file even when a non-empty output already exists.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=180,
        help="Network timeout in seconds. Default: 180",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the planned download without writing files.",
    )
    return parser.parse_args()


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_xlsx(path: Path) -> str:
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty XLSX file: {path}")

    with path.open("rb") as handle:
        if handle.read(4) != b"PK\x03\x04":
            raise ValueError(f"file does not have an XLSX ZIP signature: {path}")

    try:
        with zipfile.ZipFile(path) as archive:
            members = set(archive.namelist())
            missing = sorted(REQUIRED_XLSX_MEMBERS - members)
            if missing:
                raise ValueError(
                    f"XLSX is missing required member(s): {', '.join(missing)}"
                )
            corrupt_member = archive.testzip()
    except zipfile.BadZipFile as exc:
        raise ValueError(f"invalid XLSX ZIP archive: {path}") from exc

    if corrupt_member is not None:
        raise ValueError(f"corrupt XLSX ZIP member: {corrupt_member}")
    return "xlsx_ok"


def download_file(url: str, output: Path, *, force: bool, timeout: int) -> str:
    if output.exists() and output.stat().st_size > 0 and not force:
        validate_xlsx(output)
        return "skipped_existing_file"

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
        delete=False,
    ) as handle:
        tmp_path = Path(handle.name)

    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=timeout) as response:
            with tmp_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        validate_xlsx(tmp_path)
        tmp_path.replace(output)
        output.chmod(0o644)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise

    return "downloaded"


def write_manifest(raw_dir: Path, row: dict[str, str]) -> Path:
    path = raw_dir / "download_manifest.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(row)
    return path


def main() -> int:
    args = parse_args()
    raw_dir = args.raw_dir.expanduser().resolve()
    output = raw_dir / FILENAME

    if args.dry_run:
        print(f"{DOWNLOAD_URL}\t{output}")
        return 0

    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        status = download_file(
            DOWNLOAD_URL,
            output,
            force=args.force,
            timeout=args.timeout,
        )
        validation_status = validate_xlsx(output)
    except (OSError, urllib.error.URLError, ValueError, RuntimeError) as exc:
        print(f"[failed] {FILENAME}: {exc}", file=sys.stderr)
        return 1

    manifest = write_manifest(
        raw_dir,
        {
            "download_date": date.today().isoformat(),
            "file": FILENAME,
            "source_url": DOWNLOAD_URL,
            "status": status,
            "size_bytes": str(output.stat().st_size),
            "sha256": sha256sum(output),
            "validation_status": validation_status,
        },
    )
    print(f"[{status}] {output} ({output.stat().st_size} bytes)")
    print(f"validation_status: {validation_status}")
    print(f"manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
