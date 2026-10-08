#!/usr/bin/env python3
"""Download riboCIRC raw data archives."""

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

DOWNLOADS = [
    {
        "file": "All_annotation-guided_ribo-circRNA.rar",
        "url": "http://www.ribocirc.com/download/All_annotation-guided_ribo-circRNA.rar",
    },
    {
        "file": "All_context-specific_ribo-circRNAs.rar",
        "url": "http://www.ribocirc.com/download/All_context-specific_ribo-circRNAs.rar",
    },
    {
        "file": "All_finite-ORF-encoded_peptides.rar",
        "url": "http://www.ribocirc.com/download/All_finite-ORF-encoded_peptides.rar",
    },
    {
        "file": "Cross-species_conserved_ribo-cirRNAs.rar",
        "url": "http://www.ribocirc.com/download/Cross-species_conserved_ribo-cirRNAs.rar",
    },
    {
        "file": "literature-reported_trans-circRNAs.rar",
        "url": "http://www.ribocirc.com/download/literature-reported_trans-circRNAs.rar",
    },
]

MANIFEST_COLUMNS = [
    "download_date",
    "file",
    "source_url",
    "status",
    "size_bytes",
    "sha256",
    "rar_validation_status",
    "rar_members",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download riboCIRC raw RAR archives.")
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help=f"Raw output directory. Default: {DEFAULT_RAW_DIR}")
    parser.add_argument("--force", action="store_true", help="Re-download existing non-empty files.")
    parser.add_argument("--timeout", type=int, default=300, help="Network timeout in seconds. Default: 300")
    parser.add_argument("--dry-run", action="store_true", help="List planned downloads without downloading.")
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


def validate_rar(path: Path) -> tuple[str, list[str]]:
    result = subprocess.run(
        ["bsdtar", "-tf", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    members = [line for line in result.stdout.splitlines() if line]
    if not members:
        raise ValueError(f"RAR archive has no readable members: {path}")
    return "ok", members


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
    for item in DOWNLOADS:
        output = raw_dir / item["file"]
        status = download_file(item["url"], output, force=args.force, timeout=args.timeout)
        validation_status, members = validate_rar(output)
        rows.append(
            {
                "download_date": date.today().isoformat(),
                "file": item["file"],
                "source_url": item["url"],
                "status": status,
                "size_bytes": str(output.stat().st_size),
                "sha256": sha256sum(output),
                "rar_validation_status": validation_status,
                "rar_members": ";".join(members),
            }
        )
        print(f"{status}: {output}")
        print(f"rar_validation_status: {validation_status}")

    write_manifest(raw_dir, rows)
    print(f"manifest: {raw_dir / 'download_manifest.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
