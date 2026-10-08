import csv
import logging
import os
import platform
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parent))

from scripts.feature_schema import (
    circ_id_for_key,
    parse_coordinate,
)


REQUIRED_COLUMNS = ("chrom", "start", "end", "strand")
mapped_output = Path(snakemake.output.mapped)
unmapped_output = Path(snakemake.output.unmapped)
mapped_output.parent.mkdir(parents=True, exist_ok=True)
unmapped_output.parent.mkdir(parents=True, exist_ok=True)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


logger = logging.getLogger("liftover_mouse_presence")
logger.setLevel(logging.INFO)
logger.propagate = False
log_formatter = logging.Formatter(
    "%(asctime)s %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
file_handler.setFormatter(log_formatter)
logger.handlers[:] = [file_handler]


def log_uncaught_exception(exc_type, exc_value, exc_traceback):
    logger.critical(
        "Unhandled exception",
        exc_info=(exc_type, exc_value, exc_traceback),
    )
    for handler in logger.handlers:
        handler.flush()


def log_signal(signum, _frame):
    signal_name = signal.Signals(signum).name
    logger.error("Received %s; stopping", signal_name)
    for handler in logger.handlers:
        handler.flush()
    if signum == signal.SIGINT:
        raise KeyboardInterrupt
    raise SystemExit(128 + signum)


sys.excepthook = log_uncaught_exception
signal.signal(signal.SIGINT, log_signal)
signal.signal(signal.SIGTERM, log_signal)
started_at = time.perf_counter()
logger.info("Starting mouse presence liftOver")
logger.info("Python %s on %s", platform.python_version(), platform.platform())
logger.info("Input table: %s", snakemake.input.table)
logger.info("Chain file: %s", snakemake.input.chain)
logger.info("Requested threads: %s", snakemake.threads)


def create_result_store(path):
    connection = sqlite3.connect(path, timeout=600)
    connection.execute("PRAGMA busy_timeout = 600000")
    connection.execute("PRAGMA journal_mode = WAL")
    connection.execute("PRAGMA synchronous = NORMAL")
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS mapping_results (
            chrom TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL,
            strand TEXT NOT NULL, status TEXT NOT NULL,
            target_chrom TEXT, target_start INTEGER, target_end INTEGER,
            target_strand TEXT, reason TEXT NOT NULL,
            PRIMARY KEY(chrom, start, end, strand)
        ) WITHOUT ROWID
        """
    )
    connection.commit()
    return connection


def mapped_summary_query():
    return """
        SELECT s.row_id, s.chrom, s.start, s.end, s.strand,
               SUM(CASE WHEN m.side='left' THEN 1 ELSE 0 END) AS left_count,
               SUM(CASE WHEN m.side='right' THEN 1 ELSE 0 END) AS right_count,
               MAX(CASE WHEN m.side='left' THEN m.chrom END),
               MAX(CASE WHEN m.side='left' THEN m.pos END),
               MAX(CASE WHEN m.side='left' THEN m.strand END),
               MAX(CASE WHEN m.side='right' THEN m.chrom END),
               MAX(CASE WHEN m.side='right' THEN m.pos END),
               MAX(CASE WHEN m.side='right' THEN m.strand END)
        FROM source AS s
        LEFT JOIN mapping AS m USING (row_id)
        WHERE s.is_missing=1
        GROUP BY s.row_id
        ORDER BY s.row_id
    """


def run_liftover(task):
    worker_id, source_bed, mapped_bed, rejected_bed, worker_log = task
    with worker_log.open("w", encoding="utf-8") as log:
        subprocess.run(
            [
                snakemake.params.liftover,
                str(source_bed),
                snakemake.input.chain,
                str(mapped_bed),
                str(rejected_bed),
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    return worker_id


def append_worker_logs(tasks):
    for worker_id, source_bed, _, _, worker_log in tasks:
        logger.info("liftOver worker %d log (%s)", worker_id, source_bed.name)
        if worker_log.exists():
            with worker_log.open(encoding="utf-8", errors="replace") as source:
                for line in source:
                    logger.info("[worker %d] %s", worker_id, line.rstrip("\n"))


with tempfile.TemporaryDirectory(prefix="mouse-liftover.", dir=mapped_output.parent) as work:
    work = Path(work)
    result_store_path = work / "mapping-results.sqlite"
    logger.info("Temporary workspace: %s", work)
    work_db = sqlite3.connect(work / "mouse.sqlite")
    work_db.execute("PRAGMA journal_mode = OFF")
    work_db.execute("PRAGMA synchronous = OFF")
    work_db.execute(
        "CREATE TABLE source (row_id INTEGER PRIMARY KEY, chrom TEXT, start INTEGER, "
        "end INTEGER, strand TEXT, mouse_circ_id TEXT, is_missing INTEGER)"
    )
    result_store = create_result_store(result_store_path)
    result_store.close()
    logger.info("Initialized temporary working databases")

    stage_started = time.perf_counter()
    logger.info("Loading mouse presence table")
    with open(snakemake.input.table, newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Mouse database presence table has no header")
        missing_columns = [c for c in REQUIRED_COLUMNS if c not in reader.fieldnames]
        if missing_columns:
            raise ValueError("Mouse database presence table is missing: " + ", ".join(missing_columns))
        source_batch = []
        row_count = 0
        missing_strand_count = 0
        for row_id, row in enumerate(reader):
            label = f"mouse presence line {row_id + 2}"
            if row["strand"].strip().upper() == "NA":
                normalized = parse_coordinate(
                    row["chrom"], row["start"], row["end"], "+", label
                )
                key = (*normalized[:3], "NA")
                missing_strand_count += 1
            else:
                key = parse_coordinate(
                    row["chrom"], row["start"], row["end"], row["strand"], label
                )
            source_batch.append((row_id, *key, circ_id_for_key(key), 1))
            row_count += 1
            if len(source_batch) >= 100_000:
                work_db.executemany("INSERT INTO source VALUES (?, ?, ?, ?, ?, ?, ?)", source_batch)
                source_batch.clear()
                logger.info("Loaded %d mouse presence rows", row_count)
        if source_batch:
            work_db.executemany("INSERT INTO source VALUES (?, ?, ?, ?, ?, ?, ?)", source_batch)
    logger.info(
        "Loaded %d mouse presence rows (%d with missing strand) in %.2f seconds",
        row_count,
        missing_strand_count,
        time.perf_counter() - stage_started,
    )
    if row_count == 0:
        raise ValueError("Mouse database presence table contains no records")
    stage_started = time.perf_counter()
    logger.info("Preparing records for liftOver")
    work_db.execute("ATTACH DATABASE ? AS result_db", (str(result_store_path),))
    work_db.execute(
        """
        INSERT OR IGNORE INTO result_db.mapping_results
            (chrom, start, end, strand, status, target_chrom,
             target_start, target_end, target_strand, reason)
        SELECT chrom, start, end, strand, 'unmapped', NULL, NULL, NULL,
               NULL, 'missing_source_strand'
        FROM source
        WHERE strand='NA'
        """
    )
    work_db.execute(
        """
        UPDATE source
        SET is_missing = CASE WHEN strand='NA' THEN 0 ELSE 1 END
        """
    )
    work_db.commit()
    missing_count = work_db.execute(
        "SELECT COUNT(*) FROM source WHERE is_missing=1"
    ).fetchone()[0]
    logger.info(
        "Prepared records in %.2f seconds: %d skipped, %d eligible",
        time.perf_counter() - stage_started,
        missing_strand_count,
        missing_count,
    )

    worker_count = (
        min(max(1, int(snakemake.threads)), missing_count)
        if missing_count
        else 0
    )
    tasks = []
    for worker_id in range(worker_count):
        tasks.append((
            worker_id,
            work / f"mouse.boundaries.mm10.part_{worker_id:03d}.bed",
            work / f"mouse.boundaries.hg38.part_{worker_id:03d}.bed",
            work / f"mouse.boundaries.unmapped.part_{worker_id:03d}.bed",
            work / f"liftOver.part_{worker_id:03d}.log",
        ))

    bed_handles = [task[1].open("w", encoding="utf-8") for task in tasks]
    stage_started = time.perf_counter()
    logger.info("Writing boundary BED files for %d cache misses", missing_count)
    try:
        for record_index, (row_id, chrom, start, end, strand) in enumerate(
            work_db.execute(
                "SELECT row_id,chrom,start,end,strand FROM source "
                "WHERE is_missing=1 ORDER BY row_id"
            )
        ):
            bed = bed_handles[record_index % worker_count]
            bed.write(f"{chrom}\t{start - 1}\t{start}\t{row_id}|left\t0\t{strand}\n")
            bed.write(f"{chrom}\t{end - 1}\t{end}\t{row_id}|right\t0\t{strand}\n")
            if (record_index + 1) % 100_000 == 0:
                logger.info("Wrote boundaries for %d records", record_index + 1)
    finally:
        for bed in bed_handles:
            bed.close()
    logger.info("Boundary BED files written in %.2f seconds", time.perf_counter() - stage_started)

    work_db.execute(
        "CREATE TABLE mapping (row_id INTEGER, side TEXT, chrom TEXT, pos INTEGER, "
        "strand TEXT, UNIQUE(row_id,side,chrom,pos,strand))"
    )
    if missing_count:
        stage_started = time.perf_counter()
        logger.info("Running liftOver with %d workers", worker_count)
        try:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                completed_workers = list(executor.map(run_liftover, tasks))
            if len(completed_workers) != worker_count:
                raise RuntimeError("Not all liftOver workers completed")
        except Exception:
            append_worker_logs(tasks)
            raise
        append_worker_logs(tasks)
        logger.info("liftOver workers completed in %.2f seconds", time.perf_counter() - stage_started)

        stage_started = time.perf_counter()
        logger.info("Importing mapped liftOver boundaries")
        batch = []
        mapped_boundary_count = 0
        for worker_id, _, mapped_bed, _, _ in tasks:
            with mapped_bed.open(encoding="utf-8") as handle:
                for line_number, line in enumerate(handle, start=1):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 6:
                        raise ValueError(
                            f"Mapped BED worker {worker_id} line {line_number} "
                            "has fewer than six columns"
                        )
                    name = fields[3].rsplit("|", 1)
                    if len(name) != 2 or name[1] not in {"left", "right"}:
                        raise ValueError(
                            f"Mapped BED worker {worker_id} line {line_number} "
                            "has invalid name"
                        )
                    try:
                        row_id = int(name[0])
                        start0 = int(fields[1])
                        end = int(fields[2])
                    except ValueError as exc:
                        raise ValueError(
                            f"Mapped BED worker {worker_id} line {line_number} "
                            "has invalid coordinates"
                        ) from exc
                    if end - start0 != 1:
                        raise ValueError(
                            f"Mapped BED worker {worker_id} line {line_number} "
                            "is not one base"
                        )
                    batch.append((row_id, name[1], fields[0], end, fields[5]))
                    mapped_boundary_count += 1
                    if len(batch) >= 100_000:
                        work_db.executemany(
                            "INSERT OR IGNORE INTO mapping VALUES (?, ?, ?, ?, ?)",
                            batch,
                        )
                        batch.clear()
                        logger.info("Imported %d mapped boundaries", mapped_boundary_count)
        if batch:
            work_db.executemany(
                "INSERT OR IGNORE INTO mapping VALUES (?, ?, ?, ?, ?)", batch
            )
        logger.info(
            "Imported %d mapped boundaries in %.2f seconds",
            mapped_boundary_count,
            time.perf_counter() - stage_started,
        )
    else:
        logger.info("No eligible records; skipping liftOver")
    work_db.commit()

    stage_started = time.perf_counter()
    logger.info("Classifying liftOver results")
    work_db.execute("BEGIN IMMEDIATE")
    classified_count = 0
    for row in work_db.execute(mapped_summary_query()):
        row_id, chrom, start, end, strand, nl, nr, lchrom, lpos, lstrand, rchrom, rpos, rstrand = row
        status = "mapped"
        reason = ""
        target = (None, None, None, None)
        if nl != 1 or nr != 1:
            status = "unmapped"
            reason = "unmapped_endpoint" if nl == 0 or nr == 0 else "ambiguous_endpoint"
        elif lchrom != rchrom:
            status = "unmapped"
            reason = "different_target_chromosomes"
        elif lstrand not in {"+", "-"} or lstrand != rstrand:
            status = "unmapped"
            reason = "inconsistent_target_strand"
        else:
            target_start, target_end = sorted((lpos, rpos))
            target = (lchrom, target_start, target_end, lstrand)
        work_db.execute(
            "INSERT OR IGNORE INTO result_db.mapping_results VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (chrom, start, end, strand, status, *target, reason),
        )
        classified_count += 1
        if classified_count % 100_000 == 0:
            logger.info("Classified %d liftOver records", classified_count)
    work_db.commit()
    logger.info(
        "Classified %d liftOver records in %.2f seconds",
        classified_count,
        time.perf_counter() - stage_started,
    )

    stage_started = time.perf_counter()
    logger.info("Writing mapped and unmapped output tables")
    mapped_handle = tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=mapped_output.parent,
        prefix=f".{mapped_output.name}.", suffix=".tmp", delete=False,
    )
    unmapped_handle = tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=unmapped_output.parent,
        prefix=f".{unmapped_output.name}.", suffix=".tmp", delete=False,
    )
    mapped_tmp = Path(mapped_handle.name)
    unmapped_tmp = Path(unmapped_handle.name)
    try:
        mapped_writer = csv.writer(mapped_handle, delimiter="\t", lineterminator="\n")
        unmapped_writer = csv.writer(unmapped_handle, delimiter="\t", lineterminator="\n")
        mapped_writer.writerow([
            "chrom", "start", "end", "strand", "original_chrom", "original_start",
            "original_end", "original_strand", "mouse_circ_id",
        ])
        unmapped_writer.writerow(["chrom", "start", "end", "strand", "mouse_circ_id", "reason"])
        query = """
            SELECT s.chrom,s.start,s.end,s.strand,s.mouse_circ_id,
                   c.status,c.target_chrom,c.target_start,c.target_end,
                   c.target_strand,c.reason
            FROM source AS s
            LEFT JOIN result_db.mapping_results AS c
              ON c.chrom=s.chrom AND c.start=s.start
             AND c.end=s.end AND c.strand=s.strand
            ORDER BY s.row_id
        """
        output_count = 0
        mapped_count = 0
        unmapped_count = 0
        for result_row in work_db.execute(query):
            chrom, start, end, strand, mouse_id, status, target_chrom, target_start, target_end, target_strand, reason = result_row
            if status is None:
                raise ValueError(f"Missing mouse liftOver result for {result_row[:5]}")
            if status == "mapped":
                mapped_writer.writerow([
                    target_chrom, target_start, target_end, target_strand,
                    chrom, start, end, strand, mouse_id,
                ])
                mapped_count += 1
            else:
                unmapped_writer.writerow([chrom, start, end, strand, mouse_id, reason])
                unmapped_count += 1
            output_count += 1
            if output_count % 100_000 == 0:
                logger.info("Wrote %d output records", output_count)
        for handle in (mapped_handle, unmapped_handle):
            handle.flush()
            os.fsync(handle.fileno())
            handle.close()
        os.replace(mapped_tmp, mapped_output)
        os.replace(unmapped_tmp, unmapped_output)
        logger.info(
            "Wrote %d mapped and %d unmapped records in %.2f seconds",
            mapped_count,
            unmapped_count,
            time.perf_counter() - stage_started,
        )
    except Exception:
        mapped_handle.close()
        unmapped_handle.close()
        mapped_tmp.unlink(missing_ok=True)
        unmapped_tmp.unlink(missing_ok=True)
        raise
    finally:
        work_db.close()

logger.info("Mouse presence liftOver completed in %.2f seconds", time.perf_counter() - started_at)
for handler in logger.handlers:
    handler.flush()
