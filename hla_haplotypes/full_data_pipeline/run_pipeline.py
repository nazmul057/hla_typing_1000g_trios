"""HLA typing pipeline driver for 1000 Genomes Project trios.

Streams HLA-relevant reads from remote CRAMs and runs SpecHLA inside a pinned
Docker image. Designed to be run from a single command:

    python run_pipeline.py

On startup, the pipeline verifies its setup (working directories, reference
genome, contig lists). If anything is missing, setup is performed before
sample processing begins.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

# --- Configuration ----------------------------------------------------------
# Edit ROOT_DIR_HOST to point at your pipeline working directory. Everything
# the pipeline reads or writes lives under this path. Relative paths are fine;
# they are resolved to absolute form before being passed to Docker.
ROOT_DIR_HOST = Path("hla_haplotypes/full_data_pipeline").resolve()

DOCKER_IMAGE = "spechla_pipeline:v1"

REFERENCE_FASTA_NAME = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
REFERENCE_FASTA_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "GRCh38_reference_genome/" + REFERENCE_FASTA_NAME
)

# --- Directory layout (host-side) -------------------------------------------
# All subdirectories the pipeline reads or writes on the host machine. The
# `_host` suffix marks these as paths in the host's filesystem (not visible
# to containers unless explicitly mounted).
REFERENCE_DIR_HOST = ROOT_DIR_HOST / "reference"
INPUT_DIR_HOST     = ROOT_DIR_HOST / "input"
OUTPUTS_DIR_HOST   = ROOT_DIR_HOST / "outputs"
MARKERS_DIR_HOST   = ROOT_DIR_HOST / "work" / "markers"
LOGS_DIR_HOST      = ROOT_DIR_HOST / "work" / "logs"

# Tuple form used by setup() to create the layout in one pass.
SUBDIRS_HOST = (
    REFERENCE_DIR_HOST,
    INPUT_DIR_HOST,
    OUTPUTS_DIR_HOST,
    MARKERS_DIR_HOST,
    LOGS_DIR_HOST,
)

# --- Mount paths (container-side) -------------------------------------------
# Paths inside the per-sample container. The `_ctr` suffix marks these as
# paths the container sees, not host paths. Changing any of these requires
# updating the corresponding -v mount in process_sample.
REFERENCE_DIR_CTR    = "/data/reference"
OUTPUT_DIR_CTR       = "/data/output"
REFERENCE_FASTA_CTR  = f"{REFERENCE_DIR_CTR}/{REFERENCE_FASTA_NAME}"

# Contig list filenames written into INPUT_DIR_HOST. These are derived once
# from the reference .fai and reused by every sample.
HLA_CONTIGS_FILE = "hla_contigs.txt"
CHR6_ALT_CONTIGS_FILE = "chr6_alt_contigs.txt"

# Sanity-check ranges for the GRCh38 1000 Genomes analysis-set reference.
# Expected: ~525 HLA decoy contigs, ~16 chr6 alt contigs. Anything far outside
# these ranges suggests the wrong reference was downloaded.
EXPECTED_HLA_CONTIG_COUNT_RANGE = (500, 600)
EXPECTED_CHR6_ALT_CONTIG_COUNT_RANGE = (10, 30)

# Threads passed to samtools (-@) and SpecHLA (-j) inside one container.
# Kept modest for now; will be tuned alongside worker count when we add
# parallel sample execution.
THREADS_PER_SAMPLE = 4


# --- Setup ------------------------------------------------------------------
def setup() -> None:
    """Prepare the pipeline working environment under ROOT_DIR_HOST.

    Verifies Docker is reachable, creates the standard subdirectories,
    downloads and indexes the reference FASTA if missing, and derives the
    HLA / chr6-alt contig lists from the reference index. Idempotent: safe
    to call on every pipeline run.
    """
    print(f"[setup] Root: {ROOT_DIR_HOST}")
    _check_docker_running()
    _create_subdirs()
    _ensure_reference()
    _ensure_contig_lists()
    print("[setup] Done.")


def _check_docker_running() -> None:
    """Bail out early with a clear message if the Docker daemon isn't reachable.

    Catches the common "forgot to start Docker Desktop" mistake before it
    surfaces as a confusing per-step container failure.
    """
    try:
        result = subprocess.run(
            ["docker", "info"], capture_output=True, text=True, timeout=10,
        )
    except FileNotFoundError:
        print("[setup] ERROR: docker CLI not found on PATH.", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print("[setup] ERROR: `docker info` timed out; is the daemon hung?", file=sys.stderr)
        sys.exit(1)
    if result.returncode != 0:
        print("[setup] ERROR: Docker is not reachable. Is Docker Desktop running?",
              file=sys.stderr)
        sys.exit(1)


def _create_subdirs() -> None:
    """Create the standard subdirectory tree if any pieces are missing."""
    for subdir_host in SUBDIRS_HOST:
        subdir_host.mkdir(parents=True, exist_ok=True)


def _ensure_reference() -> None:
    """Download and index the reference FASTA if either is absent.

    Both download and indexing run inside the Docker image, so the host needs
    nothing beyond Docker and Python. The reference directory is mounted rw
    during setup; at pipeline runtime it is mounted ro by the per-sample
    workers.
    """
    fasta_path_host = REFERENCE_DIR_HOST / REFERENCE_FASTA_NAME
    fai_path_host = REFERENCE_DIR_HOST / (REFERENCE_FASTA_NAME + ".fai")

    if fasta_path_host.exists() and fai_path_host.exists():
        print(f"[setup] Reference already present: {fasta_path_host.name} (+ .fai)")
        return

    mount = f"{REFERENCE_DIR_HOST.as_posix()}:{REFERENCE_DIR_CTR}"

    if not fasta_path_host.exists():
        print(f"[setup] Downloading reference FASTA (~3 GB) to {fasta_path_host}")
        # `wget -c` resumes a partial download instead of restarting from zero
        # if a previous attempt was interrupted.
        download_cmd = [
            "docker", "run", "--rm",
            "-v", mount,
            "-w", REFERENCE_DIR_CTR,
            DOCKER_IMAGE,
            "wget", "-c", REFERENCE_FASTA_URL,
        ]
        _run_or_die(download_cmd, "reference download failed")

    print(f"[setup] Building FASTA index: {fai_path_host.name}")
    faidx_cmd = [
        "docker", "run", "--rm",
        "-v", mount,
        DOCKER_IMAGE,
        "samtools", "faidx", REFERENCE_FASTA_CTR,
    ]
    _run_or_die(faidx_cmd, "samtools faidx failed")


def _ensure_contig_lists() -> None:
    """Derive HLA decoy and chr6 alt contig lists from the reference .fai.

    The .fai is a tab-separated index with one line per sequence in the FASTA;
    the first column is the contig name. Filtering it gives the same lists
    that the pipeline doc derives from a CRAM header, with no network call:

        cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \\
            | awk '/^HLA-/'         > hla_contigs.txt
        cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \\
            | awk '/^chr6_.*_alt$/' > chr6_alt_contigs.txt

    This works because every CRAM in the high-coverage 1kGP cohort was aligned
    against this exact reference, so the CRAM @SQ contigs and the FASTA
    contigs are identical by construction.
    """
    hla_path_host = INPUT_DIR_HOST / HLA_CONTIGS_FILE
    chr6_alt_path_host = INPUT_DIR_HOST / CHR6_ALT_CONTIGS_FILE

    if hla_path_host.exists() and chr6_alt_path_host.exists():
        print(f"[setup] Contig lists already present in {INPUT_DIR_HOST}")
        return

    fai_path_host = REFERENCE_DIR_HOST / (REFERENCE_FASTA_NAME + ".fai")
    if not fai_path_host.exists():
        print(f"[setup] ERROR: expected {fai_path_host} to exist by now", file=sys.stderr)
        sys.exit(1)

    hla_contigs: list[str] = []
    chr6_alt_contigs: list[str] = []
    with fai_path_host.open() as f:
        for line in f:
            name = line.split("\t", 1)[0]
            if name.startswith("HLA-"):
                hla_contigs.append(name)
            elif name.startswith("chr6_") and name.endswith("_alt"):
                chr6_alt_contigs.append(name)

    _check_count("HLA decoy", len(hla_contigs), EXPECTED_HLA_CONTIG_COUNT_RANGE)
    _check_count("chr6 alt", len(chr6_alt_contigs), EXPECTED_CHR6_ALT_CONTIG_COUNT_RANGE)

    hla_path_host.write_text("\n".join(hla_contigs) + "\n")
    chr6_alt_path_host.write_text("\n".join(chr6_alt_contigs) + "\n")
    print(f"[setup] Wrote {len(hla_contigs)} HLA contigs to {hla_path_host.name}")
    print(f"[setup] Wrote {len(chr6_alt_contigs)} chr6 alt contigs to {chr6_alt_path_host.name}")


def _check_count(label: str, count: int, expected_range: tuple[int, int]) -> None:
    """Warn loudly if a contig count is outside its expected range.

    A wildly off count usually means the wrong reference was downloaded.
    Non-fatal: we still write the file, but the user should investigate.
    """
    low, high = expected_range
    if not (low <= count <= high):
        print(
            f"[setup] WARNING: {label} contig count is {count}, "
            f"expected {low}-{high}. Wrong reference?",
            file=sys.stderr,
        )


def _read_contig_list(path: Path) -> list[str]:
    """Read a contig list file written by setup, one contig per line."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def _run_or_die(cmd: list[str], error_message: str) -> None:
    """Run a subprocess, streaming output to the console; exit on failure.

    Used during setup, where failures are fatal to the whole pipeline run.
    Per-sample failures use _run_step instead, which logs and raises.
    """
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[setup] ERROR: {error_message} (exit {exc.returncode})", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"[setup] ERROR: command not found: {cmd[0]}", file=sys.stderr)
        print("[setup] Is Docker installed and on PATH?", file=sys.stderr)
        sys.exit(1)


# --- Per-sample processing --------------------------------------------------
def process_sample(
    sample_id: str,
    cram_url: str,
    population: str,
    chr6_alt_contigs: list[str],
    hla_contigs: list[str],
) -> None:
    """Run the full 8-step HLA typing pipeline for one sample.

    The contig lists are loaded once at startup and passed in by the caller;
    they do not change between samples in a single pipeline run.

    The eight steps (slice, sort, index, extract, sanity-check, type, rename
    placeholder, cleanup) each run in their own short-lived container. Every
    container's stdout+stderr is appended to work/logs/<sample_id>.log along
    with timestamps; high-level progress lines are echoed to the console.

    Concurrency-safe via three marker files in work/markers/:
        <sample_id>.running.json  — created exclusively before processing
        <sample_id>.done.json     — written on success (full slice metadata)
        <sample_id>.failed.json   — written on failure (error info)

    On entry: if .done.json exists, the sample is skipped. Stale markers
    from previous runs (.failed.json, .running.json) are NOT handled here;
    run_pipeline's pre-flight pass wipes them before any worker starts.
    Calling process_sample directly without that pre-flight on a sample
    with stale markers will produce a "another worker is processing" skip
    or a duplicate output write.

    After SpecHLA exits, hla.result.g.group.txt is checked to confirm at
    least one allele was actually called; an all-dashes result is treated
    as a silent SpecHLA failure and converted into an explicit .failed.json.
    """
    if population not in ("Asian", "Caucasian", "Black", "Unknown"):
        raise ValueError(
            f"population must be one of Asian/Caucasian/Black/Unknown, got {population!r}"
        )

    done_path_host = MARKERS_DIR_HOST / f"{sample_id}.done.json"
    failed_path_host = MARKERS_DIR_HOST / f"{sample_id}.failed.json"
    running_path_host = MARKERS_DIR_HOST / f"{sample_id}.running.json"
    sample_out_dir_host = OUTPUTS_DIR_HOST / sample_id
    log_path_host = LOGS_DIR_HOST / f"{sample_id}.log"

    if done_path_host.exists():
        print(f"[{sample_id}] Already done; skipping.")
        return

    # Claim the sample with an exclusive create. Within a single pipeline run
    # this protects against two workers racing on the same sample id (e.g. a
    # duplicate entry in the input list). Stale markers from previous runs
    # have already been wiped during run_pipeline's pre-flight pass.
    sample_out_dir_host.mkdir(parents=True, exist_ok=True)
    try:
        with running_path_host.open("x") as f:
            json.dump(
                {"started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")},
                f,
            )
    except FileExistsError:
        print(f"[{sample_id}] Another worker is processing this sample; skipping.")
        return

    # Container sees only the reference (ro) and its own per-sample output
    # (rw). The contig lists are read host-side and spliced into the slice
    # command's argv directly, so no input dir mount is needed.
    mounts = [
        "-v", f"{REFERENCE_DIR_HOST.as_posix()}:{REFERENCE_DIR_CTR}:ro",
        "-v", f"{sample_out_dir_host.as_posix()}:{OUTPUT_DIR_CTR}",
    ]

    # Filenames inside OUTPUT_DIR_CTR, kept consistent with the pipeline doc.
    bam_unsorted_name = f"{sample_id}.hla_full.bam"
    bam_sorted_name = f"{sample_id}.hla_full.coordsorted.bam"
    bam_sorted_idx_name = f"{sample_id}.hla_full.coordsorted.bam.bai"
    fq1_relpath = f"hla_reads/{sample_id}_extract_1.fq.gz"
    fq2_relpath = f"hla_reads/{sample_id}_extract_2.fq.gz"

    started = datetime.now(timezone.utc)
    durations: dict[str, float] = {}
    print(f"[{sample_id}] Start (population={population})")

    try:
        with log_path_host.open("a") as log:
            _log_header(log, sample_id, cram_url, population, DOCKER_IMAGE)

            # Step 1 (was Step 3 in the bash doc): slice HLA-relevant reads
            # from the remote CRAM. Contig names are passed as individual
            # argv elements rather than via shell expansion; this avoids
            # platform-dependent shell quoting issues (the reason the
            # previous bash -c version failed on Windows). samtools is the
            # container's command directly — no shell in the loop.
            slice_argv = [
                "samtools", "view", "-b",
                "-@", str(THREADS_PER_SAMPLE),
                "-T", REFERENCE_FASTA_CTR,
                "-o", f"{OUTPUT_DIR_CTR}/{bam_unsorted_name}",
                cram_url,
                "chr6:25000000-35000000",
                *chr6_alt_contigs,
                *hla_contigs,
                "*",
            ]
            _run_step("slice", sample_id, log, durations, mounts, slice_argv)

            # Step 2: coordinate-sort.
            _run_step(
                "sort", sample_id, log, durations, mounts,
                ["samtools", "sort", "-@", str(THREADS_PER_SAMPLE),
                 "-o", f"{OUTPUT_DIR_CTR}/{bam_sorted_name}",
                 f"{OUTPUT_DIR_CTR}/{bam_unsorted_name}"],
            )

            # Step 3: index the sorted BAM.
            _run_step(
                "index", sample_id, log, durations, mounts,
                ["samtools", "index", f"{OUTPUT_DIR_CTR}/{bam_sorted_name}"],
            )

            # Step 4: SpecHLA's HLA read extractor produces paired FASTQs.
            _run_step(
                "extract_hla_reads", sample_id, log, durations, mounts,
                ["spechla-extract-hla-reads",
                 "-s", sample_id,
                 "-b", f"{OUTPUT_DIR_CTR}/{bam_sorted_name}",
                 "-r", "hg38",
                 "-o", f"{OUTPUT_DIR_CTR}/hla_reads/"],
            )

            # Step 5: sanity-check the FASTQs are non-empty before paying for
            # typing.
            _sanity_check_fastqs(sample_id, sample_out_dir_host, log)

            # Step 6: SpecHLA typing.
            _run_step(
                "spechla_type", sample_id, log, durations, mounts,
                ["spechla",
                 "-n", sample_id,
                 "-1", f"{OUTPUT_DIR_CTR}/{fq1_relpath}",
                 "-2", f"{OUTPUT_DIR_CTR}/{fq2_relpath}",
                 "-p", population,
                 "-o", f"{OUTPUT_DIR_CTR}/spechla_output/",
                 "-j", str(THREADS_PER_SAMPLE)],
            )

            # SpecHLA can exit 0 even when typing has produced nothing useful
            # (e.g., the realignment stage failed silently and every gene came
            # back empty). Read the result file and verify it has at least one
            # real allele call before declaring success.
            _validate_hla_result(sample_id, sample_out_dir_host, log)

            # Step 7: capture sizes, then cleanup intermediates. Keep FASTQs
            # and SpecHLA results; remove the BAMs (large, regenerable).
            bam_size = (sample_out_dir_host / bam_unsorted_name).stat().st_size
            fq1_size = (sample_out_dir_host / fq1_relpath).stat().st_size
            fq2_size = (sample_out_dir_host / fq2_relpath).stat().st_size

            for name in (bam_unsorted_name, bam_sorted_name, bam_sorted_idx_name):
                (sample_out_dir_host / name).unlink(missing_ok=True)

        finished = datetime.now(timezone.utc)
        metadata = {
            "sample_id": sample_id,
            "cram_url": cram_url,
            "population": population,
            "image": DOCKER_IMAGE,
            "started_utc": started.isoformat(timespec="seconds"),
            "finished_utc": finished.isoformat(timespec="seconds"),
            "duration_seconds": round((finished - started).total_seconds(), 2),
            "step_durations_seconds": durations,
            "bam_size_bytes": bam_size,
            "fastq_1_size_bytes": fq1_size,
            "fastq_2_size_bytes": fq2_size,
        }
        done_path_host.write_text(json.dumps(metadata, indent=2) + "\n")
        print(f"[{sample_id}] Done in {metadata['duration_seconds']}s "
              f"(BAM {bam_size / 1e6:.0f} MB)")

    except BaseException as exc:
        # Catch BaseException so KeyboardInterrupt also writes a .failed.json
        # before re-raising. Any unhandled exit path other than a hard kill
        # will record the failure for the user.
        finished = datetime.now(timezone.utc)
        failure_info = {
            "sample_id": sample_id,
            "cram_url": cram_url,
            "population": population,
            "image": DOCKER_IMAGE,
            "started_utc": started.isoformat(timespec="seconds"),
            "finished_utc": finished.isoformat(timespec="seconds"),
            "step_durations_seconds": durations,
            "error_type": type(exc).__name__,
            "error_message": str(exc),
        }
        failed_path_host.write_text(json.dumps(failure_info, indent=2) + "\n")
        print(f"[{sample_id}] FAILED ({type(exc).__name__}); see {log_path_host}")
        raise
    finally:
        # Remove the .running.json lock unconditionally. If we got here at
        # all, the Python process is still alive and the sample is no longer
        # in progress regardless of success or failure.
        running_path_host.unlink(missing_ok=True)


def _run_step(
    step_name: str,
    sample_id: str,
    log,
    durations: dict[str, float],
    mounts: list[str],
    container_cmd: list[str],
) -> None:
    """Run one pipeline step in a fresh container; tee output to the log file.

    Records the step's wall-clock duration in `durations`. Raises
    RuntimeError on non-zero exit, leaving partial outputs on disk for
    inspection.
    """
    started = time.monotonic()
    print(f"[{sample_id}] {step_name}...")
    log.write(f"\n=== {step_name} @ {datetime.now(timezone.utc).isoformat(timespec='seconds')} ===\n")
    log.flush()

    docker_cmd = [
        "docker", "run", "--rm", *mounts,
        "-w", OUTPUT_DIR_CTR,
        DOCKER_IMAGE, *container_cmd,
    ]
    log.write("$ " + " ".join(docker_cmd) + "\n")
    log.flush()

    result = subprocess.run(docker_cmd, stdout=log, stderr=subprocess.STDOUT)
    duration = round(time.monotonic() - started, 2)
    durations[step_name] = duration
    log.write(f"--- {step_name} exit {result.returncode} in {duration}s ---\n")
    log.flush()

    if result.returncode != 0:
        raise RuntimeError(
            f"step '{step_name}' failed (exit {result.returncode})"
        )


def _wipe_sample_state(
    sample_id: str,
    failed_path_host: Path,
    running_path_host: Path,
    sample_out_dir_host: Path,
    log_path_host: Path,
) -> None:
    """Remove all on-disk state for a sample so it can be retried clean.

    Deletes the .failed.json and .running.json markers (whichever exist), the
    per-sample output directory, and the per-sample log file. Called by
    run_pipeline's pre-flight pass when a previous run left the sample in a
    non-done state. The .done.json marker, if present, is never touched here
    — that case is handled separately as a skip.
    """
    failed_path_host.unlink(missing_ok=True)
    running_path_host.unlink(missing_ok=True)
    log_path_host.unlink(missing_ok=True)
    if sample_out_dir_host.exists():
        shutil.rmtree(sample_out_dir_host)


def _sanity_check_fastqs(sample_id: str, sample_out_dir_host: Path, log) -> None:
    """Verify both paired FASTQs exist and are non-empty before typing.

    Empty FASTQs from a botched slice would silently produce garbage HLA
    types; failing here is much cheaper than running SpecHLA on nothing.
    """
    hla_reads_dir_host = sample_out_dir_host / "hla_reads"
    fq1_path_host = hla_reads_dir_host / f"{sample_id}_extract_1.fq.gz"
    fq2_path_host = hla_reads_dir_host / f"{sample_id}_extract_2.fq.gz"
    log.write(f"\n=== sanity_check_fastqs ===\n")
    for fq_host in (fq1_path_host, fq2_path_host):
        if not fq_host.exists():
            raise RuntimeError(f"[{sample_id}] expected FASTQ missing: {fq_host}")
        size = fq_host.stat().st_size
        log.write(f"  {fq_host.name}: {size} bytes\n")
        if size == 0:
            raise RuntimeError(f"[{sample_id}] FASTQ is empty: {fq_host}")
    log.flush()


def _validate_hla_result(sample_id: str, sample_out_dir_host: Path, log) -> None:
    """Confirm SpecHLA's G-group result file contains at least one real allele.

    We validate against hla.result.g.group.txt rather than hla.result.txt
    because the former is written by SpecHLA's typing/phasing stage and is
    populated reliably even when the later annotation stage has trouble
    (which happens on samples with unusual DRB1 structural variation or
    other hard-to-annotate cases). The G-group file is also the right
    resolution for both population genetics and clinical HLA work.

    SpecHLA can exit 0 even when typing has produced nothing useful (e.g.,
    the realignment stage failed silently and every gene came back empty).
    This check reads the data row and raises if every allele cell is a
    dash. A few dashes are normal (some genes may have insufficient coverage
    in a given sample); all-dashes is the unambiguous failure signal.
    """
    result_path_host = sample_out_dir_host / "spechla_output" / sample_id / "hla.result.g.group.txt"
    log.write(f"\n=== validate_hla_result ===\n")
    if not result_path_host.exists():
        raise RuntimeError(
            f"[{sample_id}] hla.result.g.group.txt missing: {result_path_host}"
        )

    # The file has a comment line, a header line, and one data row whose
    # first cell is the sample id followed by allele columns.
    data_row = None
    for line in result_path_host.read_text().splitlines():
        if line.startswith(sample_id + "\t"):
            data_row = line
            break
    if data_row is None:
        raise RuntimeError(f"[{sample_id}] no data row in {result_path_host.name}")

    allele_cells = data_row.split("\t")[1:]
    log.write(f"  alleles: {allele_cells}\n")
    if all(cell.strip() == "-" for cell in allele_cells):
        raise RuntimeError(
            f"[{sample_id}] every allele in hla.result.g.group.txt is '-'; "
            f"SpecHLA likely failed silently"
        )
    log.flush()


def _log_header(log, sample_id: str, cram_url: str, population: str, image: str) -> None:
    """Write a banner to the per-sample log file marking the start of a run."""
    log.write("\n" + "=" * 72 + "\n")
    log.write(f"Sample:     {sample_id}\n")
    log.write(f"CRAM URL:   {cram_url}\n")
    log.write(f"Population: {population}\n")
    log.write(f"Image:      {image}\n")
    log.write(f"Started:    {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
    log.write("=" * 72 + "\n")
    log.flush()


# --- Parallel driver --------------------------------------------------------
def run_pipeline(
    samples: list[tuple[str, str, str]],
    chr6_alt_contigs: list[str],
    hla_contigs: list[str],
    max_workers: int = 4,
) -> None:
    """Process a list of samples in parallel via a thread pool.

    Each entry in `samples` is a (sample_id, cram_url, population) tuple.

    Pre-flight: walks the input list once, classifying each sample as
    already-done (skip), to-retry (wipe stale markers/outputs/log, then
    process), or fresh (process). Only to-retry and fresh samples are
    submitted to the worker pool.

    Per-sample failures are logged but do not stop the pool — other samples
    keep running. On Ctrl-C, in-flight samples are allowed to finish (their
    own try/finally writes the appropriate marker); no new samples start.

    A summary block is appended to work/logs/_driver.log at the end,
    including the pre-flight breakdown and the post-run succeeded/failed
    counts.
    """
    started = datetime.now(timezone.utc)

    # --- Pre-flight pass --------------------------------------------------
    already_done: list[str] = []
    to_process: list[tuple[str, str, str]] = []
    retried_count = 0

    for sample_id, cram_url, population in samples:
        done_path_host = MARKERS_DIR_HOST / f"{sample_id}.done.json"
        failed_path_host = MARKERS_DIR_HOST / f"{sample_id}.failed.json"
        running_path_host = MARKERS_DIR_HOST / f"{sample_id}.running.json"

        if done_path_host.exists():
            already_done.append(sample_id)
            continue

        if failed_path_host.exists() or running_path_host.exists():
            _wipe_sample_state(
                sample_id,
                failed_path_host,
                running_path_host,
                OUTPUTS_DIR_HOST / sample_id,
                LOGS_DIR_HOST / f"{sample_id}.log",
            )
            retried_count += 1

        to_process.append((sample_id, cram_url, population))

    fresh_count = len(to_process) - retried_count
    print(
        f"[driver] Input: {len(samples)} samples | "
        f"already done: {len(already_done)} | "
        f"to retry: {retried_count} | fresh: {fresh_count}"
    )
    print(f"[driver] Processing {len(to_process)} on {max_workers} workers")

    # --- Worker pool ------------------------------------------------------
    succeeded: list[str] = []
    failed: list[tuple[str, str]] = []  # (sample_id, error message)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                process_sample,
                sample_id, cram_url, population,
                chr6_alt_contigs, hla_contigs,
            ): sample_id
            for sample_id, cram_url, population in to_process
        }
        try:
            for future in as_completed(futures):
                sample_id = futures[future]
                try:
                    future.result()
                    succeeded.append(sample_id)
                except BaseException as exc:
                    failed.append((sample_id, f"{type(exc).__name__}: {exc}"))
        except KeyboardInterrupt:
            print("[driver] Ctrl-C received; waiting for in-flight samples to finish...")
            executor.shutdown(wait=True, cancel_futures=True)
            raise

    finished = datetime.now(timezone.utc)
    duration = round((finished - started).total_seconds(), 1)
    print(f"[driver] Done in {duration}s: "
          f"{len(succeeded)} succeeded, {len(failed)} failed")

    _append_driver_log(
        started, finished, max_workers,
        input_total=len(samples),
        already_done=len(already_done),
        retried=retried_count,
        fresh=fresh_count,
        succeeded=succeeded,
        failed=failed,
    )


def _append_driver_log(
    started: datetime,
    finished: datetime,
    max_workers: int,
    input_total: int,
    already_done: int,
    retried: int,
    fresh: int,
    succeeded: list[str],
    failed: list[tuple[str, str]],
) -> None:
    """Append a summary block for this run to work/logs/_driver.log."""
    driver_log_host = LOGS_DIR_HOST / "_driver.log"
    duration = round((finished - started).total_seconds(), 1)
    with driver_log_host.open("a") as f:
        f.write("\n" + "=" * 72 + "\n")
        f.write(f"Started:     {started.isoformat(timespec='seconds')}\n")
        f.write(f"Finished:    {finished.isoformat(timespec='seconds')}\n")
        f.write(f"Duration:    {duration}s\n")
        f.write(f"Max workers: {max_workers}\n")
        f.write("\n")
        f.write(f"Input total: {input_total}\n")
        f.write(f"Already done (skipped): {already_done}\n")
        f.write(f"Retried:     {retried}\n")
        f.write(f"Fresh:       {fresh}\n")
        f.write(f"Processed:   {retried + fresh}\n")
        f.write("\n")
        f.write(f"Succeeded:   {len(succeeded)}\n")
        f.write(f"Failed:      {len(failed)}\n")
        for sample_id, message in failed:
            f.write(f"  - {sample_id}: {message}\n")
        f.write("=" * 72 + "\n")


# --- Entry point ------------------------------------------------------------
def main() -> None:
    setup()

    # Contig lists are constant for the whole run; load once after setup
    # has guaranteed the files exist, then pass to every process_sample call.
    chr6_alt_contigs = _read_contig_list(INPUT_DIR_HOST / CHR6_ALT_CONTIGS_FILE)
    hla_contigs = _read_contig_list(INPUT_DIR_HOST / HLA_CONTIGS_FILE)

    # Smoke test list. Replace with the JSON-driven loop once the loader
    # is built.
    samples = [
        ("HG01882", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242189/HG01882.final.cram", "Black"),
        ("HG01883", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242190/HG01883.final.cram", "Black"),
        ("HG01888", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3988942/HG01888.final.cram", "Black"),

        ("NA12878", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram", "Caucasian"),

        ("HG00403", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241665/HG00403.final.cram", "Asian"),
        ("HG00404", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241666/HG00404.final.cram", "Asian"),

        ("HG00406", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241667/HG00406.final.cram", "Asian"),
        ("HG00407", "https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241668/HG00407.final.cram", "Asian")
    ]
    run_pipeline(samples, chr6_alt_contigs, hla_contigs, max_workers=4)


if __name__ == "__main__":
    main()

