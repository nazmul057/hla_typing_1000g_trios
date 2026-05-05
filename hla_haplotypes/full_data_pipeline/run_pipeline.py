"""HLA typing pipeline driver for 1000 Genomes Project trios.

Streams HLA-relevant reads from remote CRAMs and runs SpecHLA inside a pinned
Docker image. Designed to be run from a single command:

    python run_pipeline.py

On startup, the pipeline verifies its setup (working directories, reference
genome, contig lists). If anything is missing, setup is performed before
sample processing begins.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# --- Configuration ----------------------------------------------------------
# Edit ROOT_DIR to point at your pipeline working directory. Everything the
# pipeline reads or writes lives under this path. Relative paths are fine;
# they are resolved to absolute form before being passed to Docker.
ROOT_DIR = Path("hla_haplotypes/full_data_pipeline").resolve()

DOCKER_IMAGE = "spechla_pipeline:v1"

REFERENCE_FASTA_NAME = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
REFERENCE_FASTA_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "GRCh38_reference_genome/" + REFERENCE_FASTA_NAME
)

# Subdirectory layout under ROOT_DIR. See architecture doc for the full picture.
SUBDIRS = ("reference", "input", "outputs", "work/markers", "work/logs")

# Contig list filenames written into input/. These are derived once from the
# reference .fai and reused by every sample.
HLA_CONTIGS_FILE = "hla_contigs.txt"
CHR6_ALT_CONTIGS_FILE = "chr6_alt_contigs.txt"

# Sanity-check ranges for the GRCh38 1000 Genomes analysis-set reference.
# Expected: ~525 HLA decoy contigs, ~16 chr6 alt contigs. Anything far outside
# these ranges suggests the wrong reference was downloaded.
EXPECTED_HLA_CONTIG_COUNT_RANGE = (500, 600)
EXPECTED_CHR6_ALT_CONTIG_COUNT_RANGE = (10, 30)


# --- Setup ------------------------------------------------------------------
def setup(root_dir: Path, image: str) -> None:
    """Prepare the pipeline working environment under root_dir.

    Creates the standard subdirectories, downloads and indexes the reference
    FASTA if missing, and derives the HLA / chr6-alt contig lists from the
    reference index. Idempotent: safe to call on every pipeline run.
    """
    print(f"[setup] Root: {root_dir}")
    _create_subdirs(root_dir)
    _ensure_reference(root_dir, image)
    _ensure_contig_lists(root_dir)
    print("[setup] Done.")


def _create_subdirs(root_dir: Path) -> None:
    """Create the standard subdirectory tree if any pieces are missing."""
    for subdir in SUBDIRS:
        (root_dir / subdir).mkdir(parents=True, exist_ok=True)


def _ensure_reference(root_dir: Path, image: str) -> None:
    """Download and index the reference FASTA if either is absent.

    Both download and indexing run inside the Docker image, so the host needs
    nothing beyond Docker and Python. The reference directory is mounted rw
    during setup; at pipeline runtime it is mounted ro by the per-sample
    workers.
    """
    reference_dir = root_dir / "reference"
    fasta_path = reference_dir / REFERENCE_FASTA_NAME
    fai_path = reference_dir / (REFERENCE_FASTA_NAME + ".fai")

    if fasta_path.exists() and fai_path.exists():
        print(f"[setup] Reference already present: {fasta_path.name} (+ .fai)")
        return

    mount = f"{reference_dir.as_posix()}:/data/reference"

    if not fasta_path.exists():
        print(f"[setup] Downloading reference FASTA (~3 GB) to {fasta_path}")
        # `wget -c` resumes a partial download instead of restarting from zero
        # if a previous attempt was interrupted.
        download_cmd = [
            "docker", "run", "--rm",
            "-v", mount,
            "-w", "/data/reference",
            image,
            "wget", "-c", REFERENCE_FASTA_URL,
        ]
        _run_or_die(download_cmd, "reference download failed")

    print(f"[setup] Building FASTA index: {fai_path.name}")
    faidx_cmd = [
        "docker", "run", "--rm",
        "-v", mount,
        image,
        "samtools", "faidx", f"/data/reference/{REFERENCE_FASTA_NAME}",
    ]
    _run_or_die(faidx_cmd, "samtools faidx failed")


def _ensure_contig_lists(root_dir: Path) -> None:
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
    input_dir = root_dir / "input"
    hla_path = input_dir / HLA_CONTIGS_FILE
    chr6_alt_path = input_dir / CHR6_ALT_CONTIGS_FILE

    if hla_path.exists() and chr6_alt_path.exists():
        print(f"[setup] Contig lists already present in {input_dir}")
        return

    fai_path = root_dir / "reference" / (REFERENCE_FASTA_NAME + ".fai")
    if not fai_path.exists():
        print(f"[setup] ERROR: expected {fai_path} to exist by now", file=sys.stderr)
        sys.exit(1)

    hla_contigs: list[str] = []
    chr6_alt_contigs: list[str] = []
    with fai_path.open() as f:
        for line in f:
            name = line.split("\t", 1)[0]
            if name.startswith("HLA-"):
                hla_contigs.append(name)
            elif name.startswith("chr6_") and name.endswith("_alt"):
                chr6_alt_contigs.append(name)

    _check_count("HLA decoy", len(hla_contigs), EXPECTED_HLA_CONTIG_COUNT_RANGE)
    _check_count("chr6 alt", len(chr6_alt_contigs), EXPECTED_CHR6_ALT_CONTIG_COUNT_RANGE)

    hla_path.write_text("\n".join(hla_contigs) + "\n")
    chr6_alt_path.write_text("\n".join(chr6_alt_contigs) + "\n")
    print(f"[setup] Wrote {len(hla_contigs)} HLA contigs to {hla_path.name}")
    print(f"[setup] Wrote {len(chr6_alt_contigs)} chr6 alt contigs to {chr6_alt_path.name}")


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


def _run_or_die(cmd: list[str], error_message: str) -> None:
    """Run a subprocess, streaming output to the console; exit on failure."""
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[setup] ERROR: {error_message} (exit {exc.returncode})", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"[setup] ERROR: command not found: {cmd[0]}", file=sys.stderr)
        print("[setup] Is Docker installed and on PATH?", file=sys.stderr)
        sys.exit(1)


# --- Entry point ------------------------------------------------------------
def main() -> None:
    setup(ROOT_DIR, DOCKER_IMAGE)


if __name__ == "__main__":
    main()