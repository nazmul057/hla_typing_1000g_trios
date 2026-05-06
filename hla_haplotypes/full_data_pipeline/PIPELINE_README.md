# HLA Typing Pipeline for Remote 1kGP CRAMs

A Python-driven pipeline that performs full-resolution HLA typing on samples
from public archives (1000 Genomes Project) by streaming HLA-relevant reads
directly from remote CRAM files. The pipeline orchestrates SpecHLA inside a
pinned Docker image and is designed for batch processing of trio-organized
cohorts (~600 trios, ~1,800 individuals).

The pipeline avoids downloading 15–25 GB CRAMs per sample, instead
transferring ~200 MB of HLA-relevant data per sample via HTTP range requests
on the indexed remote CRAM. Per-sample wall time is dominated by network
I/O, so the driver runs many samples concurrently against the same remote
archive.

---

## Table of contents

- [Prerequisites](#prerequisites)
- [Repository layout](#repository-layout)
- [Building the Docker image](#building-the-docker-image)
- [Configuring the pipeline](#configuring-the-pipeline)
- [Running the pipeline](#running-the-pipeline)
- [Output structure](#output-structure)
- [Re-runs and retries](#re-runs-and-retries)
- [Pipeline internals](#pipeline-internals)
- [Design decisions and rationale](#design-decisions-and-rationale)
- [Issues encountered and resolutions](#issues-encountered-and-resolutions)
- [Known limitations](#known-limitations)

---

## Prerequisites

- **Docker** — Docker Desktop on Windows/macOS, or a Docker Engine
  installation on Linux. The pipeline launches one container per sample
  step; nothing else is installed on the host.
- **Python 3.10+** — standard library only. No third-party packages.
- **Disk** — ~3 GB for the reference FASTA, plus ~250 MB transient per
  in-flight sample (BAM intermediates, deleted after typing succeeds).
- **Network** — outbound HTTPS to `ftp.sra.ebi.ac.uk` (CRAM source) and
  `ftp.1000genomes.ebi.ac.uk` (reference FASTA, fetched once during setup).

Tested on Windows 11 with Docker Desktop and on Debian 12. Production target
is a Debian server.

---

## Repository layout

Everything is in a single working directory:

```
<project_root>/
├── toolkit.dockerfile           # Docker image definition
├── run_pipeline.py              # Pipeline driver (single file)
├── README.md                    # This file
├── reference/                   # Created by setup; reference FASTA + .fai
├── input/                       # Created by setup; contig lists, families JSON
├── outputs/<sample_id>/         # Created per sample; FASTQs + SpecHLA results
├── work/
│   ├── markers/<sample_id>.{done,failed,running}.json
│   └── logs/
│       ├── _driver.log          # One block per pipeline run
│       └── <sample_id>.log      # Per-sample full output of all 8 steps
└── test/                        # Optional; ad-hoc test scripts/data
```

`reference/`, `input/`, `outputs/`, and `work/` are created automatically on
first run. The only files that must exist before the first run are
`toolkit.dockerfile` and `run_pipeline.py`.

---

## Building the Docker image

From the project directory:

```bash
docker build -f toolkit.dockerfile -t spechla_pipeline:v1 .
```

This produces an image tagged `spechla_pipeline:v1` containing:

- Debian Bookworm base
- Miniconda + a `spechla_env` conda environment with `spechla=1.0.10` from
  bioconda (which brings samtools, bwa, freebayes, and other dependencies)
- `tini` as PID 1 (for clean signal propagation when stopping containers)
- `wget`, `curl`, `bash`, and a few utilities used during setup

The Dockerfile contains build-time verifications (correct SpecHLA version,
samtools with HTTPS support, SpecHLA database/scripts present) so a
successful build implies a runnable image.

The image is the only container artifact the pipeline uses. The Python
driver invokes `docker run --rm` per step.

---

## Configuring the pipeline

Open `run_pipeline.py` and edit the path constants at the top of the file:

```python
ROOT_DIR_HOST = Path("hla_haplotypes/full_data_pipeline").resolve()
DOCKER_IMAGE = "spechla_pipeline:v1"
THREADS_PER_SAMPLE = 4
```

`ROOT_DIR_HOST` is the only constant you usually need to change. It must
resolve to an absolute path (Docker rejects relative paths in `-v` mounts).

All other paths are derived from `ROOT_DIR_HOST`. Paths visible to the host
machine carry the `_host` suffix throughout the code; paths visible inside
the container carry `_ctr`. This convention makes the Docker mount
boundary explicit at every use site.

To run on a different cohort, edit the smoke-test list inside `main()`:

```python
samples = [
    ("HG01882", "https://ftp.sra.ebi.ac.uk/.../HG01882.final.cram", "Black"),
    ("HG01888", "https://ftp.sra.ebi.ac.uk/.../HG01888.final.cram", "Black"),
]
```

For the full 1kGP cohort, replace this with code that reads
`families_1000G_trios.json` and yields one tuple per trio member. The JSON
loader is left as an exercise for the user; the `(sample_id, cram_url,
population)` tuple shape is the contract `run_pipeline` consumes.

`population` must be one of `Asian`, `Caucasian`, `Black`, `Unknown` (these
are SpecHLA's `-p` flag values).

---

## Running the pipeline

From the project directory:

```bash
python run_pipeline.py
```

On a fresh machine, the first invocation does several things in order:

1. Verifies Docker is reachable (`docker info`).
2. Creates the directory layout under `ROOT_DIR_HOST`.
3. Downloads `GRCh38_full_analysis_set_plus_decoy_hla.fa` (~3 GB) from
   the 1kGP EBI mirror, into a container with `wget -c` (resumes
   interrupted downloads).
4. Indexes the FASTA with `samtools faidx`.
5. Reads the `.fai` and writes the HLA-decoy and chr6-alt contig lists
   to `input/`. Sanity-checks counts (~525 HLA, ~16 chr6 alt) and warns
   on mismatch.

Steps 1–5 are idempotent: subsequent runs detect existing artifacts and
skip them in milliseconds.

Then the pipeline processes the sample list in parallel. Default is 4
concurrent workers; raise or lower via the `max_workers` argument to
`run_pipeline()`.

A run prints:

- Setup status lines (`[setup] ...`)
- A driver-level scope summary (`[driver] Input: 8 samples | already done:
  6 | to retry: 0 | fresh: 2`)
- One progress line per step per sample, prefixed with `[sample_id]`
- Per-sample completion summaries with duration and BAM size
- A final driver-level totals line

Detailed output (full container stdout/stderr per step) is written to
`work/logs/<sample_id>.log`. The console only carries the high-level
narrative.

---

## Output structure

Per successfully-completed sample:

```
outputs/<sample_id>/
├── hla_reads/
│   ├── <sample_id>_extract_1.fq.gz       # Paired-end HLA-only FASTQ (~1 MB)
│   └── <sample_id>_extract_2.fq.gz
└── spechla_output/<sample_id>/
    ├── hla.result.g.group.txt            # G-group calls (canonical output)
    ├── hla.result.txt                    # 4-field calls (may be incomplete)
    ├── hla.result.details.txt            # Per-allele scoring (may be incomplete)
    └── ... (intermediate phasing/variant files)

work/markers/<sample_id>.done.json        # Run metadata (see below)
work/logs/<sample_id>.log                 # Full per-step log
```

The `.done.json` marker contains a structured summary of the run:

```json
{
  "sample_id": "HG01882",
  "cram_url": "https://...",
  "population": "Black",
  "image": "spechla_pipeline:v1",
  "started_utc": "2026-05-05T14:26:42+00:00",
  "finished_utc": "2026-05-05T14:49:27+00:00",
  "duration_seconds": 1365.17,
  "step_durations_seconds": {
    "slice": 1229.79,
    "sort": 7.57,
    "index": 6.16,
    "extract_hla_reads": 10.26,
    "spechla_type": 111.35
  },
  "bam_size_bytes": 221232859,
  "fastq_1_size_bytes": 1051248,
  "fastq_2_size_bytes": 1180014
}
```

The intermediate BAMs (full and coord-sorted) are deleted after typing
succeeds. They're large (~200 MB each), regenerable from the CRAM, and
not needed for downstream HLA analysis.

**Use `hla.result.g.group.txt` as the canonical output**, not
`hla.result.txt`. See [Issues encountered](#issues-encountered-and-resolutions)
for why.

---

## Re-runs and retries

Each pipeline run is a self-contained pass over the input list. Behavior
per sample:

| Marker present | Behavior |
|----------------|----------|
| `.done.json` | Skip — already succeeded. |
| `.failed.json` or `.running.json` | Wipe per-sample state, process fresh. |
| None | Process fresh. |

The pre-flight pass (at the top of `run_pipeline()`) walks the input list
once, classifies each sample, and wipes stale state for any sample being
retried. "Wipe" means: delete the marker files, delete `outputs/<sample_id>/`,
delete `work/logs/<sample_id>.log`. Only then are workers launched.

This means: investigate failures by reading `work/logs/<sample_id>.log` and
`work/markers/<sample_id>.failed.json`, fix the underlying issue (or the
input data), then re-run the pipeline normally. No manual marker deletion,
no force flags. The driver-level `_driver.log` accumulates a history block
per run, so you can see across runs how many samples were attempted, how
many succeeded, and which ones repeatedly failed.

A `.running.json` left from a hard-killed previous run (kill -9, OS crash)
is treated as a stale lock and wiped during pre-flight, just like a
`.failed.json`. The exclusive-create on `.running.json` during a normal run
still protects against intra-run races between workers.

---

## Pipeline internals

### Per-sample workflow (eight steps inside the worker)

Each step is its own short-lived `docker run --rm` invocation:

1. **slice** — `samtools view` streams HLA-relevant reads from the remote
   CRAM. The query is the chr6 MHC region (10 Mb), all chr6 alt contigs,
   all HLA decoy contigs, and unmapped reads. Output: `<sample>.hla_full.bam`.
2. **sort** — coordinate-sort the BAM.
3. **index** — `.bai` for the sorted BAM.
4. **extract_hla_reads** — `spechla-extract-hla-reads` filters the BAM to
   HLA-gene-specific reads and produces paired FASTQs.
5. **sanity_check_fastqs** — host-side check that both FASTQs exist and
   are non-empty before paying for typing. Failed-fast on bad slices.
6. **spechla_type** — `spechla` performs full-resolution typing.
7. **validate_hla_result** — host-side check that
   `hla.result.g.group.txt` has at least one non-dash allele call.
   Catches silent SpecHLA failures.
8. **cleanup** — record BAM/FASTQ sizes in metadata; delete the BAM
   intermediates.

Each step's stdout/stderr is appended to `work/logs/<sample_id>.log` along
with a banner, the exact docker command run, and a duration line. Per-step
durations are also recorded in the `.done.json` marker.

### Concurrency model

- **Host-side `ThreadPoolExecutor`** — one Python thread per concurrent
  sample. Each thread runs `process_sample`, which itself launches eight
  short-lived containers in sequence.
- **Default 4 workers.** The bottleneck is network throughput from the
  remote CRAM source, not local CPU. EBI rate-limits per IP at some
  unspecified threshold; 4 has been observed to work well, 6+ may
  trigger throttling. Tune by trial.
- **Pre-flight runs serially**, before any worker starts. Per-sample wipes
  happen before any container is launched. This eliminates the entire
  class of "what if two workers race on the same sample dir" question.
- **Each worker mounts only its own per-sample output directory** as
  read/write. The reference directory is mounted read-only. There is no
  shared writable state inside the container layer.
- **Marker files use atomic exclusive create** (`open("x")`) for
  `.running.json` as in-run protection against duplicate sample IDs in
  the input list.
- **`Ctrl-C`** triggers a clean shutdown: no new samples are submitted,
  in-flight samples finish naturally (their try/finally writes the
  appropriate marker), then the script exits.
- **Per-sample failures don't stop the pool** — other samples keep
  running, the failure is recorded in `_driver.log`'s summary, and you
  re-run after investigating.

---

## Design decisions and rationale

### One Python file, no extra dependencies

The driver is a single `run_pipeline.py` using only Python's standard
library. No `Click`, no `tenacity`, no `docker-py`. Reasons: minimal
deployment surface (Python + Docker is enough), every dependency is one
the team already understands, and the file is short enough to read
end-to-end in one sitting.

### Eight `docker run` calls per sample, not one

The validated pipeline is eight sequential samtools/SpecHLA invocations.
We could have wrapped them in a single `bash -c "step1 && step2 && ..."`
container but chose not to. Each step in its own container gives:

- Per-step exit codes and durations, surfaced to Python directly.
- Per-step error messages without parsing a combined stdout.
- No nested shell quoting (the source of a real bug we hit, see below).

The cost — ~50ms × 8 ≈ 400ms of container startup overhead per sample —
is invisible against multi-minute network slices.

### Contig lists derived from `.fai`, not from CRAM headers

The pipeline doc we started from derives the HLA-decoy and chr6-alt
contig lists by reading the `@SQ` headers of a remote CRAM. We
derive them from the locally-indexed reference FASTA's `.fai` instead.
Both produce identical lists for the 1kGP cohort because every CRAM was
aligned against this exact reference. The `.fai` approach removes a
network call from setup and removes a sample-URL dependency from
`setup()`, which keeps the setup function self-contained.

### Container-side path constants vs host-side path constants

The code uses two clearly distinguished namespaces:

- `*_HOST` constants — paths in the host filesystem (e.g.,
  `OUTPUTS_DIR_HOST`).
- `*_CTR` constants — paths visible inside containers (e.g.,
  `OUTPUT_DIR_CTR = "/data/output"`).

Local variables and function parameters inherit the same suffix
convention. Mixing host and container paths is the kind of bug that's
silent until it isn't, and visible suffixes make every use site
self-checking. For a Docker-driven pipeline, this convention is worth
the modest ceremony.

### No auto-retry classification

We considered classifying failures (retriable network errors vs.
permanent errors) and auto-retrying only the retriable ones. We chose
the simpler model: every re-run retries every non-succeeded sample, and
investigation is the user's job. Heuristic classifiers on opaque error
strings are brittle. The user already knows whether they fixed something
between runs; the system shouldn't second-guess.

### Validation against `hla.result.g.group.txt`

SpecHLA produces several output files, not all reliably populated. We
validate against the G-group file because (a) it's the typing-stage
output, populated before annotation issues that affect other files;
(b) it's the right resolution for both clinical and population genetics
work; (c) it survived sample-specific edge cases that left other output
files incomplete. See the issues section below.

---

## Issues encountered and resolutions

This pipeline went through several iterations. The substantive issues
are documented here for two reasons: they're the kind of thing that
matters more than the happy path, and they document the reasoning
behind specific choices in the code.

### 1. Silent typing failures from mangled contig names

**Symptom.** Early test runs produced `hla.result.txt` files with every
allele as `-`. The slice step appeared to succeed and the BAMs were the
expected size, but no HLA was getting typed.

**Cause.** The slice step originally used `bash -c "..."` with shell
substitutions like `$(cat /data/input/hla_contigs.txt | tr '\n' ' ')` to
expand the contig file into a positional argument list. On Windows host
+ Linux container, the argv-to-command-line marshalling collapsed the
quoting of the bash arg, so the `tr '\n' ' '` invocation received the
literal two characters `\` and `n` rather than a newline character. The
contig names arrived at samtools with embedded newlines; samtools warned
and silently skipped them. Only the unmapped reads and the primary
chr6 region made it through, missing the alt-haplotype reads SpecHLA
needs to type the more divergent HLA alleles.

**Fix.** Read the contig files in Python, splice the contig names
directly into the `subprocess.run` argv list. No shell, no quoting.
Drop the `bash -c` wrapper for the slice step and the `/data/input`
mount along with it (the container no longer needs to read those files).

### 2. SpecHLA scripts failing under `dash`

**Symptom.** After fix #1, slice produced legitimate BAMs but SpecHLA
still produced empty results. Logs contained the line:

```
.../run.assembly.realign.sh: 15: source: not found
```

followed downstream by `pysam.libcbcf.VariantFile` errors complaining
about missing VCFs.

**Cause.** SpecHLA's shell scripts use `source` (a bash builtin). On
Debian, `/bin/sh` is `dash`, a strict POSIX shell that doesn't
recognize `source`. Some part of SpecHLA's pipeline invoked the script
as `sh script.sh` (which discards the shebang and uses `dash`), so the
script's `#!/bin/bash` was bypassed. The first `source` call failed
silently from samtools' perspective, the realignment stage produced no
VCF, and downstream stages failed reading a file that was never written.

**Fix.** Add `ln -sf /bin/bash /bin/sh` to the Dockerfile. This makes
`/bin/sh` resolve to bash for everything in the container, which is
the convention bioinformatics images typically follow. One Dockerfile
line, no patches to SpecHLA's source.

### 3. The all-dashes results problem (the most important one)

**Symptom.** After fixes #1 and #2, most samples typed correctly. But
some samples — particularly ones with hard-to-resolve DRB1 due to
structural variation — produced `hla.result.txt` files with all-dash
rows even though the run "succeeded" (SpecHLA exited 0). Worse, the
silently-bad results would have been written as `.done.json` markers
and presented as valid output for downstream analysis.

**Investigation.** Reading the per-sample log carefully showed that
SpecHLA's typing/phasing stage actually completed successfully —
`Mean depth` and `Phasing of HLA_X is done!` lines for every gene.
The all-dashes appeared in `hla.result.txt` because of failures in
SpecHLA's *annotation* stage, which runs after typing and tries to
attach detail metadata. Annotation crashed partway through (a
`samtools mpileup -t` invocation that's incompatible with newer
samtools), leaving `hla.result.txt` and `hla.result.details.txt`
incomplete.

But: SpecHLA also produces `hla.result.g.group.txt`, written by the
typing stage *before* annotation runs. That file had complete, correct
G-group calls for all eight HLA genes — including DRB1, even on the
samples where annotation had failed. For HG01882, the G-group calls
matched the published 1kGP reference exactly.

**Fix.** Two changes:

(a) The validation step now reads `hla.result.g.group.txt` rather
than `hla.result.txt`. If every cell is a dash, the run is treated
as a silent SpecHLA failure and converted to a `.failed.json`. If
at least one allele is called, the run is accepted.

(b) The README and downstream consumers should treat
`hla.result.g.group.txt` as the canonical output. The 4-field file
(`hla.result.txt`) is sometimes incomplete and shouldn't be relied
on. The G-group resolution is also the right one for clinical and
population HLA work — IPD-IMGT/HLA defines G-groups as alleles that
encode identical antigen-binding protein sequences, which is the
biologically meaningful unit.

This is the single most important caveat to know about the pipeline's
output.

### 4. Docker daemon preflight

Early runs sometimes failed at the slice step with "failed to connect
to the docker API at npipe://..." — the symptom of forgetting to start
Docker Desktop. The error appeared inside a per-sample log file, ~3
minutes into the run, looking like a step failure.

**Fix.** Setup now runs `docker info` as its first action. If the
daemon isn't reachable, the script exits immediately with a clear
"Is Docker Desktop running?" message. Mistake caught in 1 second
instead of after queueing several minutes of work.

### 5. Slice step bandwidth and tuning

Network throughput on the slice step varied between ~150 KB/s and ~2
MB/s depending on time of day and (probably) EBI's load. At ~200 MB
per slice, this gives per-sample slice times between 100 seconds and
20+ minutes. The CPU-bound steps (sort, index, extract, type) are 2–3
minutes total per sample, so the slice dominates.

This shape — long network-bound stages followed by short CPU-bound
stages — is what makes parallelism so effective. Four samples slicing
concurrently saturate the network pipe; their CPU-bound stages
overlap nicely (one sample's sort happens while another's slice is
still running). Empirically, 4 workers cuts total wall time roughly
3–4× compared to serial. Pushing workers higher risks EBI throttling
without proportional gain.

---

## Known limitations

- **Docker is required.** The pipeline doesn't have a host-only fallback.
  Running on a system without Docker would mean reimplementing every step
  to use host-installed samtools/SpecHLA.
- **No CLI yet.** Sample list is hardcoded in `main()`. The
  `families_1000G_trios.json` loader is a planned addition.
- **`max_workers` is a constant in code**, not a CLI flag. Edit
  `run_pipeline.py`'s `main()` to change it.
- **No streaming of multi-machine work.** The marker model is local
  filesystem only; running two driver processes against the same
  `ROOT_DIR_HOST` is not supported and not protected against. NFS-shared
  output directories are also not supported.
- **Annotation incompleteness for some samples.** As described in issue
  #3, `hla.result.details.txt` may be missing entries for one or more
  genes on samples with unusual HLA structural variation. The G-group
  calls in `hla.result.g.group.txt` remain valid; the affected output
  is the per-allele scoring metadata. If you need scoring confidence
  for every gene in every sample, you'll need to address SpecHLA's
  annotation issues separately (e.g., by patching the `mpileup -t`
  invocation in SpecHLA's annotation script or pinning compatible
  samtools/bcftools versions).
- **Hard kills (kill -9, power loss) leave a stale `.running.json`.**
  The next pipeline run wipes it during pre-flight as part of the retry
  logic, so this is not a manual-cleanup case in practice. But if you
  somehow have two driver processes against the same `ROOT_DIR_HOST`
  (which we don't support anyway), the wipe would be wrong.

---

## License and acknowledgments

The pipeline orchestrates [SpecHLA](https://github.com/deepomicslab/SpecHLA)
(Wang et al.), which is the actual HLA typer. This repository contains the
driver and Docker recipe; cite SpecHLA appropriately when publishing
results derived from this pipeline.

The 1000 Genomes Project high-coverage cohort used here is described in
Byrska-Bishop et al., Cell 2022.

