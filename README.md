# SpecHLA in Docker

A reproducible Docker setup for running [SpecHLA](https://github.com/deepomicslab/SpecHLA) HLA typing on paired-end reads, without installing it on the host.

## What's in this folder

- `Dockerfile` — builds an image with SpecHLA 1.0.10 installed via bioconda.
- `spec_hla_ctr.py` — config-in-code Python wrapper that runs SpecHLA in a disposable container against a mounted data directory.

## Setup attributes

| Attribute        | Value                                                       |
|------------------|-------------------------------------------------------------|
| Base image       | `debian:bookworm` (Debian 12, full — not slim)              |
| SpecHLA version  | `1.0.10` (pinned via `spechla=1.0.10` in conda install)     |
| Conda            | Latest Miniconda installed at `/opt/miniconda3`             |
| Conda env name   | `spechla`                                                   |
| Channels         | `bioconda`, `conda-forge`                                   |
| Working dir      | `/data` (host directory bind-mounted here)                  |
| Entrypoint / CMD | None — caller passes the full command                       |

A non-slim Debian base is used because SpecHLA's internal shell scripts depend on common Linux utilities (`less`, etc.) that are stripped from `*-slim` variants. The image also installs a generous bundle (`less`, `procps`, `bash`, `perl`, `file`, `vim-tiny`) on top of the conda dependencies.

## Build the image (once)

From the directory containing the Dockerfile:

```bash
docker build -t spechla:1.0.10 .
```

If the Dockerfile is renamed (e.g., `Dockerfile.spechla`):

```bash
docker build -f Dockerfile.spechla -t spechla:1.0.10 .
```

The trailing `.` is the build context (current directory). Build takes ~10–30 minutes depending on bandwidth and the bioconda solver — most of that is the conda layer.

### Verify the build

```bash
docker run --rm spechla:1.0.10 conda run -n spechla conda list spechla
```

Should print `spechla 1.0.10 ...`.

## Run via the Python script

Edit the `CONFIG` block at the top of `spec_hla_ctr.py`:

```python
DATA_DIR = Path("ops/test_data")               # mounted at /data inside container
SAMPLE = "HG00096"
READ1 = "in/HG00096/HG00096_HLA.R1.fastq.gz"   # relative to DATA_DIR
READ2 = "in/HG00096/HG00096_HLA.R2.fastq.gz"   # relative to DATA_DIR
OUTPUT_SUBDIR = "out"                          # relative to DATA_DIR
THREADS = 5
IMAGE = "spechla:1.0.10"
EXTRA_ARGS = []                                # e.g. ["-p", "Asian"] or ["-u", "1"]
```

Then run:

```bash
uv run spec_hla_ctr.py
# or:
python spec_hla_ctr.py
```

The script prints the underlying `docker run` command for transparency, then executes it. Output lands in `<DATA_DIR>/<OUTPUT_SUBDIR>/<SAMPLE>/`; the main result is `hla.result.txt`.

## Run directly with Docker (no Python)

If you'd rather skip the script:

```bash
docker run --rm \
    -v /absolute/path/to/data:/data \
    spechla:1.0.10 \
    conda run --no-capture-output -n spechla \
        spechla -j 5 -n HG00096 \
                -1 /data/in/HG00096/HG00096_HLA.R1.fastq.gz \
                -2 /data/in/HG00096/HG00096_HLA.R2.fastq.gz \
                -o /data/out
```

Notes on the flags:
- `--rm` deletes the container after it exits (each run is disposable).
- `-v <host>:/data` bind-mounts your data directory; paths inside the container are then `/data/...`.
- `conda run --no-capture-output -n spechla` activates the env without buffering stdout/stderr (so progress streams live).

On Linux/macOS, add `--user $(id -u):$(id -g)` so output files aren't owned by `root` on the host. (The Python script does this automatically.)

## Common SpecHLA options

Pass these via `EXTRA_ARGS` (script) or appended to `spechla …` (direct Docker):

| Args              | Purpose                                                                        |
|-------------------|--------------------------------------------------------------------------------|
| `-u 1`            | Exon-only typing (for WES / RNA-seq).                                          |
| `-p Asian`        | Population-aware annotation (also: `Black`, `Caucasian`, `Unknown`, `nonuse`). |
| `-j N`            | Threads. Set via `THREADS` in the script.                                      |

Long-read and hybrid modes use different SpecHLA entrypoints (`spechla-long-read`, `-t`/`-e` flags); call those directly via `docker run` rather than through the script.

## Output layout

```
<DATA_DIR>/<OUTPUT_SUBDIR>/<SAMPLE>/
├── hla.result.txt                # main result (4-field resolution)
├── hla.result.details.txt        # per-allele detail
├── hla.result.g.group.txt        # G-group resolution
├── hla.allele.{1,2}.HLA_*.fasta  # reconstructed sequences
├── HLA_*_freq.txt                # haplotype frequencies (~0.5:0.5 = confident)
└── HLA_*.rephase.vcf.gz          # phased VCFs
```

## Known harmless log lines

These appear during normal runs and can be ignored:

- `run.assembly.realign.sh: 15: source: not found` — non-fatal; SpecHLA continues past it.
- `SyntaxWarning: invalid escape sequence` in `phase_variants.py` — Python 3.12 deprecation warnings in SpecHLA's source; not errors.
- `[W::bcf_hdr_check_sanity] GQ should be declared as Type=Integer` — bcftools nag about an upstream VCF header.
- Low alignment rate (e.g., ~7%) when feeding HLA-extracted reads against just the HLA database — expected.

A high "ratio of N (masked)" for a gene (>0.3) indicates low coverage on that locus and lower-confidence calls; check `HLA_<gene>_freq.txt` for balance and `Mean depth` in the log.