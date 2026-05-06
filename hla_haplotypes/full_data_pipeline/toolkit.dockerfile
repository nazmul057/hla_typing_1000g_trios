# SpecHLA pipeline image for HLA typing on remote 1kGP CRAMs.
#
# Build:
#   docker build -f toolkit.dockerfile -t spechla_pipeline:v1 .
#
# Test:
#    docker run --rm spechla_pipeline:v1 samtools --version | head -3
#    docker run --rm spechla_pipeline:v1 spechla -h 2>&1 | head -5
#    docker run --rm spechla_pipeline:v1 which spechla-extract-hla-reads
#
#    Windows:
#    docker run --rm spechla_pipeline:v1 bash -c "samtools --version | head -3"
#    docker run --rm spechla_pipeline:v1 bash -c "spechla -h 2>&1 | head -5"
#    docker run --rm spechla_pipeline:v1 which spechla-extract-hla-reads

#    :: Verification 1: SpecHLA is the pinned version.
#    docker run --rm spechla_pipeline:v1 bash -c "conda run -n spechla_env conda list spechla | grep -E '^spechla\s+1\.0\.10\b' || (echo 'SpecHLA version mismatch' && exit 1)"
#    :: Verification 2: samtools came in and supports HTTPS (libcurl).
#     docker run --rm spechla_pipeline:v1 bash -c "samtools --version | head -3 && samtools --version | grep -qi 'libcurl' || (echo 'samtools missing libcurl support' && exit 1)"
#     :: Verification 3: SpecHLA entrypoints + database/scripts directories.
#     docker run --rm spechla_pipeline:v1 bash -c "which spechla && which spechla-extract-hla-reads && test -d \"\${CONDA_PREFIX}/share/spechla/db\" && test -d \"\${CONDA_PREFIX}/share/spechla/script\""

#
# Usage (per-sample, from driver):
#   docker run --rm \
#     -v /host/.../reference:/data/reference:ro \
#     -v /host/.../input:/data/input:ro \
#     -v /host/.../outputs/<sample>:/data/output \
#     spechla_pipeline:v1 \
#     bash -c "<the pipeline steps for one sample>"

FROM debian:bookworm

# CONDA_PREFIX is intentionally NOT set here — it would cause `conda tos accept`
# and `conda create` to fail by trying to operate inside an env that doesn't
# exist yet. It's set later, after `conda create`, for runtime use.
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/miniconda3/envs/spechla_env/bin:/opt/miniconda3/bin:$PATH

# --- System dependencies ----------------------------------------------------
# tini: clean signal handling so `docker stop` actually kills mid-slice samtools.
# curl + ca-certificates: HTTPS access to EBI for remote CRAM streaming.
# bash + coreutils: SpecHLA's shell scripts assume a real bash environment.
# procps: ps/top for debugging long-running containers if needed.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        wget git curl bzip2 ca-certificates \
        gcc g++ make \
        libgl1 libegl1 \
        bash coreutils procps \
        less perl file vim-tiny \
        tini \
 && rm -rf /var/lib/apt/lists/* \
 && ln -sf /bin/bash /bin/sh

# --- Miniconda --------------------------------------------------------------
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh \
 && bash /tmp/miniconda.sh -b -p /opt/miniconda3 \
 && rm /tmp/miniconda.sh

# Accept Anaconda terms of service for the default channels (required for
# non-interactive installs on recent conda versions).
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# --- SpecHLA + samtools (pinned) --------------------------------------------
# spechla=1.0.10 brings samtools, bwa, and other dependencies via bioconda.
RUN conda create -n spechla_env -c bioconda -c conda-forge spechla=1.0.10 -y \
 && conda clean -afy

# --- Declare the conda env as active for runtime ----------------------------
# Set CONDA_PREFIX now, AFTER `conda create`, so SpecHLA's wrapper script can
# locate its bundled database and scripts at runtime. Setting this earlier
# would break conda's own commands during the build.
ENV CONDA_PREFIX=/opt/miniconda3/envs/spechla_env

# --- Verification: fail the build if anything is wrong ----------------------
# 1. Confirm SpecHLA is the pinned version.
RUN conda run -n spechla_env conda list spechla \
        | grep -E '^spechla\s+1\.0\.10\b' \
        || (echo "SpecHLA version mismatch" && exit 1)

# 2. Confirm samtools came in and supports HTTPS (libcurl).
RUN samtools --version | head -3 \
 && samtools --version | grep -qi 'libcurl' \
        || (echo "samtools missing libcurl support" && exit 1)

# 3. Confirm SpecHLA's entrypoints are on PATH and the wrapper can locate
#    its bundled database/scripts via CONDA_PREFIX.
RUN which spechla \
 && which spechla-extract-hla-reads \
 && test -d "${CONDA_PREFIX}/share/spechla/db" \
 && test -d "${CONDA_PREFIX}/share/spechla/script"

# --- Working directory ------------------------------------------------------
# This is where the container "stands" when commands run. The driver passes
# absolute paths anyway, so this is mostly for interactive debugging.
WORKDIR /data/output

# --- Entrypoint -------------------------------------------------------------
# tini as PID 1 ensures SIGTERM/SIGINT propagate cleanly to samtools/spechla.
# Without this, `docker stop` waits 10 seconds then SIGKILLs, which can leave
# partial files on disk.
ENTRYPOINT ["/usr/bin/tini", "--"]

# No CMD: caller supplies the full command per `docker run` invocation.

