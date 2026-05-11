# SpecHLA pipeline image for HLA typing on remote 1kGP CRAMs.
#
# Build:
#   docker build -f toolkit.dockerfile -t hla_toolkit:v1 .

FROM debian:bookworm@sha256:8a8cd02c5912770b4980228a54d4aff9e4f986f1eb2525d2d371dec5232cefcc

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/miniconda3/envs/spechla_env/bin:/opt/miniconda3/bin:$PATH

# --- System dependencies ----------------------------------------------------
# Unpinned by design: these don't affect SpecHLA's output. The conda env is
# what's locked down. If apt ever drifts in a way that breaks the build,
# that's a build problem, not a reproducibility problem.
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

# --- Miniconda (pinned, sha256-verified) ------------------------------------
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-py312_26.3.2-2-Linux-x86_64.sh \
        -O /tmp/miniconda.sh \
 && echo "32673413a39a21ae3997c9b38236e9df15c9fcef930b510487c64fe259e03f95  /tmp/miniconda.sh" \
        | sha256sum -c - \
 && bash /tmp/miniconda.sh -b -p /opt/miniconda3 \
 && rm /tmp/miniconda.sh

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# --- SpecHLA env (pinned from lockfile) -------------------------------------
# spechla_env.yml is a `conda env export --no-builds` snapshot. Every
# scientifically-relevant package version is fixed here.
COPY spechla_env.yml /tmp/spechla_env.yml
RUN conda env create -n spechla_env -f /tmp/spechla_env.yml \
 && conda clean -afy

ENV CONDA_PREFIX=/opt/miniconda3/envs/spechla_env

# --- Verification: fail the build if anything is wrong ----------------------
RUN conda run -n spechla_env conda list spechla \
        | grep -E '^spechla\s+1\.0\.10\b' \
        || (echo "SpecHLA version mismatch" && exit 1)

RUN samtools --version | head -3 \
 && { samtools --version | grep -qi 'libcurl' \
        || { echo "samtools missing libcurl support"; exit 1; }; }

RUN which spechla \
 && which spechla-extract-hla-reads \
 && test -d "${CONDA_PREFIX}/share/spechla/db" \
 && test -d "${CONDA_PREFIX}/share/spechla/script"

WORKDIR /data/output

ENTRYPOINT ["/usr/bin/tini", "--"]

