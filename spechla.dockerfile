# SpecHLA 1.0.10 in a container.
# Build:  docker build -t spechla:1.0.10 .
# Build: docker build -f spechla.dockerfile -t dckr_spechla:v5 .

FROM debian:bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/miniconda3/bin:$PATH

# --- System dependencies ----------------------------------------------------
# Includes a generous bundle of utilities that SpecHLA's shell scripts assume
# are present on a typical Linux system (less, procps, bash, perl, file, vim-tiny).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        wget git curl bzip2 ca-certificates \
        gcc g++ make \
        libgl1 libegl1 \
        less procps bash perl file vim-tiny \
 && rm -rf /var/lib/apt/lists/*

# --- Miniconda --------------------------------------------------------------
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
 && bash /tmp/miniconda.sh -b -p /opt/miniconda3 \
 && rm /tmp/miniconda.sh

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# --- SpecHLA (pinned) -------------------------------------------------------
# Version pinned to 1.0.10 for reproducibility.
RUN conda create -n spechla_env -c bioconda -c conda-forge spechla=1.0.10 -y \
 && conda clean -afy

# Verify the pinned version is what landed in the env.
RUN conda run -n spechla_env conda list spechla | grep -E '^spechla\s+1\.0\.10\b'

WORKDIR /hla_data

# No ENTRYPOINT / CMD: the caller (Python script) supplies the full command,
# wrapping it in `conda run -n spechla_env ...` as needed.
