#!/usr/bin/env python3
"""
Run SpecHLA inside the spechla:1.0.10 Docker image.

Edit the CONFIG block below for your run, then:
    uv run spec_hla_ctr.py
or:
    python spec_hla_ctr.py

The host directory DATA_DIR is bind-mounted at /hla_data inside the container.
READ1, READ2, and OUTPUT_SUBDIR are interpreted relative to DATA_DIR.
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# CONFIG — edit these for your run
# ---------------------------------------------------------------------------

# Host directory containing input FASTQs; mounted into the container at /hla_data.
# Outputs land here too.
HOST_DATA_DIR = Path("ops/test_data")

'''
SAMPLE = "HG00733"
READ1 = "HG00733.final_extract_1.fq.gz"   # relative to HOST_DATA_DIR
READ2 = "HG00733.final_extract_2.fq.gz"   # relative to HOST_DATA_DIR
OUTPUT_SUBDIR = "output"                  # relative to HOST_DATA_DIR
'''
SAMPLE = "HG00096"
READ1 = "in/HG00096/HG00096_HLA.R1.fastq.gz"   # relative to HOST_DATA_DIR
READ2 = "in/HG00096/HG00096_HLA.R2.fastq.gz"   # relative to HOST_DATA_DIR
OUTPUT_SUBDIR = "out"                  # relative to HOST_DATA_DIR

THREADS = 5

CTR_MOUNT_POINT = "/hla_data"

# Image to run.
# IMAGE = "spechla:1.0.10"
IMAGE = "dckr_spechla:v5"

# Extra args appended to the spechla command, e.g. ["-p", "Asian"] or ["-u", "1"].
EXTRA_ARGS: list[str] = []

# ---------------------------------------------------------------------------


def main() -> int:
    host_data_dir = HOST_DATA_DIR.resolve()
    if not host_data_dir.is_dir():
        sys.exit(f"HOST_DATA_DIR does not exist: {host_data_dir}")

    for label, name in (("READ1", READ1), ("READ2", READ2)):
        if not (host_data_dir / name).is_file():
            sys.exit(f"{label} not found: {host_data_dir / name}")

    (host_data_dir / OUTPUT_SUBDIR).mkdir(parents=True, exist_ok=True)

    # The command we want to run inside the container.
    spechla_cmd = [
        "conda", "run", "--no-capture-output", "-n", "spechla_env",
        "spechla",
        "-j", str(THREADS),
        "-n", SAMPLE,
        "-1", f"{CTR_MOUNT_POINT}/{READ1}",
        "-2", f"{CTR_MOUNT_POINT}/{READ2}",
        "-o", f"{CTR_MOUNT_POINT}/{OUTPUT_SUBDIR}",
        *EXTRA_ARGS,
    ]

    # Run as the host user so outputs aren't owned by root (Linux/macOS).
    user_flag = [f"--user={os.getuid()}:{os.getgid()}"] if hasattr(os, "getuid") else []

    docker_cmd = [
        "docker", "run", "--rm",
        *user_flag,
        "-v", f"{host_data_dir}:{CTR_MOUNT_POINT}",
        IMAGE,
        *spechla_cmd,
    ]

    print("$", " ".join(shlex.quote(c) for c in docker_cmd), flush=True)
    result = subprocess.run(docker_cmd)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())