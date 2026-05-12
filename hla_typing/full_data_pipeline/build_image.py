#!/usr/bin/env python3
"""Build the Docker image for the HLA toolkit pipeline."""
import os
import subprocess
import sys
from pathlib import Path
from dotenv import load_dotenv

def main() -> int:
    load_dotenv()

    required = [
        "HLA_TOOLKIT_BUILD_IMAGE_NAME",
        "HLA_TOOLKIT_IMAGE_TAG",
        "HLA_TOOLKIT_DOCKERFILE_NAME",
    ]
    missing = [k for k in required if not os.environ.get(k)]
    if missing:
        print(f"ERROR: missing in .env: {', '.join(missing)}", file=sys.stderr)
        return 1

    image = f"{os.environ['HLA_TOOLKIT_BUILD_IMAGE_NAME']}:{os.environ['HLA_TOOLKIT_IMAGE_TAG']}"
    dockerfile_dir = Path(__file__).parent
    dockerfile_path = dockerfile_dir / os.environ["HLA_TOOLKIT_DOCKERFILE_NAME"]

    print(f"Building {image}...")
    result = subprocess.run(
        ["docker", "build", "--no-cache", "-f", str(dockerfile_path), "-t", image, str(dockerfile_dir)],
    )
    if result.returncode != 0:
        print(f"ERROR: build failed for {image}", file=sys.stderr)
        return result.returncode

    print(f"Successfully built {image}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

