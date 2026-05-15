#!/usr/bin/env python3
"""
Collect HLA typing outputs from a SpecHLA pipeline directory and write
two CSVs: one from hla.result.txt files, one from hla.result.g.group.txt files.

Expected layout:
    <root>/
        <sample_id>/
            spechla_output/
                <sample_id>/
                    hla.result.txt
                    hla.result.g.group.txt

The outer and inner <sample_id> directory names must match.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

try:
    from dotenv import load_dotenv
except ImportError:
    print(
        "[error] python-dotenv is required; install with: pip install python-dotenv",
        file=sys.stderr,
    )
    sys.exit(1)

# Expected column order in both result files
HLA_COLUMNS = [
    "HLA_A_1", "HLA_A_2",
    "HLA_B_1", "HLA_B_2",
    "HLA_C_1", "HLA_C_2",
    "HLA_DPA1_1", "HLA_DPA1_2",
    "HLA_DPB1_1", "HLA_DPB1_2",
    "HLA_DQA1_1", "HLA_DQA1_2",
    "HLA_DQB1_1", "HLA_DQB1_2",
    "HLA_DRB1_1", "HLA_DRB1_2",
]

RESULT_FILE = "hla.result.txt"
GGROUP_FILE = "hla.result.g.group.txt"


class SampleMismatchError(Exception):
    """Raised when a file's contents don't line up with its sample folder."""


def parse_hla_file(file_path: Path, sample_id: str) -> dict:
    """
    Parse an HLA result file and return a row dict for the CSV.

    Empty row (sample id + empty cells) is returned when:
      - the file is missing
      - the file is empty or contains only comments/blank lines
      - the header is present but there is no data row

    SampleMismatchError is raised when:
      - the header and data row have different column counts
      - the Sample cell inside the file disagrees with the folder name
      - any allele value doesn't start with the expected gene prefix
        (e.g. an 'A*...' value appearing in an HLA_B_* column)
    """
    # Start with an empty row; we fill it in if the file has data.
    row = {"Sample": sample_id}
    for col in HLA_COLUMNS:
        row[col] = ""

    if not file_path.is_file():
        print(f"  [warn] missing file: {file_path}", file=sys.stderr)
        return row

    # Keep only non-blank, non-comment lines.
    lines = []
    for line in file_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            lines.append(stripped)

    # Need at least a header line and a data line.
    if len(lines) < 2:
        return row

    headers = [h.strip() for h in lines[0].split("\t")]
    values = [v.strip() for v in lines[1].split("\t")]

    if len(headers) != len(values):
        raise SampleMismatchError(
            f"column count mismatch in {file_path}: header has "
            f"{len(headers)} fields, data row has {len(values)}"
        )

    file_map = dict(zip(headers, values))

    # If headers contained duplicates, dict() collapsed them and file_map
    # is now shorter than headers — meaning some columns were silently lost.
    if len(file_map) != len(headers):
        raise SampleMismatchError(
            f"duplicate column names in {file_path}: {headers}"
        )

    # The Sample cell in the file must agree with the folder name.
    file_sample = file_map.get("Sample", "")
    if file_sample and file_sample != sample_id:
        raise SampleMismatchError(
            f"sample id in file ({file_sample}) does not match "
            f"folder ({sample_id}) for {file_path}"
        )

    # Copy over each HLA column, validating the gene prefix.
    for col in HLA_COLUMNS:
        value = file_map.get(col, "")
        row[col] = value

        # Expected gene prefix: HLA_A_1 -> "A", HLA_DPB1_2 -> "DPB1".
        prefix = col[len("HLA_"):].rsplit("_", 1)[0]

        # Each slash-separated allele must start with that prefix.
        # Empty cells and "-" are valid "no value" markers.
        for part in value.split("/"):
            part = part.strip()
            if part == "" or part == "-":
                continue
            if part.split("*", 1)[0] != prefix:
                raise SampleMismatchError(
                    f"value '{value}' in column {col} of {file_path} does "
                    f"not start with expected gene prefix '{prefix}*'"
                )

    return row


def discover_samples(outputs_dir: Path, markers_dir: Path) -> list:
    """
    List samples to process based on '<sample_id>.done.json' files in
    markers_dir. Files ending in anything else (e.g. '.failed.json') are
    skipped. Returns a sorted list of (sample_id, sample_inner_dir) tuples.
    Raises ValueError if a sample is marked done but its output folder
    does not exist.
    """
    if not outputs_dir.is_dir():
        raise FileNotFoundError(f"outputs directory does not exist: {outputs_dir}")
    if not markers_dir.is_dir():
        raise FileNotFoundError(f"markers directory does not exist: {markers_dir}")

    samples = []
    for marker in sorted(markers_dir.iterdir()):
        if not marker.is_file() or not marker.name.endswith(".done.json"):
            if marker.is_file():
                print(f"[info] skipping non-done marker: {marker.name}", file=sys.stderr)
            continue

        sample_id = marker.name[: -len(".done.json")]
        inner = outputs_dir / sample_id / "spechla_output" / sample_id
        if not inner.is_dir():
            raise ValueError(
                f"sample {sample_id} is marked done but expected folder "
                f"{inner} does not exist"
            )

        samples.append((sample_id, inner))

    return samples


def write_csv(rows: list, out_path: Path) -> None:
    fieldnames = ["Sample"] + HLA_COLUMNS
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Collect SpecHLA outputs into CSVs. Reads outputs/ and "
            "work/markers/ from $ROOT_DIR_FULL_PIPELINE."
        )
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory where the two CSVs are written. Defaults to the "
             "directory containing this script.",
    )
    args = parser.parse_args()

    load_dotenv()
    pipeline_root = os.environ.get("ROOT_DIR_FULL_PIPELINE")
    if not pipeline_root:
        print("[error] ROOT_DIR_FULL_PIPELINE is not set", file=sys.stderr)
        return 1
    pipeline_root = Path(pipeline_root)
    outputs_dir = pipeline_root / "outputs"
    markers_dir = pipeline_root / "work" / "markers"

    try:
        samples = discover_samples(outputs_dir, markers_dir)
    except (FileNotFoundError, ValueError) as e:
        print(f"[error] {e}", file=sys.stderr)
        return 1

    if not samples:
        print("[error] no .done.json markers found", file=sys.stderr)
        return 1

    result_rows = []
    ggroup_rows = []

    try:
        for sample_id, inner_dir in samples:
            print(f"Processing {sample_id}")
            result_rows.append(
                parse_hla_file(inner_dir / RESULT_FILE, sample_id)
            )
            ggroup_rows.append(
                parse_hla_file(inner_dir / GGROUP_FILE, sample_id)
            )
    except SampleMismatchError as e:
        print(f"[error] {e}", file=sys.stderr)
        return 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    result_csv = args.out_dir / "hla_results.csv"
    ggroup_csv = args.out_dir / "hla_results_g_group.csv"
    write_csv(result_rows, result_csv)
    write_csv(ggroup_rows, ggroup_csv)

    print(f"\nWrote {len(result_rows)} rows to {result_csv}")
    print(f"Wrote {len(ggroup_rows)} rows to {ggroup_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
    