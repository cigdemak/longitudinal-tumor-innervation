#!/usr/bin/env python3
"""Compare regenerated PNG files with the reference manuscript exports."""

import argparse
import csv
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    args = parser.parse_args()
    with (ROOT / "manifests/reference_png_sha256.csv").open(newline="") as handle:
        references = list(csv.DictReader(handle))
    failures = []
    for row in references:
        path = args.output_dir / "figures" / row["filename"]
        if not path.exists():
            failures.append(f"MISSING: {row['filename']}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            failures.append(f"DIFFERENT: {row['filename']}")
    for message in failures:
        print(message)
    print(f"{len(references) - len(failures)}/{len(references)} PNG exports match the reference bytes.")
    if failures:
        print("This is an exact image-file check, not a comparison of statistical results. Rendering differences also change checksums.")
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
