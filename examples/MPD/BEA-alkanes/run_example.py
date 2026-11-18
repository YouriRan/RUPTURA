#!/usr/bin/env python3
"""Generate, run, and compare the BEA alkane MPD validation example."""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=REPOSITORY / "build" / "ruptura")
    arguments = parser.parse_args()
    binary = arguments.binary.expanduser().resolve()
    if not binary.is_file():
        parser.error(f"Ruptura binary not found: {binary}")

    subprocess.run([sys.executable, ROOT / "generate_particle_distribution.py"], check=True)
    for directory in (ROOT / "mpd", ROOT / "langmuir"):
        print(f"Running {directory.name} example", flush=True)
        subprocess.run([binary], cwd=directory, check=True)
    subprocess.run([sys.executable, ROOT / "compare_results.py"], check=True)


if __name__ == "__main__":
    main()
