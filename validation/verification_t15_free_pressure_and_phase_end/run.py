#!/usr/bin/env python3
"""Run every simulation of simulations/index.json, 8 at a time.

    python run.py                 # build/ruptura
    python run.py --binary PATH   # another Ruptura binary, e.g. a frozen copy while build/ruptura is rebuilt
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUPTURA_BINARY = ROOT.parent.parent / "build" / "ruptura"


def run(directory: Path, binary: Path) -> str:
    started = time.perf_counter()
    with (directory / "run.log").open("w", encoding="utf-8") as log:
        code = subprocess.run([str(binary)], cwd=directory, stdout=log, stderr=subprocess.STDOUT).returncode
    return f"{'ok' if code == 0 else f'failed ({code})'} {time.perf_counter() - started:7.1f} s  {directory.relative_to(ROOT)}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path, default=RUPTURA_BINARY, help="Ruptura binary (default: build/ruptura)")
    arguments = parser.parse_args()
    directories = [ROOT / entry["directory"] for entry in json.loads((ROOT / "simulations" / "index.json").read_text())]
    with ThreadPoolExecutor(max_workers=8) as pool:
        for line in pool.map(lambda directory: run(directory, arguments.binary.resolve()), directories):
            print(line, flush=True)


if __name__ == "__main__":
    main()
