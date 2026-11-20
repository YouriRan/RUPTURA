#!/usr/bin/env python3
"""Run the full-breakthrough benchmark and record the wall time of every simulation.

    python run.py                 # run every simulation in manifest.csv, 8 at a time
    python run.py --jobs 4        # fewer concurrent runs, for less contention in the timings
    python run.py --match methods # only directories containing this substring
    python run.py --match a b     # only directories containing any of the substrings
    python run.py --binary PATH   # another Ruptura binary, e.g. a frozen copy while build/ruptura is rebuilt
    python run.py --repeats 3     # time every run three times: elapsed_seconds is the median, elapsed_all all

Wall times and return codes go to results.csv; each run writes run.log in its directory.  Rows of
simulations that are not rerun are kept.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SIMULATION_DIR = ROOT / "simulations"
MANIFEST_PATH = ROOT / "manifest.csv"
RESULTS_PATH = ROOT / "results.csv"
RUPTURA_BINARY = ROOT.parents[2] / "build" / "ruptura"


def run_simulation(row: dict, binary: Path, repeats: int) -> dict:
    directory = SIMULATION_DIR / row["directory"]
    elapsed = []
    for _ in range(repeats):
        started = time.perf_counter()
        with (directory / "run.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run([str(binary)], cwd=directory, stdout=log, stderr=subprocess.STDOUT)
        elapsed.append(time.perf_counter() - started)
        if completed.returncode != 0:
            break
    return {
        "directory": row["directory"],
        "return_code": completed.returncode,
        "elapsed_seconds": f"{statistics.median(elapsed):.6f}",
        "elapsed_all": " ".join(f"{seconds:.6f}" for seconds in elapsed),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--match", nargs="+", help="run only directories containing any of these substrings")
    parser.add_argument("--binary", type=Path, default=RUPTURA_BINARY, help="Ruptura binary (default: build/ruptura)")
    parser.add_argument("--repeats", type=int, default=1, help="time every simulation this many times")
    arguments = parser.parse_args()
    if arguments.repeats < 1:
        raise SystemExit("--repeats must be at least 1")
    if not arguments.binary.is_file():
        raise SystemExit(f"Ruptura binary not found at {arguments.binary}")

    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for row in rows if not arguments.match or any(match in row["directory"] for match in arguments.match)]
    if not selected:
        raise SystemExit("no matching simulations; run generate.py first")
    # Longest runs first, so that the pool does not end with one long run on its own.
    selected.sort(key=lambda row: float(row["time_step"]))

    preserved = []
    if RESULTS_PATH.is_file():
        chosen = {row["directory"] for row in selected}
        with RESULTS_PATH.open(newline="", encoding="utf-8") as stream:
            preserved = [row for row in csv.DictReader(stream) if row["directory"] not in chosen]

    results = []
    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = [executor.submit(run_simulation, row, arguments.binary, arguments.repeats) for row in selected]
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            status = "ok" if result["return_code"] == 0 else f"FAILED ({result['return_code']})"
            print(f"[{count}/{len(selected)}] {status} {float(result['elapsed_seconds']):8.1f} s  {result['directory']}",
                  flush=True)

    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("directory", "return_code", "elapsed_seconds", "elapsed_all"))
        writer.writeheader()
        writer.writerows(preserved)
        writer.writerows(results)
    failures = sum(int(result["return_code"]) != 0 for result in results)
    print(f"\n{len(results) - failures}/{len(results)} simulations completed")


if __name__ == "__main__":
    main()
