#!/usr/bin/env python3
"""Run the convection-scheme study and record the wall and CPU time of every simulation.

    python run.py                    # run every simulation in manifest.csv, 8 at a time
    python run.py --jobs 4           # fewer concurrent runs, for less contention in the timings
    python run.py --match cost       # only directories containing this substring
    python run.py --binary PATH      # another Ruptura binary, e.g. a frozen copy while build/ruptura is rebuilt

Wall times, CPU times (user + system of the run, from /usr/bin/time), and return codes go to results.csv; each run
writes run.log in its directory.  Rows of simulations that are not rerun are kept.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SIMULATION_DIR = ROOT / "simulations"
MANIFEST_PATH = ROOT / "manifest.csv"
RESULTS_PATH = ROOT / "results.csv"
RUPTURA_BINARY = ROOT.parents[2] / "build" / "ruptura"
FIELDS = ("directory", "return_code", "elapsed_seconds", "cpu_seconds")


def run_simulation(row: dict, binary: Path) -> dict:
    directory = SIMULATION_DIR / row["directory"]
    for path in directory.glob("*.data"):
        path.unlink()
    started = time.perf_counter()
    with (directory / "run.log").open("w", encoding="utf-8") as log:
        completed = subprocess.run(["/usr/bin/time", "-p", str(binary)], cwd=directory, stdout=log,
                                   stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    cpu = 0.0
    for line in (directory / "run.log").read_text(encoding="utf-8").splitlines():
        words = line.split()
        if len(words) == 2 and words[0] in ("user", "sys"):
            cpu += float(words[1])
    return {"directory": row["directory"], "return_code": completed.returncode,
            "elapsed_seconds": f"{elapsed:.6f}", "cpu_seconds": f"{cpu:.2f}"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--match", nargs="+", help="run only directories containing any of these substrings")
    parser.add_argument("--binary", type=Path, default=RUPTURA_BINARY, help="Ruptura binary (default: build/ruptura)")
    arguments = parser.parse_args()
    if not arguments.binary.is_file():
        raise SystemExit(f"Ruptura binary not found at {arguments.binary}")

    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for row in rows if not arguments.match or any(match in row["directory"] for match in arguments.match)]
    if not selected:
        raise SystemExit("no matching simulations; run generate.py first")
    # Longest runs first (references, then the finest grids), so that the pool does not end with one long run.
    selected.sort(key=lambda row: (row["study"] != "reference", -int(row["grid_points"]), row["integrator"] != "rk3"))

    preserved = []
    if RESULTS_PATH.is_file():
        chosen = {row["directory"] for row in selected}
        with RESULTS_PATH.open(newline="", encoding="utf-8") as stream:
            preserved = [row for row in csv.DictReader(stream) if row["directory"] not in chosen]

    results = []
    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = [executor.submit(run_simulation, row, arguments.binary) for row in selected]
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            status = "ok" if result["return_code"] == 0 else f"FAILED ({result['return_code']})"
            print(f"[{count}/{len(selected)}] {status} {float(result['elapsed_seconds']):8.1f} s  {result['directory']}",
                  flush=True)

    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(preserved)
        writer.writerows(results)
    failures = sum(int(result["return_code"]) != 0 for result in results)
    print(f"{len(results) - failures} succeeded, {failures} failed; results in {RESULTS_PATH}")


if __name__ == "__main__":
    main()
