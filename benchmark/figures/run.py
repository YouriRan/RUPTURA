#!/usr/bin/env python3
"""Run the generated figure simulations and record repeated wall timings."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


FIGURE_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = FIGURE_ROOT.parents[1]
SIMULATION_ROOT = FIGURE_ROOT / "simulations"
MANIFEST_PATH = FIGURE_ROOT / "manifest.csv"
RESULTS_PATH = FIGURE_ROOT / "results.csv"
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "ruptura"


@dataclass(frozen=True)
class RepeatedResult:
    directory: str
    rows: list[dict]


def selected_rows(suite: str | None, match: str | None) -> list[dict]:
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if suite:
        rows = [row for row in rows if row["suite"] == suite]
    if match:
        rows = [row for row in rows if match in row["directory"]]
    return rows


def run_repeated(row: dict, executable: Path, repeats: int) -> RepeatedResult:
    directory = SIMULATION_ROOT / row["directory"]
    results = []
    for repetition in range(repeats):
        started = time.perf_counter()
        with (directory / f"run-{repetition}.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(executable)], cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        results.append(
            {
                "suite": row["suite"],
                "directory": row["directory"],
                "repetition": repetition,
                "return_code": completed.returncode,
                "elapsed_seconds": f"{time.perf_counter() - started:.9f}",
            }
        )
    return RepeatedResult(row["directory"], results)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("cvode-rk3", "methods"))
    parser.add_argument("--match", help="run only directories containing this substring")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--executable", type=Path, default=DEFAULT_EXECUTABLE)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1 or args.repeats < 1:
        raise SystemExit("--jobs and --repeats must be at least one")
    rows = selected_rows(args.suite, args.match)
    if not rows:
        raise SystemExit("No matching simulations; run generate.py first")
    if args.dry_run:
        for row in rows:
            print(row["directory"])
        print(f"Selected {len(rows)} simulations x {args.repeats} repetitions")
        return 0

    executable = args.executable.expanduser().resolve()
    selected_directories = {row["directory"] for row in rows}
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        current_directories = {row["directory"] for row in csv.DictReader(stream)}
    preserved = []
    if RESULTS_PATH.is_file():
        with RESULTS_PATH.open(newline="", encoding="utf-8") as stream:
            preserved = [
                row
                for row in csv.DictReader(stream)
                if row["directory"] in current_directories and row["directory"] not in selected_directories
            ]

    failures = 0
    completed_rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_repeated, row, executable, args.repeats): row for row in rows}
        for count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            completed_rows.extend(result.rows)
            case_failures = sum(int(row["return_code"]) != 0 for row in result.rows)
            failures += case_failures
            status = "ok" if not case_failures else f"{case_failures} failed"
            print(f"[{count}/{len(rows)}] {status}: {result.directory}", flush=True)

    with RESULTS_PATH.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = ("suite", "directory", "repetition", "return_code", "elapsed_seconds")
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(preserved)
        writer.writerows(completed_rows)

    print(f"Finished {len(rows)} simulations x {args.repeats} repetitions with {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
