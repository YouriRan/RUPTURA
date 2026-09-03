#!/usr/bin/env python3
"""Run the generated Ruptura benchmarks with eight concurrent processes."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


BENCHMARK_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = BENCHMARK_ROOT.parent
DEFAULT_EXECUTABLE = REPOSITORY_ROOT / "build" / "ruptura"
MANIFEST_PATH = BENCHMARK_ROOT / "manifest.csv"


@dataclass(frozen=True)
class RunResult:
    directory: Path
    return_code: int
    elapsed_seconds: float


def simulation_directories(case: str | None, include_sirk3: bool) -> list[Path]:
    if not MANIFEST_PATH.is_file():
        return []

    with MANIFEST_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    selected = []
    for row in rows:
        if not include_sirk3 and row["integrator"] == "sirk3":
            continue
        if case is not None and row["base_case"] != case:
            continue
        selected.append(BENCHMARK_ROOT / row["directory"])
    return sorted(selected)


def run_simulation(directory: Path, executable: Path) -> RunResult:
    started = time.perf_counter()
    log_path = directory / "run.log"

    try:
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(executable)],
                cwd=directory,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        return_code = completed.returncode
    except OSError as error:
        log_path.write_text(f"Could not start {executable}: {error}\n", encoding="utf-8")
        return_code = 127

    return RunResult(directory, return_code, time.perf_counter() - started)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="maximum number of simultaneous Ruptura processes (default: 8)",
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=DEFAULT_EXECUTABLE,
        help=f"Ruptura executable (default: {DEFAULT_EXECUTABLE})",
    )
    parser.add_argument(
        "--case",
        choices=("bea-alkanes-C7",),
        help="run only one base case",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list the selected simulations without running them",
    )
    parser.add_argument(
        "--include-sirk3",
        action="store_true",
        help="include the currently non-runnable SIRK3 configurations",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")

    executable = args.executable.expanduser().resolve()
    directories = simulation_directories(args.case, args.include_sirk3)
    if not directories:
        raise SystemExit("No benchmark simulations found; run benchmark/generate.py first")

    if args.dry_run:
        for directory in directories:
            print(directory.relative_to(BENCHMARK_ROOT))
        print(f"Selected {len(directories)} simulations")
        return 0

    if not executable.is_file():
        raise SystemExit(f"Ruptura executable not found: {executable}")

    results_path = BENCHMARK_ROOT / "results.csv"
    selected_directories = {
        directory.relative_to(BENCHMARK_ROOT).as_posix() for directory in directories
    }
    preserved_rows = []
    if args.case is not None and results_path.is_file():
        with results_path.open(newline="", encoding="utf-8") as result_stream:
            preserved_rows = [
                row
                for row in csv.DictReader(result_stream)
                if row["directory"] not in selected_directories
            ]

    failures = 0
    with results_path.open("w", newline="", encoding="utf-8") as result_stream:
        writer = csv.DictWriter(
            result_stream,
            fieldnames=("directory", "return_code", "elapsed_seconds"),
        )
        writer.writeheader()
        writer.writerows(preserved_rows)

        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(run_simulation, directory, executable): directory
                for directory in directories
            }
            for completed_count, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                relative_directory = result.directory.relative_to(BENCHMARK_ROOT)
                status = "ok" if result.return_code == 0 else "failed"
                failures += result.return_code != 0

                writer.writerow(
                    {
                        "directory": relative_directory.as_posix(),
                        "return_code": result.return_code,
                        "elapsed_seconds": f"{result.elapsed_seconds:.6f}",
                    }
                )
                result_stream.flush()
                print(
                    f"[{completed_count}/{len(directories)}] {status}: "
                    f"{relative_directory} ({result.elapsed_seconds:.2f} s)"
                )

    print(f"Finished {len(directories)} simulations with {failures} failures")
    print(f"Summary: {results_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
