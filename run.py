#!/usr/bin/env python3
"""Rerun every simulation in the repository with the current code, 8 at a time.

    python run.py                       # build build/ruptura, then rerun everything
    nohup python run.py > rerun.out 2>&1 &
    python run.py --dry-run             # list what would run, longest first, with the expected core time
    python run.py --resume              # after an interruption: skip what already finished in rerun_results.csv
    python run.py --only liu marx2      # only simulations under these paths
    python run.py --skip-benchmarks     # leave out the timed benchmark studies

Every directory with a simulation.json is run in place, outside .worktrees, build*, and EXCLUDED. The binary is
built first and frozen as build/ruptura-rerun, so a rebuild of build/ruptura during the run does not mix binaries.
Each input is validated first; invalid inputs are reported and skipped.

- A directory with a run.sh (the Dantas runs) is run through that script, which prepares its restart files and
  calls build/ruptura itself; its own output goes to run.sh.log.
- Any other simulation writes its output to the log its directory already uses (ruptura.log or run.log), so the
  per-project run.py scripts still find their results.
- The timed benchmark studies (TIMED_STUDIES) are run afterwards through their own run.py with the frozen binary,
  so that their wall times are measured without other load and recorded in their results.csv.

Longest runs start first (as estimated from the "Total time" of the previous log). A run longer than
--timeout-hours is killed. Every finished run is appended to rerun_results.csv.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUILD_BINARY = ROOT / "build" / "ruptura"
FROZEN_BINARY = ROOT / "build" / "ruptura-rerun"
RESULTS_PATH = ROOT / "rerun_results.csv"
# Marked for deletion, and the run records of Ruptura Lab.
EXCLUDED = ("stampi-bombelli/_to_delete", "simulations")
# Timed studies with their own run.py, and its option for the binary.
TIMED_STUDIES = {
    "benchmark/breakthrough": "--binary",
    "benchmark/cvode-max-step": "--binary",
    "benchmark/convection": "--binary",
    "benchmark/figures": "--executable",
}
# Expected core time of a simulation without a previous log, for the run order only.
DEFAULT_SECONDS = 600.0
FIELDS = ("directory", "status", "return_code", "elapsed_seconds", "finished_at")
STAMP = "%Y-%m-%d %H:%M:%S"


def run_job(directory: Path, timeout: float | None, running: set, lock: threading.Lock) -> dict:
    """Run one simulation directory in its own process group, so that a timeout also stops run.sh's children."""
    if (directory / "run.sh").is_file():
        command, log_name = ["sh", "run.sh"], "run.sh.log"
    else:
        command, log_name = [str(FROZEN_BINARY)], "ruptura.log" if (directory / "ruptura.log").is_file() else "run.log"
    started = time.perf_counter()
    with (directory / log_name).open("w") as log:
        process = subprocess.Popen(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        with lock:
            running.add(process)
        try:
            return_code = process.wait(timeout=timeout)
            status = "ok" if return_code == 0 else "failed"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            return_code = process.wait()
            status = "timeout"
            log.write(f"\nKilled by run.py after {timeout / 3600:g} h; the simulation is unfinished.\n")
        finally:
            with lock:
                running.discard(process)
    return {"directory": str(directory.relative_to(ROOT)), "status": status, "return_code": return_code,
            "elapsed_seconds": f"{time.perf_counter() - started:.1f}",
            "finished_at": datetime.now().strftime(STAMP)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-j", "--workers", type=int, default=8)
    parser.add_argument("--only", nargs="+", default=[], help="run only directories under these paths")
    parser.add_argument("--timeout-hours", type=float, default=16.0, help="kill a run after this long (0: never)")
    parser.add_argument("--resume", action="store_true", help="skip directories that finished in rerun_results.csv")
    parser.add_argument("--skip-benchmarks", action="store_true", help="leave out the timed benchmark studies")
    parser.add_argument("--no-build", action="store_true", help="use build/ruptura as it is")
    parser.add_argument("--dry-run", action="store_true", help="list the runs, longest first, and stop")
    arguments = parser.parse_args()
    only = [Path(path).as_posix().rstrip("/") for path in arguments.only]

    if not arguments.no_build and not arguments.dry_run:
        print(f"{datetime.now().strftime(STAMP)}  building build/ruptura", flush=True)
        subprocess.run(["cmake", "--build", "build", "--target", "ruptura"], cwd=ROOT, check=True)
    if not BUILD_BINARY.is_file():
        raise SystemExit(f"Ruptura binary not found at {BUILD_BINARY}; build it first")
    if not arguments.dry_run:
        shutil.copy2(BUILD_BINARY, FROZEN_BINARY)

    # Simulation directories, pruning the worktrees, build trees, and the excluded and timed paths.
    directories = []
    skipped_prefixes = EXCLUDED + tuple(TIMED_STUDIES)
    for current, subdirectories, files in os.walk(ROOT):
        relative = Path(current).relative_to(ROOT).as_posix()
        subdirectories[:] = sorted(
            name for name in subdirectories
            if not name.startswith((".", "build", "__pycache__", "node_modules"))
            and (Path(relative) / name).as_posix() not in skipped_prefixes)
        if "simulation.json" in files and relative != "." and (
                not only or any(relative == path or relative.startswith(path + "/") for path in only)):
            directories.append(Path(current))
        if "run.sh" in files:
            subdirectories[:] = []  # run.sh runs its nested simulations (the Dantas preconditioning) itself
    binary = FROZEN_BINARY if not arguments.dry_run else BUILD_BINARY

    # Validate every input with the binary, in parallel; invalid inputs are skipped.
    def validate(directory: Path) -> tuple[Path, str]:
        result = subprocess.run([str(binary), "--validate", "simulation.json"], cwd=directory, capture_output=True,
                                text=True)
        output = (result.stdout + result.stderr).strip()
        if result.returncode == 0 and "Input is valid" in output:
            return directory, ""
        return directory, (output.splitlines() or [f"validation failed with return code {result.returncode}"])[-1]

    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        validation = list(executor.map(validate, directories))
    invalid = [(directory, message) for directory, message in validation if message]
    valid = [directory for directory, message in validation if not message]

    finished = set()
    if arguments.resume and RESULTS_PATH.is_file():
        with RESULTS_PATH.open(newline="") as stream:
            finished = {row["directory"] for row in csv.DictReader(stream) if row["status"] in ("ok", "timeout")}
    pending = [directory for directory in valid if str(directory.relative_to(ROOT)) not in finished]

    # Longest first, from the "Total time" of the previous log, so that the pool does not end with one long run.
    total_time = re.compile(r"Total time:\s*([0-9.eE+-]+)\s*s")
    expected = {}
    for directory in pending:
        times = [float(match) for name in ("run.log", "ruptura.log", "nohup.out") if (directory / name).is_file()
                 for match in total_time.findall((directory / name).read_text(errors="ignore")[-4000:])]
        expected[directory] = times[-1] if times else DEFAULT_SECONDS
    pending.sort(key=lambda directory: expected[directory], reverse=True)
    studies = [study for study in TIMED_STUDIES if not arguments.skip_benchmarks
               and (not only or any(study == path or study.startswith(path + "/") for path in only))
               and not (arguments.resume and study in finished)]

    print(f"{datetime.now().strftime(STAMP)}  {len(directories)} simulations found: {len(pending)} to run, "
          f"{len(invalid)} invalid, {len(valid) - len(pending)} already finished; expected core time "
          f"{sum(expected.values()) / 3600:.1f} h on {arguments.workers} workers, then {len(studies)} timed "
          f"benchmark studies", flush=True)
    for directory, message in invalid:
        print(f"  invalid, skipped: {directory.relative_to(ROOT)}: {message}", flush=True)
    if arguments.dry_run:
        for directory in pending:
            print(f"  {expected[directory] / 3600:8.2f} h  {directory.relative_to(ROOT)}")
        for study in studies:
            print(f"  benchmark study {study} (its own run.py)")
        return

    if not arguments.resume or not RESULTS_PATH.is_file():
        with RESULTS_PATH.open("w", newline="") as stream:
            csv.DictWriter(stream, fieldnames=FIELDS).writeheader()
    timeout = arguments.timeout_hours * 3600.0 or None
    # A kill (SIGTERM, e.g. of a nohup run) stops the runs like Ctrl-C does.
    signal.signal(signal.SIGTERM, signal.default_int_handler)
    running, lock = set(), threading.Lock()
    results, started, interrupted = [], time.perf_counter(), False
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [executor.submit(run_job, directory, timeout, running, lock) for directory in pending]
        try:
            for count, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                results.append(result)
                with RESULTS_PATH.open("a", newline="") as stream:
                    csv.DictWriter(stream, fieldnames=FIELDS).writerow(result)
                print(f"{datetime.now().strftime(STAMP)}  [{count}/{len(pending)}] {result['status']:7s} "
                      f"{float(result['elapsed_seconds']) / 60:8.1f} min  {result['directory']}", flush=True)
        except KeyboardInterrupt:
            # Leaving the pool waits for its runs, so stop the queue and kill the running process groups, which do not
            # see the interrupt themselves.
            interrupted = True
            for future in futures:
                future.cancel()
            with lock:
                for process in running:
                    os.killpg(process.pid, signal.SIGKILL)
    if interrupted:
        raise SystemExit("interrupted; continue with python run.py --resume")

    for study in studies:
        print(f"{datetime.now().strftime(STAMP)}  timed benchmark study {study}", flush=True)
        study_started = time.perf_counter()
        completed = subprocess.run([sys.executable, "run.py", "--jobs", str(arguments.workers),
                                    TIMED_STUDIES[study], str(FROZEN_BINARY)], cwd=ROOT / study)
        result = {"directory": study, "status": "ok" if completed.returncode == 0 else "failed",
                  "return_code": completed.returncode, "finished_at": datetime.now().strftime(STAMP),
                  "elapsed_seconds": f"{time.perf_counter() - study_started:.1f}"}
        results.append(result)
        with RESULTS_PATH.open("a", newline="") as stream:
            csv.DictWriter(stream, fieldnames=FIELDS).writerow(result)

    print(f"\n{datetime.now().strftime(STAMP)}  done in {(time.perf_counter() - started) / 3600:.1f} h: "
          f"{sum(result['status'] == 'ok' for result in results)} ok, "
          f"{sum(result['status'] == 'failed' for result in results)} failed, "
          f"{sum(result['status'] == 'timeout' for result in results)} timed out, {len(invalid)} invalid")
    for result in results:
        if result["status"] != "ok":
            print(f"  {result['status']}: {result['directory']} (return code {result['return_code']})")
    print(f"Results per simulation: {RESULTS_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
