"""Run the generated RUPTURA simulations, eight at a time, and record the wall time of each.

Usage: python run.py [name ...]  (default: every generated simulation)
Set RUPTURA_BIN to use another executable than build/ruptura of the repository (other sessions rebuild that binary;
for timed runs use a frozen copy).
"""

import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).parent
RUPTURA = os.environ.get("RUPTURA_BIN", str(HERE.resolve().parents[3] / "build" / "ruptura"))


def run(directory):
    start = time.perf_counter()
    with open(directory / "run.log", "w") as log:
        result = subprocess.run([RUPTURA], cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=False)
    seconds = time.perf_counter() - start
    (directory / "wall_time.txt").write_text(f"{seconds:.3f}\n")
    return directory.name, result.returncode, seconds


if __name__ == "__main__":
    directories = sorted(d for d in HERE.iterdir() if (d / "simulation.json").is_file()
                         and (len(sys.argv) == 1 or d.name in sys.argv[1:]))
    print(f"{len(directories)} simulations with {RUPTURA}", flush=True)
    with ThreadPoolExecutor(max_workers=8) as executor:
        for future in as_completed([executor.submit(run, directory) for directory in directories]):
            name, returncode, seconds = future.result()
            print(f"{'done' if returncode == 0 else f'FAILED ({returncode})'}: {name} ({seconds:.1f} s)", flush=True)
