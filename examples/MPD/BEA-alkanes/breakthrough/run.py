#!/usr/bin/env python3

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import sys
import glob

def run_ruptura(workdir):
    print(f"Running in {workdir}")

    log_path = f"{workdir}/ruptura.log"
    with open(log_path, "w") as log:
        result = subprocess.run(
            ["../../../../../build/ruptura"],
            cwd=workdir,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    return {
        "dir": str(workdir),
        "returncode": result.returncode,
        "log": str(log_path),
    }


dirs = glob.glob("*_r1")

with ThreadPoolExecutor(max_workers=8) as executor:
    futures = {executor.submit(run_ruptura, d): d for d in dirs}

    for future in as_completed(futures):
        result = future.result()

        print(f"\n=== {result['dir']} ===")
        print(f"exit code: {result['returncode']}")
        print(f"log: {result['log']}")
