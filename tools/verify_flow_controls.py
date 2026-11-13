#!/usr/bin/env python3
"""Validate inputs and execute the short flow-control integration examples."""
import argparse
import json
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, default=ROOT / "build/ruptura")
    parser.add_argument("--output", type=Path, default=ROOT / "verification/flow_update")
    parser.add_argument("--validate-all", action="store_true")
    args = parser.parse_args()
    binary, output = args.executable.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    examples = sorted((ROOT / "examples/new_examples/integrated/flow_controls").glob("*/simulation.json"))
    validation = sorted((ROOT / "examples/new_examples/integrated").rglob("simulation.json")) if args.validate_all else examples
    checks, runs = [], []
    for path in validation:
        start = time.monotonic()
        result = subprocess.run([str(binary), str(path), "--validate"], capture_output=True, text=True, timeout=60)
        checks.append({"input": str(path.relative_to(ROOT)), "returncode": result.returncode,
                       "seconds": time.monotonic() - start})
        if result.returncode:
            checks[-1]["error"] = (result.stderr + result.stdout)[-3000:]
        (output / "input_validation.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(f"Input validation: {sum(c['returncode']==0 for c in checks)}/{len(checks)}", flush=True)
    for path in examples:
        folder = output / "runs" / path.parent.name
        folder.mkdir(parents=True, exist_ok=True)
        data = json.loads(path.read_text())
        data["Output"]["Directory"] = str(folder / "results")
        data["Output"]["Plots"] = False
        case = folder / "simulation.json"
        case.write_text(json.dumps(data, indent=2) + "\n")
        start = time.monotonic()
        result = subprocess.run([str(binary), str(case), "--no-plots"], capture_output=True, text=True, timeout=60)
        (folder / "run.log").write_text(result.stdout + result.stderr)
        record = {"example": path.parent.name, "returncode": result.returncode,
                  "seconds": time.monotonic() - start}
        summary = folder / "results/summary.json"
        if summary.exists():
            report = json.loads(summary.read_text())
            record["summary"] = report
            (output / (path.parent.name + "_summary.json")).write_text(json.dumps(report, indent=2) + "\n")
        if result.returncode:
            record["error"] = (result.stderr + result.stdout)[-3000:]
        runs.append(record)
        (output / "complete_runs.json").write_text(json.dumps(runs, indent=2) + "\n")
        print(f"{path.parent.name}: {'passed' if result.returncode == 0 else 'FAILED'}", flush=True)
    return 0 if all(c["returncode"] == 0 for c in checks + runs) else 1

if __name__ == "__main__":
    raise SystemExit(main())
