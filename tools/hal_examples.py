#!/usr/bin/env python3
"""Catalog, validate, submit and collect the reorganized HAL examples."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import socket
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples/new_examples"
GROUPS = ("all", "tutorials", "mfi_breakthrough", "mfi_cycles", "unified",
          "chemisorption", "shared_chemisorption", "integrated")


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def group_of(relative):
    parts = relative.parts
    if parts[:3] == ("integrated", "campaigns", "co2_ch4_breakthrough"):
        return "mfi_breakthrough"
    if parts[:3] == ("integrated", "campaigns", "co2_ch4_cycles"):
        return "mfi_cycles"
    return {"integrated": "tutorials", "unified": "unified",
            "shared_chemisorption": "chemisorption"}[parts[0]]


def catalog(group="all"):
    rows, configurations = [], {}
    for path in sorted(EXAMPLES.rglob("simulation.json")):
        data = json.loads(path.read_text())
        normalized = json.dumps(data, sort_keys=True, separators=(",", ":"))
        if normalized in configurations:
            raise ValueError("Duplicate configurations: {} and {}".format(
                configurations[normalized], path.relative_to(ROOT)))
        configurations[normalized] = str(path.relative_to(ROOT))
        relative = path.relative_to(EXAMPLES)
        case_group = group_of(relative)
        if group == "shared_chemisorption":
            wanted = case_group == "chemisorption"
        elif group == "integrated":
            wanted = relative.parts[0] == "integrated"
        else:
            wanted = group == "all" or case_group == group
        if not wanted:
            continue
        css = data.get("CyclicSteadyState", {})
        rows.append({
            "case": relative.parent.as_posix(),
            "input": path.relative_to(ROOT).as_posix(),
            "group": case_group,
            "simulation_type": data["SimulationType"],
            "fluid": data.get("FluidPhase", "Gas"),
            "spatial_discretization": data.get("SpatialDiscretization", "Legacy"),
            "convection_scheme": data.get("ConvectionScheme", "legacy default"),
            "integrator": data.get("BreakthroughIntegrator", "legacy default"),
            "css_requested": bool(css.get("Enabled", False)),
            "plots_requested": bool(data.get("Output", {}).get("Plots", False)),
        })
    if not rows:
        raise ValueError("No simulation inputs found for " + group)
    return rows


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def refresh():
    rows = catalog()
    (EXAMPLES / "manifest.txt").write_text(
        "\n".join(row["input"] for row in rows) + "\n")
    write_json(EXAMPLES / "cases.json", rows)
    write_csv(EXAMPLES / "cases.csv", rows)
    print(json.dumps({"total": len(rows),
                      "groups": dict(Counter(row["group"] for row in rows)),
                      "simulation_types": dict(Counter(row["simulation_type"] for row in rows)),
                      "convection_schemes": dict(Counter(row["convection_scheme"] for row in rows))},
                     indent=2))


def resolve_engine(value):
    candidate = Path(value).expanduser()
    if candidate.is_file():
        candidate = candidate.resolve()
    else:
        found = shutil.which(value)
        if found:
            candidate = Path(found).resolve()
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise ValueError("Build RUPTURA first, or set --executable/RUPTURA_EXECUTABLE: " + value)
    return candidate


def validate(rows, engine, jobs):
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    def check(row):
        path = ROOT / row["input"]
        try:
            result = subprocess.run([str(engine), str(path), "--validate"],
                                    cwd=path.parent, env=environment,
                                    capture_output=True, text=True, timeout=60)
            return {"input": row["input"], "exit_code": result.returncode,
                    "error": result.stderr[-4000:]}
        except subprocess.TimeoutExpired:
            return {"input": row["input"], "exit_code": 124,
                    "error": "Input validation exceeded 60 seconds"}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(check, rows))


def submission_command(args, run_dir, engine, count):
    array = "0-{}".format(count - 1)
    if args.max_parallel:
        array += "%" + str(args.max_parallel)
    command = ["sbatch", "--parsable", "--partition=" + args.partition,
               "--nodes=1", "--ntasks=1", "--cpus-per-task=1", "--time=" + args.time,
               "--array=" + array, "--job-name=ruptura-" + args.group,
               "--output=" + str(run_dir / "logs/%A_%a.out"),
               "--error=" + str(run_dir / "logs/%A_%a.err")]
    if args.exclude:
        command.append("--exclude=" + args.exclude)
    command.extend([str(ROOT / "tools/run_new_examples_hal.sh"),
                    str(run_dir / "manifest.json")])
    return command


def submit(args):
    rows = catalog(args.group)
    run_dir = (args.run_dir or ROOT / "hal_runs" /
               (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" +
                args.group + "_" + uuid.uuid4().hex[:8])).expanduser().resolve()
    engine_value = args.executable
    if args.dry_run:
        command = submission_command(args, run_dir, engine_value, len(rows))
        print("Cases: {} | {}".format(len(rows), dict(Counter(r["group"] for r in rows))))
        print("Dry run: no validation, output creation or submission.")
        print(shlex.join(command))
        return 0
    engine = resolve_engine(engine_value)
    if not args.validate_only and not shutil.which("sbatch"):
        raise ValueError("sbatch is unavailable; use --dry-run or --validate-only locally")
    if run_dir.exists() and any(run_dir.iterdir()):
        raise ValueError("Run directory is not empty; choose a new --run-dir: " + str(run_dir))
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    engine_options = list(args.engine_options)
    if engine_options[:1] == ["--"]:
        engine_options.pop(0)
    if args.no_plots:
        engine_options.append("--no-plots")
    for row in rows:
        row["input_sha256"] = digest(ROOT / row["input"])
    manifest = {"schema_version": 1, "created_utc": utc_now(), "project_root": str(ROOT),
                "run_dir": str(run_dir), "selection": args.group, "case_count": len(rows),
                "engine": str(engine), "engine_sha256": digest(engine),
                "engine_size_bytes": engine.stat().st_size,
                "engine_mtime_ns": engine.stat().st_mtime_ns,
                "engine_options": engine_options, "cases": rows}
    write_json(run_dir / "manifest.json", manifest)
    if not args.skip_validation:
        results = validate(rows, engine, args.validation_jobs)
        write_json(run_dir / "validation.json", results)
        failures = [row for row in results if row["exit_code"] != 0]
        print("Input validation: {}/{} passed".format(len(rows)-len(failures), len(rows)), flush=True)
        if failures:
            for failure in failures:
                print(failure["input"] + ": " + failure["error"], file=sys.stderr)
            return 1
    if args.validate_only:
        print("Report: " + str(run_dir / "validation.json"))
        return 0
    command = submission_command(args, run_dir, engine, len(rows))
    result = subprocess.run(command, capture_output=True, text=True)
    write_json(run_dir / "submission.json", {"command": command, "exit_code": result.returncode,
                                           "stdout": result.stdout, "stderr": result.stderr})
    if result.returncode:
        print(result.stderr, file=sys.stderr)
        return result.returncode
    print("Submitted {} cases: {}".format(len(rows), result.stdout.strip()))
    print("Run directory: " + str(run_dir))
    print("Collect with: python3 tools/collect_new_examples_hal.py " + shlex.quote(str(run_dir)))
    return 0


def stage_input(source, destination):
    """Keep physical settings; resolve local assets and isolate execution outputs."""
    data = json.loads(source.read_text())
    assets = []
    def visit(value):
        if isinstance(value, dict):
            for key, item in list(value.items()):
                if isinstance(item, str) and "file" in key.lower():
                    reference = Path(item)
                    source_asset = (source.parent / reference).resolve()
                    if source_asset.is_file():
                        try:
                            relative = source_asset.relative_to(source.parent.resolve())
                        except ValueError:
                            value[key] = str(source_asset)
                        else:
                            target = destination / relative
                            target.parent.mkdir(parents=True, exist_ok=True)
                            if source_asset != target:
                                shutil.copy2(source_asset, target)
                            value[key] = relative.as_posix()
                        assets.append({"source": str(source_asset), "sha256": digest(source_asset)})
                else:
                    visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)
    visit(data)
    if "Output" in data:
        data["Output"]["Directory"] = "results"
    staged = destination / "simulation.json"
    write_json(staged, data)
    return staged, assets


def run_case(args):
    manifest_path = args.manifest.expanduser().resolve()
    manifest = json.loads(manifest_path.read_text())
    if not 0 <= args.index < len(manifest["cases"]):
        raise ValueError("Array index is outside this manifest")
    row = manifest["cases"][args.index]
    project = Path(manifest["project_root"])
    source = project / row["input"]
    directory = Path(manifest["run_dir"]) / "cases" / row["case"]
    directory.mkdir(parents=True, exist_ok=True)
    receipt = directory / "execution.json"
    if receipt.exists():
        raise ValueError("This task already has an execution record; use a new run directory")
    started = time.monotonic()
    record = dict(row, status="running", started_utc=utc_now(), hostname=socket.gethostname(),
                  slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                  slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"),
                  slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                  engine_sha256=manifest["engine_sha256"])
    write_json(receipt, record)
    process = None
    caught_signal = []
    def forward(signum, frame):
        caught_signal.append(signum)
        if process is not None and process.poll() is None:
            process.send_signal(signum)
    old_handlers = {s: signal.signal(s, forward) for s in (signal.SIGTERM, signal.SIGINT)}
    try:
        if digest(source) != row["input_sha256"]:
            raise ValueError("Source input changed after submission: " + row["input"])
        engine = resolve_engine(manifest["engine"])
        if (engine.stat().st_size != manifest["engine_size_bytes"] or
                engine.stat().st_mtime_ns != manifest["engine_mtime_ns"]):
            raise ValueError("The executable changed after submission; use the recorded build")
        staged, assets = stage_input(source, directory)
        record["effective_input_sha256"] = digest(staged)
        record["assets"] = assets
        command = [str(engine), str(staged)] + manifest["engine_options"]
        if not args.local:
            command = ["srun", "--ntasks=1", "--cpus-per-task=1"] + command
        record["command"] = command
        write_json(receipt, record)
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                           MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
        environment.setdefault("RUPTURA_PLOT_BACKEND", "matplotlib")
        print("Case {}: {}".format(args.index, row["input"]), flush=True)
        process = subprocess.Popen(command, cwd=directory, env=environment)
        code = process.wait()
        if code < 0:
            code = 128 - code
        record["status"] = ("interrupted" if caught_signal else
                            "completed" if code == 0 else "failed")
        record["exit_code"] = code
    except Exception as error:
        record.update(status="failed", exit_code=1, error=str(error))
        print(str(error), file=sys.stderr)
    finally:
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)
        record.update(finished_utc=utc_now(), elapsed_seconds=round(time.monotonic()-started, 6))
        write_json(receipt, record)
    return record["exit_code"]


def collect(run_dir):
    run_dir = run_dir.expanduser().resolve()
    manifest = json.loads((run_dir / "manifest.json").read_text())
    records = []
    for case in manifest["cases"]:
        row = dict(case)
        directory = run_dir / "cases" / case["case"]
        receipt = directory / "execution.json"
        summary_path = directory / "results/summary.json"
        row["status"] = "not_started"
        if receipt.exists():
            execution = json.loads(receipt.read_text())
            row.update({k: execution.get(k) for k in
                        ("exit_code", "elapsed_seconds", "hostname", "slurm_job_id")})
            row["status"] = execution.get("status", "unknown")
            if row["status"] == "running":
                row["status"] = "running_or_interrupted"
            if execution.get("error"):
                row["error"] = execution["error"]
            if row["status"] == "completed":
                if summary_path.exists():
                    summary = json.loads(summary_path.read_text())
                    row["cycles"] = summary.get("cycles")
                    row["converged"] = summary.get("converged")
                    if case["css_requested"]:
                        row["status"] = ("css_reached" if summary.get("converged") is True
                                         else "completed_not_css" if summary.get("converged") is False
                                         else "completed_css_unknown")
                    else:
                        row["status"] = "completed_fixed_cycles" if case["simulation_type"] == "SwingAdsorption" else "completed"
                    for key, value in summary.get("final_metrics", {}).items():
                        row["metric_" + key] = (json.dumps(value, separators=(",", ":"))
                                                if isinstance(value, (dict, list)) else value)
                else:
                    row["status"] = ("completed_legacy" if case["spatial_discretization"] == "Legacy"
                                     else "completed_missing_summary")
        plot_directory = directory / "results/plots"
        plot_manifest = plot_directory / "plot_manifest.json"
        files = []
        if plot_manifest.is_file():
            try:
                names = json.loads(plot_manifest.read_text()).get("files", [])
                files = [plot_directory / Path(name).name for name in names
                         if (plot_directory / Path(name).name).is_file()]
            except (ValueError, OSError):
                row["plot_error"] = "The plot manifest could not be read"
        row["plot_png_count"] = sum(path.suffix.lower() == ".png" for path in files)
        row["plot_pdf_count"] = sum(path.suffix.lower() == ".pdf" for path in files)
        row["plot_svg_count"] = sum(path.suffix.lower() == ".svg" for path in files)
        disabled = "--no-plots" in manifest.get("engine_options", [])
        row["plot_status"] = ("png_pdf" if row["plot_png_count"] and row["plot_pdf_count"] else
                              "svg" if row["plot_svg_count"] else
                              "disabled" if disabled or case.get("plots_requested") is False else "missing")
        records.append(row)
    counts = dict(Counter(row["status"] for row in records))
    write_json(run_dir / "results_summary.json",
               {"selection": manifest["selection"], "total": len(records),
                "status_counts": counts, "cases": records})
    write_csv(run_dir / "results_summary.csv", records)
    print(json.dumps({"total": len(records), "status_counts": counts,
                      "json": str(run_dir / "results_summary.json"),
                      "csv": str(run_dir / "results_summary.csv")}, indent=2))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("refresh", help="Rebuild the portable example catalogs")
    listing = sub.add_parser("list")
    listing.add_argument("group", choices=GROUPS, nargs="?", default="all")
    submission = sub.add_parser("submit")
    submission.add_argument("group", choices=GROUPS, nargs="?", default="all")
    submission.add_argument("--executable", default=os.environ.get("RUPTURA_EXECUTABLE", str(ROOT / "build/ruptura")))
    submission.add_argument("--run-dir", type=Path)
    submission.add_argument("--dry-run", action="store_true")
    submission.add_argument("--validate-only", action="store_true")
    submission.add_argument("--skip-validation", action="store_true")
    submission.add_argument("--validation-jobs", type=int, default=4)
    submission.add_argument("--partition", default=os.environ.get("RUPTURA_HAL_PARTITION", "serial"))
    submission.add_argument("--time", default=os.environ.get("RUPTURA_HAL_TIME", "168:00:00"))
    submission.add_argument("--exclude", default=os.environ.get("RUPTURA_HAL_EXCLUDE", "c171"))
    submission.add_argument("--max-parallel", type=int)
    submission.add_argument("--no-plots", action="store_true")
    running = sub.add_parser("run")
    running.add_argument("manifest", type=Path)
    running.add_argument("--index", type=int, required=True)
    running.add_argument("--local", action="store_true", help="Execute without srun for local checks")
    collecting = sub.add_parser("collect")
    collecting.add_argument("run_dir", type=Path)
    arguments = list(sys.argv[1:] if argv is None else argv)
    engine_options = []
    if "--" in arguments:
        separator = arguments.index("--")
        engine_options = arguments[separator + 1:]
        arguments = arguments[:separator]
    args = parser.parse_args(arguments)
    if args.command == "submit":
        args.engine_options = engine_options
        if args.validation_jobs < 1 or (args.max_parallel is not None and args.max_parallel < 1):
            parser.error("Job counts must be positive")
        if args.validate_only and args.skip_validation:
            parser.error("--validate-only cannot be combined with --skip-validation")
    elif engine_options:
        parser.error("Engine options after -- are supported only for submit")
    try:
        if args.command == "refresh":
            refresh()
            return 0
        if args.command == "list":
            rows = catalog(args.group)
            print("\n".join(row["input"] for row in rows))
            print("Total: " + str(len(rows)), file=sys.stderr)
            return 0
        if args.command == "submit":
            return submit(args)
        if args.command == "run":
            return run_case(args)
        return collect(args.run_dir)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as error:
        print("Error: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

