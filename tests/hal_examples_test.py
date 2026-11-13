"""Exercise HAL submission, staging, failure reporting and CSS classification."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_tools():
    spec = importlib.util.spec_from_file_location("hal_examples_under_test", ROOT / "tools/hal_examples.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_example_catalog_has_unique_weno3_inputs_and_matched_pairs():
    module = load_tools()
    rows = module.catalog()
    assert len(rows) == 296
    assert {row["convection_scheme"] for row in rows} == {"WENO3"}
    assert all(row["plots_requested"] for row in rows)
    assert sum(row["simulation_type"] == "Breakthrough" for row in rows) == 92
    for group in ("unified", "shared_chemisorption"):
        for fdm in (module.EXAMPLES / group / "fdm").rglob("simulation.json"):
            fvm = module.EXAMPLES / group / "fvm" / fdm.relative_to(module.EXAMPLES / group / "fdm")
            a, b = json.loads(fdm.read_text()), json.loads(fvm.read_text())
            assert a.pop("SpatialDiscretization") == "FDM"
            assert b.pop("SpatialDiscretization") == "FVM"
            assert a == b


def test_dry_run_recognizes_options_after_group(tmp_path):
    result = subprocess.run([str(ROOT / "tools/submit_new_examples_hal.sh"), "all",
                             "--dry-run", "--run-dir", str(tmp_path / "not_created")],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--array=0-295" in result.stdout
    assert "%" not in result.stdout.split("--array=", 1)[1].split()[0]
    assert not (tmp_path / "not_created").exists()


def test_generators_preserve_layout_and_weno3_without_recreating_duplicates(tmp_path):
    for relative in ["marx2/ternary/T15/simulation.json"] + [
            "ruptura/webapp/presets/fvm-{}-{}.json".format(fluid, mode)
            for fluid in ("gas", "liquid") for mode in ("breakthrough", "cycle")]:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, target)
    for name in ("unified", "integrated", "flow_control", "chemisorption"):
        path = ROOT / ("tools/generate_{}_examples.py".format(name))
        spec = importlib.util.spec_from_file_location("generator_" + name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.ROOT = tmp_path
        if name == "integrated":
            module.OUT = tmp_path / "examples/new_examples/integrated"
        if name == "flow_control":
            module.TARGET = tmp_path / "examples/new_examples/integrated/flow_controls"
        module.main()
    inputs = list((tmp_path / "examples/new_examples").rglob("simulation.json"))
    assert len(inputs) == 86  # Generated cases; migrated donor cases are separate.
    assert all(json.loads(p.read_text())["ConvectionScheme"] == "WENO3" for p in inputs)
    assert all(json.loads(p.read_text())["Output"]["Plots"] is True for p in inputs)
    assert not (tmp_path / "examples/new_examples/integrated/breakthrough/gas_fixed_bed").exists()
    for old in ("integrated", "unified", "shared_chemisorption"):
        assert not (tmp_path / "examples" / old).exists()


@pytest.fixture
def project(tmp_path):
    root = tmp_path / "project with spaces"
    (root / "tools").mkdir(parents=True)
    for name in ("hal_examples.py", "submit_new_examples_hal.sh",
                 "run_new_examples_hal.sh", "collect_new_examples_hal.py"):
        shutil.copy2(ROOT / "tools" / name, root / "tools" / name)
    binary = root / "engine with spaces"
    binary.write_text("#!" + sys.executable + "\n" + """
import json, pathlib, sys
if '--validate' in sys.argv:
    raise SystemExit(0)
path = pathlib.Path(sys.argv[1])
data = json.loads(path.read_text())
asset = path.parent / data['MPDSettings']['FileName']
assert asset.read_text() == 'asset data'
if data.get('TestFailure'):
    raise SystemExit(7)
out = path.parent / 'results'
out.mkdir()
(out/'summary.json').write_text(json.dumps({
    'cycles': 2, 'converged': data.get('TestConverged', False),
    'final_metrics': {'mass_balance_error': 1e-10, 'energy_balance_error': 2e-9}}))
(path.parent/'arguments.json').write_text(json.dumps(sys.argv[2:]))
""")
    binary.chmod(0o755)
    bin_dir = tmp_path / "mock bin"
    bin_dir.mkdir()
    sbatch = bin_dir / "sbatch"
    sbatch.write_text("#!" + sys.executable + "\n" + """
import json, os, pathlib, sys
pathlib.Path(os.environ['SBATCH_CAPTURE']).write_text(json.dumps(sys.argv[1:]))
print('12345')
""")
    sbatch.chmod(0o755)
    srun = bin_dir / "srun"
    srun.write_text("#!" + sys.executable + "\n" + """
import os, sys
args = [value for value in sys.argv[1:] if not value.startswith('--ntasks=') and not value.startswith('--cpus-per-task=')]
os.execv(args[0], args)
""")
    srun.chmod(0o755)
    case = root / "examples/new_examples/unified/fvm/example"
    case.mkdir(parents=True)
    (case / "distribution.data").write_text("asset data")
    data = {"SimulationType": "Breakthrough", "SpatialDiscretization": "FVM",
            "ConvectionScheme": "WENO3", "MPDSettings": {"FileName": "distribution.data"},
            "Output": {"Directory": "results", "Plots": True}}
    (case / "simulation.json").write_text(json.dumps(data))
    environment = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ["PATH"],
                       SBATCH_CAPTURE=str(tmp_path / "sbatch_arguments.json"),
                       SLURM_ARRAY_TASK_ID="0", SLURM_JOB_ID="12345",
                       SLURM_ARRAY_JOB_ID="12345")
    return root, case, data, binary, environment


def submit(project, *, no_plots=True):
    root, case, data, engine, env = project
    run_dir = root / "HAL output"
    arguments = [str(root / "tools/submit_new_examples_hal.sh"), "all",
                 "--executable", str(engine), "--run-dir", str(run_dir)]
    if no_plots:
        arguments.append("--no-plots")
    result = subprocess.run(arguments, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return run_dir


def test_submission_uses_frozen_array_and_keeps_plot_choice(project):
    run_dir = submit(project, no_plots=False)
    root, case, data, engine, env = project
    command = json.loads(Path(env["SBATCH_CAPTURE"]).read_text())
    assert "--array=0-0" in command
    assert command[-2] == str(root / "tools/run_new_examples_hal.sh")
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["engine_options"] == []
    assert manifest["case_count"] == 1
    assert manifest["engine"] == str(engine)
    assert json.loads((run_dir / "validation.json").read_text())[0]["exit_code"] == 0


@pytest.mark.parametrize("kind,css,converged,status", [
    ("Breakthrough", False, False, "completed"),
    ("SwingAdsorption", False, False, "completed_fixed_cycles"),
    ("SwingAdsorption", True, False, "completed_not_css"),
    ("SwingAdsorption", True, True, "css_reached"),
])
def test_worker_stages_assets_and_collector_distinguishes_css(project, kind, css, converged, status):
    root, case, data, engine, env = project
    data.update(SimulationType=kind, CyclicSteadyState={"Enabled": css}, TestConverged=converged)
    (case / "simulation.json").write_text(json.dumps(data))
    original = (case / "simulation.json").read_bytes()
    run_dir = submit(project)
    # SLURM executes a copy in a spool directory, away from the project tree.
    spool = root.parent / "slurm spool"
    spool.mkdir()
    copied_worker = spool / "slurm_script"
    shutil.copy2(root / "tools/run_new_examples_hal.sh", copied_worker)
    run = subprocess.run([str(copied_worker), str(run_dir / "manifest.json")],
                         cwd=spool, env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert (case / "simulation.json").read_bytes() == original
    assert not (case / "results").exists()
    staged = run_dir / "cases/unified/fvm/example"
    assert (staged / "distribution.data").read_text() == "asset data"
    assert json.loads((staged / "arguments.json").read_text()) == ["--no-plots"]
    result = subprocess.run([sys.executable, str(root / "tools/collect_new_examples_hal.py"), str(run_dir)],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    report = json.loads((run_dir / "results_summary.json").read_text())
    assert report["cases"][0]["status"] == status
    assert report["cases"][0]["metric_mass_balance_error"] == 1e-10
    assert (run_dir / "results_summary.csv").exists()


def test_worker_records_engine_failure(project):
    root, case, data, engine, env = project
    data["TestFailure"] = True
    (case / "simulation.json").write_text(json.dumps(data))
    run_dir = submit(project)
    run = subprocess.run([str(root / "tools/run_new_examples_hal.sh"), str(run_dir / "manifest.json")],
                         env=env, capture_output=True, text=True)
    assert run.returncode == 7
    receipt = json.loads((run_dir / "cases/unified/fvm/example/execution.json").read_text())
    assert receipt["status"] == "failed" and receipt["exit_code"] == 7


def test_worker_rejects_input_changed_after_submission(project):
    root, case, data, engine, env = project
    run_dir = submit(project)
    data["ConvectionScheme"] = "TVD"
    (case / "simulation.json").write_text(json.dumps(data))
    run = subprocess.run([str(root / "tools/run_new_examples_hal.sh"), str(run_dir / "manifest.json")],
                         env=env, capture_output=True, text=True)
    assert run.returncode != 0
    assert "changed after submission" in run.stderr
    assert not (run_dir / "cases/unified/fvm/example/results").exists()


def test_collector_recognizes_replotted_files_after_optout(project):
    root, case, data, engine, env = project
    run_dir = submit(project, no_plots=True)
    result = subprocess.run([str(root / "tools/run_new_examples_hal.sh"), str(run_dir / "manifest.json")],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    plots = run_dir / "cases/unified/fvm/example/results/plots"
    plots.mkdir()
    # File-accounting fixtures emulate a later, separate postprocessing step.
    (plots / "figure.png").write_bytes(b"PNG accounting fixture")
    (plots / "figure.pdf").write_bytes(b"PDF accounting fixture")
    (plots / "plot_manifest.json").write_text(json.dumps({"files": ["figure.png", "figure.pdf"]}))
    result = subprocess.run([sys.executable, str(root / "tools/collect_new_examples_hal.py"), str(run_dir)],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    row = json.loads((run_dir / "results_summary.json").read_text())["cases"][0]
    assert row["plot_status"] == "png_pdf"
    assert row["plot_png_count"] == row["plot_pdf_count"] == 1

