"""Verify headless automatic figures and replotting on the legacy CLI path."""
import hashlib
import json
import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = Path(os.environ.get("RUPTURA_CLI_EXECUTABLE", ROOT / "build/ruptura"))
pytestmark = pytest.mark.skipif(not CLI.is_file(), reason="Build the native CLI first")


def legacy_input(directory):
    path = ROOT / "examples/new_examples/integrated/validation/marx_T15_original/simulation.json"
    document = json.loads(path.read_text())
    document.update(NumberOfGridPoints=6, SimulationTime=.02, WriteInterval=.01)
    document["Output"] = {"Plots": True, "Directory": "results"}
    target = directory / "simulation.json"
    target.write_text(json.dumps(document))
    return target


def run(path, *options, env=None):
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", RUPTURA_PLOT_BACKEND="matplotlib")
    environment.update(env or {})
    return subprocess.run([str(CLI), str(path), *options], cwd=path.parent,
                          env=environment, capture_output=True, text=True, timeout=60)


def figures(directory):
    plots = directory / "results/plots"
    report = json.loads((plots / "plot_manifest.json").read_text())
    png = {Path(name).stem for name in report["files"] if name.endswith(".png")}
    pdf = {Path(name).stem for name in report["files"] if name.endswith(".pdf")}
    assert len(png) >= 5 and png == pdf
    assert all((plots / name).stat().st_size > 100 for name in report["files"])
    assert report["legacy"] is True
    assert report["final_time_min"] == pytest.approx(.02 / 60)
    return report


@pytest.mark.parametrize("skip_initial_plots", [False, True])
def test_legacy_automatic_figures_and_plot_only_preserve_raw_results(tmp_path, skip_initial_plots):
    path = legacy_input(tmp_path)
    result = run(path, *(["--no-plots"] if skip_initial_plots else []))
    assert result.returncode == 0, result.stderr + result.stdout[-1000:]
    assert (tmp_path / "results/plots/plot_manifest.json").exists() is not skip_initial_plots
    if not skip_initial_plots:
        figures(tmp_path)
    raw = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.glob("*.data")}
    assert "column.data" in raw and any(name.startswith("component_") for name in raw)
    result = run(path, "--plot-only")
    assert result.returncode == 0, result.stderr
    figures(tmp_path)
    assert raw == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.glob("*.data")}


def test_legacy_plotting_failure_is_reported(tmp_path):
    path = legacy_input(tmp_path)
    result = run(path, env={"RUPTURA_PYTHON": "/bin/false"})
    assert result.returncode != 0
    assert "Legacy plotting failed" in result.stderr
    assert (tmp_path / "results/plots/plot_log.txt").exists()

