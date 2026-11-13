"""Execute the public configuration through the CLI and Python entry points.

RUPTURA_CLI_EXECUTABLE=/path/to/ruptura python -m pytest tests/shared_column_test.py
"""
from copy import deepcopy
import io
import json
import os
from pathlib import Path
import subprocess
import zipfile

import numpy as np
import pytest

from ruptura.process import ProcessSimulation
from ruptura.webapp.cluster import cluster_bundle
from ruptura.webapp.config import ISOTHERMS, KINETICS, presets
from ruptura.webapp.cycle_templates import templates

ROOT = Path(__file__).resolve().parents[1]
CLI = Path(os.environ.get("RUPTURA_CLI_EXECUTABLE", ROOT / "build/ruptura"))
pytestmark = pytest.mark.skipif(not CLI.is_file(), reason="Build the native CLI first")


def small():
    d = json.loads((ROOT / "ruptura/webapp/presets/fvm-gas-breakthrough.json").read_text())
    d.update(NumberOfGridPoints=4, SimulationTime=.003, WriteInterval=.001)
    d["Solver"] = {"MaximumTimeStep": .0005, "FirstTimeStep": 1e-5,
                   "RelativeTolerance": 1e-6, "AbsoluteTolerance": 1e-9}
    d["Output"] = {"Directory": "results", "Plots": False, "MainFormatFiles": False,
                   "Checkpoint": False}
    return d


def execute(config, directory):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "simulation.json"
    path.write_text(json.dumps(config, allow_nan=False))
    run = subprocess.run([str(CLI), str(path), "--no-plots"],
                         capture_output=True, text=True, timeout=45)
    assert run.returncode == 0, run.stderr + run.stdout[-2000:]
    with np.load(directory / "results/states.npz", allow_pickle=False) as f:
        arrays = {k: f[k].copy() for k in f.files}
    summary = json.loads((directory / "results/summary.json").read_text())
    assert np.isfinite(arrays["state"]).all()
    return arrays, summary


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
@pytest.mark.parametrize("scheme", ["Upwind", "TVD", "WENO3", "WENO5"])
@pytest.mark.parametrize("solver", ["BDF1", "BDF", "CVODE", "SSPRK3", "SemiImplicitSSPRK3"])
def test_identical_breakthrough_and_cycle_execution(tmp_path, spatial, scheme, solver):
    d = small()
    d.update(SpatialDiscretization=spatial, ConvectionScheme=scheme, BreakthroughIntegrator=solver)
    bt, bts = execute(d, tmp_path / "bt")
    cycle = deepcopy(d)
    cycle["SimulationType"] = "SwingAdsorption"
    cycle["SwingAdsorptionPhases"] = [{"Name": "Adsorption", "Duration": cycle.pop("SimulationTime")}]
    cy, cys = execute(cycle, tmp_path / "cycle")
    np.testing.assert_array_equal(bt["time"], cy["time"])
    np.testing.assert_allclose(bt["state"], cy["state"], atol=1e-11, rtol=1e-11)
    assert abs(bts["final_metrics"]["mass_balance_error"]) < 1e-6
    assert abs(cys["final_metrics"]["mass_balance_error"]) < 1e-6


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
@pytest.mark.parametrize("liquid", [False, True])
@pytest.mark.parametrize("shape,power", [("Sphere", 2), ("CylinderRadial", 1)])
def test_film_and_no_film_hsdm_trajectories_conserve_inventory(tmp_path, spatial, liquid, shape, power):
    d = small()
    d.update(SpatialDiscretization=spatial, SimulationType="SwingAdsorption", NumberOfCycles=1)
    duration = d.pop("SimulationTime")
    d["InitialConditions"] = {"Loading": "zero", **({} if liquid else {"Type": "FeedPressurized"})}
    if liquid:
        d["FluidPhase"] = "Liquid"
        d["MixturePredictionMethod"] = "SPI"
        d["Components"][0]["InitialLiquidPhaseConcentration"] = .5
        d["Components"][0]["PhysisorptionSites"][0]["Parameters"] = [2., .4]
    d["SwingAdsorptionPhases"] = [{"Name": "Closed relaxation", "Duration": duration,
        "Columns": {"A": {"LeftPort": {"Mode": "closed"}, "RightPort": {"Mode": "closed"}}}}]
    d["ParticleTransport"] = {"Model": "HSDM", "Geometry": shape, "Radius": .0005, "RadialCells": 4,
        "SurfaceDiffusivity": {"A": 1e-10, "B": 0.}, "FilmCoefficient": {"A": 1e-5, "B": 0.}}
    means = []
    for film in [False, True]:
        d["ParticleTransport"]["FilmResistance"] = film
        arrays, summary = execute(d, tmp_path / str(film))
        y = arrays["state"]
        count, nz, nc = len(y), 5, 2
        mean = y[:, nc*nz:2*nc*nz].reshape(count, nc, nz)
        shells = y[:, 2*nc*nz:6*nc*nz].reshape(count, 4, nc, nz)
        weights = np.diff(np.linspace(0., 1., 5) ** (power + 1))
        np.testing.assert_allclose(mean, np.einsum("r,trcn->tcn", weights, shells), rtol=1e-7, atol=2e-10)
        assert abs(summary["final_metrics"]["mass_balance_error"]) < 1e-6
        means.append(mean[-1, 0].mean())
    assert means[0] > means[1] > 0., "Finite film resistance must reduce uptake in this matched test"


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
@pytest.mark.parametrize("model", ["BET", "Redlich-Peterson", "Unilan", "OBrien&Myers", "Quadratic", "Temkin", "Bingel&Walton"])
def test_original_main_isotherms_on_both_axial_methods(tmp_path, spatial, model):
    d = small()
    d.update(SpatialDiscretization=spatial, MixturePredictionMethod="SPI")
    d["Components"][0]["PhysisorptionSites"] = [
        {"Type": model, "Parameters": [field[2] for field in ISOTHERMS[model]]}]
    execute(d, tmp_path / "isotherm")


@pytest.mark.parametrize("name", list(templates()))
def test_gui_templates_and_cluster_json_are_the_same_configuration(tmp_path, name):
    d = templates()[name]
    original = deepcopy(d)
    archive = cluster_bundle(d)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        assert json.loads(z.read("simulation.json")) == original
        for script in ["run.sh", "submit_hal.sh"]:
            assert subprocess.run(["bash", "-n"], input=z.read(script), capture_output=True).returncode == 0
        path = tmp_path / "simulation.json"
        path.write_bytes(z.read("simulation.json"))
    run = subprocess.run([str(CLI), str(path), "--validate"], capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr
    assert d == original


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
def test_python_and_cli_execute_identical_configuration(tmp_path, spatial):
    pytest.importorskip("ruptura._ruptura")
    d = small()
    d.update(SpatialDiscretization=spatial, BreakthroughIntegrator="CVODE")
    arrays, _ = execute(d, tmp_path / "cli")
    simulation = ProcessSimulation(d, directory=tmp_path / "python")
    assert simulation.validate()["valid"]
    simulation.run(plots=False)
    np.testing.assert_allclose(simulation.results()["state"], arrays["state"], atol=1e-11, rtol=1e-11)
    assert d["Output"]["Directory"] == "results"


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
@pytest.mark.parametrize("model", ["FirstOrder", "PseudoNth", "Avrami", "General", "Elovich"])
def test_main_chemisorption_kinetics_execute_and_balance(tmp_path, spatial, model):
    d = small()
    d.update(SpatialDiscretization=spatial, MixturePredictionMethod="SPI")
    d["InitialConditions"] = {"Type": "FeedPressurized", "Loading": "zero"}
    parameters = {key: value for key, _, value in KINETICS[model]}
    d["Components"][0]["ChemisorptionSites"] = [{"Type": "Langmuir", "Parameters": [1., 1e-5],
        "RateEquation": model, "RateEquationParameters": parameters, "HeatOfAdsorption": 42000., "MaximumLoading": 1.}]
    _, summary = execute(d, tmp_path / "chemical")
    assert abs(summary["final_metrics"]["mass_balance_error"]) < 1e-6


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
def test_transported_ph_on_both_methods(tmp_path, spatial):
    d = small()
    d.update(SpatialDiscretization=spatial, FluidPhase="Liquid", MixturePredictionMethod="SPI",
             pHMode="HPlus", pHComponent=0, pKw=14.)
    for component in d["Components"]:
        component.pop("CarrierGas", None)
    d["Components"][0]["InitialLiquidPhaseConcentration"] = .5
    d["Components"][0]["PhysisorptionSites"] = [{"Type": "pH-Langmuir", "Parameters": [2., .4, 7.]}]
    d["Components"][1]["PhysisorptionSites"] = [{"Type": "Henry", "Parameters": [0.]}]
    execute(d, tmp_path / "ph")


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
@pytest.mark.parametrize("phase", ["Physisorbed", "BulkConcentration"])
def test_breakthrough_reactions_are_conservative_with_production(tmp_path, spatial, phase):
    d = small()
    d.update(SpatialDiscretization=spatial, MixturePredictionMethod="SPI")
    d["InitialConditions"] = {"Type": "FeedPressurized", "Loading": "equilibrium"}
    d["Reactions"] = [{"Phase": phase, "Style": "GeneralPowerLaw", "Reactants": [0], "Products": [1],
                      "Stoichiometry": {"Reactants": [1.], "Products": [1.]},
                      "Kinetics": {"forwardRateCoefficient": .001, "forwardActivationEnergy": 0.,
                                   "equilibriumConstant": 1., "gibbsFreeEnergy": 0.}}]
    _, summary = execute(d, tmp_path / "reaction")
    assert abs(summary["final_metrics"]["mass_balance_error"]) < 1e-6


@pytest.mark.parametrize("spatial", ["Legacy", "FDM", "FVM"])
def test_cycles_reject_reaction_networks_on_every_entry_path(tmp_path, spatial):
    d = presets()["swing"]
    if spatial != "Legacy":
        d.update(SpatialDiscretization=spatial, BreakthroughIntegrator="BDF")
    d["Reactions"] = [{"Phase": "Physisorbed", "Style": "GeneralPowerLaw", "Reactants": [1], "Products": [0],
        "Stoichiometry": {"Reactants": [1.], "Products": [1.]},
        "Kinetics": {"forwardRateCoefficient": .001, "forwardActivationEnergy": 0.,
                     "equilibriumConstant": 1., "gibbsFreeEnergy": 0.}}]
    path = tmp_path / "simulation.json"
    path.write_text(json.dumps(d))
    run = subprocess.run([str(CLI), str(path), "--validate"], capture_output=True, text=True, timeout=30)
    assert run.returncode != 0
    assert "reaction networks are excluded from cycles" in run.stderr.lower()


@pytest.mark.parametrize("spatial", ["FDM", "FVM"])
def test_cycle_column_overrides_cannot_enable_reaction_networks(tmp_path, spatial):
    d = small()
    d.update(SpatialDiscretization=spatial, SimulationType="SwingAdsorption")
    d["SwingAdsorptionPhases"] = [{"Name": "Adsorption", "Duration": d.pop("SimulationTime")}]
    d["Columns"] = {"A": {"Reactions": [{"Phase": "Physisorbed", "Style": "GeneralPowerLaw",
        "Reactants": [0], "Products": [1], "Stoichiometry": {"Reactants": [1.], "Products": [1.]},
        "Kinetics": {"forwardRateCoefficient": .001, "forwardActivationEnergy": 0.,
                     "equilibriumConstant": 1., "gibbsFreeEnergy": 0.}}]}}
    path = tmp_path / "simulation.json"
    path.write_text(json.dumps(d))
    run = subprocess.run([str(CLI), str(path), "--validate"], capture_output=True, text=True, timeout=30)
    assert run.returncode != 0
    assert "reaction networks are excluded from cycles" in run.stderr.lower()
