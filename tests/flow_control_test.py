import copy
import json
from pathlib import Path

import pytest
import ruptura
from ruptura.process import ProcessSimulation
from ruptura.webapp.config import SCHEMA, presets

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((ROOT / "examples/new_examples/integrated/flow_controls").glob("*/simulation.json"))


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.parent.name)
def test_flow_examples_validate_through_python(path, tmp_path):
    data = json.loads(path.read_text())
    data["Output"]["Directory"] = str(tmp_path / "results")
    before = copy.deepcopy(data)
    assert ProcessSimulation(data, directory=tmp_path).validate()["valid"]
    assert data == before


def test_original_main_sphericity_and_removed_boundary(tmp_path):
    data = presets()["breakthrough"]
    data["Geometry"]["ParticleSphericity"] = 0.7
    path = tmp_path / "simulation.json"
    path.write_text(json.dumps(data))
    assert ruptura.InputReader(str(path)).particle_sphericity == 0.7
    assert not hasattr(ruptura.Column.BoundaryCondition, "INLET_PRESSURE_" + "INLET_VELOCITY")
    data["BoundaryCondition"] = "InletPressure" + "InletVelocity"
    path.write_text(json.dumps(data))
    with pytest.raises((ValueError, RuntimeError)):
        ruptura.InputReader(str(path))


def test_lab_boundary_menu_and_thermal_fields():
    condition = next(f for f in SCHEMA["column"] if f["key"] == "BoundaryCondition")
    assert set(condition["options"]) == {
        "InletVelocityOutletPressure", "InletPressureOutletPressure", "FixedVelocity", "FixedPressureGradient"
    }
    assert {"ThermalModel", "FluidSolidThermalEquilibrium"} <= {f["key"] for f in SCHEMA["thermal"]}


def test_phase_neutral_thermal_alias_conflict(tmp_path):
    data = json.loads((ROOT / "examples/new_examples/integrated/flow_controls/liquid_common_temperature_coolant/simulation.json").read_text())
    data["Output"]["Directory"] = str(tmp_path / "results")
    data["GasSolidThermalEquilibrium"] = False
    with pytest.raises((ValueError, RuntimeError), match="conflicts"):
        ProcessSimulation(data, directory=tmp_path).validate()
