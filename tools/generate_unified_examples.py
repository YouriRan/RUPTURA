#!/usr/bin/env python3
"""Generate matched FDM/FVM examples using the same public simulation schema."""
from copy import deepcopy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ruptura.webapp.cycle_templates import templates


MPD_DATA = "# Synthetic single-site probabilities; A then B, last index fastest.\n" \
    "0.16666666666666666\n0.33333333333333331\n0.5\n0\n"


def mpd_breakthrough():
    """Synthetic competitive Langmuir limit, not a fitted material model."""
    return {
        "SimulationType": "Breakthrough", "FluidPhase": "Gas",
        "MixturePredictionMethod": "MPD", "Temperature": 300., "ColumnLength": .02,
        "NumberOfGridPoints": 10, "ParticleDensity": 1000., "InletPressure": 2.5e5,
        "OutletPressure": 2.5e5, "ColumnEntranceVelocity": .05, "BoundaryCondition": "FixedVelocity",
        "Geometry": {"Type": "PackedBed", "InternalDiameter": .01, "OuterDiameter": .012,
                     "ColumnVoidFraction": .4, "ParticleDiameter": .001},
        "Components": [
            {"Name": "A", "GasPhaseMolFraction": .2, "MolecularWeight": .044,
             "GasHeatCapacity": 37., "MassTransferCoefficient": .2,
             "AxialDispersionCoefficient": 1e-6, "HeatOfAdsorption": 20000.},
            {"Name": "He", "GasPhaseMolFraction": .2, "MolecularWeight": .004,
             "GasHeatCapacity": 20.8, "CarrierGas": True},
            {"Name": "B", "GasPhaseMolFraction": .6, "MolecularWeight": .028,
             "GasHeatCapacity": 29., "MassTransferCoefficient": .1,
             "AxialDispersionCoefficient": 1e-6, "HeatOfAdsorption": 18000.}],
        "MPDSettings": {"FileName": "distribution.data", "ReferenceTemperature": 300.,
                        "ReferenceFugacity": 1e5, "ReferenceFrameworkMass": 1./6.02214076e23,
                        "ComponentBounds": [{"Component": name, "NMin": 0, "NMax": 1, "DeltaN": 1}
                                            for name in ["A", "B"]]},
        "InitialConditions": {"Type": "CarrierGas", "Loading": "zero"},
        "BreakthroughIntegrator": "CVODE", "ConvectionScheme": "WENO3",
        "SimulationTime": 60., "WriteInterval": .5,
        "Solver": {"MaximumTimeStep": .1, "FirstTimeStep": 1e-5,
                   "RelativeTolerance": 1e-7, "AbsoluteTolerance": 1e-10},
        "Output": {"Directory": "results", "Plots": True, "Checkpoint": True,
                   "MainFormatFiles": True}}


def main():
    destination = ROOT / "examples/new_examples/unified"
    cases = {}
    for fluid in ["gas", "liquid"]:
        for mode in ["breakthrough", "cycle"]:
            cases[f"{fluid}_{mode}"] = json.loads(
                (ROOT / f"ruptura/webapp/presets/fvm-{fluid}-{mode}.json").read_text())
        for shape in ["Sphere", "CylinderRadial"]:
            for film in [False, True]:
                d = deepcopy(cases[f"{fluid}_breakthrough"])
                d["ParticleTransport"] = {"Model": "HSDM", "Geometry": shape, "Radius": .0005,
                    "RadialCells": 12, "SurfaceDiffusivity": {"A": 1e-10, "B": 0.},
                    "FilmResistance": film, "FilmCoefficient": {"A": 1e-5, "B": 0.}}
                d["CVODELinearSolver"] = "SPGMR"
                cases[f"{fluid}_hsdm_{shape.lower()}_{'film' if film else 'no_film'}"] = d
    for name, document in templates().items():
        cases[name] = document
    cases["gas_mpd_breakthrough"] = mpd_breakthrough()
    paths = []
    for spatial in ["FDM", "FVM"]:
        for name, original in cases.items():
            d = deepcopy(original)
            d["SpatialDiscretization"] = spatial
            d["ConvectionScheme"] = "WENO3"
            d.setdefault("Output", {})["Plots"] = True
            path = destination / spatial.lower() / name / "simulation.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(d, indent=2, allow_nan=False) + "\n")
            if name == "gas_mpd_breakthrough":
                (path.parent / "distribution.data").write_text(MPD_DATA)
            paths.append(path.relative_to(ROOT).as_posix())
    (destination / "manifest.txt").write_text("\n".join(paths) + "\n")
    print(f"Generated {len(paths)} matched examples in {destination.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
