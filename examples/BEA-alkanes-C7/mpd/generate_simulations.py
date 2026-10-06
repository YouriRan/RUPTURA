#!/usr/bin/env python3
"""Write the four Ruptura inputs of this example from ``fitted_parameters.json``.

Keeping the inputs generated means the MPD and the Langmuir/IAST runs can never
drift apart: both are emitted from the same fitted saturation capacities and
affinities, and the MPD framework mass is the one the particle distribution was
built with.

Column, feed and integrator settings follow
``examples/BEA-alkanes-C7/breakthrough``.  Two of them are reduced for cost:
``NumberOfGridPoints`` (100 -> 20) and the output ``TimeStep`` (0.001 -> 0.105 s),
matching the MPD benchmark in ``examples/MPD/BEA-alkanes/breakthrough``.  An MPD
equilibrium evaluation walks every populated macrostate, so its cost scales with
(grid points) x (integration steps).
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARAMETERS = json.loads((HERE / "fitted_parameters.json").read_text(encoding="utf-8"))

COMPONENTS = list(PARAMETERS["components"])
TEMPERATURE = PARAMETERS["referenceTemperature"]
FRAMEWORK_MASS = PARAMETERS["frameworkMass"]
TOTAL_SITE_COUNT = PARAMETERS["totalSiteCount"]
DELTA_N = 2
REFERENCE_FUGACITY = 1.0e5
MOLECULAR_WEIGHT = 0.100198  # kg/mol, C7H16

PRESSURE_START = 1.0e3
PRESSURE_END = 1.0e7
NUMBER_OF_PRESSURE_POINTS = 64

MPD_SETTINGS = {
    "ReferenceTemperature": TEMPERATURE,
    "ReferenceFugacity": REFERENCE_FUGACITY,
    "ReferenceFrameworkMass": FRAMEWORK_MASS,
    "ComponentBounds": [
        {"Component": component, "NMin": 0, "NMax": TOTAL_SITE_COUNT, "DeltaN": DELTA_N}
        for component in COMPONENTS
    ],
}


def langmuir_sites(component: str) -> list[dict]:
    return [
        {"Type": "Langmuir", "Parameters": [family["saturationCapacity"], family["affinities"][component]]}
        for family in PARAMETERS["siteFamilies"]
    ]


def write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(HERE)}")


def mixture_prediction(method: str, display_name: str, *, mpd: bool) -> dict:
    payload = {
        "SimulationType": "MixturePrediction",
        "MixturePredictionMethod": method,
        "DisplayName": display_name,
        "Temperature": TEMPERATURE,
        "PressureStart": PRESSURE_START,
        "PressureEnd": PRESSURE_END,
        "NumberOfPressurePoints": NUMBER_OF_PRESSURE_POINTS,
        "PressureScale": "Log",
        "Components": [
            {"Name": component, "GasPhaseMolFraction": 0.5, "MolecularWeight": MOLECULAR_WEIGHT}
            for component in COMPONENTS
        ],
    }
    if mpd:
        payload["MPDSettings"] = {"FileName": "../../particle_distribution.data", **MPD_SETTINGS}
    else:
        for entry, component in zip(payload["Components"], COMPONENTS):
            entry["PhysisorptionSites"] = langmuir_sites(component)
            entry["HeatOfAdsorption"] = PARAMETERS["heatOfAdsorption"][component]
            entry["referenceTemperature"] = TEMPERATURE
    return payload


def breakthrough(method: str, display_name: str, *, mpd: bool) -> dict:
    carrier = {
        "Name": "Helium",
        "GasPhaseMolFraction": 0.98,
        "MolecularWeight": 0.0040026,
        "CarrierGas": True,
    }
    adsorbates = []
    for component in COMPONENTS:
        entry = {
            "Name": component,
            "GasPhaseMolFraction": 0.01,
            "MassTransferCoefficient": 0.06,
            "AxialDispersionCoefficient": 0.0,
            "MolecularWeight": MOLECULAR_WEIGHT,
        }
        if not mpd:
            entry["PhysisorptionSites"] = langmuir_sites(component)
            entry["HeatOfAdsorption"] = PARAMETERS["heatOfAdsorption"][component]
            entry["referenceTemperature"] = TEMPERATURE
        adsorbates.append(entry)

    payload = {
        "SimulationType": "Breakthrough",
        "MixturePredictionMethod": method,
        "DisplayName": display_name,
        "Temperature": TEMPERATURE,
        "ColumnVoidFraction": 0.4,
        "ParticleDensity": 1508.52,
        "PressureGradient": -1000.0,
        "ColumnEntranceVelocity": 0.019,
        "ColumnLength": 0.1,
        "InletPressure": 120000.0,
        "NumberOfTimeSteps": "auto",
        "BreakthroughIntegrator": "CVODE",
        "BoundaryCondition": "FixedPressureInletVelocity",
        "PrintEvery": 5000,
        "WriteEvery": 100,
        "DynamicViscosity": 3e-05,
        "ParticleDiameter": 0.001,
        "TimeStep": 0.10526315789473686,
        "NumberOfGridPoints": 20,
        "Components": [carrier, *adsorbates],
        "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": 0.4, "ParticleDiameter": 0.001},
    }
    if mpd:
        payload["MPDSettings"] = {"FileName": "../../particle_distribution.data", **MPD_SETTINGS}
    return payload


def main() -> None:
    label = " / ".join(COMPONENTS)
    write(HERE / "mixture" / "mpd" / "simulation.json",
          mixture_prediction("MPD", f"BEA {label} - macrostate particle distribution", mpd=True))
    write(HERE / "mixture" / "iast" / "simulation.json",
          mixture_prediction("IAST", f"BEA {label} - dual-site Langmuir with IAST", mpd=False))
    write(HERE / "breakthrough" / "mpd" / "simulation.json",
          breakthrough("MPD", f"BEA {label} breakthrough - macrostate particle distribution", mpd=True))
    write(HERE / "breakthrough" / "iast" / "simulation.json",
          breakthrough("IAST", f"BEA {label} breakthrough - dual-site Langmuir with IAST", mpd=False))


if __name__ == "__main__":
    main()
