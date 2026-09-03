#!/usr/bin/env python3
"""Generate the simulation matrices used by the publication-style figures."""

from __future__ import annotations

import csv
import itertools
import json
import math
from copy import deepcopy
from pathlib import Path


FIGURE_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = FIGURE_ROOT.parents[1]
SIMULATION_ROOT = FIGURE_ROOT / "simulations"

BASE_CASES = {
    "bea-alkanes-C7": REPOSITORY_ROOT
    / "examples"
    / "BEA-alkanes-C7"
    / "breakthrough"
    / "simulation.json",
    "marx-403030_298_1e6_2e-1": REPOSITORY_ROOT
    / "marx"
    / "403030_298_1e6_2e-1"
    / "simulation.json",
}

DEFAULT_NCOMP = {"bea-alkanes-C7": 7, "marx-403030_298_1e6_2e-1": 3}
NCOMP_SWEEP = {"bea-alkanes-C7": (2, 4, 6), "marx-403030_298_1e6_2e-1": (1, 2, 3)}
GRID_SWEEP = (10, 20, 100)
BEA_TIME_STEPS = (
    ("1e-5", 1e-5),
    ("3.16e-5", 3.16e-5),
    ("1e-4", 1e-4),
    ("3.16e-4", 3.16e-4),
    ("1e-3", 1e-3),
    ("3.16e-3", 3.16e-3),
    ("1e-2", 1e-2),
)
MARX_TIME_STEPS = (
    ("1e-6", 1e-6),
    ("3.16e-6", 3.16e-6),
    ("1e-5", 1e-5),
    ("3.16e-5", 3.16e-5),
    ("1e-4", 1e-4),
    ("3.16e-4", 3.16e-4),
    ("1e-3", 1e-3),
)
INTEGRATORS = {"rk3": "RungeKutta3", "cvode": "CVODE"}
METHODS = ("IAST", "SIAST", "SEI", "SPI")
HORIZON_SECONDS = 0.1
METHOD_TIME_STEP = 5e-4


CHEMISORPTION_SITE = {
    "Type": "General",
    "Parameters": {
        "maximumLoading": 0.55,
        "heatOfChemisorption": 52000.0,
        "adsorptionRateCoefficient": 2.5e-5,
        "adsorptionActivationEnergy": 0.0,
        "desorptionRateCoefficient": 1e-6,
        "desorptionActivationEnergy": 0.0,
        "poreConcentrationOrder": 1,
        "capacityOrder": 1,
        "desorptionOrder": 1,
        "filmMassTransferCoefficient": 0.0015,
        "poreDiffusivity": 8e-11,
        "usePoreSurfaceTransport": True,
        "Isotherm": {"Type": "Langmuir", "Parameters": [0.55, 1.8e-5]},
    },
}


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def make_current_geometry(simulation: dict) -> None:
    geometry = simulation.get("Geometry")
    if isinstance(geometry, dict) and geometry.get("Type") == "HollowTube":
        geometry["Type"] = "PackedBed"


def select_adsorbing_components(simulation: dict, ncomp: int) -> None:
    components = simulation["Components"]
    carrier = next(component for component in components if component.get("CarrierGas", False))
    adsorbates = [component for component in components if component is not carrier][:ncomp]
    if len(adsorbates) != ncomp:
        raise ValueError(f"Requested {ncomp} adsorbates but only found {len(adsorbates)}")

    adsorbate_fraction = sum(float(component.get("GasPhaseMolFraction", 0.0)) for component in adsorbates)
    carrier["GasPhaseMolFraction"] = 1.0 - adsorbate_fraction
    simulation["Components"] = [carrier, *adsorbates]


def set_run_controls(simulation: dict, time_step: float, write_final_state: bool) -> int:
    steps = max(2, math.ceil(HORIZON_SECONDS / time_step))
    simulation["TimeStep"] = time_step
    simulation["NumberOfTimeSteps"] = steps
    simulation["PrintEvery"] = steps + 1
    simulation["WriteEvery"] = steps - 1 if write_final_state else steps + 1
    return steps


def use_single_langmuir_sites(simulation: dict) -> None:
    """Project Marx Sips sites to the common model supported by all four methods."""
    for component in simulation["Components"]:
        for site in component.get("PhysisorptionSites", []):
            if site.get("Type") == "Sips":
                site["Type"] = "Langmuir"
                site["Parameters"] = site["Parameters"][:2]


def write_simulation(relative_directory: Path, simulation: dict) -> None:
    directory = SIMULATION_ROOT / relative_directory
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "simulation.json").open("w", encoding="utf-8") as stream:
        json.dump(simulation, stream, indent=2)
        stream.write("\n")


def generate_integrator_rows() -> list[dict]:
    rows = []
    for case, source in BASE_CASES.items():
        base = load_json(source)
        default_ncomp = DEFAULT_NCOMP[case]
        time_steps = MARX_TIME_STEPS if case.startswith("marx-") else BEA_TIME_STEPS
        physical_settings = {(grid, default_ncomp) for grid in GRID_SWEEP}
        physical_settings.update((100, ncomp) for ncomp in NCOMP_SWEEP[case])

        for (grid_points, ncomp), (dt_label, time_step), (integrator, ruptura_integrator) in itertools.product(
            sorted(physical_settings), time_steps, INTEGRATORS.items()
        ):
            simulation = deepcopy(base)
            make_current_geometry(simulation)
            select_adsorbing_components(simulation, ncomp)
            steps = set_run_controls(simulation, time_step, write_final_state=True)
            simulation["NumberOfGridPoints"] = grid_points
            simulation["BreakthroughIntegrator"] = ruptura_integrator

            relative_directory = (
                Path("cvode-rk3")
                / case
                / f"grid-{grid_points}"
                / f"ncomp-{ncomp}"
                / f"dt-{dt_label}"
                / integrator
            )
            write_simulation(relative_directory, simulation)
            rows.append(
                {
                    "suite": "cvode-rk3",
                    "directory": relative_directory.as_posix(),
                    "case": case,
                    "grid_points": grid_points,
                    "ncomp": ncomp,
                    "time_step": time_step,
                    "integrator": integrator,
                    "method": simulation.get("MixturePredictionMethod", "IAST"),
                    "energy_balance": str(bool(simulation.get("energyBalance", False))).lower(),
                    "chemisorption_transport": "false",
                    "number_of_time_steps": steps,
                    "simulated_duration_seconds": steps * time_step,
                }
            )
    return rows


def generate_method_rows() -> list[dict]:
    rows = []
    case = "marx-403030_298_1e6_2e-1"
    base = load_json(BASE_CASES[case])
    for method, (integrator, ruptura_integrator), energy_balance, chemisorption_transport in itertools.product(
        METHODS, INTEGRATORS.items(), (False, True), (False, True)
    ):
        simulation = deepcopy(base)
        make_current_geometry(simulation)
        select_adsorbing_components(simulation, DEFAULT_NCOMP[case])
        use_single_langmuir_sites(simulation)
        steps = set_run_controls(simulation, METHOD_TIME_STEP, write_final_state=False)
        simulation["NumberOfGridPoints"] = 100
        simulation["BreakthroughIntegrator"] = ruptura_integrator
        simulation["MixturePredictionMethod"] = method
        simulation["energyBalance"] = energy_balance
        if chemisorption_transport:
            simulation["Components"][1]["ChemisorptionSites"] = [deepcopy(CHEMISORPTION_SITE)]
        else:
            for component in simulation["Components"]:
                component.pop("ChemisorptionSites", None)

        relative_directory = (
            Path("methods")
            / case
            / f"energy-{'on' if energy_balance else 'off'}"
            / f"chemisorption-transport-{'on' if chemisorption_transport else 'off'}"
            / method
            / integrator
        )
        write_simulation(relative_directory, simulation)
        rows.append(
            {
                "suite": "methods",
                "directory": relative_directory.as_posix(),
                "case": case,
                "grid_points": 100,
                "ncomp": DEFAULT_NCOMP[case],
                "time_step": METHOD_TIME_STEP,
                "integrator": integrator,
                "method": method,
                "energy_balance": str(energy_balance).lower(),
                "chemisorption_transport": str(chemisorption_transport).lower(),
                "number_of_time_steps": steps,
                "simulated_duration_seconds": steps * METHOD_TIME_STEP,
            }
        )
    return rows


def main() -> None:
    rows = [*generate_integrator_rows(), *generate_method_rows()]
    manifest_path = FIGURE_ROOT / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Generated {len(rows)} figure simulations below {SIMULATION_ROOT}")


if __name__ == "__main__":
    main()
