#!/usr/bin/env python3
"""Generate the Ruptura integrator benchmark simulation matrix."""

from __future__ import annotations

import csv
import itertools
import json
import math
from copy import deepcopy
from pathlib import Path


BENCHMARK_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = BENCHMARK_ROOT.parent

BASE_CASES = {
    "bea-alkanes-C7": REPOSITORY_ROOT
    / "examples"
    / "BEA-alkanes-C7"
    / "breakthrough"
    / "simulation.json",
}

GRID_POINTS = (10, 20, 50, 100, 200)
TIME_STEPS = (
    ("1e-5", 1e-5),
    ("3.16e-5", 3.16e-5),
    ("1e-4", 1e-4),
    ("3.16e-4", 3.16e-4),
    ("1e-3", 1e-3),
    ("3.16e-3", 3.16e-3),
    ("1e-2", 1e-2),
)
DISPERSION_COEFFICIENTS = (
    ("0", 0.0),
    ("1e-6", 1e-6),
    ("1e-5", 1e-5),
)
INTEGRATORS = {
    "rk3": "RungeKutta3",
    "cvode": "CVODE",
    "sirk3": "SIRK3",
}
BENCHMARK_DURATION_SECONDS = 0.1


def load_base_case(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        settings = json.load(stream)

    components = settings.get("Components")
    if not isinstance(components, list) or not components:
        raise ValueError(f"Base case has no components: {path}")
    return settings


def generate() -> int:
    rows = []

    for case_name, source_path in BASE_CASES.items():
        base = load_base_case(source_path)

        combinations = itertools.product(
            GRID_POINTS,
            TIME_STEPS,
            DISPERSION_COEFFICIENTS,
            INTEGRATORS.items(),
        )
        for grid_points, (dt_label, dt), (disp_label, dispersion), (
            integrator_label,
            integrator,
        ) in combinations:
            simulation = deepcopy(base)
            number_of_time_steps = max(1, math.ceil(BENCHMARK_DURATION_SECONDS / dt))
            simulation["NumberOfGridPoints"] = grid_points
            simulation["TimeStep"] = dt
            simulation["NumberOfTimeSteps"] = number_of_time_steps
            simulation["BreakthroughIntegrator"] = integrator

            for component in simulation["Components"]:
                component["AxialDispersionCoefficient"] = dispersion

            relative_directory = (
                Path(case_name)
                / f"grid-{grid_points}"
                / f"dt-{dt_label}"
                / f"dispersion-{disp_label}"
                / integrator_label
            )
            output_directory = BENCHMARK_ROOT / relative_directory
            output_directory.mkdir(parents=True, exist_ok=True)

            simulation_path = output_directory / "simulation.json"
            with simulation_path.open("w", encoding="utf-8") as stream:
                json.dump(simulation, stream, indent=2)
                stream.write("\n")

            rows.append(
                {
                    "directory": relative_directory.as_posix(),
                    "base_case": case_name,
                    "source": source_path.relative_to(REPOSITORY_ROOT).as_posix(),
                    "grid_points": grid_points,
                    "time_step": dt_label,
                    "axial_dispersion_coefficient": disp_label,
                    "integrator": integrator_label,
                    "ruptura_integrator": integrator,
                    "number_of_time_steps": number_of_time_steps,
                    "simulated_duration_seconds": number_of_time_steps * dt,
                }
            )

    manifest_path = BENCHMARK_ROOT / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


if __name__ == "__main__":
    count = generate()
    print(f"Generated {count} simulations below {BENCHMARK_ROOT}")
