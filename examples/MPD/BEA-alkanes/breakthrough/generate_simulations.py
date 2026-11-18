#!/usr/bin/env python3
"""Materialize the 24 BEA MPD/Langmuir breakthrough benchmark inputs."""

from __future__ import annotations

from copy import deepcopy
from itertools import product
from pathlib import Path
import csv
import json
import math


OUTPUT_DIRECTORY = Path(__file__).resolve().parent
EXAMPLE_DIRECTORY = OUTPUT_DIRECTORY.parent

MODELS = ("MPD", "Langmuir")
INTEGRATORS = ("CVODE", "RungeKutta3")
GRID_POINTS = (20, 50, 100)
VELOCITIES = (0.019, 0.19)
COMPONENTS = ("nC7", "C6m2", "C6m3", "C5m2m2", "C5m2m3", "C5m2m4", "C5m3m3")

TARGET_TAU = 900.0
CFL = 0.4
NUMBER_OF_SNAPSHOTS = 180
MOLECULAR_WEIGHTS = {"Helium": 0.0040026, **{name: 0.100198 for name in COMPONENTS}}


def read_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


MPD_MIXTURE_INPUT = read_json(EXAMPLE_DIRECTORY / "mpd" / "simulation.json")
LANGMUIR_MIXTURE_INPUT = read_json(EXAMPLE_DIRECTORY / "langmuir" / "simulation.json")


def velocity_tag(velocity: float) -> str:
    return f"{velocity:g}".replace(".", "p")


def run_name(model: str, integrator: str, grid_points: int, velocity: float) -> str:
    return f"{model.lower()}_{integrator.lower()}_n{grid_points}_v{velocity_tag(velocity)}_r1"


def make_simulation(model: str, integrator: str, grid_points: int, velocity: float) -> dict:
    source = MPD_MIXTURE_INPUT if model == "MPD" else LANGMUIR_MIXTURE_INPUT
    components = deepcopy(source["Components"])
    for component in components:
        component["MolecularWeight"] = MOLECULAR_WEIGHTS[component["Name"]]
        if not component.get("CarrierGas", False):
            component["MassTransferCoefficient"] = 0.06
            component["AxialDispersionCoefficient"] = 0.0

    column_length = 0.1
    time_step = CFL * column_length / (velocity * grid_points)
    requested_steps = math.ceil(TARGET_TAU * column_length / (velocity * time_step))
    write_every = max(1, math.ceil(requested_steps / NUMBER_OF_SNAPSHOTS))
    number_of_steps = write_every * NUMBER_OF_SNAPSHOTS + 1

    simulation = {
        "SimulationType": "Breakthrough",
        "MixturePredictionMethod": source["MixturePredictionMethod"],
        "DisplayName": f"BEA C7 {model} breakthrough benchmark",
        "Temperature": 552.0,
        "ParticleDensity": 1508.52,
        "InletPressure": 120000.0,
        "PressureGradient": -1000.0,
        "ColumnEntranceVelocity": velocity,
        "ColumnLength": column_length,
        "SimulationTime": number_of_steps * time_step,
        "WriteInterval": write_every * time_step,
        "BreakthroughIntegrator": integrator,
        "BoundaryCondition": "FixedPressureInletVelocity",
        "DynamicViscosity": 3.0e-5,
        # CVODE chooses its own internal steps and takes no TimeStep.
        **({} if integrator == "CVODE" else {"TimeStep": time_step}),
        "NumberOfGridPoints": grid_points,
        "energyBalance": False,
        "Components": components,
        "Geometry": {
            "Type": "PackedBed",
            "ColumnVoidFraction": 0.4,
            "ParticleDiameter": 0.001,
        },
    }
    if model == "MPD":
        simulation["MPDSettings"] = deepcopy(source["MPDSettings"])
        simulation["MPDSettings"]["FileName"] = "../../mpd/particle_distribution.data"
    return simulation


def main() -> None:
    manifest = []
    for model, integrator, grid_points, velocity in product(MODELS, INTEGRATORS, GRID_POINTS, VELOCITIES):
        name = run_name(model, integrator, grid_points, velocity)
        directory = OUTPUT_DIRECTORY / name
        directory.mkdir(parents=True, exist_ok=True)
        simulation = make_simulation(model, integrator, grid_points, velocity)
        with (directory / "simulation.json").open("w", encoding="utf-8") as stream:
            json.dump(simulation, stream, indent=2)
            stream.write("\n")
        manifest.append(
            {
                "directory": name,
                "model": model,
                "integrator": integrator,
                "grid_points": grid_points,
                "velocity_m_s": velocity,
                "time_step_s": simulation.get("TimeStep", ""),
                "simulation_time_s": simulation["SimulationTime"],
                "write_interval_s": simulation["WriteInterval"],
            }
        )

    with (OUTPUT_DIRECTORY / "manifest.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=manifest[0].keys())
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Wrote {len(manifest)} simulations beneath {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
