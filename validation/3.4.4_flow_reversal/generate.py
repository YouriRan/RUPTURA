#!/usr/bin/env python3
"""Generate the flow-reversal study: co-current and counter-current regeneration against equilibrium theory.

A trace adsorbate A (mole fraction 1e-3 in helium) with a favourable Langmuir isotherm, b p_F = 4, is fed to a clean
column until its shock front has crossed half the column, and the column is then purged with helium, either
co-currently (FlowDirection Forward, the purge enters at z = 0) or counter-currently (FlowDirection Reverse, the purge
enters at z = L). Both runs are two-phase SwingAdsorption runs. In the trace limit the velocity is constant, which
FixedVelocity imposes, the pressure is uniform, and with fast mass transfer and no axial dispersion the column follows
the equilibrium theory of Rhee, Aris, and Amundson (method of characteristics), against which analysis.ipynb compares
the outlet histories and the column profiles (README.md).

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI: Pa, K, s, m, mol/kg.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SIMULATION_DIR = ROOT / "simulations"

TEMPERATURE = 298.15  # K
PRESSURE = 1.0e5  # Pa
VELOCITY = 0.05  # m/s, interstitial
LENGTH = 0.1  # m
VOID_FRACTION = 0.4
PARTICLE_DENSITY = 1000.0  # kg/m3
PARTICLE_DIAMETER = 2.0e-3  # m
VISCOSITY = 2.0e-5  # Pa s
FEED_MOLE_FRACTION = 1.0e-3
SATURATION_LOADING = 1.0  # mol/kg
AFFINITY = 0.04  # 1/Pa, b p_F = 4 at the feed partial pressure of 100 Pa
MASS_TRANSFER_COEFFICIENT = 1.0  # 1/s, local equilibrium on the time scale of the fronts (about 1e4 s)

LOADING_DURATION = 30_000.0  # s, the shock front reaches z = 0.504 L
PURGE_DURATION = 320_000.0  # s, beyond the end of the co-current desorption tail (297 500 s)
WRITE_INTERVAL = 500.0  # s

GRID_POINTS = (25, 50, 100, 200, 400)
SCHEMES = ("WENO3", "Upwind")
DIRECTIONS = {"co-current": "Forward", "counter-current": "Reverse"}


def simulation(grid_points: int, scheme: str, direction: str) -> dict:
    return {
        "SimulationType": "SwingAdsorption",
        "DisplayName": f"Trace Langmuir adsorbate, {direction} helium purge, {scheme}, {grid_points} cells",
        "BreakthroughIntegrator": "ROS34PW2",
        "ConvectionScheme": scheme,
        "InletCondition": "Danckwerts",
        "BoundaryCondition": "FixedVelocity",
        "Temperature": TEMPERATURE,
        "InletPressure": PRESSURE,
        "ColumnEntranceVelocity": VELOCITY,
        "ColumnLength": LENGTH,
        "ColumnVoidFraction": VOID_FRACTION,
        "ParticleDensity": PARTICLE_DENSITY,
        "ParticleDiameter": PARTICLE_DIAMETER,
        "DynamicViscosity": VISCOSITY,
        "NumberOfGridPoints": grid_points,
        "WriteInterval": WRITE_INTERVAL,
        "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": VOID_FRACTION, "ParticleDiameter": PARTICLE_DIAMETER},
        "Components": [
            {"Name": "Helium", "GasPhaseMolFraction": 1.0 - FEED_MOLE_FRACTION, "MolecularWeight": 0.0040026,
             "CarrierGas": True},
            {"Name": "A", "GasPhaseMolFraction": FEED_MOLE_FRACTION, "MolecularWeight": 0.044,
             "MassTransferCoefficient": MASS_TRANSFER_COEFFICIENT, "AxialDispersionCoefficient": 0.0,
             "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [SATURATION_LOADING, AFFINITY]}]},
        ],
        "SwingAdsorptionPhases": [
            {"Name": "Loading", "Duration": LOADING_DURATION},
            {"Name": "Purge", "Duration": PURGE_DURATION, "FlowDirection": DIRECTIONS[direction],
             "GasPhaseMolFractions": {"Helium": 1.0}},
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true",
                        help="delete simulations/ first, discarding any Ruptura output already there")
    arguments = parser.parse_args()

    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)
    SIMULATION_DIR.mkdir(parents=True, exist_ok=True)

    index = []
    for direction in DIRECTIONS:
        for scheme in SCHEMES:
            for grid_points in GRID_POINTS:
                directory = SIMULATION_DIR / f"{direction}_purge_{scheme}_N{grid_points}"
                directory.mkdir(parents=True, exist_ok=True)
                data = simulation(grid_points, scheme, direction)
                (directory / "simulation.json").write_text(json.dumps(data, indent=4) + "\n")
                index.append({"directory": str(directory.relative_to(ROOT)), "direction": direction,
                              "scheme": scheme, "grid_points": grid_points})

    parameters = {
        "temperature_K": TEMPERATURE, "pressure_Pa": PRESSURE, "velocity_m_per_s": VELOCITY, "length_m": LENGTH,
        "void_fraction": VOID_FRACTION, "particle_density_kg_per_m3": PARTICLE_DENSITY,
        "feed_mole_fraction": FEED_MOLE_FRACTION, "saturation_loading_mol_per_kg": SATURATION_LOADING,
        "affinity_per_Pa": AFFINITY, "loading_duration_s": LOADING_DURATION, "purge_duration_s": PURGE_DURATION,
    }
    (SIMULATION_DIR / "index.json").write_text(json.dumps({"parameters": parameters, "runs": index}, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
