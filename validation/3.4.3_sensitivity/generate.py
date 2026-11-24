#!/usr/bin/env python3
"""Generate the parameter-sensitivity study on the Marx et al. (2014) ternary experiment T15.

T15 (CO2/N2/H2 = 40/30/30, 298 K, 1e6 Pa, 1e-5 m3/s) is the strongest-sorbing
experiment at 1e6 Pa, is validated in the paper (Section 3.1), and is the same
experiment as the Marx case of the integrator benchmark (benchmark case
marx-403030_298_1e6_2e-1) and of the boundary-condition study (3.4.1).

Three kinds of runs are written, each a copy of marx2/ternary/T15/simulation.json
with one parameter changed:

  * reference           the unchanged base input;
  * local               each model parameter multiplied by 0.9 and 1.1, for the
                        normalized sensitivity coefficients;
  * sweep               wide one-at-a-time sweeps of the CO2 mass-transfer
                        coefficient, the axial dispersion coefficient, and the
                        number of grid points.  The dispersion and grid levels
                        are those of the integrator benchmark (benchmark/generate.py).

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

Reference:
    D. Marx, L. Joss, N. Casas, J. Schell, M. Mazzotti, Prediction of
    non-isothermal ternary gas-phase breakthrough experiments based on binary
    data, Adsorption 20 (2014) 493-510.  DOI 10.1007/s10450-013-9593-5.

All quantities are SI: Pa, K, m, s, m2/s, mol/kg, J/mol, W/(m2 K), J/(kg K), kg/m3.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parent.parent
BASE_INPUT = REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json"
SIMULATION_DIR = ROOT / "simulations"

# Grid and axial-dispersion levels of the integrator benchmark, loaded from benchmark/generate.py so that both
# studies stay on the same levels.  Two further levels, 1e-4 and 1e-3 m2/s, bracket the numerical dispersion of
# the first-order upwind convection term, v dz / 2 = 1.6e-4 m2/s at 100 grid points.  The benchmark level 1e-6 m2/s
# is left out: it lies between 0 and the T15 value 4.3e-6 m2/s, all far below the numerical dispersion.
_benchmark_spec = importlib.util.spec_from_file_location("benchmark_generate", REPOSITORY / "benchmark" / "generate.py")
_benchmark = importlib.util.module_from_spec(_benchmark_spec)
_benchmark_spec.loader.exec_module(_benchmark)
GRID_POINTS = _benchmark.GRID_POINTS
DISPERSION_COEFFICIENTS = (
    *(level for level in _benchmark.DISPERSION_COEFFICIENTS if level[0] != "1e-6"),
    ("1e-4", 1e-4),
    ("1e-3", 1e-3),
)

# ROS34PW2 (from the T15 input) adapts its internal step and ends a step at every output (README.md).
SIMULATION_TIME = 2000.0  # s: complete CO2 breakthrough and the outlet temperature peak
WRITE_INTERVAL = 2.0  # s, one output block every 2 s

LOCAL_MULTIPLIERS = (0.9, 1.1)
NUMERICAL_DISPERSION_CHECK = 8.7e-5  # m2/s, every component, with twice the grid points (see main)
MASS_TRANSFER_MULTIPLIERS = (0.25, 0.5, 2.0, 4.0)


def co2(simulation: dict) -> dict:
    return next(component for component in simulation["Components"] if component["Name"] == "CO2")


# Parameter name -> (container, key) of the scalar that is scaled.  Axial dispersion is scaled for every component
# and handled separately.  The three Sips parameters of CO2 are temperature-law objects (VantHoff q_sat and affinity,
# Arctan exponent; marx2/conversions.md, section 7), whose Value is scaled.  Every law carries its own energy, so
# HeatOfAdsorption sets only the heat released in the energy balance.
SCALED_PARAMETERS = {
    "saturation_loading_CO2": lambda simulation: (co2(simulation)["PhysisorptionSites"][0]["Parameters"][0], "Value"),
    "affinity_CO2": lambda simulation: (co2(simulation)["PhysisorptionSites"][0]["Parameters"][1], "Value"),
    "sips_exponent_CO2": lambda simulation: (co2(simulation)["PhysisorptionSites"][0]["Parameters"][2], "Value"),
    "heat_of_adsorption_CO2": lambda simulation: (co2(simulation), "HeatOfAdsorption"),
    "mass_transfer_coefficient_CO2": lambda simulation: (co2(simulation), "MassTransferCoefficient"),
    "particle_density": lambda simulation: (simulation, "ParticleDensity"),
    "heat_transfer_gas_solid": lambda simulation: (simulation, "heatTransferGasSolid"),
    "heat_transfer_gas_wall": lambda simulation: (simulation, "heatTransferGasWall"),
    "heat_transfer_wall_external": lambda simulation: (simulation, "heatTransferWallExternal"),
    "heat_capacity_solid": lambda simulation: (simulation, "heatCapacitySolid"),
    "heat_capacity_wall": lambda simulation: (simulation, "heatCapacityWall"),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete simulations/ first, discarding any Ruptura output already there",
    )
    arguments = parser.parse_args()

    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)
    SIMULATION_DIR.mkdir(parents=True, exist_ok=True)

    base = json.loads(BASE_INPUT.read_text())
    base["SimulationTime"] = SIMULATION_TIME
    base["WriteInterval"] = WRITE_INTERVAL
    base["WriteOutletOnly"] = True  # the analysis reads the outlet only

    index = []

    def write(relative_directory: str, simulation: dict, parameter: str, kind: str, multiplier, value) -> None:
        directory = SIMULATION_DIR / relative_directory
        directory.mkdir(parents=True, exist_ok=True)
        simulation["DisplayName"] = f"Marx T15 sensitivity {relative_directory}"
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
        index.append(
            {
                "directory": str(directory.relative_to(ROOT)),
                "parameter": parameter,
                "kind": kind,
                "multiplier": multiplier,
                "value": value,
            }
        )

    write("reference", copy.deepcopy(base), "reference", "reference", 1.0, None)

    multipliers = {name: LOCAL_MULTIPLIERS for name in SCALED_PARAMETERS}
    multipliers["mass_transfer_coefficient_CO2"] = tuple(sorted((*LOCAL_MULTIPLIERS, *MASS_TRANSFER_MULTIPLIERS)))
    for parameter, accessor in SCALED_PARAMETERS.items():
        for multiplier in multipliers[parameter]:
            simulation = copy.deepcopy(base)
            container, key = accessor(simulation)
            container[key] *= multiplier
            kind = "local" if multiplier in LOCAL_MULTIPLIERS else "sweep"
            write(f"{parameter}/x{multiplier:.2f}", simulation, parameter, kind, multiplier, container[key])

    for multiplier in LOCAL_MULTIPLIERS:
        simulation = copy.deepcopy(base)
        for component in simulation["Components"]:
            component["AxialDispersionCoefficient"] *= multiplier
        write(f"axial_dispersion_coefficient/x{multiplier:.2f}", simulation, "axial_dispersion_coefficient", "local",
              multiplier, None)

    # As in the integrator benchmark, a dispersion level is applied to every component, including the carrier gas.
    for label, dispersion in DISPERSION_COEFFICIENTS:
        simulation = copy.deepcopy(base)
        for component in simulation["Components"]:
            component["AxialDispersionCoefficient"] = dispersion
        write(f"axial_dispersion_coefficient/{label}", simulation, "axial_dispersion_coefficient", "sweep", None,
              dispersion)

    for grid_points in GRID_POINTS:
        if grid_points == base["NumberOfGridPoints"]:
            continue
        simulation = copy.deepcopy(base)
        simulation["NumberOfGridPoints"] = grid_points
        write(f"grid_points/N{grid_points}", simulation, "grid_points", "sweep", None, grid_points)

    # The first-order upwind convection term adds a numerical dispersion of about v dz / 2.  Halving dz and adding
    # the difference, 0.0275 * (0.012 - 0.006) / 2 = 8.2e-5 m2/s, plus the T15 value, as physical dispersion should
    # reproduce the reference run.
    simulation = copy.deepcopy(base)
    simulation["NumberOfGridPoints"] = 2 * base["NumberOfGridPoints"]
    for component in simulation["Components"]:
        component["AxialDispersionCoefficient"] = NUMERICAL_DISPERSION_CHECK
    write("numerical_dispersion_check/N200_D8.7e-5", simulation, "numerical_dispersion_check", "check", None,
          NUMERICAL_DISPERSION_CHECK)

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
