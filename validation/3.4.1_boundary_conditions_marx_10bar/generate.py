#!/usr/bin/env python3
"""Generate the boundary-condition check on the Marx et al. (2014) ternary runs at 1e6 Pa.

Each base experiment is simulated once per gas-phase pressure-velocity boundary
mode.  The modes are parameterized so that a column filled with the feed, i.e.
after complete breakthrough, has the same inlet pressure, outlet pressure, and
inlet velocity in every mode.  Differences during breakthrough therefore come
only from how each mode couples pressure and velocity.  conversions.md derives
the pressure drop that sets the outlet pressure and the pressure gradient.

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

The base inputs are read from ../3.3.1_energy_balance_marx_ternary/marx2/ternary/<run>/simulation.json.  Only
the boundary-condition keys, the simulated duration, the output interval, the
display name, and the integration are changed: ROS34PW2 at its default
tolerances instead of CVODE, which stalls in the settings without continuity of
the molar flux, and the Danckwerts inlet.  Re-running this script rewrites every
simulation.json but leaves existing Ruptura output in place.

Reference:
    D. Marx, L. Joss, N. Casas, J. Schell, M. Mazzotti, Prediction of
    non-isothermal ternary gas-phase breakthrough experiments based on binary
    data, Adsorption 20 (2014) 493-510.  DOI 10.1007/s10450-013-9593-5.

All quantities are SI: Pa, K, m, m/s, m3/s, kg/m3, kg/mol, Pa s.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parent.parent
MARX_DIR = REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2"
SIMULATION_DIR = ROOT / "simulations"

GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1, the value in src/utils.h

# The four experiments at 1e6 Pa and the lowest flow rate, 1e-5 m3/s: both feeds
# at both temperatures.  T21 is included because its marx2 run had to switch to
# InletVelocityOutletPressure and was stopped at 1174 s.
BASE_RUNS = ["T3", "T9", "T15", "T21"]

# 4000 s covers the complete breakthrough of all three components in every base run.
SIMULATION_TIME = 4000.0  # s

# One output block every 4 s instead of every 1 s keeps the 20 runs near 2 GB.
WRITE_INTERVAL = 4.0  # s

# JSON name of each gas-phase boundary mode and its label in the SI table of the paper.
BOUNDARY_MODES = {
    "InletVelocityOutletPressure": "BC1",
    "InletPressureOutletPressure": "BC2",
    "InletPressureInletVelocity": "BC3",
    "FixedPressureInletVelocity": "BC4",
    "FixedVelocity": "BC5",
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

    index = []
    for run in BASE_RUNS:
        base = json.loads((MARX_DIR / "ternary" / run / "simulation.json").read_text())

        # Feed of the adsorbing components; the helium carrier only fills the initial column.
        adsorbates = [component for component in base["Components"] if not component.get("CarrierGas", False)]
        fraction_sum = sum(component["GasPhaseMolFraction"] for component in adsorbates)
        feed = {component["Name"]: component["GasPhaseMolFraction"] / fraction_sum for component in adsorbates}
        molar_mass = sum(feed[component["Name"]] * component["MolecularWeight"] for component in adsorbates)

        geometry = base["Geometry"]
        void_fraction = geometry["ColumnVoidFraction"]
        particle_diameter = geometry["ParticleDiameter"]
        solid_to_fluid = (1.0 - void_fraction) / void_fraction
        inlet_pressure = base["InletPressure"]
        feed_temperature = base["InfluxTemperature"]
        velocity = base["ColumnEntranceVelocity"]

        # Ergun gradient -dP/dz of a feed-filled column, with the coefficients of src/geometry.cpp.
        feed_density = inlet_pressure * molar_mass / (GAS_CONSTANT * feed_temperature)
        pressure_gradient = (
            150.0 * solid_to_fluid**2 / particle_diameter**2 * base["DynamicViscosity"] * velocity
            + 1.75 * solid_to_fluid / particle_diameter * feed_density * velocity**2
        )
        outlet_pressure = inlet_pressure - pressure_gradient * base["ColumnLength"]

        # The fluid fraction of the Marx inputs is the total porosity, with which marx2 converts the volumetric feed
        # flow to ColumnEntranceVelocity (../3.3.1_energy_balance_marx_ternary/marx2/conversions.md, section 6).
        flow_rate = velocity * math.pi * geometry["InternalDiameter"] ** 2 / 4.0 * void_fraction
        composition = "-".join(f"{100.0 * fraction:.0f}" for fraction in feed.values())
        label = (f"{run}_{composition}_{feed_temperature:.0f}K_{inlet_pressure:.0e}Pa_{flow_rate:.0e}m3s").replace(
            "e+0", "e"
        ).replace("e-0", "e-")

        for mode, paper_label in BOUNDARY_MODES.items():
            simulation = copy.deepcopy(base)
            for key in ("InletPressure", "OutletPressure", "PressureGradient"):
                simulation.pop(key, None)
            for key in [key for key in simulation if key.startswith("CVODE")]:
                simulation.pop(key)
            simulation["BreakthroughIntegrator"] = "ROS34PW2"
            simulation["InletCondition"] = "Danckwerts"
            simulation["BoundaryCondition"] = mode
            simulation["DisplayName"] = f"Marx {run} {mode}"
            simulation["SimulationTime"] = SIMULATION_TIME
            simulation["WriteInterval"] = WRITE_INTERVAL
            if mode != "InletVelocityOutletPressure":
                simulation["InletPressure"] = inlet_pressure
            if mode in ("InletPressureOutletPressure", "InletVelocityOutletPressure"):
                simulation["OutletPressure"] = outlet_pressure
            if mode == "FixedPressureInletVelocity":
                simulation["PressureGradient"] = -pressure_gradient

            directory = SIMULATION_DIR / label / mode
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
            index.append(
                {
                    "directory": str(directory.relative_to(ROOT)),
                    "run": run,
                    "boundary_condition": mode,
                    "paper_label": paper_label,
                    "feed_mole_fractions": feed,
                    "feed_temperature_K": feed_temperature,
                    "inlet_pressure_Pa": inlet_pressure,
                    "outlet_pressure_Pa": outlet_pressure,
                    "pressure_gradient_Pa_per_m": -pressure_gradient,
                    "inlet_velocity_m_per_s": velocity,
                    "flow_rate_m3_per_s": flow_rate,
                    "reference_data": str((MARX_DIR / "reference_digitized" / run).relative_to(REPOSITORY)),
                }
            )

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
