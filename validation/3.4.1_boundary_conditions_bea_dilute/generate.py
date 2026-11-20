#!/usr/bin/env python3
"""Generate the dilute control case of the boundary-condition check: BEA with seven C7 isomers in helium.

The Marx et al. (2014) runs in ../3.4.1_boundary_conditions_marx_10bar take up a large share of the
feed at the adsorption front.  This control case is the opposite limit: seven alkanes at 1% each in
helium, so sorption changes the molar flow by a few percent at most.  It is the BEA base case of the
integrator benchmark (benchmark/generate.py), examples/BEA-alkanes-C7/breakthrough/simulation.json.

The example runs FixedPressureInletVelocity with PressureGradient = -1000 Pa/m.  That operating point is
kept: the inlet velocity of every other mode is the Ergun velocity of the feed at this gradient, and the
outlet pressure is P_in - G L, so that all five modes coincide once the column is filled with the feed,
as in the Marx case (../3.4.1_boundary_conditions_marx_10bar/conversions.md).

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

Only the boundary keys, the integration and output settings, and the display name differ from the
example.  All quantities are SI: Pa, K, m, m/s, kg/m3, kg/mol, Pa s.
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
BASE_INPUT = REPOSITORY / "examples" / "BEA-alkanes-C7" / "breakthrough" / "simulation.json"
SIMULATION_DIR = ROOT / "simulations"

GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1, the value in src/utils.h

# ROS34PW2 at its default tolerances, as in ../3.4.1_boundary_conditions_marx_10bar, adapts its step and ends one
# on every output.  Without physical axial dispersion, the Danckwerts and the Dirichlet inlet coincide here.
SIMULATION_TIME = 2000.0  # s: complete breakthrough of all seven isomers
WRITE_INTERVAL = 2.0  # s, one output block every 2 s

# JSON name of each gas-phase boundary mode and its label in ../3.4.1_boundary_conditions_marx_10bar.
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

    base = json.loads(BASE_INPUT.read_text())

    # Feed including the helium carrier, which here is part of the feed.
    fraction_sum = sum(component["GasPhaseMolFraction"] for component in base["Components"])
    molar_mass = sum(
        component["GasPhaseMolFraction"] / fraction_sum * component["MolecularWeight"] for component in base["Components"]
    )

    void_fraction = base["Geometry"]["ColumnVoidFraction"]
    particle_diameter = base["Geometry"]["ParticleDiameter"]
    solid_to_fluid = (1.0 - void_fraction) / void_fraction
    inlet_pressure = base["InletPressure"]
    temperature = base["Temperature"]
    pressure_gradient = -base["PressureGradient"]

    # Ergun velocity of the feed at the example's gradient, with the coefficients of src/geometry.cpp and the
    # stable root of src/geometry.h (../3.4.1_boundary_conditions_marx_10bar/conversions.md, section 4).
    feed_density = inlet_pressure * molar_mass / (GAS_CONSTANT * temperature)
    viscous = 150.0 * solid_to_fluid**2 * base["DynamicViscosity"] / particle_diameter**2
    inertial = 1.75 * solid_to_fluid * feed_density / particle_diameter
    velocity = 2.0 * pressure_gradient / (viscous + math.sqrt(viscous * viscous + 4.0 * inertial * pressure_gradient))
    outlet_pressure = inlet_pressure - pressure_gradient * base["ColumnLength"]

    label = f"BEA-C7_7x1pct-in-He_{temperature:.0f}K_{inlet_pressure:.1e}Pa".replace("e+0", "e")
    index = []
    for mode, paper_label in BOUNDARY_MODES.items():
        simulation = copy.deepcopy(base)
        for key in ("InletPressure", "OutletPressure", "PressureGradient"):
            simulation.pop(key, None)
        simulation["BoundaryCondition"] = mode
        simulation["DisplayName"] = f"BEA C7 dilute {mode}"
        simulation["ColumnEntranceVelocity"] = velocity
        simulation["BreakthroughIntegrator"] = "ROS34PW2"
        simulation["InletCondition"] = "Danckwerts"
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
                "boundary_condition": mode,
                "paper_label": paper_label,
                "feed_molar_mass_kg_per_mol": molar_mass,
                "feed_density_kg_per_m3": feed_density,
                "inlet_pressure_Pa": inlet_pressure,
                "outlet_pressure_Pa": outlet_pressure,
                "pressure_gradient_Pa_per_m": -pressure_gradient,
                "inlet_velocity_m_per_s": velocity,
            }
        )

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"feed velocity {velocity:.6f} m/s, outlet pressure {outlet_pressure:.3f} Pa")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
