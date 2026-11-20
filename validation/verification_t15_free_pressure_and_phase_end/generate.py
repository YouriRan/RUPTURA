#!/usr/bin/env python3
"""Generate the T15 checks of the verification table: the free gas pressure and the ends of phases.

Both start from experiment T15 of Marx et al. (2014) with the inputs of their model
(../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json) on 30 cells:

  * free_pressure/<closure>/<thermal>: one 1500 s phase of the T15 feed at its velocity into the column at 1e6 Pa,
    with the gas leaving to 1e6 Pa at z = L, once with the free gas pressure (MomentumBalance, ports at both ends) and
    once with the continuity closure InletVelocityOutletPressure, both with ROS34PW2, isothermal and with the energy
    balance.  free_pressure/momentum_balance/isothermal_cvode_klu repeats the isothermal free-pressure run with CVODE
    and the KLU solver, the one CVODE setting the reader accepts with MomentumBalance;
  * phase_end/<integrator>: the T15 breakthrough (energy balance, InletPressureInletVelocity) as one phase that ends
    when the outlet CO2 mole fraction reaches half its feed value, 0.2, with CVODE and ROS34PW2 at their default
    tolerances and RK3 at 1e-2 s (a Courant number of 0.007 on 30 cells).

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI.
"""

from __future__ import annotations

import argparse
import copy
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parent.parent
BASE_INPUT = REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json"
SIMULATION_DIR = ROOT / "simulations"

GRID_POINTS = 30
DURATION = 1500.0  # s, beyond the complete CO2 breakthrough of T15
PRESSURE = 1.0e6  # Pa, inlet and outlet of T15
RK3_TIME_STEP = 1.0e-2  # s
HALF_FEED_CO2 = 0.2  # outlet CO2 mole fraction that ends the phase: half the 40 % of the feed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true", help="delete simulations/ first, including any output")
    arguments = parser.parse_args()
    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)

    base = json.loads(BASE_INPUT.read_text())
    base.update(SimulationType="SwingAdsorption", NumberOfGridPoints=GRID_POINTS, WriteInterval=1.0)
    base.pop("SimulationTime", None)
    index = []

    def write(relative: str, simulation: dict, **fields) -> None:
        directory = SIMULATION_DIR / relative
        directory.mkdir(parents=True, exist_ok=True)
        simulation["DisplayName"] = f"T15 verification {relative}"
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
        index.append({"directory": str(directory.relative_to(ROOT)), **fields})

    velocity = base["ColumnEntranceVelocity"]
    for thermal, energy_balance in (("isothermal", False), ("energy_balance", True)):
        momentum = copy.deepcopy(base)
        momentum.update(BoundaryCondition="MomentumBalance", BreakthroughIntegrator="ROS34PW2", OutletPressure=PRESSURE,
                        energyBalance=energy_balance)
        momentum["SwingAdsorptionPhases"] = [{
            "Name": "Adsorption", "Duration": DURATION,
            "FeedEnd": {"Mode": "Velocity", "Velocity": velocity, "Pressure": PRESSURE},
            "ProductEnd": {"Mode": "Pressure", "Pressure": PRESSURE}}]
        write(f"free_pressure/momentum_balance/{thermal}", momentum, check="free_pressure",
              closure="MomentumBalance", thermal=thermal, integrator="ROS34PW2")
        if not energy_balance:
            klu = copy.deepcopy(momentum)
            klu.update(BreakthroughIntegrator="CVODE", CVODELinearSolver="KLU")
            write(f"free_pressure/momentum_balance/{thermal}_cvode_klu", klu, check="free_pressure",
                  closure="MomentumBalance", thermal=thermal, integrator="CVODE KLU")

        closure = copy.deepcopy(base)
        closure.update(BoundaryCondition="InletVelocityOutletPressure", BreakthroughIntegrator="ROS34PW2",
                       OutletPressure=PRESSURE, energyBalance=energy_balance)
        closure["SwingAdsorptionPhases"] = [{"Name": "Adsorption", "Duration": DURATION}]
        write(f"free_pressure/continuity_closure/{thermal}", closure, check="free_pressure",
              closure="InletVelocityOutletPressure", thermal=thermal, integrator="ROS34PW2")

    for integrator in ("CVODE", "ROS34PW2", "RungeKutta3"):
        simulation = copy.deepcopy(base)
        simulation["BreakthroughIntegrator"] = integrator
        if integrator == "RungeKutta3":
            simulation["TimeStep"] = RK3_TIME_STEP
        simulation["SwingAdsorptionPhases"] = [{
            "Name": "Adsorption", "Duration": DURATION,
            "StopCriterion": {"Metric": "MoleFraction", "Component": "CO2", "End": "ProductEnd", "Operator": ">=",
                              "Value": HALF_FEED_CO2}}]
        write(f"phase_end/{integrator}", simulation, check="phase_end", integrator=integrator)

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
