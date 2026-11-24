#!/usr/bin/env python3
"""Generate the convection-scheme comparison on T15 quoted in the SI (Method of lines): the CO2 front width t90 - t10
on 40 cells with upwinding, TVD, WENO3, and WENO5, each with ROS34PW2 and CVODE at their default settings.

Base: experiment T15 of Marx et al. (2014) with the inputs of their model
(../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json), 1500 s.

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
BASE_INPUT = ROOT.parent / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json"
SIMULATION_DIR = ROOT / "simulations"

GRID_POINTS = 40
SIMULATION_TIME = 1500.0  # s, beyond the complete CO2 breakthrough
SCHEMES = ("Upwind", "TVD", "WENO3", "WENO5")
INTEGRATORS = ("ROS34PW2", "CVODE")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true", help="delete simulations/ first, including any output")
    arguments = parser.parse_args()
    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)

    base = json.loads(BASE_INPUT.read_text())
    base.update(NumberOfGridPoints=GRID_POINTS, SimulationTime=SIMULATION_TIME, WriteInterval=1.0,
                WriteOutletOnly=True)  # analysis.py reads the outlet only
    index = []
    for scheme in SCHEMES:
        for integrator in INTEGRATORS:
            simulation = copy.deepcopy(base)
            simulation.update(ConvectionScheme=scheme, BreakthroughIntegrator=integrator,
                              DisplayName=f"T15 {scheme} {integrator} {GRID_POINTS} cells")
            directory = SIMULATION_DIR / scheme / integrator
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
            index.append({"directory": str(directory.relative_to(ROOT)), "scheme": scheme, "integrator": integrator})
    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
