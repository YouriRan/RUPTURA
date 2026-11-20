#!/usr/bin/env python3
"""Generate CVODE copies of the 3.4.1 Marx runs without continuity of the molar flux.

The near-vacuum fronts of the FixedVelocity and FixedPressureInletVelocity settings
(../3.4.1_boundary_conditions_marx_10bar) stall CVODE.  This study copies those eight inputs (T3, T9, T15, T21) and
changes only the integrator: CVODE at its default settings (SPGMR, default tolerances) instead of ROS34PW2.  run.py
gives each run 900 s of wall time; how far a run gets, and its total gas concentration against the ROS34PW2 run up to
that time, is what the paper quotes (analysis.py).

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / "3.4.1_boundary_conditions_marx_10bar" / "simulations"
SIMULATION_DIR = ROOT / "simulations"
SETTINGS = ("FixedVelocity", "FixedPressureInletVelocity")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true", help="delete simulations/ first, including any output")
    arguments = parser.parse_args()
    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)

    index = []
    for source in sorted(SOURCE.glob("*/*/simulation.json")):
        if source.parent.name not in SETTINGS:
            continue
        simulation = json.loads(source.read_text())
        simulation["BreakthroughIntegrator"] = "CVODE"
        for key in [key for key in simulation if key.startswith(("ROS34PW2", "CVODE"))]:
            simulation.pop(key)
        simulation["DisplayName"] += " CVODE"
        directory = SIMULATION_DIR / source.parent.relative_to(SOURCE)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
        index.append({"directory": str(directory.relative_to(ROOT)), "reference": str(source.parent.relative_to(ROOT.parent))})

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
