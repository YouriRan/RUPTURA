#!/usr/bin/env python3
"""Generate the validation cases of the paper with the ROS34PW2 integrator and the Danckwerts inlet.

Every breakthrough and swing-adsorption input of the studies in STUDIES is copied to

    simulations/<study>/<case>/simulation.json

with "BreakthroughIntegrator": "ROS34PW2" and "InletCondition": "Danckwerts" (the default since 2026-11-03, written
out so that the copies do not depend on it). The settings of the previous integrators are removed: TimeStep and
NumberOfInitTimeSteps (RK3; ROS34PW2 rejects them) and every CVODE key. ROS34PW2 therefore runs at its default
tolerances (relative 1e-6; absolute 1e-7 mol/m^3, 1e-8 mol/kg, 1e-5 K) without a step cap. Nothing else changes, and
the source inputs are only read. The studies and case paths are those of ../inlet_condition_dirichlet_vs_danckwerts
where they overlap, so that compare.py finds the CVODE runs with the Danckwerts inlet there. The integrator benchmark
(benchmark/breakthrough) and the convection study compare integrators and schemes themselves; they are left out.

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI, as in the source inputs.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parent.parent
SIMULATION_DIR = ROOT / "simulations"

# Study name -> (directory the case paths are taken relative to, directories searched for simulation.json).
STUDIES = {
    "marx2_ternary": (REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary", [REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary"]),
    "stampi-bombelli": (
        REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli",
        [REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli" / "dl_from_literature", REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli" / "dl_from_stampibombelli"],
    ),
    "stadler": (REPOSITORY / "validation" / "3.3.5_reactions_stadler" / "stadler", [REPOSITORY / "validation" / "3.3.5_reactions_stadler" / "stadler"]),
    # Elfving is left out until its water model is decided: ROS34PW2 crawls on the stored inputs (TODO B-EX-9, B-EX-16).
    "liu": (REPOSITORY / "validation" / "3.3.4_humid_dac_liu" / "liu" / "simulations", [REPOSITORY / "validation" / "3.3.4_humid_dac_liu" / "liu" / "simulations"]),
    "3.4.1_boundary_conditions_marx_10bar": (
        REPOSITORY / "validation" / "3.4.1_boundary_conditions_marx_10bar" / "simulations",
        [REPOSITORY / "validation" / "3.4.1_boundary_conditions_marx_10bar" / "simulations"],
    ),
    "3.4.1_boundary_conditions_bea_dilute": (
        REPOSITORY / "validation" / "3.4.1_boundary_conditions_bea_dilute" / "simulations",
        [REPOSITORY / "validation" / "3.4.1_boundary_conditions_bea_dilute" / "simulations"],
    ),
    "3.4.2_swing_adsorption": (
        REPOSITORY / "validation" / "3.4.2_swing_adsorption" / "simulations",
        [REPOSITORY / "validation" / "3.4.2_swing_adsorption" / "simulations"],
    ),
    "3.4.3_sensitivity": (
        REPOSITORY / "validation" / "3.4.3_sensitivity" / "simulations",
        [REPOSITORY / "validation" / "3.4.3_sensitivity" / "simulations"],
    ),
}

SIMULATION_TYPES = {"Breakthrough", "SwingAdsorption"}
REMOVED_KEYS = {"breakthroughintegrator", "inletcondition", "timestep", "numberofinittimesteps"}


def source_inputs():
    """Yields (study, case path relative to the study directory, source simulation.json)."""
    for study, (base, roots) in STUDIES.items():
        for root in roots:
            for path in sorted(root.rglob("simulation.json")):
                # Archived and scratch directories (_legacy*, _to_delete, legacy) are not part of the studies.
                if any(part.startswith(("_", "legacy")) for part in path.relative_to(base).parts[:-1]):
                    continue
                yield study, path.parent.relative_to(base), path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true", help="delete simulations/ first, including any output")
    arguments = parser.parse_args()

    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)

    index = []
    for study, case, source in source_inputs():
        data = json.loads(source.read_text())
        keys = {key.lower(): key for key in data}
        if data.get(keys.get("simulationtype", "SimulationType"), "Breakthrough") not in SIMULATION_TYPES:
            continue
        integrator = data.get(keys.get("breakthroughintegrator", "BreakthroughIntegrator"), "RungeKutta3")
        data = {
            key: value
            for key, value in data.items()
            if key.lower() not in REMOVED_KEYS and not key.lower().startswith("cvode")
        }
        directory = SIMULATION_DIR / study / case
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "simulation.json").write_text(
            json.dumps({**data, "BreakthroughIntegrator": "ROS34PW2", "InletCondition": "Danckwerts"}, indent=4) + "\n"
        )
        index.append(
            {
                "study": study,
                "case": case.as_posix(),
                "directory": directory.relative_to(ROOT).as_posix(),
                "source": source.parent.relative_to(REPOSITORY).as_posix(),
                "source_integrator": integrator,
            }
        )

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} cases to {SIMULATION_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
