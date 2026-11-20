#!/usr/bin/env python3
"""Generate the validation cases of the paper with a Dirichlet and with a Danckwerts inlet.

Every breakthrough and swing-adsorption input of the studies in STUDIES is copied twice, once with
"InletCondition": "Dirichlet" (the default: the inlet face holds the feed state) and once with "Danckwerts" (the
total convective-dispersive inflow equals the feed flux), into

    simulations/<inlet>/<study>/<case>/simulation.json

Nothing else in the inputs changes, except that CVODE inputs without a CVODELinearSolver get "Dense", the default
with which their stored results were computed (the default is SPGMR since 2026-11-02, TODO V-GEN-4; with SPGMR the
steep 5.6% Stampi-Bombelli fronts stall). Both inlets therefore run with the same binary and the same settings, and
the Dirichlet copies also show how far the current binary reproduces the stored results. The source inputs are only
read. The BEA cases (3.4.1 dilute control, integrator benchmark) have no physical axial dispersion, so both inlets
are identical there; they are left out.

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
    "3.4.1_boundary_conditions_marx_10bar": (
        REPOSITORY / "validation" / "3.4.1_boundary_conditions_marx_10bar" / "simulations",
        [REPOSITORY / "validation" / "3.4.1_boundary_conditions_marx_10bar" / "simulations"],
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

INLET_CONDITIONS = {"dirichlet": "Dirichlet", "danckwerts": "Danckwerts"}
SIMULATION_TYPES = {"Breakthrough", "SwingAdsorption"}
FORMER_DEFAULT_LINEAR_SOLVER = "Dense"  # CVODELinearSolver default until 2026-11-02


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
        data.pop(keys.get("inletcondition", "InletCondition"), None)
        if data.get(keys.get("breakthroughintegrator", "BreakthroughIntegrator")) == "CVODE" and "cvodelinearsolver" not in keys:
            data["CVODELinearSolver"] = FORMER_DEFAULT_LINEAR_SOLVER
        for inlet, value in INLET_CONDITIONS.items():
            directory = SIMULATION_DIR / inlet / study / case
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "simulation.json").write_text(json.dumps({**data, "InletCondition": value}, indent=4) + "\n")
            index.append(
                {
                    "study": study,
                    "case": case.as_posix(),
                    "inlet": value,
                    "directory": directory.relative_to(ROOT).as_posix(),
                    "source": source.parent.relative_to(REPOSITORY).as_posix(),
                }
            )

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    cases = len(index) // len(INLET_CONDITIONS)
    print(f"wrote {cases} cases x {len(INLET_CONDITIONS)} inlet conditions to {SIMULATION_DIR.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
