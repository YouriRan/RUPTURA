#!/usr/bin/env python3
"""Materialize the simulation inputs used by presentation.ipynb.

mixture/iast, mixture/mpd     mixture isotherms, dual-site Langmuir + IAST versus the multinomial MPD
variance/mpd_<name>_<sign>    MPD mixture isotherms with y_<name> shifted by a relative 1e-3 against the
                              carrier gas, so only f_<name> changes; used for d<N_i>/d ln f_i
breakthrough/iast, breakthrough/mpd   the standard settings of examples/BEA-alkanes-C7/breakthrough, with the
                              matrix-free SPGMR linear solver instead of CVODE's dense finite-difference Jacobian
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import json


OUTPUT_DIRECTORY = Path(__file__).resolve().parent
EXAMPLE_DIRECTORY = OUTPUT_DIRECTORY.parent
STANDARD_BREAKTHROUGH_INPUT = EXAMPLE_DIRECTORY.parents[1] / "BEA-alkanes-C7" / "breakthrough" / "simulation.json"

# Every simulation directory sits two levels below this one.
MPD_FILE_NAME = "../../../mpd/particle_distribution.data"
RELATIVE_PERTURBATION = 1.0e-3

LANGMUIR_MIXTURE = json.loads((EXAMPLE_DIRECTORY / "langmuir" / "simulation.json").read_text(encoding="utf-8"))
MPD_MIXTURE = json.loads((EXAMPLE_DIRECTORY / "mpd" / "simulation.json").read_text(encoding="utf-8"))
STANDARD_BREAKTHROUGH = json.loads(STANDARD_BREAKTHROUGH_INPUT.read_text(encoding="utf-8"))
# Numerics only: the dense Jacobian costs one right-hand side per state variable (~1700 here) per setup.
STANDARD_BREAKTHROUGH["CVODELinearSolver"] = "SPGMR"

SITES = {component["Name"]: component["PhysisorptionSites"]
         for component in LANGMUIR_MIXTURE["Components"] if not component.get("CarrierGas", False)}


def main() -> None:
    simulations = {}

    iast_mixture = deepcopy(LANGMUIR_MIXTURE)
    iast_mixture["MixturePredictionMethod"] = "IAST"
    iast_mixture["DisplayName"] = "BEA C7 alkanes — dual-site Langmuir + IAST"
    simulations["mixture/iast"] = iast_mixture

    mpd_mixture = deepcopy(MPD_MIXTURE)
    mpd_mixture["MPDSettings"]["FileName"] = MPD_FILE_NAME
    simulations["mixture/mpd"] = mpd_mixture

    # The carrier gas does not adsorb, so moving the shift onto it keeps every other f_j = y_j p fixed.
    for name in SITES:
        for sign, label in ((-1.0, "minus"), (1.0, "plus")):
            perturbed = deepcopy(mpd_mixture)
            target = next(component for component in perturbed["Components"] if component["Name"] == name)
            carrier = next(component for component in perturbed["Components"] if component.get("CarrierGas", False))
            shift = sign * RELATIVE_PERTURBATION * target["GasPhaseMolFraction"]
            target["GasPhaseMolFraction"] = round(target["GasPhaseMolFraction"] + shift, 12)
            carrier["GasPhaseMolFraction"] = round(carrier["GasPhaseMolFraction"] - shift, 12)
            perturbed["DisplayName"] = f"BEA C7 alkanes — multinomial MPD, y_{name} {label} {RELATIVE_PERTURBATION:g}"
            simulations[f"variance/mpd_{name}_{label}"] = perturbed

    iast_breakthrough = deepcopy(STANDARD_BREAKTHROUGH)
    iast_breakthrough["MixturePredictionMethod"] = "IAST"
    iast_breakthrough["DisplayName"] = "BEA C7 dual-site Langmuir + IAST breakthrough"
    for component in iast_breakthrough["Components"]:
        if component["Name"] in SITES:
            component["PhysisorptionSites"] = deepcopy(SITES[component["Name"]])
    simulations["breakthrough/iast"] = iast_breakthrough

    mpd_breakthrough = deepcopy(STANDARD_BREAKTHROUGH)
    mpd_breakthrough["MixturePredictionMethod"] = "MPD"
    mpd_breakthrough["DisplayName"] = "BEA C7 multinomial MPD breakthrough"
    mpd_breakthrough["MPDSettings"] = deepcopy(mpd_mixture["MPDSettings"])
    for component in mpd_breakthrough["Components"]:
        component.pop("PhysisorptionSites", None)
    simulations["breakthrough/mpd"] = mpd_breakthrough

    for name, simulation in simulations.items():
        directory = OUTPUT_DIRECTORY / name
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(simulations)} simulations beneath {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
