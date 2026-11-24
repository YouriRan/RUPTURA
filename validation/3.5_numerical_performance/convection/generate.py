#!/usr/bin/env python3
"""Generate the convection-scheme study: first-order upwind against WENO3 (ConvectionScheme) on the three
breakthrough cases of ../breakthrough.

Base cases (as in ../breakthrough/generate.py):

  * bea-alkanes-C7: seven C7 isomers at 1% in helium on BEA, isothermal, 1200 s, no axial dispersion
    (../../../examples/BEA-alkanes-C7/breakthrough/simulation.json); convection-limited;
  * marx-T15: experiment T15 of Marx et al. (2014) with energy balance, 1500 s (../../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15);
  * monolith-400ppm-fitted-dl: the Stampi-Bombelli et al. (2024) monolith at 400 ppm CO2 with the fitted axial
    dispersion coefficient, 3000 s, N = 50 (../../3.3.2_monoliths_stampi_bombelli/stampi-bombelli/dl_from_stampibombelli/monolith_400ppm_1mmol);
    dispersion-limited.

Studies (scheme = Upwind or WENO3; the grid study also TVD and WENO5):

  * cost       input grid; RK3 at the step of the integrator benchmark (BEA 1e-2 s, T15 0.125 s, monolith 1e-3 s)
               and CVODE (SPGMR) at the default tolerances; two repeats each;
  * stability  RK3 at steps around the stability limits (BEA 1.25e-2, 1.5e-2, 1.75e-2 s, Courant numbers 1.23,
               1.48, 1.72; T15 0.125 and 0.15 s; monolith 1.5e-3 and 2e-3 s);
  * grid       CVODE (SPGMR) at a relative tolerance of 1e-6 with N = 25, 50, 100, 200, for Upwind, TVD, WENO3, and
               WENO5;
  * reference  CVODE, WENO3, N = 400, relative tolerance 1e-6.

All runs write the column every 2 s (WriteInterval).  The CVODE absolute tolerances scale with the relative
tolerance from their defaults (1e-5 mol/m3, 1e-6 mol/kg, 1e-3 K at 1e-4), as in ../breakthrough.

    python generate.py            # write simulations/ and manifest.csv
    python generate.py --clean    # delete simulations/ first

All quantities are SI.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
SIMULATION_DIR = ROOT / "simulations"

SCHEMES = ("Upwind", "WENO3")
GRID_SCHEMES = ("Upwind", "TVD", "WENO3", "WENO5")  # the grid study also runs TVD and WENO5
CASES = {
    "bea-alkanes-C7": {
        "input": REPOSITORY / "examples" / "BEA-alkanes-C7" / "breakthrough" / "simulation.json",
        "horizon": 1200.0,
        "rk3_time_step": 1.0e-2,
        "rk3_stability_time_steps": (1.25e-2, 1.5e-2, 1.75e-2),
    },
    "marx-T15": {
        "input": REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json",
        "horizon": 1500.0,
        "rk3_time_step": 0.125,  # the Marx inputs (2026-11-20); 2.5e-2 s before
        "rk3_stability_time_steps": (0.125, 0.15),  # (3e-2, 5e-2) s before
    },
    "monolith-400ppm-fitted-dl": {
        "input": REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli" / "dl_from_stampibombelli" / "monolith_400ppm_1mmol" / "simulation.json",
        "horizon": 3000.0,
        "rk3_time_step": 1.0e-3,
        "rk3_stability_time_steps": (1.5e-3, 2.0e-3),
    },
}
REPEATS = 2
GRID_POINTS = (25, 50, 100, 200)
REFERENCE_GRID_POINTS = 400
GRID_RELATIVE_TOLERANCE = 1.0e-6
DEFAULT_RELATIVE_TOLERANCE = 1.0e-4
DEFAULT_ABSOLUTE_TOLERANCES = {
    "CVODEAbsoluteToleranceConcentration": 1.0e-5,  # mol/m3
    "CVODEAbsoluteToleranceLoading": 1.0e-6,  # mol/kg
    "CVODEAbsoluteToleranceTemperature": 1.0e-3,  # K
}
WRITE_INTERVAL = 2.0  # s
# Time keys of the base inputs (current and legacy), replaced by configure().
TIME_KEYS = ("TimeStep", "SimulationTime", "WriteInterval", "NumberOfTimeSteps", "WriteEvery", "PrintEvery",
             "NumberOfInitTimeSteps")


def configure(simulation: dict, scheme: str, integrator: str, horizon: float, time_step: float | None = None,
              relative_tolerance: float | None = None, grid_points: int | None = None) -> None:
    for key in TIME_KEYS:
        simulation.pop(key, None)
    simulation["ConvectionScheme"] = scheme
    simulation["BreakthroughIntegrator"] = integrator
    simulation["SimulationTime"] = horizon
    simulation["WriteInterval"] = WRITE_INTERVAL
    simulation["WriteOutletOnly"] = True  # the analysis reads the outlet only (since 2026-11-20)
    # The runs of 2026-10-31 used the Dirichlet inlet, the default until 2026-11-03, as ../breakthrough does; the base
    # inputs carry Danckwerts since then.
    simulation["InletCondition"] = "Dirichlet"
    if integrator == "RungeKutta3":
        simulation["TimeStep"] = time_step
    if grid_points is not None:
        simulation["NumberOfGridPoints"] = grid_points
    if integrator == "CVODE":
        simulation["CVODELinearSolver"] = "SPGMR"
        factor = relative_tolerance / DEFAULT_RELATIVE_TOLERANCE
        simulation["CVODERelativeTolerance"] = relative_tolerance
        for key, value in DEFAULT_ABSOLUTE_TOLERANCES.items():
            simulation[key] = value * factor


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true", help="delete simulations/ first, including any output")
    arguments = parser.parse_args()
    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)

    rows = []

    def write(relative: str, simulation: dict, **fields) -> None:
        directory = SIMULATION_DIR / relative
        directory.mkdir(parents=True, exist_ok=True)
        simulation["DisplayName"] = f"convection {relative}"
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=2) + "\n")
        rows.append({"directory": relative, "grid_points": simulation.get("NumberOfGridPoints", 100), **fields})

    for case, settings in CASES.items():
        base = json.loads(settings["input"].read_text())
        horizon = settings["horizon"]
        for scheme in SCHEMES:
            label = scheme.lower()
            for repeat in range(1, REPEATS + 1):
                simulation = deepcopy(base)
                configure(simulation, scheme, "RungeKutta3", horizon, time_step=settings["rk3_time_step"])
                write(f"cost/{case}/{label}/rk3_dt-{settings['rk3_time_step']:g}_repeat-{repeat}", simulation,
                      study="cost", case=case, scheme=scheme, integrator="rk3", time_step=settings["rk3_time_step"],
                      relative_tolerance="", repeat=repeat)
                simulation = deepcopy(base)
                configure(simulation, scheme, "CVODE", horizon, relative_tolerance=DEFAULT_RELATIVE_TOLERANCE)
                write(f"cost/{case}/{label}/cvode_rtol-1e-04_repeat-{repeat}", simulation, study="cost", case=case,
                      scheme=scheme, integrator="cvode", time_step="", relative_tolerance=DEFAULT_RELATIVE_TOLERANCE,
                      repeat=repeat)
            for time_step in settings["rk3_stability_time_steps"]:
                simulation = deepcopy(base)
                configure(simulation, scheme, "RungeKutta3", horizon, time_step=time_step)
                write(f"stability/{case}/{label}/rk3_dt-{time_step:g}", simulation, study="stability", case=case,
                      scheme=scheme, integrator="rk3", time_step=time_step, relative_tolerance="", repeat=1)
        for scheme in GRID_SCHEMES:
            label = scheme.lower()
            for grid_points in GRID_POINTS:
                simulation = deepcopy(base)
                configure(simulation, scheme, "CVODE", horizon, relative_tolerance=GRID_RELATIVE_TOLERANCE,
                          grid_points=grid_points)
                write(f"grid/{case}/{label}/cvode_n-{grid_points}", simulation, study="grid", case=case,
                      scheme=scheme, integrator="cvode", time_step="", relative_tolerance=GRID_RELATIVE_TOLERANCE,
                      repeat=1)
        simulation = deepcopy(base)
        configure(simulation, "WENO3", "CVODE", horizon, relative_tolerance=GRID_RELATIVE_TOLERANCE,
                  grid_points=REFERENCE_GRID_POINTS)
        write(f"reference/{case}/weno3/cvode_n-{REFERENCE_GRID_POINTS}", simulation, study="reference", case=case,
              scheme="WENO3", integrator="cvode", time_step="", relative_tolerance=GRID_RELATIVE_TOLERANCE, repeat=1)

    with (ROOT / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} runs written to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
