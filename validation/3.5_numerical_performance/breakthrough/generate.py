#!/usr/bin/env python3
"""Generate the full-breakthrough integrator benchmark (RK3 against CVODE).

The 0.1 s matrices in ../figures time single integration steps; this study integrates complete
breakthroughs, so that cost and accuracy refer to the same, physically complete result.  Base cases:

  * bea-alkanes-C7: the BEA base case of the integrator benchmark, isothermal, seven C7 isomers at 1% in
    helium (../../../examples/BEA-alkanes-C7/breakthrough/simulation.json), 1200 s; convection-limited;
  * marx-T15: the validated Marx et al. (2014) experiment T15 with energy balance
    (../../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json), the same experiment as the benchmark case
    marx-403030_298_1e6_2e-1, 1500 s;
  * monolith-400ppm-fitted-dl: the Stampi-Bombelli et al. (2024) monolith at 400 ppm CO2 and 1 mmol/s with
    the authors' fitted axial dispersion coefficient, 2.6e-3 m2/s
    (../../3.3.2_monoliths_stampi_bombelli/stampi-bombelli/dl_from_stampibombelli/monolith_400ppm_1mmol/simulation.json), isothermal, 3000 s;
    diffusion-limited: v dz / D = 0.04, so the explicit stability limit 0.63 dz^2 / D = 1.7e-3 s lies far
    below dz / v = 0.062 s;
  * monolith-400ppm-fitted-dl-danckwerts: the same monolith with the Danckwerts inlet of the validation runs (the
    Dirichlet inlet of the other cases admits a dispersive inflow comparable to the feed here).

All three cases share the RK3 time steps 1e-4, 3.125e-4, 1e-3, 3.125e-3, and 1e-2 s (logspace(-4, -2, 5), with 3.125
instead of 3.16 so that every step divides the 2 s output interval and all runs sample the outlet at the same
times).  Steps above the stability limit of a case are run as well and fail (the monolith at 3.125e-3 and
1e-2 s).

Studies:

  * precision   RK3 over its stable time-step range, CVODE over a range of tolerances with the SPGMR and with
                the Dense linear solver, and ROS34PW2 over the same tolerances (added 2026-11-03), for the
                work-precision diagram; reference: CVODE (SPGMR) with a relative tolerance of 1e-8; and, for a fair
                comparison with ROS34PW2, CVODE without its step cap: unpreconditioned (cvode-nocap), with the
                analytical sorption Jacobian as SPGMR preconditioner (cvode-precond), and with it in the KLU direct
                solver (cvode-klu, BEA only: KLU rejects the energy balance and chemisorption);
  * components  BEA with 2, 4, and 6 adsorbing components (the component sweep of ../figures);
  * methods     Marx T15 with single-site Langmuir isotherms and IAST, SIAST, SEI, and SPI (the method
                matrix of ../figures; SEI requires single-site Langmuir isotherms);
  * stability   one RK3 run per base case at the smallest time step found unstable, just above the largest
                step of the precision study (BEA 1.5e-2 s, T15 0.15 s, monolith 2e-3 s).  BEA runs to the end
                with an unbounded solution; the other two stop with "negative outlet pressure".

CVODE uses the SPGMR linear solver, except for the Dense series of the precision study (finite-difference
Jacobian, the default of CVODELinearSolver).  All runs record the column every 2 s (WriteInterval), so every run
samples the outlet at the same times; the linear interpolation of t50 between samples taken at different
times differs by up to 3e-4 (relative) on the sharp H2 front of marx-T15.  The CVODE absolute tolerances are scaled with the relative tolerance from their defaults
(1e-6 mol/m3, 1e-7 mol/kg, 1e-4 K at a relative tolerance of 1e-5).  The CVODE step cap is left at its default,
4 dz/v, except in the cvode-nocap, cvode-precond, and cvode-klu series.

    python generate.py            # write simulations/ and manifest.csv
    python generate.py --clean    # delete simulations/ first

All quantities are SI.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import shutil
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[2]
SIMULATION_DIR = ROOT / "simulations"

# Helpers of the 0.1 s figure matrices, so that both use the same component subsets and Langmuir projection.
_figures_spec = importlib.util.spec_from_file_location("figures_generate", ROOT.parent / "figures" / "generate.py")
_figures = importlib.util.module_from_spec(_figures_spec)
_figures_spec.loader.exec_module(_figures)

RK3_TIME_STEPS = (1.0e-2, 3.125e-3, 1.0e-3, 3.125e-4, 1.0e-4)  # s, shared by all cases

CASES = {
    "bea-alkanes-C7": {
        "input": REPOSITORY / "examples" / "BEA-alkanes-C7" / "breakthrough" / "simulation.json",
        "horizon": 1200.0,  # s; every isomer is within 1% of its feed value after 1000 s
        "rk3_scaling_time_step": 1.0e-2,  # largest stable step in RK3_TIME_STEPS, a Courant number of 0.98
        "rk3_unstable_time_step": 1.5e-2,  # 1.25e-2 s is stable
        "klu": True,  # isothermal, no chemisorption: the only case CVODELinearSolver KLU accepts
    },
    "marx-T15": {
        "input": REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json",
        "horizon": 1500.0,  # s; complete CO2 breakthrough
        # The Marx inputs (2026-11-20): 0.125 s is stable, 0.15 s is not; before, 3e-2 s was stable and 5e-2 s was not.
        "rk3_scaling_time_step": 0.125,  # largest stable step that divides 2 s
        "rk3_unstable_time_step": 0.15,
    },
    "monolith-400ppm-fitted-dl": {
        "input": REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli" / "dl_from_stampibombelli" / "monolith_400ppm_1mmol" / "simulation.json",
        "horizon": 3000.0,  # s; t50 = 1284 s, 95% of the feed at 2860 s
        "rk3_unstable_time_step": 2.0e-3,  # 1.5e-3 s is stable
    },
    "monolith-400ppm-fitted-dl-danckwerts": {
        "input": REPOSITORY / "validation" / "3.3.2_monoliths_stampi_bombelli" / "stampi-bombelli" / "dl_from_stampibombelli" / "monolith_400ppm_1mmol" / "simulation.json",
        "horizon": 3000.0,  # s
        "rk3_unstable_time_step": 2.0e-3,  # the interior dispersion limit of the Dirichlet case; the inlet does not move it
        "inlet_condition": "Danckwerts",
    },
}
# CVODE without the 4 dz/v cap, which ROS34PW2 does not have, alone and with the analytical sorption Jacobian.
CVODE_FAIR_SERIES = {
    "cvode-nocap": {"CVODELinearSolver": "SPGMR", "CVODEMaximumTimeStep": 0.0},
    "cvode-precond": {"CVODELinearSolver": "SPGMR", "CVODEPreconditioner": "SorptionJacobian",
                      "CVODEMaximumTimeStep": 0.0},
    "cvode-klu": {"CVODELinearSolver": "KLU", "CVODEMaximumTimeStep": 0.0},
}
RELATIVE_TOLERANCES = (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6)
REFERENCE_RELATIVE_TOLERANCE = 1.0e-8
DEFAULT_RELATIVE_TOLERANCE = 1.0e-5
# The absolute tolerances of CVODE and ROS34PW2 at the default relative tolerance (keys without the integrator prefix).
DEFAULT_ABSOLUTE_TOLERANCES = {
    "AbsoluteToleranceConcentration": 1.0e-6,  # mol/m3
    "AbsoluteToleranceLoading": 1.0e-7,  # mol/kg
    "AbsoluteToleranceTemperature": 1.0e-4,  # K
}
CVODE_OUTPUT_INTERVAL = 2.0  # s; equal to WRITE_INTERVAL; listed as the CVODE time_step in the manifest
WRITE_INTERVAL = 2.0  # s
BEA_COMPONENT_COUNTS = (2, 4, 6)
METHODS = ("IAST", "SIAST", "SEI", "SPI")


def configure(simulation: dict, integrator: str, horizon: float, time_step: float, relative_tolerance=None,
              linear_solver: str = "SPGMR", inlet_condition: str = "Dirichlet") -> None:
    """Set the integrator and times; time_step is the RK3 step and is not used by the adaptive CVODE and ROS34PW2,
    which pick their own and get the relative tolerance with the absolute tolerances scaled from their defaults."""
    simulation["BreakthroughIntegrator"] = integrator
    simulation["SimulationTime"] = horizon
    simulation["WriteInterval"] = WRITE_INTERVAL
    # The batches of 2026-10-28 and 2026-11-03 ran with the Dirichlet inlet, the default until 2026-11-03.
    simulation["InletCondition"] = inlet_condition
    if integrator == "RungeKutta3":
        simulation["TimeStep"] = time_step
        return
    simulation.pop("TimeStep", None)
    if integrator == "CVODE":
        simulation["CVODELinearSolver"] = linear_solver
    if relative_tolerance is not None:
        factor = relative_tolerance / DEFAULT_RELATIVE_TOLERANCE
        simulation[f"{integrator}RelativeTolerance"] = relative_tolerance
        for key, value in DEFAULT_ABSOLUTE_TOLERANCES.items():
            simulation[f"{integrator}{key}"] = value * factor


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
        simulation["DisplayName"] = f"benchmark {relative}"
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=2) + "\n")
        rows.append({"directory": relative, **fields})

    for case, settings in CASES.items():
        base = json.loads(settings["input"].read_text())
        horizon = settings["horizon"]
        inlet = settings.get("inlet_condition", "Dirichlet")
        for time_step in RK3_TIME_STEPS:
            simulation = deepcopy(base)
            configure(simulation, "RungeKutta3", horizon, time_step, inlet_condition=inlet)
            write(f"precision/{case}/rk3/dt-{time_step:g}", simulation, study="precision", case=case,
                  integrator="rk3", time_step=time_step, relative_tolerance="", ncomp="", method="")
        for tolerance in (*RELATIVE_TOLERANCES, REFERENCE_RELATIVE_TOLERANCE):
            simulation = deepcopy(base)
            configure(simulation, "CVODE", horizon, CVODE_OUTPUT_INTERVAL, tolerance, inlet_condition=inlet)
            kind = "reference" if tolerance == REFERENCE_RELATIVE_TOLERANCE else "cvode"
            write(f"precision/{case}/{kind}/rtol-{tolerance:g}", simulation, study="precision", case=case,
                  integrator=kind, time_step=CVODE_OUTPUT_INTERVAL, relative_tolerance=tolerance, ncomp="", method="")
        for tolerance in RELATIVE_TOLERANCES:
            simulation = deepcopy(base)
            configure(simulation, "CVODE", horizon, CVODE_OUTPUT_INTERVAL, tolerance, linear_solver="Dense",
                      inlet_condition=inlet)
            write(f"precision/{case}/cvode-dense/rtol-{tolerance:g}", simulation, study="precision", case=case,
                  integrator="cvode-dense", time_step=CVODE_OUTPUT_INTERVAL, relative_tolerance=tolerance, ncomp="",
                  method="")
        for series, keys in CVODE_FAIR_SERIES.items():
            if series == "cvode-klu" and not settings.get("klu"):
                continue
            for tolerance in RELATIVE_TOLERANCES:
                simulation = deepcopy(base)
                configure(simulation, "CVODE", horizon, CVODE_OUTPUT_INTERVAL, tolerance, inlet_condition=inlet)
                simulation.update(keys)
                write(f"precision/{case}/{series}/rtol-{tolerance:g}", simulation, study="precision", case=case,
                      integrator=series, time_step=CVODE_OUTPUT_INTERVAL, relative_tolerance=tolerance, ncomp="",
                      method="")
        for tolerance in RELATIVE_TOLERANCES:
            simulation = deepcopy(base)
            configure(simulation, "ROS34PW2", horizon, CVODE_OUTPUT_INTERVAL, tolerance, inlet_condition=inlet)
            write(f"precision/{case}/ros34pw2/rtol-{tolerance:g}", simulation, study="precision", case=case,
                  integrator="ros34pw2", time_step=CVODE_OUTPUT_INTERVAL, relative_tolerance=tolerance, ncomp="",
                  method="")
        simulation = deepcopy(base)
        time_step = settings["rk3_unstable_time_step"]
        configure(simulation, "RungeKutta3", horizon, time_step, inlet_condition=inlet)
        write(f"stability/{case}/rk3/dt-{time_step:g}", simulation, study="stability", case=case,
              integrator="rk3", time_step=time_step, relative_tolerance="", ncomp="", method="")

    bea = CASES["bea-alkanes-C7"]
    for ncomp in BEA_COMPONENT_COUNTS:
        for integrator, time_step in (("RungeKutta3", bea["rk3_scaling_time_step"]), ("CVODE", CVODE_OUTPUT_INTERVAL)):
            simulation = json.loads(bea["input"].read_text())
            _figures.select_adsorbing_components(simulation, ncomp)
            configure(simulation, integrator, bea["horizon"], time_step)
            label = "rk3" if integrator == "RungeKutta3" else "cvode"
            write(f"components/bea-alkanes-C7/ncomp-{ncomp}/{label}", simulation, study="components",
                  case="bea-alkanes-C7", integrator=label, time_step=time_step, relative_tolerance="",
                  ncomp=ncomp, method="IAST")

    marx = CASES["marx-T15"]
    for method in METHODS:
        for integrator, time_step in (("RungeKutta3", marx["rk3_scaling_time_step"]), ("CVODE", CVODE_OUTPUT_INTERVAL)):
            simulation = json.loads(marx["input"].read_text())
            _figures.use_single_langmuir_sites(simulation)
            simulation["MixturePredictionMethod"] = method
            configure(simulation, integrator, marx["horizon"], time_step)
            label = "rk3" if integrator == "RungeKutta3" else "cvode"
            write(f"methods/marx-T15-langmuir/{method}/{label}", simulation, study="methods", case="marx-T15",
                  integrator=label, time_step=time_step, relative_tolerance="", ncomp=3, method=method)

    with (ROOT / "manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} simulations below {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
