#!/usr/bin/env python3
"""Generate the temperature-swing cycle study on the Marx et al. (2014) system.

The column, adsorbent, isotherms, and heat-transfer parameters are those of the
validated ternary experiment T15 (CO2/N2/H2 = 40/30/30, 1e6 Pa, 1e-5 m3/s; paper
Section 3.1).  Each run repeats a two-step temperature-swing cycle with the
SwingAdsorption phase driver:

    Adsorption     feed and surroundings at 298 K, duration t_half
    Regeneration   feed and surroundings at 373 K, duration t_half

The feed composition, pressure, and velocity are the same in both steps, because
the phase driver changes only the temperature and the inlet pressure.  With the
energy balance on, a phase temperature sets the inlet-gas and external
temperatures; the column temperatures follow from the energy balances.  The first
adsorption step starts from the helium-filled column of the experiment.  Every
run is one half-cycle duration, repeated until cyclic steady state.

Two MixturePrediction runs give the IAST CO2 loading of the feed at 298 K and
373 K, i.e. the equilibrium working capacity that bounds the cyclic one.

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI: Pa, K, s, m, mol/kg.
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
BASE_INPUT = REPOSITORY / "validation" / "3.3.1_energy_balance_marx_ternary" / "marx2" / "ternary" / "T15" / "simulation.json"
SIMULATION_DIR = ROOT / "simulations"

ADSORPTION_TEMPERATURE = 298.0
REGENERATION_TEMPERATURE = 373.0
HALF_CYCLE_DURATIONS = (1000, 2000, 4000, 8000, 16000)  # s

# The wall dominates the thermal inertia of this column (relaxation time about 5000 s), so the mean column
# temperature needs about 60000 s to settle.  Every run covers at least that time and at least five cycles.
MINIMUM_SIMULATED_TIME = 60_000.0  # s
MINIMUM_NUMBER_OF_CYCLES = 5

# ROS34PW2 (from the T15 input) adapts its internal step and ends a step at every output.
WRITE_INTERVAL = 20.0  # s, one output block every 20 s


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
    index = []

    for half_cycle in HALF_CYCLE_DURATIONS:
        number_of_cycles = max(MINIMUM_NUMBER_OF_CYCLES, math.ceil(MINIMUM_SIMULATED_TIME / (2.0 * half_cycle)))
        phases = [
            {"Name": "Adsorption", "Temperature": ADSORPTION_TEMPERATURE, "Duration": float(half_cycle)},
            {"Name": "Regeneration", "Temperature": REGENERATION_TEMPERATURE, "Duration": float(half_cycle)},
        ] * number_of_cycles

        simulation = copy.deepcopy(base)
        simulation["SimulationType"] = "SwingAdsorption"
        simulation["DisplayName"] = f"Marx T15 TSA, half cycle {half_cycle} s"
        # The phase durations set the simulated time of a swing run.
        simulation.pop("SimulationTime", None)
        simulation["WriteInterval"] = WRITE_INTERVAL
        simulation["SwingAdsorptionPhases"] = phases

        directory = SIMULATION_DIR / f"T15_TSA_{ADSORPTION_TEMPERATURE:.0f}K-{REGENERATION_TEMPERATURE:.0f}K_half-cycle-{half_cycle}s"
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
        index.append(
            {
                "directory": str(directory.relative_to(ROOT)),
                "kind": "cycle",
                "half_cycle_s": half_cycle,
                "number_of_cycles": number_of_cycles,
                "adsorption_temperature_K": ADSORPTION_TEMPERATURE,
                "regeneration_temperature_K": REGENERATION_TEMPERATURE,
                "output_interval_s": WRITE_INTERVAL,
            }
        )

    # Equilibrium loading of the feed (adsorbing components only, fractions normalized) at the two step temperatures,
    # with the mixture method of the input.
    adsorbates = [copy.deepcopy(component) for component in base["Components"] if not component.get("CarrierGas", False)]
    fraction_sum = sum(component["GasPhaseMolFraction"] for component in adsorbates)
    for component in adsorbates:
        component["GasPhaseMolFraction"] /= fraction_sum
    for temperature in (ADSORPTION_TEMPERATURE, REGENERATION_TEMPERATURE):
        simulation = {
            "SimulationType": "MixturePrediction",
            "DisplayName": f"Marx T15 feed, {base['MixturePredictionMethod']} at {temperature:.0f} K",
            "Temperature": temperature,
            "MixturePredictionMethod": base["MixturePredictionMethod"],
            "PressureScale": "Linear",
            "PressureStart": base["InletPressure"],
            "PressureEnd": base["InletPressure"],
            "NumberOfPressurePoints": 1,
            "ColumnVoidFraction": base["ColumnVoidFraction"],
            "ParticleDensity": base["ParticleDensity"],
            "Components": adsorbates,
        }
        directory = SIMULATION_DIR / f"T15_feed_equilibrium_{temperature:.0f}K"
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=4) + "\n")
        index.append({"directory": str(directory.relative_to(ROOT)), "kind": "equilibrium", "temperature_K": temperature})

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
