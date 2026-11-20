#!/usr/bin/env python3
"""Gas pressure against the Ergun pressure, and its effect on the CO2 inventory, for the TSA runs (SI remark).

The `InletPressureInletVelocity` closure carries the thermal expansion of the gas, so the gas pressure c_tot R T_g
departs from the Ergun pressure P that the equilibrium uses only by the integration error.  For every run this prints the
range of c_tot R T_g / P - 1, the gas share of the CO2 in the column, and the relative change of the CO2 inventory
at the step ends if the gas is put at the Ergun pressure.

    python gas_inventory.py

All quantities are SI.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1, the value in src/utils.h
HALF_CYCLES = (1000, 2000, 4000, 8000, 16000)  # s


def read_columns(path: Path, headings: list[str], rows: int) -> np.ndarray:
    """Read the columns whose '# column N:' headings start with the given texts."""
    numbers = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("# column "):
                number, description = line[9:].split(":", 1)
                for heading in headings:
                    if description.strip().casefold().startswith(heading.casefold()):
                        numbers.setdefault(heading, int(number) - 1)
    columns = [numbers[heading] for heading in headings]
    table = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=columns).to_numpy()
    if rows:
        table = table[:rows]
    # pandas returns the columns in file order; put them in the order of the headings.
    return table[:, np.argsort(np.argsort(columns))]


def main() -> None:
    for half_cycle in HALF_CYCLES:
        directory = ROOT / "simulations" / f"T15_TSA_298K-373K_half-cycle-{half_cycle}s"
        simulation = json.loads((directory / "simulation.json").read_text())
        nodes = simulation["NumberOfGridPoints"] + 1
        void_fraction = simulation["ColumnVoidFraction"]
        density = simulation["ParticleDensity"]

        co2 = read_columns(directory / "component_1_CO2.data", ["Time, t", "Concentration, c_i", "Physisorption, q_phy"], 0)
        outputs = len(co2) // nodes
        co2 = co2[: outputs * nodes].reshape(outputs, nodes, 3)
        time = 60.0 * co2[:, 0, 0]
        total_concentration = np.zeros((outputs, nodes))
        for index, component in enumerate(simulation["Components"]):
            path = directory / f"component_{index}_{component['Name']}.data"
            total_concentration += read_columns(path, ["Concentration, c_i"], outputs * nodes).reshape(outputs, nodes)
        column = read_columns(directory / "column.data", ["Total pressure", "Gas temperature,"], outputs * nodes)
        column = column.reshape(outputs, nodes, 2)

        # Nodes 1..N; node 0 is the inlet node, held at the feed state.
        deviation = total_concentration[:, 1:] * GAS_CONSTANT * column[:, 1:, 1] / column[:, 1:, 0] - 1.0
        gas = void_fraction * co2[:, 1:, 1]
        solid = (1.0 - void_fraction) * density * co2[:, 1:, 2]
        share = gas.sum(axis=1) / (gas + solid).sum(axis=1)
        step_ends = [int(np.argmin(np.abs(time - k * half_cycle))) for k in range(1, int(round(time[-1] / half_cycle)) + 1)]
        inventory = (gas + solid).sum(axis=1)[step_ends]
        corrected = (gas / (1.0 + deviation) + solid).sum(axis=1)[step_ends]
        change = corrected / inventory - 1.0
        print(f"t_half {half_cycle:5d} s: c_tot R T_g / P - 1 from {deviation.min():+.1e} to {deviation.max():+.1e}; "
              f"gas share {share[1:].min():.3f}-{share[1:].max():.3f}; "
              f"inventory change at the step ends {change.min():+.1e} to {change.max():+.1e}")


if __name__ == "__main__":
    main()
