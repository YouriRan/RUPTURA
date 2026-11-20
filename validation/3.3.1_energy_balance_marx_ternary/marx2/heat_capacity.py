#!/usr/bin/env python3
"""Set the gas heat capacity of the Marx inputs per unit mass of the feed.

RUPTURA's gas energy balance uses the volumetric heat capacity rho_g C_p,g, with rho_g = sum_j M_j c_j the
mass density, so `heatCapacityGas` is a specific heat capacity in J/(kg K) (SI, input table).  Marx et al.
(2014), Table 4, give the molar heat capacity of the gas, 29.51 J/(mol K), which the inputs used directly;
the gas volumetric heat capacity was then smaller than c_tot C_p by the factor 1/M_mean (about 38 for the
40/30/30 feed).  The specific value of the feed is

    C_p,g = C_p,g^mol / M_mean,    M_mean = sum_j y_j,in M_j,

with the feed mole fractions and molecular weights of each input (see conversions.md).

    python heat_capacity.py            # print the table and write heatCapacityGas into base.json, reference_settings_hassan_chung/T*, binary/B*
    python heat_capacity.py --dry-run  # print the table only
    python heat_capacity.py ../../../benchmark/marx-403030_298_1e6_2e-1  # any simulation.json, or directories holding them

Only the value of heatCapacityGas is changed; the rest of each file keeps its formatting.  All quantities are SI.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent

MOLAR_HEAT_CAPACITY = 29.51  # J/(mol K), Marx et al. (2014), Table 4

DEFAULT_INPUTS = [
    ROOT / "base.json",
    *(ROOT / "reference_settings_hassan_chung" / f"T{index}" / "simulation.json" for index in range(1, 25)),
    *(ROOT / "binary" / f"B{index}" / "simulation.json" for index in range(1, 37)),
]


def mean_molar_mass(simulation: dict) -> float:
    """Mean molar mass [kg/mol] of the feed."""
    components = [component for component in simulation["Components"] if component["GasPhaseMolFraction"] > 0.0]
    total = sum(component["GasPhaseMolFraction"] for component in components)
    return sum(component["GasPhaseMolFraction"] * component["MolecularWeight"] for component in components) / total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="print the table without changing the inputs")
    parser.add_argument("inputs", nargs="*", type=Path,
                        help="simulation.json files or directories searched for them (default: base.json, reference_settings_hassan_chung/T1-T24, B1-B36)")
    arguments = parser.parse_args()

    paths = [path for given in arguments.inputs
             for path in ([given] if given.is_file() else sorted(given.rglob("simulation.json")))] or DEFAULT_INPUTS

    print(f"{'M_mean/(kg/mol)':>16s} {'C_p,g/(J/kg/K)':>15s} {'previous':>9s}  input")
    for path in paths:
        text = path.read_text(encoding="utf-8")
        simulation = json.loads(text)
        molar_mass = mean_molar_mass(simulation)
        heat_capacity = round(MOLAR_HEAT_CAPACITY / molar_mass, 1)
        previous = simulation["heatCapacityGas"]
        print(f"{molar_mass:16.6f} {heat_capacity:15.1f} {previous:9.2f}  {path}")
        if not arguments.dry_run:
            text = re.sub(r'("heatCapacityGas":\s*)[-+0-9.eE]+', rf"\g<1>{heat_capacity}", text, count=1)
            path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
