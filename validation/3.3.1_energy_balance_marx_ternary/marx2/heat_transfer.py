#!/usr/bin/env python3
"""Set the gas-wall heat-transfer coefficient of the ternary runs T1-T24 from the Marx et al. correlation.

Marx et al. (2014), Eq. 3, with eta_1 = 41.13 and eta_2 = 0.32 from Casas et al. (2012):

    Nu = h_L d_i / K_L = eta_1 Re_0^eta_2 exp(-6 d_p / d_i),    Re_0 = rho_0 u d_p / mu.

Re_0 is the initial Reynolds number: rho_0 is the density of the helium that fills the column before the
feed is switched on, at the pressure and temperature of the experiment, and u is the superficial feed
velocity (Casas, Diss. ETH 20681, 2012, Sect. 6.4.2).  Marx et al. use one representative K_L for all
experiments but do not report it; K_L = 0.0400 W/(m K) reproduces the correlation curve of their Fig. 6,
h_L = 32.0 Re_0^0.32 W/(m2 K) (digitized; see conversions.md).

    python heat_transfer.py            # print the table and write heatTransferGasWall into ternary/T*/simulation.json
    python heat_transfer.py --dry-run  # print the table only

Only the key heatTransferGasWall is changed.  All quantities are SI.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1, the value in src/utils.h
HELIUM_MOLAR_MASS = 0.004  # kg/mol, as in the inputs

# Marx et al. (2014), Table 4 and Eq. 3.
ETA_1 = 41.13
ETA_2 = 0.32
VISCOSITY = 1.9e-5  # Pa s
TOTAL_POROSITY = 0.742096  # 0.403 + (1 - 0.403) * 0.568; the inputs use v = Q / (A eps_t)
AXIAL_THERMAL_CONDUCTIVITY = 0.0400  # W/(m K), from the correlation curve of Marx et al., Fig. 6


def heat_transfer_coefficient(simulation: dict) -> tuple[float, float]:
    """Initial Reynolds number [-] and bed-wall heat-transfer coefficient [W/(m2 K)] of one run."""
    geometry = simulation["Geometry"]
    inner_diameter = geometry["InternalDiameter"]
    particle_diameter = geometry["ParticleDiameter"]
    superficial_velocity = simulation["ColumnEntranceVelocity"] * TOTAL_POROSITY
    helium_density = (
        simulation["InletPressure"] * HELIUM_MOLAR_MASS / (GAS_CONSTANT * simulation["InfluxTemperature"])
    )
    reynolds = helium_density * superficial_velocity * particle_diameter / VISCOSITY
    nusselt = ETA_1 * reynolds**ETA_2 * math.exp(-6.0 * particle_diameter / inner_diameter)
    return reynolds, nusselt * AXIAL_THERMAL_CONDUCTIVITY / inner_diameter


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="print the table without changing the inputs")
    arguments = parser.parse_args()

    print(f"{'run':4s} {'P_in/Pa':>10s} {'T_in/K':>7s} {'u/(m/s)':>9s} {'Re_0':>7s} {'h_L/(W/m2/K)':>13s} {'previous':>9s}")
    for index in range(1, 25):
        path = ROOT / "ternary" / f"T{index}" / "simulation.json"
        simulation = json.loads(path.read_text())
        reynolds, coefficient = heat_transfer_coefficient(simulation)
        previous = simulation["heatTransferGasWall"]
        print(
            f"T{index:<3d} {simulation['InletPressure']:10.0f} {simulation['InfluxTemperature']:7.0f} "
            f"{simulation['ColumnEntranceVelocity'] * TOTAL_POROSITY:9.5f} {reynolds:7.3f} {coefficient:13.2f} {previous:9.2f}"
        )
        if not arguments.dry_run:
            simulation["heatTransferGasWall"] = round(coefficient, 2)
            path.write_text(json.dumps(simulation, indent=4))


if __name__ == "__main__":
    main()
