#!/usr/bin/env python3
"""Write the inputs of the ternary experiments T1-T24 of Marx et al. (2014) with the model inputs of that paper.

Marx, Joss, Casas, Schell, Mazzotti, Adsorption 20 (2014) 493-510 (s10450-013-9593-5-1.pdf in this directory):

  * runs: Table 2 (feed, temperature, pressure, volumetric flow at column conditions);
  * equilibrium: the temperature-dependent Sips isotherm of Schell et al. (2011), Table 5, combined as their
    multicomponent Sips equation (Sect. 4.2.1), i.e. RUPTURA's SCI:
        q_sat = a exp(-b/RT),  K = A exp(-B/RT),  s = alpha atan[beta (T - T_ref)] + s_ref,
    written as uncentred VantHoff laws (Energy -b and -B) and, for the exponent 1/nu = s, the Arctan law;
  * heats of adsorption Delta H, mass-transfer coefficients k_i, and h_w: Tables 5 and 6;
  * column (Table 4): their mass and energy balances hold the gas in the total porosity eps_t with the bulk density
    rho_b, so the fluid fraction is eps_t and the particle density rho_b/(1 - eps_t); the adsorbent per unit volume,
    the solid and gas heat capacities, and the bed-to-wall coupling are then theirs (conversions.md, section 6);
  * gas heat capacity: their constant molar C_g, 29.51 J/(mol K), for every component (GasHeatCapacity);
  * bed-wall coefficient h_L: their Eq. 3 (heat_transfer.py);
  * ambient: the experimental temperature; the column is heated by an electrical heater controlled at the outer wall.

Marx et al. assume gas and adsorbent at one temperature.  RUPTURA resolves both; the gas-solid coefficient follows
the correlation of Wakao, Kaguei and Funazkri, Chem. Eng. Sci. 34 (1979) 325, Nu = 2 + 1.1 Pr^(1/3) Re^0.6, with
the feed properties of each run.  Not in RUPTURA: the heat capacity of the adsorbed phase in the bed energy balance.

    python generate.py            # write ternary/T1-T24/simulation.json (outputs of earlier runs are kept)
    python generate.py --dry-run  # print the run table only

All quantities are SI.
"""

from __future__ import annotations

import argparse
import json
import math
from copy import deepcopy
from pathlib import Path

from heat_transfer import GAS_CONSTANT, TOTAL_POROSITY, heat_transfer_coefficient

ROOT = Path(__file__).resolve().parent

BED_POROSITY = 0.403  # Table 4
BULK_DENSITY = 502.0  # kg/m3, Table 4
MOLAR_GAS_HEAT_CAPACITY = 29.51  # J/(mol K), C_g of Marx et al.
GAS_THERMAL_CONDUCTIVITY = 0.028  # W/(m K), molecular, as gasThermalConductivity of the inputs
WALL_EXTERNAL = 5.0  # W/(m2 K), h_w of Table 6

# Table 5 (Schell et al. 2011): a in mol/kg, b and B in J/mol, A in 1/Pa, Delta H in J/mol.
ISOTHERMS = {
    "CO2": {"a": 1.38, "b": -5628.0, "A": 1.68e-2 * 1.0e-6, "B": -9159.0, "s_ref": 0.83,
            "alpha": 0.072, "beta": 0.106, "T_ref": 329.0, "heat": 26000.0},
    "N2": {"a": 2.82, "b": -1706.0, "A": 1.74e-3 * 1.0e-6, "B": -12661.0, "s_ref": 0.86, "heat": 15600.0},
    "H2": {"a": 6.66, "b": 0.0, "A": 6.97e-4 * 1.0e-6, "B": -9826.0, "s_ref": 0.96, "heat": 9800.0},
}
# Table 6: k_CO2 by feed and temperature; k_N2 = 0.3 and k_H2 = 1.0 1/s throughout.
MASS_TRANSFER_CO2 = {("15/80/5", 298.0): 0.04, ("15/80/5", 318.0): 0.07,
                     ("40/30/30", 298.0): 0.11, ("40/30/30", 318.0): 0.15}
MASS_TRANSFER = {"N2": 0.3, "H2": 1.0}
FEEDS = {"15/80/5": {"CO2": 0.15, "N2": 0.80, "H2": 0.05}, "40/30/30": {"CO2": 0.40, "N2": 0.30, "H2": 0.30}}
SIMULATION_TIME = {1.0e-5: 4000.0, 2.0e-5: 2000.0, 3.0e-5: 1400.0}  # s, by flow rate
# Table 2: per feed, temperature, the (pressure / Pa, flow / m3 s-1) of the six runs, numbered T1-T24 in this order.
CONDITIONS = [(1.0e5, 1.0e-5), (5.0e5, 1.0e-5), (1.0e6, 1.0e-5), (1.0e6, 2.0e-5), (1.0e6, 3.0e-5), (2.0e6, 1.0e-5)]
RUNS = {
    f"T{6 * block + index + 1}": (feed, temperature, *conditions)
    for block, (feed, temperature) in enumerate(
        [("15/80/5", 298.0), ("15/80/5", 318.0), ("40/30/30", 298.0), ("40/30/30", 318.0)])
    for index, conditions in enumerate(CONDITIONS)
}


def sips_site(name: str) -> dict:
    """RUPTURA Sips site q_sat (k p)^(1/nu) / (1 + (k p)^(1/nu)) with the laws of Table 5."""
    p = ISOTHERMS[name]
    saturation = {"Value": p["a"], "Law": "VantHoff", "Energy": -p["b"]} if p["b"] else p["a"]
    affinity = {"Value": p["A"], "Law": "VantHoff", "Energy": -p["B"]}
    exponent = (
        {"Value": 1.0 / p["s_ref"], "Law": "Arctan", "Alpha": p["alpha"], "Beta": p["beta"],
         "ReferenceTemperature": p["T_ref"]}
        if "alpha" in p else 1.0 / p["s_ref"]
    )
    return {"Type": "Sips", "Parameters": [saturation, affinity, exponent]}


def gas_solid_heat_transfer(pressure: float, temperature: float, superficial_velocity: float, feed: dict,
                            molar_mass: dict, viscosity: float, particle_diameter: float) -> float:
    """Wakao-Kaguei-Funazkri coefficient [W/(m2 K)] per unit external particle area, with the feed properties."""
    mean_molar_mass = sum(fraction * molar_mass[name] for name, fraction in feed.items())
    density = pressure * mean_molar_mass / (GAS_CONSTANT * temperature)
    reynolds = density * superficial_velocity * particle_diameter / viscosity
    prandtl = MOLAR_GAS_HEAT_CAPACITY / mean_molar_mass * viscosity / GAS_THERMAL_CONDUCTIVITY
    nusselt = 2.0 + 1.1 * prandtl ** (1.0 / 3.0) * reynolds**0.6
    return nusselt * GAS_THERMAL_CONDUCTIVITY / particle_diameter


def simulation(run: str, base: dict) -> dict:
    feed_name, temperature, pressure, flow = RUNS[run]
    feed = FEEDS[feed_name]
    s = deepcopy(base)
    for key in ("heatCapacityGas", "OutletPressure"):
        s.pop(key, None)
    geometry = s["Geometry"]
    area = math.pi * geometry["InternalDiameter"] ** 2 / 4.0
    s.update({
        "DisplayName": f"Marx et al. (2014) {run}",
        "MixturePredictionMethod": "SCI",
        "BoundaryCondition": "InletPressureInletVelocity",
        "SimulationTime": SIMULATION_TIME[flow],
        "WriteInterval": 1.0,
        "Temperature": temperature,
        "InfluxTemperature": temperature,
        "InletPressure": pressure,
        "ColumnVoidFraction": TOTAL_POROSITY,
        "ParticleDensity": round(BULK_DENSITY / (1.0 - TOTAL_POROSITY), 2),
        "ColumnEntranceVelocity": flow / (area * TOTAL_POROSITY),
        "heatTransferWallExternal": WALL_EXTERNAL,
        "gasThermalConductivity": GAS_THERMAL_CONDUCTIVITY,
    })
    geometry["ColumnVoidFraction"] = TOTAL_POROSITY
    molar_mass = {c["Name"]: c["MolecularWeight"] for c in s["Components"]}
    for component in s["Components"]:
        name = component["Name"]
        component.pop("ReferenceTemperature", None)  # uncentred VantHoff laws, the form of Table 5
        component["GasHeatCapacity"] = MOLAR_GAS_HEAT_CAPACITY
        if component.get("CarrierGas"):
            continue
        component["GasPhaseMolFraction"] = round(100.0 * feed[name], 6)  # percent, as create.ipynb reads them
        component["MassTransferCoefficient"] = MASS_TRANSFER.get(name, MASS_TRANSFER_CO2[(feed_name, temperature)])
        component["HeatOfAdsorption"] = ISOTHERMS[name]["heat"]
        component["PhysisorptionSites"] = [sips_site(name)]
    s["heatTransferGasWall"] = round(heat_transfer_coefficient(s)[1], 2)
    # RUPTURA's particles fill 1 - eps_t of the column; this keeps h a per unit volume at its packed-bed value.
    h_fs = gas_solid_heat_transfer(pressure, temperature, flow / area, feed, molar_mass, s["DynamicViscosity"],
                                   geometry["ParticleDiameter"])
    s["heatTransferGasSolid"] = round(h_fs * (1.0 - BED_POROSITY) / (1.0 - TOTAL_POROSITY), 2)
    return s


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="print the run table without writing inputs")
    arguments = parser.parse_args()
    base = json.loads((ROOT / "base.json").read_text())

    print(f"{'run':4s} {'feed':9s} {'T/K':>5s} {'P/Pa':>9s} {'Q/(m3/s)':>9s} {'h_fw':>7s} {'h_fs':>7s} {'k_CO2':>6s}")
    for run in RUNS:
        s = simulation(run, base)
        feed_name, temperature, pressure, flow = RUNS[run]
        print(f"{run:4s} {feed_name:9s} {temperature:5.0f} {pressure:9.0f} {flow:9.1e} {s['heatTransferGasWall']:7.2f} "
              f"{s['heatTransferGasSolid']:7.2f} {s['Components'][1]['MassTransferCoefficient']:6.2f}")
        if not arguments.dry_run:
            path = ROOT / "ternary" / run / "simulation.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(s, indent=4) + "\n")


if __name__ == "__main__":
    main()
