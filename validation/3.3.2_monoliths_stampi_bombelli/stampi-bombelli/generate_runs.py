#!/usr/bin/env python3
"""Write one RUPTURA `simulation.json` per breakthrough experiment of Table 5.

All physical parameters are taken from Stampi-Bombelli, Storione, Grossmann and
Mazzotti, Ind. Eng. Chem. Res. 63 (2024) 11637-11653:

  * Toth isotherm of the amine-functionalised gamma-alumina pellets (Table 1);
    the monolith isotherm is 3.5 % of the pellet isotherm (Section 2.1).
  * Packed-bed and monolith dimensions (Tables 2 and 3).
  * Feed flow rates and interstitial velocities (Table 5).
  * Mass-transfer coefficients k1 and k2 and the surface-site fraction eta
    (Tables 7 and 8).  The packed bed uses the pseudo-first-order model
    (k = k1); the monolith uses the dual-kinetic model, mapped onto RUPTURA as
    a Toth physisorption site of capacity eta*q_sat relaxing with k1 plus a
    first-order "chemisorption" site of capacity (1-eta)*q_sat relaxing with k2.
  * Axial dispersion from the literature correlations of Table 4.  The fitted
    D_L of Table 7 (6.95u + 0.02 m^2/s on the packed bed) is three orders of
    magnitude above the literature value and drives the 1-D dispersed plug-flow
    model into the well-mixed limit, so it is recorded here but not used.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

R = 8.31446261815324          # J/(mol K)
P0 = 101325.0                 # Pa

# --------------------------------------------------------------- Toth isotherm
NS0 = 1.23                    # mol/kg
B0 = 4839.0 / 1000.0          # 4839 kPa^-1 expressed in Pa^-1
T_EXPONENT0 = 0.25
ALPHA = 0.11
CHI = 0.0
DELTA_H = 70.0e3              # J/mol
T_REFERENCE = 298.0           # K
MONOLITH_UPTAKE_RATIO = 0.035  # q*_monolith / q*_pellet, Section 2.1

MOLAR_MASS = {"N2": 0.0280134, "CO2": 0.0440095}
VISCOSITY = 1.8e-5            # Pa s


def toth_exponent(temperature):
    return T_EXPONENT0 + ALPHA * (1.0 - T_REFERENCE / temperature)


def saturation_capacity(temperature):
    return NS0 * math.exp(CHI * (1.0 - temperature / T_REFERENCE))


# ------------------------------------------------------------------ contactors
PACKED_BED = dict(
    length=0.325, inner_diameter=0.0337, outer_diameter=0.055,
    particle_diameter=0.003, particle_density=1044.0,
    sorbent_mass=0.187, temperature=297.15,
)
PACKED_BED["cross_section"] = math.pi / 4.0 * PACKED_BED["inner_diameter"] ** 2
PACKED_BED["void_fraction"] = 1.0 - PACKED_BED["sorbent_mass"] / (
    PACKED_BED["particle_density"] * PACKED_BED["cross_section"] * PACKED_BED["length"]
)
PACKED_BED["open_area"] = PACKED_BED["void_fraction"] * PACKED_BED["cross_section"]

MONOLITH = dict(
    length=0.13, width=0.029, channel=0.002, wall_thickness=0.0004,
    number_of_channels=144, sorbent_mass=0.05782, temperature=295.15,
)
MONOLITH["open_area"] = MONOLITH["number_of_channels"] * MONOLITH["channel"] ** 2
MONOLITH["frontal_area"] = MONOLITH["width"] ** 2
MONOLITH["outer_diameter"] = math.sqrt(4.0 * MONOLITH["frontal_area"] / math.pi)
MONOLITH["void_fraction"] = MONOLITH["open_area"] / MONOLITH["frontal_area"]
MONOLITH["solid_to_fluid"] = (
    MONOLITH["frontal_area"] - MONOLITH["open_area"]) / MONOLITH["open_area"]
# Chosen so that solid volume x particle density reproduces the 57.82 g monolith.
MONOLITH["particle_density"] = MONOLITH["sorbent_mass"] / (
    MONOLITH["open_area"] * MONOLITH["length"] * MONOLITH["solid_to_fluid"])

# Wall model for the non-isothermal packed-bed runs.  The paper reports
# h_wall = 10 W/(m^2 K) and h_ambient = 26 W/(m^2 K) but no wall dimensions,
# density or heat capacity, so the gas-wall coefficient and the wall thickness
# below were chosen to reproduce the measured 43 C peak of Figure 4b.
THERMAL = dict(
    gasThermalConductivity=0.026, wallThermalConductivity=16.0,
    heatTransferGasSolid=20.0, heatTransferGasWall=200.0,
    heatTransferWallExternal=26.0, heatCapacityGas=1040.0,
    heatCapacitySolid=784.0, heatCapacityWall=500.0, wallDensity=8000.0,
)

# ------------------------------------------------------- transport coefficients
# Table 8: k1 and k2 per contactor and feed concentration; eta from Table 7.
KINETICS = {
    ("packed_bed", "400ppm"): dict(k1=0.000088, k2=0.000082, eta=1.00),
    ("packed_bed", "5.6pct"): dict(k1=0.0063, k2=0.0010, eta=1.00),
    ("monolith", "400ppm"): dict(k1=0.037, k2=0.0011, eta=0.75),
    ("monolith", "5.6pct"): dict(k1=2.7, k2=0.0011, eta=0.75),
}

# Table 7 lists two axial dispersion coefficients per contactor: the value the
# authors fitted to their own breakthrough curves, and the value the literature
# correlations of Table 4 predict.  Both are run here, into one directory each,
# so the effect of the choice can be read straight off the figures.
DISPERSION = {
    "dl_from_stampibombelli": {
        "packed_bed": lambda u: 6.95 * u + 0.02,
        "monolith": lambda u: 1.22 * u * u + 0.0004,
    },
    "dl_from_literature": {
        "packed_bed": lambda u: 0.0015 * u + 0.00003,
        "monolith": lambda u: 0.0013 * u * u + 0.0000162,
    },
}

# ------------------------------------------------------------ experimental runs
# (contactor, CO2 mole fraction, total feed flow in mmol/s) from Table 5,
# together with the integration settings used for that run.
# Integration settings, shared by every case.  ROS34PW2 chooses its own steps at
# its default tolerances.  The runs span 24 000 s with a profile every 12 s.
# ROS34PW2 integrates every case in a few seconds, including the
# dispersion-dominated fitted-D_L runs on which CVODE needed per-case settings
# (TODO B-ITG-11); the two integrators agree to 0.1 % in t50 (2026-11-03).
INTEGRATOR = "ROS34PW2"
INTEGRATION = {
    "InletCondition": "Danckwerts",
    "NumberOfGridPoints": 50,
    "SimulationTime": 24000.0,
    "WriteInterval": 12.0,
}

# Cases whose breakthrough does not fit in the 24 000 s that the settings above
# span, with the simulated time in s they need instead.  Only the packed bed at
# 400 ppm is affected: Table 6 puts its t90 at 164 255 s.
TIME_OVERRIDE = {
    "packed_bed_400ppm_2mmol": 200000.0,
}

# (contactor, CO2 mole fraction in the feed, total feed flow in mmol/s),
# from Table 5.
RUNS = [
    ("packed_bed", 400.0e-6, 2.00),
    ("packed_bed", 0.056, 2.00),
    ("packed_bed", 0.056, 1.00),
    ("packed_bed", 0.056, 0.50),
    ("monolith", 400.0e-6, 3.00),
    ("monolith", 400.0e-6, 2.00),
    ("monolith", 400.0e-6, 1.00),
    ("monolith", 400.0e-6, 0.18),
    ("monolith", 0.056, 0.18),
    ("monolith", 0.056, 0.06),
    ("monolith", 0.056, 0.02),
]




def concentration_tag(mole_fraction):
    return "400ppm" if mole_fraction < 1.0e-3 else "5.6pct"


def flow_tag(flow_mmol_s):
    return f"{flow_mmol_s:g}".replace(".", "p")


def case_name(contactor, mole_fraction, flow_mmol_s):
    return f"{contactor}_{concentration_tag(mole_fraction)}_{flow_tag(flow_mmol_s)}mmol"


def interstitial_velocity(contactor, flow_mmol_s):
    data = PACKED_BED if contactor == "packed_bed" else MONOLITH
    volumetric = flow_mmol_s * 1.0e-3 * R * data["temperature"] / P0
    return volumetric / data["open_area"]


def build(contactor, mole_fraction, flow_mmol_s, dispersion_model):
    data = PACKED_BED if contactor == "packed_bed" else MONOLITH
    concentration = concentration_tag(mole_fraction)
    kinetics = KINETICS[(contactor, concentration)]
    temperature = data["temperature"]
    velocity = interstitial_velocity(contactor, flow_mmol_s)
    dispersion = DISPERSION[dispersion_model][contactor](velocity)

    exponent = toth_exponent(temperature)
    capacity = saturation_capacity(temperature)
    if contactor == "monolith":
        capacity *= MONOLITH_UPTAKE_RATIO
    eta = kinetics["eta"]

    carbon_dioxide = {
        "Name": "CO2",
        "GasPhaseMolFraction": mole_fraction,
        "MassTransferCoefficient": kinetics["k1"],
        "AxialDispersionCoefficient": dispersion,
        "MolecularWeight": MOLAR_MASS["CO2"],
        "HeatOfAdsorption": DELTA_H,
        "ReferenceTemperature": T_REFERENCE,
        # The affinity b0 at T_reference follows a VantHoff law whose Energy is the HeatOfAdsorption.
        "PhysisorptionSites": [
            {"Type": "Toth", "Parameters": [eta * capacity, {"Value": B0, "Law": "VantHoff"}, exponent]}
        ],
    }
    if eta < 1.0:
        carbon_dioxide["ChemisorptionSites"] = [{
            "Type": "Toth",
            "Parameters": [(1.0 - eta) * capacity, {"Value": B0, "Law": "VantHoff"}, exponent],
            "RateEquationParameters": {"MassTransferCoefficient": kinetics["k2"]},
            "HeatOfAdsorption": DELTA_H,
            "MaximumLoading": (1.0 - eta) * capacity,
        }]

    if contactor == "packed_bed":
        geometry = {
            "Type": "PackedBed",
            "ColumnVoidFraction": data["void_fraction"],
            "ParticleDiameter": data["particle_diameter"],
            "InternalDiameter": data["inner_diameter"],
            "OuterDiameter": data["outer_diameter"],
        }
        particle_diameter = data["particle_diameter"]
    else:
        geometry = {
            "Type": "Monolith",
            "ChannelShape": "square",
            "InternalChannelDimension": data["channel"],
            "NumberOfChannels": data["number_of_channels"],
            "OuterDiameter": data["outer_diameter"],
            "WashcoatThickness": 0.0,
        }
        particle_diameter = data["wall_thickness"]

    simulation = {
        "SimulationType": "Breakthrough",
        "DisplayName": (f"Stampi-Bombelli 2024 - {contactor.replace('_', ' ')} - "
                        f"{'400 ppm' if concentration == '400ppm' else '5.6%'} - "
                        f"{flow_mmol_s:g} mmol/s"),
        "Temperature": temperature,
        "InfluxTemperature": temperature,
        "ColumnVoidFraction": data["void_fraction"],
        "DynamicViscosity": VISCOSITY,
        "ParticleDiameter": particle_diameter,
        "ParticleDensity": data["particle_density"],
        "InletPressure": P0,
        "PressureGradient": 0.0,
        "ColumnEntranceVelocity": velocity,
        "ColumnLength": data["length"],
        "MixturePredictionMethod": "SPI",
        "BreakthroughIntegrator": INTEGRATOR,
        "BoundaryCondition": "FixedVelocity",
        **INTEGRATION,
        "SimulationTime": TIME_OVERRIDE.get(
            case_name(contactor, mole_fraction, flow_mmol_s),
            INTEGRATION["SimulationTime"]),
        "Components": [
            {"Name": "N2", "GasPhaseMolFraction": 1.0 - mole_fraction,
             "CarrierGas": True, "MolecularWeight": MOLAR_MASS["N2"]},
            carbon_dioxide,
        ],
        "Geometry": geometry,
    }

    # The 5.6 % packed-bed runs release enough heat to shift the isotherm, so the
    # energy balance is switched on there.  Every other run stays within 1 K of
    # ambient in the experiments (Figure 4a,c,d) and is integrated isothermally.
    if contactor == "packed_bed" and concentration == "5.6pct":
        simulation["energyBalance"] = True
        simulation.update(THERMAL)

    return simulation


def main():
    for dispersion_model in DISPERSION:
        root = ROOT / dispersion_model
        root.mkdir(exist_ok=True)
        print(dispersion_model)
        for contactor, mole_fraction, flow in RUNS:
            name = case_name(contactor, mole_fraction, flow)
            if name == "packed_bed_400ppm_2mmol": continue
            directory = root / name
            directory.mkdir(exist_ok=True)
            simulation = build(contactor, mole_fraction, flow, dispersion_model)
            with (directory / "simulation.json").open("w") as stream:
                json.dump(simulation, stream, indent=2)
                stream.write("\n")
            dispersion = simulation["Components"][1]["AxialDispersionCoefficient"]
            peclet = (simulation["ColumnEntranceVelocity"] * simulation["ColumnLength"]
                      / dispersion)
            print(f"  {name:34s} u={simulation['ColumnEntranceVelocity']:9.6f} m/s  "
                  f"D_L={dispersion:.3e} m2/s  Pe={peclet:8.3f}  "
                  f"t_end={simulation['SimulationTime']:8.0f} s")


if __name__ == "__main__":
    main()
