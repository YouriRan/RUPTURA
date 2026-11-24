#!/usr/bin/env python3
"""Write one RUPTURA `simulation.json` per Stadler et al. (2022) experiment.

All physical parameters are taken from Stadler, Bender and Pfeifer, Front. Chem.
Eng. 4 (2022) 1000064 (fceng-04-1000064.pdf), unless marked otherwise:

  * Freundlich CO2 isotherm on K-HDMO at 523 K, dry (Eq. 11, Table 3).
  * Sorption kinetics (Eqs. 6-10, Table 3): k_ads for adsorption and the
    coverage-dependent k_des = k1 exp[beta (q/q_max)/RT] for desorption, with
    q_max = q*(20 bar).  The bulk-to-pore LDF coefficient of Eq. 13 is four orders of
    magnitude faster than either, so the pore gas follows the bulk gas and the laws
    act on the bulk partial pressure.
  * WGS rate (Eq. 4) with Moe's equilibrium constant (Eq. 5), as a catalyst-weighted
    reaction on the bulk concentrations (the authors' Fig. 4 model is reproduced to
    1e-4 by a plug-flow reactor without pores).
  * Laboratory reactor: two slits of 2 mm x 50 mm x 100 mm (Section 2.3); bed void
    fraction 0.4 and particle radius 100 um (Table 2); sorbent masses 9.70 g
    (sorbent only) and 9.35 g with 0.51 g catalyst in the first quarter (Section 2.5);
    0.50 g catalyst in the first quarter for the activity tests (Section 2.4).
  * Axial dispersion D_ax = D_mix sqrt(eps_b) + u r_p (Eq. 12), with u the
    superficial velocity and D_mix the Fuller binary diffusivity in N2 (Poling et
    al., The Properties of Gases and Liquids, 5th ed., Eq. 11-4.4).  RUPTURA holds
    one D_ax per component, so the cycles use the adsorption-step value.
  * Feed flows: 2000 mL/min (STP) for adsorption and the activity tests, 1000 mL/min
    for regeneration, from Table 2 of Stadler et al., Energies 14 (2021) 355; the
    2022 paper does not state them.  STP is taken as 273.15 K and 101325 Pa.

Every conversion is written out in conversions.md.

Cases written:

  figure4_wgs_conversion/<T>K_SG<S/G>
      Steady-state CO conversion over the 25 mm catalyst section (no sorbent) for
      225, 250 and 275 C and S/G = 0.5-4.5.
  figure4_grid_study_548K_SG4p5/N<n>
      The fastest Figure 4 case with the physical dispersion switched off, on
      successively finer grids, to show the first-order numerical dispersion of the
      upwind convection converging to the plug-flow result.
  figure5_isotherm_CO2_523K
      MixturePrediction of the pure-CO2 Freundlich isotherm over 0-800 mmHg.
  figure6_sorbent_only_fresh_bed_8bar_523K
      One sorbent-only adsorption step on a clean bed, CO2 on a physisorption site
      (adsorption kinetics only).
  figure6_8a_sorbent_only_cycles_8bar-1bar_523K
      Five sorbent-only cycles (Figures 6 and 8A): 15 min adsorption at 8 bar and
      40 min regeneration at 1 bar with H2O/N2 = 0.4/0.6.
  figure7_8b_sewgs_cycles_8bar-1bar_523K
      Five SE-WGS cycles (Figures 7 and 8B) on a two-section multibed column: catalyst
      and sorbent in the first quarter, sorbent downstream.

    python generate.py            # (re)write every case directory
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

R = 8.31446261815324          # J/(mol K)
T_STP = 273.15                # K
P_STP = 101325.0              # Pa
MMHG = 133.322387415          # Pa

# ------------------------------------------------------------- isotherm, Table 3
# q = k_Fr p^(1/n_Fr) with k_Fr in mmol/g/bar^(1/n_Fr).  RUPTURA's Freundlich site
# is q = a p^(1/nu) with p in Pa and q in mol/kg: a = k_Fr (1e5)^(-1/n_Fr), nu = n_Fr.
FREUNDLICH_K = 0.177          # mmol/g bar^(-1/n) = mol/kg bar^(-1/n)
FREUNDLICH_INVERSE_N = 0.220  # 1/n_Fr [-]
FREUNDLICH = [FREUNDLICH_K * 1.0e5 ** (-FREUNDLICH_INVERSE_N), 1.0 / FREUNDLICH_INVERSE_N]

# ------------------------------------------------------ sorption kinetics, Table 3
K_ADSORPTION = 0.1125         # 1/s, Eq. 6
K1_DESORPTION = 5.38e-4       # 1/s, Eq. 10
BETA_DESORPTION = 74000.0     # J/mol, Eq. 10
Q_MAX = FREUNDLICH_K * 20.0 ** FREUNDLICH_INVERSE_N  # mol/kg, q*(20 bar), Section 3.3

# ----------------------------------------------------------- WGS, Eqs. 4-5, Table 3
K_INFINITY = 2.01e5           # mol bar^-2 g_cat^-1 h^-1
ACTIVATION_ENERGY = 51845.0   # J/mol
MOE_CONSTANT = math.exp(-4.33)  # K_eq = exp(4577.8/T - 4.33) = K0 exp(Q_r/(R T)), a VantHoff law
MOE_HEAT = 4577.8 * R           # Q_r = -dH_r in J/mol, the Energy of the VantHoff law of K_eq

# ------------------------------------------------------------ laboratory reactor
TEMPERATURE = 523.15          # K, 250 C
NUMBER_OF_SLITS = 2
SLIT_GAP = 0.002              # m
SLIT_WIDTH = 0.050            # m
SLIT_LENGTH = 0.100           # m, flow direction
CATALYST_LENGTH = 0.025       # m, first quarter of the slits
VOID_FRACTION = 0.4           # Table 2
PARTICLE_RADIUS = 100.0e-6    # m, Table 2
CROSS_SECTION = NUMBER_OF_SLITS * SLIT_GAP * SLIT_WIDTH
# The slit cross-section is rectangular; the diameter of a tube with the same area
# is passed to the PackedBed geometry, which uses it only for the wall model.  The
# runs are isothermal, so it has no effect.
EQUIVALENT_DIAMETER = math.sqrt(4.0 * CROSS_SECTION / math.pi)
GEOMETRY = {
    "Type": "PackedBed",
    "ColumnVoidFraction": VOID_FRACTION,
    "ParticleDiameter": 2.0 * PARTICLE_RADIUS,
    "InternalDiameter": EQUIVALENT_DIAMETER,
}

# ----------------------------------------------------------------- components
MOLAR_MASS = {"N2": 0.0280134, "CO": 0.0280101, "CO2": 0.0440095, "H2": 0.00201588, "H2O": 0.01801528}
# Fuller atomic diffusion volumes (Poling et al., Table 11-1).
DIFFUSION_VOLUME = {"N2": 18.5, "CO": 18.0, "CO2": 26.9, "H2": 6.12, "H2O": 13.1}
SPECIES = ("CO", "CO2", "H2", "H2O")  # N2 is the carrier

# ------------------------------------------------------------------ operation
ADSORPTION_PRESSURE = 8.0e5   # Pa
ADSORPTION_FLOW_STP = 2000.0  # mL/min, Stadler et al. (2021), Table 2
ADSORPTION_TIME = 900.0       # s, 15 min
REGENERATION_PRESSURE = 1.0e5  # Pa
REGENERATION_FLOW_STP = 1000.0  # mL/min, Stadler et al. (2021), Table 2
REGENERATION_TIME = 2400.0    # s, 40 min
REGENERATION_FEED = {"N2": 0.6, "H2O": 0.4}
NUMBER_OF_CYCLES = 5
SORBENT_ONLY_MASS = 9.70e-3   # kg, Section 2.5
SORBENT_ONLY_FEED = {"N2": 0.65, "CO": 0.05, "CO2": 0.05, "H2": 0.05, "H2O": 0.20}  # Figure 6
SEWGS_SORBENT_MASS = 9.35e-3  # kg, Section 2.5
SEWGS_CATALYST_MASS = 0.51e-3  # kg, in the first quarter
SEWGS_FEED = {"N2": 0.75, "CO": 0.05, "H2O": 0.20}  # Figure 7
ACTIVITY_CATALYST_MASS = 0.50e-3  # kg, Section 2.4, diluted with SiC in the first quarter
ACTIVITY_CO_FRACTION = 0.10
ACTIVITY_TEMPERATURES = (498.15, 523.15, 548.15)  # K, 225, 250, 275 C
ACTIVITY_STEAM_TO_GAS = [0.5 + 0.25 * i for i in range(17)]  # 0.5-4.5
GRID_STUDY_POINTS = (25, 50, 100, 200, 400, 800, 1600)

# Fresh-bed step: ROS34PW2 adapts its internal step; the solution is written every 0.5 s.
# The gas residence time in the bed is about 1 s.
FRESH_BED_INTEGRATION = {
    "BreakthroughIntegrator": "ROS34PW2",
    "InletCondition": "Danckwerts",
    "NumberOfGridPoints": 100,
    "SimulationTime": ADSORPTION_TIME,
    "WriteInterval": 0.5,
    "WriteOutletOnly": True,
}
# Cycles: 40 cells put the catalyst/sorbent boundary on a node (2.5 mm cells); the
# outlet is written every 2 s, which keeps five cycles at a few hundred MB.
CYCLE_INTEGRATION = {
    "BreakthroughIntegrator": "ROS34PW2",
    "InletCondition": "Danckwerts",
    "NumberOfGridPoints": 40,
    "WriteInterval": 2.0,
}
# Activity tests: 25 mm column, gas residence time about 0.25 s; 3 s reach steady state. The fresh bed and the
# activity tests are analysed at the outlet only; the cycles keep the column for the bed-averaged loading.
ACTIVITY_INTEGRATION = {
    "BreakthroughIntegrator": "ROS34PW2",
    "InletCondition": "Danckwerts",
    "NumberOfGridPoints": 400,
    "SimulationTime": 3.0,
    "WriteInterval": 0.5,
    "WriteOutletOnly": True,
}


def superficial_velocity(flow_stp, temperature, pressure):
    """Superficial velocity in m/s of a flow given in mL/min at STP."""
    molar_flow = P_STP * flow_stp * 1.0e-6 / 60.0 / (R * T_STP)
    return molar_flow * R * temperature / (pressure * CROSS_SECTION)


def gas_components(feed, temperature, pressure, superficial, carbon_dioxide_sites=None, dispersion=True):
    """N2 carrier plus CO, CO2, H2 and H2O with Eq. 12 dispersion.

    `carbon_dioxide_sites` holds the PhysisorptionSites/ChemisorptionSites entries of
    CO2 (and its MassTransferCoefficient); without it every species is non-adsorbing,
    for which a zero Henry site is RUPTURA's placeholder.
    """
    components = [{"Name": "N2", "GasPhaseMolFraction": feed.get("N2", 0.0), "CarrierGas": True,
                   "MolecularWeight": MOLAR_MASS["N2"]}]
    for name in SPECIES:
        reduced_mass = 2.0 / (1.0 / (1.0e3 * MOLAR_MASS[name]) + 1.0 / (1.0e3 * MOLAR_MASS["N2"]))
        diffusivity = 1.0e-4 * 0.00143 * temperature ** 1.75 / (
            pressure / 1.0e5 * math.sqrt(reduced_mass)
            * (DIFFUSION_VOLUME[name] ** (1.0 / 3.0) + DIFFUSION_VOLUME["N2"] ** (1.0 / 3.0)) ** 2)
        component = {
            "Name": name,
            "GasPhaseMolFraction": feed.get(name, 0.0),
            "AxialDispersionCoefficient": (
                diffusivity * math.sqrt(VOID_FRACTION) + superficial * PARTICLE_RADIUS if dispersion else 0.0),
            "MolecularWeight": MOLAR_MASS[name],
            "PhysisorptionSites": [{"Type": "Henry", "Parameters": [0.0]}],
        }
        if name == "CO2" and carbon_dioxide_sites is not None:
            component.update(carbon_dioxide_sites)
        components.append(component)
    return components


def water_gas_shift(temperature, catalyst_per_bed_volume, adsorbents=None):
    """CO + H2O = CO2 + H2 on the bulk concentrations, per unit fluid volume.

    Eq. 4 per gram of catalyst with pressures in bar becomes, with c in mol/m^3,
    k0 = k_inf (RT/1e5)^2 * 1e3 [g/kg] / 3600 [s/h] * rho_cat / eps_b, where rho_cat is the
    catalyst mass per bed volume (conversions.md, Section 9).
    """
    forward = (K_INFINITY * (R * temperature / 1.0e5) ** 2 * 1.0e3 / 3600.0
               * catalyst_per_bed_volume / VOID_FRACTION)
    reaction = {
        "Phase": "BulkConcentration",
        "Style": "GeneralPowerLaw",
        "Reactants": ["CO", "H2O"],
        "Products": ["CO2", "H2"],
        "ForwardOrders": [1, 1],
        "BackwardOrders": [1, 1],
        "Kinetics": {
            "ForwardRateCoefficient": {"Value": forward, "Law": "Arrhenius", "ActivationEnergy": ACTIVATION_ENERGY},
            "EquilibriumConstant": {"Value": MOE_CONSTANT, "Law": "VantHoff"},
            "HeatOfReaction": MOE_HEAT,
        },
    }
    if adsorbents is not None:
        reaction["Adsorbents"] = adsorbents
    return reaction


def write(name, simulation):
    directory = ROOT / name
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "simulation.json").open("w") as stream:
        json.dump(simulation, stream, indent=2)
        stream.write("\n")


def main():
    # CO2 on K-HDMO.  The fresh-bed step keeps the original physisorption LDF (adsorption only); the
    # cycles need the desorption branch, which lives on the chemisorption FirstOrder law.
    physisorbed_carbon_dioxide = {
        "MassTransferCoefficient": K_ADSORPTION,
        "PhysisorptionSites": [{"Type": "Freundlich", "Parameters": FREUNDLICH}],
    }
    chemisorbed_carbon_dioxide = {
        "ChemisorptionSites": [{
            "Type": "Freundlich",
            "Parameters": FREUNDLICH,
            # Caps q*, which never reaches q*(20 bar) at 8 bar or less.
            "MaximumLoading": Q_MAX,
            "RateEquation": "FirstOrder",
            "RateEquationParameters": {
                "MassTransferCoefficient": K_ADSORPTION,
                "DesorptionMassTransferCoefficient": K1_DESORPTION,
                "DesorptionCoverageEnergy": BETA_DESORPTION,
                "CoverageReferenceLoading": Q_MAX,
            },
        }],
    }

    # Figure 5: the dry CO2 isotherm at 250 C, 0-800 mmHg.
    write("figure5_isotherm_CO2_523K", {
        "SimulationType": "MixturePrediction",
        "DisplayName": "Stadler 2022 - Figure 5 - CO2 on K-HDMO at 523 K",
        "Temperature": TEMPERATURE,
        "MixturePredictionMethod": "IAST",
        "PressureScale": "Linear",
        "PressureStart": 1.0 * MMHG,
        "PressureEnd": 800.0 * MMHG,
        "NumberOfPressurePoints": 200,
        "Components": [{
            "Name": "CO2",
            "GasPhaseMolFraction": 1.0,
            "PhysisorptionSites": [{"Type": "Freundlich", "Parameters": FREUNDLICH}],
        }],
    })
    print(f"figure5_isotherm_CO2_523K  a={FREUNDLICH[0]:.6g} mol/kg/Pa^(1/nu)  nu={FREUNDLICH[1]:.6g}")

    # Figure 4: steady-state conversion over the catalyst section, which holds no sorbent.
    activity_density = ACTIVITY_CATALYST_MASS / (CROSS_SECTION * CATALYST_LENGTH)
    for temperature in ACTIVITY_TEMPERATURES:
        superficial = superficial_velocity(ADSORPTION_FLOW_STP, temperature, ADSORPTION_PRESSURE)
        for steam_to_gas in ACTIVITY_STEAM_TO_GAS:
            feed = {"CO": ACTIVITY_CO_FRACTION, "H2O": steam_to_gas * ACTIVITY_CO_FRACTION}
            feed["N2"] = 1.0 - feed["CO"] - feed["H2O"]
            name = f"figure4_wgs_conversion/{temperature:.0f}K_SG{steam_to_gas:g}".replace(".", "p")
            write(name, {
                "SimulationType": "Breakthrough",
                "DisplayName": f"Stadler 2022 - Figure 4 - {temperature:.2f} K, S/G {steam_to_gas:g}",
                "Temperature": temperature,
                "InfluxTemperature": temperature,
                "ColumnVoidFraction": VOID_FRACTION,
                "DynamicViscosity": 2.65e-5,
                "ParticleDiameter": 2.0 * PARTICLE_RADIUS,
                "ParticleDensity": 0.0,  # no sorbent in the catalyst/SiC section
                "InletPressure": ADSORPTION_PRESSURE,
                "PressureGradient": 0.0,
                "ColumnEntranceVelocity": superficial / VOID_FRACTION,
                "ColumnLength": CATALYST_LENGTH,
                "MixturePredictionMethod": "IAST",
                "BoundaryCondition": "FixedVelocity",
                **ACTIVITY_INTEGRATION,
                "Components": gas_components(feed, temperature, ADSORPTION_PRESSURE, superficial),
                "Reactions": [water_gas_shift(temperature, activity_density)],
                "Geometry": GEOMETRY,
            })
    print(f"figure4_wgs_conversion/    {len(ACTIVITY_TEMPERATURES) * len(ACTIVITY_STEAM_TO_GAS)} cases, "
          f"rho_cat={activity_density:.1f} kg/m3 bed")

    # Figure 4 grid study: fastest case (275 C, S/G = 4.5), physical dispersion off.
    temperature = ACTIVITY_TEMPERATURES[-1]
    superficial = superficial_velocity(ADSORPTION_FLOW_STP, temperature, ADSORPTION_PRESSURE)
    feed = {"CO": ACTIVITY_CO_FRACTION, "H2O": 4.5 * ACTIVITY_CO_FRACTION}
    feed["N2"] = 1.0 - feed["CO"] - feed["H2O"]
    for points in GRID_STUDY_POINTS:
        write(f"figure4_grid_study_548K_SG4p5/N{points}", {
            "SimulationType": "Breakthrough",
            "DisplayName": f"Stadler 2022 - Figure 4 grid study - 548.15 K, S/G 4.5, {points} cells, no D_ax",
            "Temperature": temperature,
            "InfluxTemperature": temperature,
            "ColumnVoidFraction": VOID_FRACTION,
            "DynamicViscosity": 2.65e-5,
            "ParticleDiameter": 2.0 * PARTICLE_RADIUS,
            "ParticleDensity": 0.0,
            "InletPressure": ADSORPTION_PRESSURE,
            "PressureGradient": 0.0,
            "ColumnEntranceVelocity": superficial / VOID_FRACTION,
            "ColumnLength": CATALYST_LENGTH,
            "MixturePredictionMethod": "IAST",
            "BoundaryCondition": "FixedVelocity",
            **ACTIVITY_INTEGRATION,
            "NumberOfGridPoints": points,
            "Components": gas_components(feed, temperature, ADSORPTION_PRESSURE, superficial, dispersion=False),
            "Reactions": [water_gas_shift(temperature, activity_density)],
            "Geometry": GEOMETRY,
        })
    print(f"figure4_grid_study_548K_SG4p5/  N = {', '.join(map(str, GRID_STUDY_POINTS))}")

    # Figure 6: sorbent-only adsorption step on a clean bed.
    adsorption_superficial = superficial_velocity(ADSORPTION_FLOW_STP, TEMPERATURE, ADSORPTION_PRESSURE)
    regeneration_superficial = superficial_velocity(REGENERATION_FLOW_STP, TEMPERATURE, REGENERATION_PRESSURE)
    sorbent_only_density = SORBENT_ONLY_MASS / ((1.0 - VOID_FRACTION) * CROSS_SECTION * SLIT_LENGTH)
    base = {
        "Temperature": TEMPERATURE,
        "InfluxTemperature": TEMPERATURE,
        "ColumnVoidFraction": VOID_FRACTION,
        # N2 at 523 K.  Only the pressure-drop model uses it, and FixedVelocity bypasses that.
        "DynamicViscosity": 2.65e-5,
        "ParticleDiameter": 2.0 * PARTICLE_RADIUS,
        "InletPressure": ADSORPTION_PRESSURE,
        "PressureGradient": 0.0,
        "ColumnEntranceVelocity": adsorption_superficial / VOID_FRACTION,
        "ColumnLength": SLIT_LENGTH,
        "MixturePredictionMethod": "IAST",
        "BoundaryCondition": "FixedVelocity",
        "Geometry": GEOMETRY,
    }
    write("figure6_sorbent_only_fresh_bed_8bar_523K", {
        "SimulationType": "Breakthrough",
        "DisplayName": "Stadler 2022 - Figure 6 feed - sorbent only - clean bed - 8 bar, 523 K",
        **base,
        "ParticleDensity": sorbent_only_density,
        **FRESH_BED_INTEGRATION,
        "Components": gas_components(SORBENT_ONLY_FEED, TEMPERATURE, ADSORPTION_PRESSURE, adsorption_superficial,
                                     physisorbed_carbon_dioxide),
    })

    # Figures 6 and 8A: sorbent-only cycles.  The first adsorption step starts from a clean bed
    # filled with N2 at 8 bar, as in the authors' model.
    phases = [
        {"Name": "Adsorption", "InletPressure": ADSORPTION_PRESSURE,
         "Duration": ADSORPTION_TIME},
        {"Name": "Regeneration", "InletPressure": REGENERATION_PRESSURE,
         "ColumnEntranceVelocity": regeneration_superficial / VOID_FRACTION,
         "GasPhaseMolFractions": REGENERATION_FEED,
         "Duration": REGENERATION_TIME},
    ]
    write("figure6_8a_sorbent_only_cycles_8bar-1bar_523K", {
        "SimulationType": "SwingAdsorption",
        "DisplayName": "Stadler 2022 - Figures 6 and 8A - sorbent-only cycles",
        **base,
        "ParticleDensity": sorbent_only_density,
        **CYCLE_INTEGRATION,
        "NumberOfCycles": NUMBER_OF_CYCLES,
        "SwingAdsorptionPhases": phases,
        "Components": gas_components(SORBENT_ONLY_FEED, TEMPERATURE, ADSORPTION_PRESSURE, adsorption_superficial,
                                     chemisorbed_carbon_dioxide),
    })

    # Figures 7 and 8B: SE-WGS cycles.  The sorbent is spread uniformly over the slits, so the catalyst
    # quarter holds a quarter of it (as in mode C of Stadler et al. 2021); only the reaction differs.
    sewgs_density = SEWGS_SORBENT_MASS / ((1.0 - VOID_FRACTION) * CROSS_SECTION * SLIT_LENGTH)
    catalyst_density = SEWGS_CATALYST_MASS / (CROSS_SECTION * CATALYST_LENGTH)
    adsorbent = {"ColumnVoidFraction": VOID_FRACTION, "ParticleDensity": sewgs_density,
                 "ParticleDiameter": 2.0 * PARTICLE_RADIUS}
    write("figure7_8b_sewgs_cycles_8bar-1bar_523K", {
        "SimulationType": "SwingAdsorption",
        "DisplayName": "Stadler 2022 - Figures 7 and 8B - SE-WGS cycles",
        **base,
        "ParticleDensity": sewgs_density,
        **CYCLE_INTEGRATION,
        "NumberOfCycles": NUMBER_OF_CYCLES,
        "SwingAdsorptionPhases": phases,
        "Components": gas_components(SEWGS_FEED, TEMPERATURE, ADSORPTION_PRESSURE, adsorption_superficial,
                                     chemisorbed_carbon_dioxide),
        "Adsorbents": [{"Name": "catalyst-zone", **adsorbent}, {"Name": "sorbent-zone", **adsorbent}],
        # The upwind convection adds a numerical dispersion v dz/2 (first order in dz, see the Figure 4 grid
        # study), which lowers the conversion; 40 cells over the catalyst quarter keep that below 0.005.
        "ColumnSections": [
            {"Adsorbent": "catalyst-zone", "Length": CATALYST_LENGTH, "NumberOfGridPoints": 40},
            {"InterfaceLength": 0.0},
            {"Adsorbent": "sorbent-zone", "Length": SLIT_LENGTH - CATALYST_LENGTH, "NumberOfGridPoints": 30},
        ],
        "Reactions": [water_gas_shift(TEMPERATURE, catalyst_density, ["catalyst-zone"])],
    })
    print(f"cycles: v_ads={adsorption_superficial / VOID_FRACTION:.5f} m/s, "
          f"v_reg={regeneration_superficial / VOID_FRACTION:.5f} m/s, rho_p={sorbent_only_density:.1f} / "
          f"{sewgs_density:.1f} kg/m3, rho_cat={catalyst_density:.1f} kg/m3 bed, q_max={Q_MAX:.4f} mol/kg")


if __name__ == "__main__":
    main()
