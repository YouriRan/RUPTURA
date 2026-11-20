#!/usr/bin/env python3
"""Generate every Ruptura simulation directory used by ``liu_paper_figures.ipynb``.

Each figure panel of Liu et al. (2024) is backed by one or more Ruptura runs.
Re-running this script rewrites ``simulations/`` from scratch; ``run.py`` then
executes every directory it produced and ``liu_paper_figures.ipynb`` reads the
resulting ``component_*.data`` files.  No figure data is stored as a CSV in this
project except the digitized paper references in ``reference_digitized/``.

Reference:
    H. Liu, H. Lin, S. Dai, D.-e. Jiang, Ind. Eng. Chem. Res. 63 (2024) 5871-5879.
    DOI 10.1021/acs.iecr.3c04535.

All quantities are SI: Pa, K, m, m/s, m3/s, mol/m3, mol/kg, s.  The unit
conversions are derived in ``conversions.md``.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SIMULATION_DIR = ROOT / "simulations"

GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1
STANDARD_TEMPERATURE = 273.15  # K
STANDARD_PRESSURE = 100_000.0  # Pa
ATMOSPHERE = 101_325.0  # Pa

# Ruptura resolves these DAC fronts to better than 1% on 30 grid points.  ROS34PW2
# integrates them at its default tolerances (since 2026-11-03; before, CVODE, which
# needed per-figure tolerances and stopped in 36 of the 50 runs, TODO B-EX-10).
GRID_POINTS = 30
OUTPUT_BLOCKS = 300


# --------------------------------------------------------------------------
# Paper equations 3-9
# --------------------------------------------------------------------------


def toth_at_temperature(
    temperature: float,
    *,
    reference_temperature: float,
    n_s0: float,
    affinity0: float,
    heat: float,
    exponent0: float,
    chi: float = 0.0,
    alpha: float = 0.0,
) -> tuple[float, float, float]:
    """Equations 4-6: Toth ``n_s`` [mol/kg], ``b`` [1/Pa] and ``t_T`` [-] at ``T``."""

    t0 = reference_temperature
    n_s = n_s0 * math.exp(chi * (1.0 - temperature / t0))
    affinity = affinity0 * math.exp(heat / (GAS_CONSTANT * t0) * (t0 / temperature - 1.0))
    exponent = exponent0 + alpha * (1.0 - t0 / temperature)
    return n_s, affinity, exponent


def gab_loading(relative_humidity: float, c_m: float, c_g: float, k_ads: float) -> float:
    """Equation 7: GAB water loading in mol/kg."""

    phi = min(max(relative_humidity, 0.0), 0.999999)
    numerator = c_m * c_g * k_ads * phi
    denominator = (1.0 - k_ads * phi) * (1.0 + (c_g - 1.0) * k_ads * phi)
    return numerator / denominator


def antoine_saturation_pressure(temperature: float) -> float:
    """Water saturation pressure in Pa, using the same Antoine form as Ruptura's GAB site.

    Ruptura evaluates its GAB isotherm against this correlation internally, so the
    notebook must use it (and not the paper's Buck equation) whenever it converts
    a Ruptura water partial pressure back to a relative humidity.
    """

    return 10.0 ** (8.07131 - 1730.63 / (temperature - 39.724)) * 133.32236842105263


def standard_flow_to_superficial_velocity(
    standard_flow_m3_s: float, area_m2: float, temperature: float
) -> float:
    """Standard m3/s at 273.15 K and 100000 Pa to superficial m/s at ``temperature``."""

    return standard_flow_m3_s * (temperature / STANDARD_TEMPERATURE) / area_m2


# --------------------------------------------------------------------------
# Ruptura input writers
# --------------------------------------------------------------------------


def latex_number(value: float, digits: int = 3) -> str:
    """A legend-ready number: plain decimal, or LaTeX scientific when long."""

    plain = f"{value:.{digits}g}"
    if "e" not in plain:
        return plain
    mantissa, exponent = plain.split("e")
    exponent = int(exponent)
    if float(mantissa) == 1.0:
        return rf"$10^{{{exponent}}}$"
    return rf"${mantissa} \times 10^{{{exponent}}}$"


def write_isotherm_case(
    name: str,
    *,
    figure: int,
    panel: str,
    label: str,
    component: str,
    temperature: float,
    pressure_start: float,
    pressure_end: float,
    points: int,
    sites: list[dict],
    scale: str = "Linear",
    relative_humidity: float | None = None,
    metadata: dict | None = None,
) -> dict:
    """Write a pure-component ``MixturePrediction`` directory and return its index entry.

    ``relative_humidity`` sets the fixed ``RelativeHumidity`` that a humidity-dependent site sees.
    """

    directory = SIMULATION_DIR / name
    directory.mkdir(parents=True, exist_ok=True)
    simulation = {
        "SimulationType": "MixturePrediction",
        "DisplayName": label,
        "Temperature": temperature,
        "MixturePredictionMethod": "SPI",
        "PressureScale": scale,
        "PressureStart": pressure_start,
        "PressureEnd": pressure_end,
        "NumberOfPressurePoints": points,
        **({} if relative_humidity is None else {"RelativeHumidity": relative_humidity}),
        "ColumnVoidFraction": 0.5,
        "ParticleDensity": 1000.0,
        "Components": [
            {
                "Name": component,
                "GasPhaseMolFraction": 1.0,
                "MolecularWeight": 0.0440095 if component == "CO2" else 0.01801528,
                "PhysisorptionSites": sites,
            }
        ],
    }
    (directory / "simulation.json").write_text(json.dumps(simulation, indent=2) + "\n")
    return {
        "directory": name,
        "kind": "isotherm",
        "figure": figure,
        "panel": panel,
        "label": label,
        "component": component,
        "temperature_K": temperature,
        **(metadata or {}),
    }


def write_breakthrough_case(
    name: str,
    *,
    figure: int,
    panel: str,
    label: str,
    temperature: float,
    total_pressure: float,
    column_length: float,
    superficial_velocity: float,
    bulk_density: float,
    void_fraction: float,
    internal_diameter: float,
    end_time: float,
    adsorbates: list[dict],
    grid_points: int = GRID_POINTS,
    metadata: dict | None = None,
) -> dict:
    """Write a ``Breakthrough`` directory and return its index entry.

    ``adsorbates`` holds one dict per adsorbing component with the keys ``name``,
    ``mole_fraction``, ``molecular_weight``, ``mass_transfer_coefficient`` and
    ``sites``, and ``water`` set for the water component whose partial pressure
    and isotherm the humidity-dependent sites follow.  The nitrogen carrier gas
    is added automatically.
    """

    directory = SIMULATION_DIR / name
    directory.mkdir(parents=True, exist_ok=True)

    carrier_fraction = 1.0 - sum(entry["mole_fraction"] for entry in adsorbates)
    components = [
        {
            "Name": "N2",
            "GasPhaseMolFraction": carrier_fraction,
            "CarrierGas": True,
            "MolecularWeight": 0.0280134,
        }
    ]
    components.extend(
        {
            "Name": entry["name"],
            **({"Water": True} if entry.get("water") else {}),
            "GasPhaseMolFraction": entry["mole_fraction"],
            "MassTransferCoefficient": entry["mass_transfer_coefficient"],
            "AxialDispersionCoefficient": 1.0e-5,
            "MolecularWeight": entry["molecular_weight"],
            "PhysisorptionSites": entry["sites"],
        }
        for entry in adsorbates
    )

    # Ruptura's mass balance uses the particle density and the void fraction; the
    # paper writes its balance with the bulk density rho_b = (1 - eps) rho_p.
    particle_density = bulk_density / (1.0 - void_fraction)
    # column.data reports the interstitial velocity, which is what
    # ColumnEntranceVelocity sets; the paper tabulates superficial velocities.
    interstitial_velocity = superficial_velocity / void_fraction

    # Whole seconds, with about OUTPUT_BLOCKS outputs over the run.
    simulation_time = float(round(end_time))
    write_interval = float(max(1, int(round(end_time)) // OUTPUT_BLOCKS))
    simulation = {
        "SimulationType": "Breakthrough",
        "BreakthroughIntegrator": "ROS34PW2",
        "InletCondition": "Danckwerts",
        "DisplayName": label,
        "Temperature": temperature,
        "InfluxTemperature": temperature,
        "ColumnVoidFraction": void_fraction,
        "ParticleDensity": particle_density,
        "ParticleDiameter": 0.0005,
        "DynamicViscosity": 1.78e-05,
        "InletPressure": total_pressure,
        "PressureGradient": 0.0,
        "ColumnEntranceVelocity": interstitial_velocity,
        "ColumnLength": column_length,
        "MixturePredictionMethod": "SPI",
        "BoundaryCondition": "FixedVelocity",
        "SimulationTime": simulation_time,
        "WriteInterval": write_interval,
        "NumberOfGridPoints": grid_points,
        "energyBalance": False,
        "Geometry": {
            "Type": "PackedBed",
            "ColumnVoidFraction": void_fraction,
            "ParticleDiameter": 0.0005,
            "InternalDiameter": internal_diameter,
            "OuterDiameter": internal_diameter * 1.2,
        },
        "Components": components,
    }
    (directory / "simulation.json").write_text(json.dumps(simulation, indent=2) + "\n")
    return {
        "directory": name,
        "kind": "breakthrough",
        "figure": figure,
        "panel": panel,
        "label": label,
        "temperature_K": temperature,
        "end_time_s": end_time,
        "bulk_density_kg_m3": bulk_density,
        "column_length_m": column_length,
        "superficial_velocity_m_s": superficial_velocity,
        "grid_points": grid_points,
        **(metadata or {}),
    }


# --------------------------------------------------------------------------
# Paper parameter tables
# --------------------------------------------------------------------------

# Table 1 flue-gas sorbent (Figure 1), reference temperature 298.15 K.
FLUE_GAS_TOTH = {
    "reference_temperature": 298.15,
    "n_s0": 0.97,
    "affinity0": 8.3e5 / ATMOSPHERE,  # 8.3e5 atm-1 -> Pa-1
    "heat": 84_000.0,
    "exponent0": 0.21,
}

# Table 2 baseline used by the Figure 2-4 parametric sweeps at 298.15 K.
BASELINE = {
    "n_s": 1.0,
    "affinity": 1.0e5 / ATMOSPHERE,  # 1e5 atm-1 -> Pa-1
    "exponent": 0.25,
    "rate_constant": 0.02,  # s-1
    "superficial_velocity": 0.008,  # m/s
    "column_length": 0.02,  # m
    "bulk_density": 500.0,  # kg/m3
    "void_fraction": 0.5,
    "temperature": 298.15,
    "co2_concentration": 1.6e-2,  # mol/m3
    "internal_diameter": 0.01,  # m
}

# Table 3 PEI/silica hollow-fiber sorbent (Figures 5-7), reference 308 K.
FIBER_TOTH = {
    "reference_temperature": 308.0,
    "n_s0": 0.81,
    "affinity0": 6.2e3 / 1000.0,  # 6.2e3 kPa-1 -> Pa-1
    "heat": 210_000.0,
    "exponent0": 0.40,
    "chi": 6.6,
    "alpha": 10.8,
}
# Equation 8: the humid enhancement exp(1.11 RH) multiplies the whole Toth loading,
# which is Ruptura's Exponential-humidity-Toth with a = 1.11 at the local RH.
HUMID_ENHANCEMENT = 1.11

FIBER_GAB = (26.0, 0.22, 0.84)  # c_m [mol/kg], c_g [-], k [-] (Figure 5a)
SWEEP_GAB = (35.0, 0.15, 0.51)  # Figure 4 water isotherm

FIBER_RADIUS = 0.002  # m
FIBER_AREA = math.pi * FIBER_RADIUS**2
FIBER_LENGTH = 0.20  # m
FIBER_VOID_FRACTION = 0.43
FIBER_CO2_MOLE_FRACTION = 380e-6
# Effective bulk densities inferred from the printed time axes; the article does
# not report the loaded module masses (see README.md).  Both are calibrated
# against Ruptura, by scaling a trial density with the ratio of the digitized to
# the simulated 50% breakthrough time:
#   Figure 6, geometric mean over all nine temperature/flow curves:  80 x 1.197
#   Figure 7, from the dry curve alone:            110 x 1.568, then refined once
# Breakthrough time is not exactly proportional to the density, because the
# linear-driving-force spreading does not scale with it, so Figure 7 took a
# second pass.  Residuals after calibration: Figure 6 geometric mean 0.993 over
# nine curves, Figure 7 dry 1.010.
# Figure 7's humid curve wants a different density than its dry curve; the dry
# one is used because it carries no humid-model assumptions.  See README.md.
FIBER_BULK_DENSITY = 95.7  # kg/m3, Figure 6 temperature/flow campaign
HUMID_FIBER_BULK_DENSITY = 160.9  # kg/m3, Figure 7 humidity campaign
FIBER_RATE_CONSTANT = {308.0: 0.040 / 60.0, 318.0: 0.048 / 60.0, 328.0: 0.100 / 60.0}
FIGURE_7_WATER_RATE_CONSTANT = 0.36 / 60.0  # s-1


def baseline_adsorbate(
    *,
    n_s: float | None = None,
    affinity: float | None = None,
    exponent: float | None = None,
    rate_constant: float | None = None,
) -> dict:
    """One CO2 adsorbate entry for the Table 2 baseline with optional overrides."""

    mole_fraction = (
        BASELINE["co2_concentration"] * GAS_CONSTANT * BASELINE["temperature"] / STANDARD_PRESSURE
    )
    return {
        "name": "CO2",
        "mole_fraction": mole_fraction,
        "molecular_weight": 0.0440095,
        "mass_transfer_coefficient": BASELINE["rate_constant"] if rate_constant is None else rate_constant,
        "sites": [
            {
                "Type": "Toth",
                "Parameters": [
                    BASELINE["n_s"] if n_s is None else n_s,
                    BASELINE["affinity"] if affinity is None else affinity,
                    BASELINE["exponent"] if exponent is None else exponent,
                ],
            }
        ],
    }


# --------------------------------------------------------------------------
# Figure 1 - flue-gas benchmark
# --------------------------------------------------------------------------


def figure_1_cases() -> list[dict]:
    entries = []
    for temperature in (298.15, 318.15):
        n_s, affinity, exponent = toth_at_temperature(temperature, **FLUE_GAS_TOTH)
        entries.append(
            write_isotherm_case(
                f"figure_1a_isotherm_{temperature:.0f}K",
                figure=1,
                panel="a",
                label=f"Liu Figure 1a - flue-gas Toth isotherm at {temperature:.0f} K",
                component="CO2",
                temperature=temperature,
                pressure_start=8.0e-5 * ATMOSPHERE,
                pressure_end=0.12 * ATMOSPHERE,
                points=220,
                scale="Log",
                sites=[{"Type": "Toth", "Parameters": [n_s, affinity, exponent]}],
                metadata={"legend": f"{temperature:.0f} K"},
            )
        )

    n_s, affinity, exponent = toth_at_temperature(298.15, **FLUE_GAS_TOTH)
    entries.append(
        write_breakthrough_case(
            "figure_1b_breakthrough_298K",
            figure=1,
            panel="b",
            label="Liu Figure 1b - dry 10132.5 Pa CO2 breakthrough at 298 K",
            temperature=298.15,
            total_pressure=ATMOSPHERE,
            # Calibrated: the article does not repeat the apparatus dimensions
            # of its reference 21, only the printed time scale.
            column_length=0.060,
            superficial_velocity=0.015,
            bulk_density=500.0,
            void_fraction=0.5,
            internal_diameter=0.01,
            end_time=1800.0,
            adsorbates=[
                {
                    "name": "CO2",
                    "mole_fraction": 0.1,
                    "molecular_weight": 0.0440095,
                    "mass_transfer_coefficient": 0.020,
                    "sites": [{"Type": "Toth", "Parameters": [n_s, affinity, exponent]}],
                }
            ],
            metadata={"legend": "298 K"},
        )
    )
    return entries


# --------------------------------------------------------------------------
# Figure 2 - kinetic and bed-operation sweeps
# --------------------------------------------------------------------------

FIGURE_2_SWEEPS = {
    "rate_constant": {
        "panel": "a",
        "values": [0.005, 0.01, 0.02, 0.04, 0.06],
        "end_time": 25_200.0,
        "unit": r"s$^{-1}$",
        "reference_key": "rate_constant_s-1",
    },
    "velocity": {
        "panel": "b",
        "values": [0.002, 0.004, 0.006, 0.008, 0.009, 0.010],
        "end_time": 114_000.0,
        "unit": r"m s$^{-1}$",
        "reference_key": "velocity_m_s",
    },
    "length": {
        "panel": "c",
        "values": [0.005, 0.010, 0.020, 0.040, 0.080],
        "end_time": 114_000.0,
        "unit": "m",
        "reference_key": "length_m",
    },
}


def figure_2_cases() -> list[dict]:
    entries = []
    for sweep, spec in FIGURE_2_SWEEPS.items():
        for value in spec["values"]:
            entries.append(
                write_breakthrough_case(
                    f"figure_2{spec['panel']}_{sweep}_{value:g}",
                    figure=2,
                    panel=spec["panel"],
                    label=f"Liu Figure 2{spec['panel']} - {sweep} {value:g}",
                    temperature=BASELINE["temperature"],
                    total_pressure=STANDARD_PRESSURE,
                    column_length=value if sweep == "length" else BASELINE["column_length"],
                    superficial_velocity=(
                        value if sweep == "velocity" else BASELINE["superficial_velocity"]
                    ),
                    bulk_density=BASELINE["bulk_density"],
                    void_fraction=BASELINE["void_fraction"],
                    internal_diameter=BASELINE["internal_diameter"],
                    end_time=spec["end_time"],
                    adsorbates=[
                        baseline_adsorbate(
                            rate_constant=value if sweep == "rate_constant" else None
                        )
                    ],
                    metadata={
                        "sweep": sweep,
                        "sweep_value": value,
                        "reference_sweep_parameter": spec["reference_key"],
                        "reference_sweep_value": value,
                        "legend": f"{latex_number(value, 4)} {spec['unit']}".strip(),
                    },
                )
            )
    return entries


# --------------------------------------------------------------------------
# Figure 3 - Toth parameter sweeps
# --------------------------------------------------------------------------

FIGURE_3_SWEEPS = {
    "saturation_loading": {
        "panel": "a",
        "values": [0.1, 0.5, 1.0, 2.0, 3.0, 4.0],
        "end_time": 132_000.0,
        "unit": r"mol kg$^{-1}$",
        "reference_key": "n_s_mol_kg",
    },
    "affinity": {
        "panel": "b",
        "values": [1.0e4, 5.0e4, 1.0e5, 2.0e5, 3.0e5],  # atm-1, converted below
        "end_time": 43_200.0,
        "unit": r"Pa$^{-1}$",
        "reference_key": "affinity_pa-1",
    },
    "exponent": {
        "panel": "c",
        "values": [0.20, 0.22, 0.25, 0.27, 0.30],
        "end_time": 36_000.0,
        "unit": "",
        "reference_key": "toth_exponent",
    },
}


def figure_3_cases() -> list[dict]:
    entries = []
    for sweep, spec in FIGURE_3_SWEEPS.items():
        for value in spec["values"]:
            affinity = value / ATMOSPHERE if sweep == "affinity" else None
            legend_value = affinity if sweep == "affinity" else value
            entries.append(
                write_breakthrough_case(
                    f"figure_3{spec['panel']}_{sweep}_{value:g}",
                    figure=3,
                    panel=spec["panel"],
                    label=f"Liu Figure 3{spec['panel']} - {sweep} {value:g}",
                    temperature=BASELINE["temperature"],
                    total_pressure=STANDARD_PRESSURE,
                    column_length=BASELINE["column_length"],
                    superficial_velocity=BASELINE["superficial_velocity"],
                    bulk_density=BASELINE["bulk_density"],
                    void_fraction=BASELINE["void_fraction"],
                    internal_diameter=BASELINE["internal_diameter"],
                    end_time=spec["end_time"],
                    adsorbates=[
                        baseline_adsorbate(
                            n_s=value if sweep == "saturation_loading" else None,
                            affinity=affinity,
                            exponent=value if sweep == "exponent" else None,
                        )
                    ],
                    metadata={
                        "sweep": sweep,
                        "sweep_value": legend_value,
                        "reference_sweep_parameter": spec["reference_key"],
                        "reference_sweep_value": legend_value,
                        "legend": f"{latex_number(legend_value)} {spec['unit']}".strip(),
                    },
                )
            )
    return entries


# --------------------------------------------------------------------------
# Figure 4 - humidity effects
# --------------------------------------------------------------------------

FIGURE_4_HUMIDITIES = [0.50, 0.40, 0.30, 0.20, 0.0]
FIGURE_4_WATER_RATE_CONSTANT = 0.03  # s-1

# Figure 4 is by far the stiffest set here.  At the feed humidity equation 8
# puts the affinity at up to 57.6 Pa-1, which against a 0.016 mol/m3 feed gives
# a partition coefficient near 1e8.  With CVODE the only affordable
# setting was 15 grid points with an absolute tolerance of 0.6% of the feed,
# about 8.5 h per case.  ROS34PW2 runs these cases in well under a second at its
# default tolerances; the grid stays at 15 points so that the discretization of
# these runs does not change.  The notebook reports the worst backward step of
# every Figure 4 curve it draws.
FIGURE_4_GRID_POINTS = 15


# Equation 8, the modified Toth law: n_s = n_s0 / (1 - gamma q_H2O) and
# b = b0 (1 + beta q_H2O), with q_H2O the GAB water loading at the local humidity.
MODIFIED_TOTH_GAMMA = 0.006  # kg/mol
MODIFIED_TOTH_BETA = 25.0  # kg/mol


def figure_4_cases() -> list[dict]:
    entries = []
    for relative_humidity in FIGURE_4_HUMIDITIES:
        name = f"figure_4_rh_{relative_humidity:.2f}".replace(".", "p")
        saturation_pressure = antoine_saturation_pressure(BASELINE["temperature"])
        water_loading = gab_loading(relative_humidity, *SWEEP_GAB)
        # Ruptura's GAB-modified Toth follows the local water loading; the dry
        # case has no water component and keeps the plain Toth site.
        site = (
            {
                "Type": "GAB-modified-Toth",
                "Parameters": [
                    BASELINE["n_s"],
                    BASELINE["affinity"],
                    BASELINE["exponent"],
                    MODIFIED_TOTH_GAMMA,
                    MODIFIED_TOTH_BETA,
                ],
            }
            if relative_humidity > 0.0
            else {"Type": "Toth", "Parameters": [BASELINE["n_s"], BASELINE["affinity"], BASELINE["exponent"]]}
        )

        adsorbates = [
            {
                "name": "CO2",
                "mole_fraction": (
                    BASELINE["co2_concentration"]
                    * GAS_CONSTANT
                    * BASELINE["temperature"]
                    / STANDARD_PRESSURE
                ),
                "molecular_weight": 0.0440095,
                "mass_transfer_coefficient": BASELINE["rate_constant"],
                "sites": [site],
            }
        ]
        if relative_humidity > 0.0:
            adsorbates.append(
                {
                    "name": "H2O",
                    "water": True,
                    "mole_fraction": relative_humidity * saturation_pressure / STANDARD_PRESSURE,
                    "molecular_weight": 0.01801528,
                    "mass_transfer_coefficient": FIGURE_4_WATER_RATE_CONSTANT,
                    "sites": [{"Type": "GAB", "Parameters": [*SWEEP_GAB]}],
                }
            )

        entries.append(
            write_breakthrough_case(
                name,
                figure=4,
                panel="all",
                label=f"Liu Figure 4 - relative humidity {relative_humidity:.2f}",
                temperature=BASELINE["temperature"],
                total_pressure=STANDARD_PRESSURE,
                column_length=BASELINE["column_length"],
                superficial_velocity=BASELINE["superficial_velocity"],
                bulk_density=BASELINE["bulk_density"],
                void_fraction=BASELINE["void_fraction"],
                internal_diameter=BASELINE["internal_diameter"],
                end_time=90_000.0,
                adsorbates=adsorbates,
                grid_points=FIGURE_4_GRID_POINTS,
                metadata={
                    "relative_humidity": relative_humidity,
                    "water_saturation_pressure_pa": saturation_pressure,
                    "feed_water_loading_mol_kg": water_loading,
                    "legend": f"{relative_humidity:.2f}",
                },
            )
        )
    return entries


# --------------------------------------------------------------------------
# Figure 5 - PEI-fiber equilibrium fits
# --------------------------------------------------------------------------

FIGURE_5_CO2_CASES = [
    (308.0, 0.0, "308 K, dry"),
    (318.0, 0.0, "318 K, dry"),
    (328.0, 0.0, "328 K, dry"),
    (308.0, 0.85, "308 K, 0.85 RH"),
]


def figure_5_cases() -> list[dict]:
    saturation_pressure = antoine_saturation_pressure(308.0)
    # The water isotherm comes from a short breakthrough into a dry fiber module:
    # every grid point reports its local water partial pressure alongside
    # Ruptura's own GAB equilibrium loading at that pressure, and the developing
    # front sweeps the whole 0 to 0.85 humidity range within minutes.  (Until
    # 2026-11-04 a MixturePrediction run could not write a GAB site, TODO B-IO-7.)
    entries = [
        write_breakthrough_case(
            "figure_5a_water_isotherm_308K",
            figure=5,
            panel="a",
            label="Liu Figure 5a - GAB water isotherm at 308 K",
            temperature=308.0,
            total_pressure=STANDARD_PRESSURE,
            column_length=FIBER_LENGTH,
            superficial_velocity=0.05 * FIBER_VOID_FRACTION,
            bulk_density=FIBER_BULK_DENSITY,
            void_fraction=FIBER_VOID_FRACTION,
            internal_diameter=2.0 * FIBER_RADIUS,
            end_time=600.0,
            adsorbates=[
                {
                    "name": "H2O",
                    "mole_fraction": 0.85 * saturation_pressure / STANDARD_PRESSURE,
                    "molecular_weight": 0.01801528,
                    "mass_transfer_coefficient": FIGURE_7_WATER_RATE_CONSTANT,
                    "sites": [{"Type": "GAB", "Parameters": [*FIBER_GAB]}],
                }
            ],
            metadata={
                "kind_override": "isotherm_from_breakthrough",
                "component": "H2O",
                "water_saturation_pressure_pa": saturation_pressure,
                "legend": "308 K",
            },
        )
    ]
    for temperature, relative_humidity, legend in FIGURE_5_CO2_CASES:
        n_s, affinity, exponent = toth_at_temperature(temperature, **FIBER_TOTH)
        # The humid isotherm is equation 8's exp(1.11 RH) at a fixed RelativeHumidity.
        site = (
            {"Type": "Exponential-humidity-Toth", "Parameters": [n_s, affinity, exponent, HUMID_ENHANCEMENT]}
            if relative_humidity > 0.0
            else {"Type": "Toth", "Parameters": [n_s, affinity, exponent]}
        )
        name = f"figure_5b_co2_isotherm_{temperature:.0f}K_rh_{relative_humidity:.2f}".replace(".", "p")
        entries.append(
            write_isotherm_case(
                name,
                figure=5,
                panel="b",
                label=f"Liu Figure 5b - CO2 Toth isotherm, {legend}",
                component="CO2",
                temperature=temperature,
                pressure_start=0.5,
                pressure_end=100.0,
                points=240,
                sites=[site],
                relative_humidity=relative_humidity if relative_humidity > 0.0 else None,
                metadata={"relative_humidity": relative_humidity, "legend": legend},
            )
        )
    return entries


# --------------------------------------------------------------------------
# Figures 6 and 7 - PEI-fiber DAC breakthrough
# --------------------------------------------------------------------------

FIGURE_6_TEMPERATURES = [308.0, 318.0, 328.0]
FIGURE_6_FLOWS = [6.6666666667e-7, 1.0e-6, 1.5e-6]  # standard m3/s
FIGURE_7_FLOW = 3.3333333333e-6  # standard m3/s
FIGURE_7_RELATIVE_HUMIDITY = 0.85
FIGURE_7_CO2_RATE_CONSTANT = 0.020 / 60.0  # s-1


def fiber_adsorbate(
    temperature: float,
    *,
    humid: bool,
    rate_constant: float,
) -> dict:
    n_s, affinity, exponent = toth_at_temperature(temperature, **FIBER_TOTH)
    # Humid: equation 8's exp(1.11 RH) at the local humidity of the water component.
    site = (
        {"Type": "Exponential-humidity-Toth", "Parameters": [n_s, affinity, exponent, HUMID_ENHANCEMENT]}
        if humid
        else {"Type": "Toth", "Parameters": [n_s, affinity, exponent]}
    )
    return {
        "name": "CO2",
        "mole_fraction": FIBER_CO2_MOLE_FRACTION,
        "molecular_weight": 0.0440095,
        "mass_transfer_coefficient": rate_constant,
        "sites": [site],
    }


def figure_6_cases() -> list[dict]:
    entries = []
    for temperature in FIGURE_6_TEMPERATURES:
        panel = "abc"[FIGURE_6_TEMPERATURES.index(temperature)]
        for flow in FIGURE_6_FLOWS:
            name = f"figure_6{panel}_{temperature:.0f}K_flow_{flow:.3g}".replace(".", "p").replace("-", "m")
            entries.append(
                write_breakthrough_case(
                    name,
                    figure=6,
                    panel=panel,
                    label=f"Liu Figure 6{panel} - {temperature:.0f} K, {flow:.3g} m3/s",
                    temperature=temperature,
                    total_pressure=STANDARD_PRESSURE,
                    column_length=FIBER_LENGTH,
                    superficial_velocity=standard_flow_to_superficial_velocity(
                        flow, FIBER_AREA, temperature
                    ),
                    bulk_density=FIBER_BULK_DENSITY,
                    void_fraction=FIBER_VOID_FRACTION,
                    internal_diameter=2.0 * FIBER_RADIUS,
                    end_time=18_000.0,
                    adsorbates=[
                        fiber_adsorbate(
                            temperature,
                            humid=False,
                            rate_constant=FIBER_RATE_CONSTANT[temperature],
                        )
                    ],
                    metadata={
                        "standard_flow_m3_s": flow,
                        "legend": rf"{latex_number(flow)} m$^3$ s$^{{-1}}$",
                    },
                )
            )
    return entries


def figure_7_cases() -> list[dict]:
    temperature = 308.0
    velocity = standard_flow_to_superficial_velocity(FIGURE_7_FLOW, FIBER_AREA, temperature)
    saturation_pressure = antoine_saturation_pressure(temperature)
    common = dict(
        figure=7,
        panel="a",
        temperature=temperature,
        total_pressure=STANDARD_PRESSURE,
        column_length=FIBER_LENGTH,
        superficial_velocity=velocity,
        bulk_density=HUMID_FIBER_BULK_DENSITY,
        void_fraction=FIBER_VOID_FRACTION,
        internal_diameter=2.0 * FIBER_RADIUS,
        end_time=18_000.0,
    )
    entries = [
        write_breakthrough_case(
            "figure_7_dry_308K",
            label="Liu Figure 7 - dry, 308 K",
            adsorbates=[
                fiber_adsorbate(
                    temperature,
                    humid=False,
                    rate_constant=FIBER_RATE_CONSTANT[temperature],
                )
            ],
            metadata={"condition": "dry", "legend": "Dry"},
            **common,
        ),
        write_breakthrough_case(
            "figure_7_rh_0p85_308K",
            label="Liu Figure 7 - 0.85 relative humidity, 308 K",
            adsorbates=[
                fiber_adsorbate(
                    temperature,
                    humid=True,
                    rate_constant=FIGURE_7_CO2_RATE_CONSTANT,
                ),
                {
                    "name": "H2O",
                    "water": True,
                    "mole_fraction": FIGURE_7_RELATIVE_HUMIDITY
                    * saturation_pressure
                    / STANDARD_PRESSURE,
                    "molecular_weight": 0.01801528,
                    "mass_transfer_coefficient": FIGURE_7_WATER_RATE_CONSTANT,
                    "sites": [{"Type": "GAB", "Parameters": [*FIBER_GAB]}],
                },
            ],
            metadata={
                "condition": "humid_85pct",
                "relative_humidity": FIGURE_7_RELATIVE_HUMIDITY,
                "water_saturation_pressure_pa": saturation_pressure,
                "legend": "0.85 RH",
            },
            **common,
        ),
    ]
    return entries


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete simulations/ first, discarding any Ruptura output already there",
    )
    arguments = parser.parse_args()

    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)
    SIMULATION_DIR.mkdir(parents=True, exist_ok=True)

    index = []
    for builder in (
        figure_1_cases,
        figure_2_cases,
        figure_3_cases,
        figure_4_cases,
        figure_5_cases,
        figure_6_cases,
        figure_7_cases,
    ):
        index.extend(builder())

    (SIMULATION_DIR / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")


if __name__ == "__main__":
    main()
