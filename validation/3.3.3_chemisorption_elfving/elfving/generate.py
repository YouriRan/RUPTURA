#!/usr/bin/env python3
"""Generate Ruptura inputs and equilibrium reproductions for the Elfving cases.

This writes two kinds of output:

- `run.py`-executable fixed-bed breakthrough cases, one directory per case,
  named after the model/temperature/humidity/flow combination it represents.
- Standalone equilibrium reproductions (CSV + SVG) evaluated directly from
  the paper's equations, independent of the Ruptura binary.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from paper_parameters import (
    GAB_PARAMETERS,
    MODELS,
    coadsorption_loading,
    effective_co2_heat,
    empty_column_water_loading,
    fit_toth,
    gab_loading,
    toth_loading,
    toth_supplement_parameters,
    water_saturation_pressure_bar,
)


ROOT = Path(__file__).resolve().parent
WATER_LDF_CAP = 0.25


# --------------------------------------------------------------------------
# Fixed-bed breakthrough cases (runs 1-13)
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Case:
    run: int
    model: str
    temperature: float
    water_fraction: float
    flow_ml_min: float
    nominal_humidity: str  # paper's nominal H2O setpoint label: "0p2", "1", or "2"
    kf1: float | None
    kf2: float | None
    water_ldf: float
    heat_conductivity: float
    wall_heat_transfer: float
    co2_fraction: float = 0.00037

    @property
    def directory_name(self) -> str:
        temperature_c = round(self.temperature - 273.15)
        return f"{self.model}_{temperature_c}c_{self.nominal_humidity}pct_{self.flow_ml_min:.0f}mlmin"


CASES = (
    Case(1, "five", 298.15, 0.0016, 500.0, "0p2", 1.37, 1.71, 0.99, 1.74, 17.78),
    Case(2, "five", 299.15, 0.0095, 500.0, "1", 1.48, 12.67, 0.38, 2.12, 10.51),
    Case(3, "five", 298.15, 0.0190, 500.0, "2", 0.92, 58.74, 0.22, 2.10, 17.34),
    Case(4, "five", 308.15, 0.0095, 500.0, "1", 1.35, 2.06, 0.26, 1.96, 19.65),
    Case(5, "five", 308.15, 0.0190, 500.0, "2", 1.78, 14.54, 0.88, 1.96, 36.25),
    Case(6, "seven_bounded", 298.15, 0.0016, 500.0, "0p2", 0.35, 2.33, 1.22, 1.96, 11.87),
    Case(7, "seven_bounded", 299.15, 0.0095, 500.0, "1", 0.32, 4.94, 0.39, 2.00, 10.13),
    Case(8, "seven_bounded", 298.15, 0.0190, 500.0, "2", 0.20, 8.53, 0.20, 2.10, 18.90),
    Case(9, "seven_bounded", 308.15, 0.0190, 500.0, "2", 0.45, 3.57, 1.60, 3.12, 33.12),
    Case(10, "five", 298.15, 0.0200, 100.0, "2", 3.81, 30.25, 0.11, 2.03, 18.13, 0.000401),
    Case(11, "seven_bounded", 298.15, 0.0200, 100.0, "2", 0.69, 5.03, 0.10, 2.10, 18.59, 0.000401),
    Case(12, "toth_ldf", 298.15, 0.0190, 500.0, "2", None, None, 0.29, 1.86, 16.67),
    Case(13, "toth_ldf", 308.15, 0.0190, 500.0, "2", None, None, 0.20, 1.89, 30.13),
)


def pressure_grid() -> list[float]:
    return [0.00005 + index * (0.0042 - 0.00005) / 79.0 for index in range(80)]


@lru_cache(maxsize=None)
def projected_co2_isotherm(model_name: str, temperature: float, p_h2o_bar: float) -> tuple[float, float, float, float]:
    model = MODELS[model_name]
    pressures = pressure_grid()
    loadings = [sum(coadsorption_loading(model, temperature, p, p_h2o_bar)) for p in pressures]
    return fit_toth(pressures, loadings)


def co2_parameters(case: Case, total_pressure_bar: float) -> tuple[float, float, float, float, float]:
    if case.model == "toth_ldf":
        qmax, affinity_bar, exponent = toth_supplement_parameters(case.temperature)
        return qmax, affinity_bar, exponent, 0.0012 if case.temperature < 303.0 else 0.0018, 112_800.0

    model = MODELS[case.model]
    p_h2o_bar = case.water_fraction * total_pressure_bar
    qmax, affinity_bar, exponent, fit_error = projected_co2_isotherm(case.model, case.temperature, p_h2o_bar)
    p_co2_bar = case.co2_fraction * total_pressure_bar
    q1, q2 = coadsorption_loading(model, case.temperature, p_co2_bar, p_h2o_bar)
    initial_rate = (
        case.kf1 * model.qm**model.t1 * p_co2_bar
        + case.kf2 * model.qm**model.t2 * p_co2_bar * p_h2o_bar
    )
    effective_ldf = initial_rate / (q1 + q2)
    heat = effective_co2_heat(model, q1, q2)
    if fit_error > 0.025:
        raise RuntimeError(f"Toth projection error too large for run {case.run}: {fit_error:.3%}")
    return qmax, affinity_bar, exponent, effective_ldf, heat


def general_chemisorption_site(
    isotherm: dict,
    maximum_loading: float,
    heat: float,
    adsorption_rate_coefficient: float,
    capacity_order: int,
    desorption_rate_coefficient: float = 0.0,
) -> dict:
    """Return a native General chemisorption site.

    Ruptura evaluates

        dq/dt = ka * c * (qeq - q)**capacity_order - kd * q.

    This is the closest current-interface approximation to the paper's free
    amine factors without introducing a shared surface-state implementation.
    """
    return {
        "Type": isotherm["Type"],
        "Parameters": isotherm["Parameters"],
        "RateEquation": "General",
        "RateEquationParameters": {
            "AdsorptionRateCoefficient": adsorption_rate_coefficient,
            "DesorptionRateCoefficient": desorption_rate_coefficient,
            "PoreConcentrationOrder": 1,
            "CapacityOrder": capacity_order,
            "DesorptionOrder": 1,
            "FilmMassTransferCoefficient": 0.0,
            "PoreDiffusivity": 0.0,
            "UsePoreSurfaceTransport": False,
        },
        "HeatOfAdsorption": heat,
        "MaximumLoading": maximum_loading,
    }


def first_order_chemisorption_site(isotherm: dict, maximum_loading: float, heat: float, rate: float) -> dict:
    """Return a native LDF chemisorption site, dq/dt = k * (qeq - q)."""
    return {
        "Type": isotherm["Type"],
        "Parameters": isotherm["Parameters"],
        "RateEquationParameters": {"MassTransferCoefficient": rate},
        "HeatOfAdsorption": heat,
        "MaximumLoading": maximum_loading,
    }


def water_chemisorption_sites(case: Case, pressure_pa: float) -> tuple[list[dict], float]:
    """Map the GAB and empty-column water capacities to linear free-site kinetics."""
    rate = min(case.water_ldf, WATER_LDF_CAP)
    pressure_bar = pressure_pa / 100_000.0
    gab_feed_loading = gab_loading(case.temperature, case.water_fraction * pressure_bar)

    # C(T) = C_0 exp(dH_C/(R T)) and K(T) = K_0 exp(dH_K/(R T)) are uncentred VantHoff laws; GAB evaluates its
    # saturation pressure at the column temperature.
    qm_mono, c0, k0, delta_h_c, delta_h_k = GAB_PARAMETERS
    gab_parameters = [qm_mono, *({"Value": value, "Law": "VantHoff", "Energy": energy} if energy else value
                                 for value, energy in ((c0, delta_h_c), (k0, delta_h_k)))]
    main_site = first_order_chemisorption_site(
        {"Type": "GAB", "Parameters": gab_parameters},
        maximum_loading=max(GAB_PARAMETERS[0], 1.10 * gab_feed_loading),
        heat=50_700.0,
        rate=rate,
    )
    empty_site = first_order_chemisorption_site(
        {"Type": "Langmuir", "Parameters": [1.1216, 67.1217 / 100_000.0]},
        maximum_loading=1.1216,
        heat=50_700.0,
        rate=rate,
    )
    return [main_site, empty_site], rate


def case_config(case: Case) -> dict:
    pressure_pa = 105_000.0
    pressure_bar = pressure_pa / 100_000.0
    area = math.pi * 0.009**2 / 4.0
    volumetric_flow = case.flow_ml_min * 1.0e-6 / 60.0
    interstitial_velocity = volumetric_flow / area / 0.375
    low_flow = case.flow_ml_min < 200.0
    axial_dispersion = 7.44e-5 if low_flow else (3.52e-4 if abs(case.temperature - 308.15) < 0.2 else 3.50e-4)
    qmax, affinity_bar, exponent, co2_ldf, co2_heat = co2_parameters(case, pressure_bar)
    co2_feed_concentration = case.co2_fraction * pressure_pa / (8.31446261815324 * case.temperature)
    if case.model in MODELS:
        model = MODELS[case.model]
        co2_feed_loading = sum(
            coadsorption_loading(
                model,
                case.temperature,
                case.co2_fraction * pressure_bar,
                case.water_fraction * pressure_bar,
            )
        )
        # A CO2 molecule consuming two amine sites can occupy at most qm/2.
        co2_capacity_loading = 0.5 * model.qm
    else:
        co2_feed_loading = toth_loading(qmax, affinity_bar, exponent, case.co2_fraction * pressure_bar)
        co2_capacity_loading = qmax

    initial_co2_rate = co2_ldf * co2_feed_loading
    co2_adsorption_rate_coefficient = initial_co2_rate / (
        co2_feed_concentration * co2_capacity_loading**2
    )
    co2_desorption_rate_coefficient = (
        initial_co2_rate * (co2_capacity_loading - co2_feed_loading) ** 2
        / (co2_capacity_loading**2 * co2_feed_loading)
    )
    water_sites, water_rate_coefficient = water_chemisorption_sites(case, pressure_pa)
    duration_seconds = 16.0 * 3600.0 if low_flow else 5.0 * 3600.0
    # RK3 heat conduction is stable below about 0.6 dz^2 rho c_p / lambda: 1e-4 s up to 2.4 W/(m K), so run 9
    # (3.12 W/(m K), limit 7.9e-5 s) takes half the step (TODO B-ITG-24).
    time_step = 5.0e-5 if case.heat_conductivity > 2.4 else 1.0e-4
    # Write cadence carried over from the previous 20 s output interval,
    # so real-time (s) cadence is unchanged by the finer TimeStep.
    write_interval_seconds = 200.0 if low_flow else 40.0
    viscosity = 1.82e-5 if case.temperature > 303.0 else 1.78e-5
    cp_gas = 1060.0 if case.temperature > 303.0 else 1050.0
    model_label = {
        "five": "5-parameter projection",
        "seven_bounded": "7-parameter projection (t1=t2=3)",
        "toth_ldf": "supplementary Toth+LDF",
    }[case.model]

    return {
        "SimulationType": "Breakthrough",
        "DisplayName": f"Elfving 2021 - Run {case.run:02d} - {model_label}",
        "Temperature": case.temperature,
        "InfluxTemperature": case.temperature,
        "ColumnVoidFraction": 0.375,
        "DynamicViscosity": viscosity,
        "ParticleDiameter": 6.0e-4,
        "ParticleDensity": 720.0,
        "InletPressure": pressure_pa,
        "PressureGradient": 0.0,
        "ColumnEntranceVelocity": interstitial_velocity,
        "ColumnLength": 0.01746,
        "MixturePredictionMethod": "SPI",
        "BreakthroughIntegrator": "RungeKutta3",
        # "CVODELinearSolver": "Dense",
        "BoundaryCondition": "FixedVelocity",
        "SimulationTime": duration_seconds,
        "WriteInterval": write_interval_seconds,
        "TimeStep": time_step,
        "NumberOfGridPoints": 15 if low_flow else 30,
        "energyBalance": True,
        "heatCapacityGas": cp_gas,
        "heatCapacitySolid": 1580.0,
        "heatCapacityWall": 500.0,
        "wallDensity": 8000.0,
        "gasThermalConductivity": case.heat_conductivity,
        "wallThermalConductivity": 16.0,
        "heatTransferGasSolid": 250.0,
        "heatTransferGasWall": case.wall_heat_transfer,
        "heatTransferWallExternal": 1000.0,
        "Geometry": {
            "Type": "PackedBed",
            "ColumnVoidFraction": 0.375,
            "ParticleDiameter": 6.0e-4,
            "InternalDiameter": 0.009,
            "OuterDiameter": 0.010,
        },
        "Components": [
            {
                "Name": "N2",
                "GasPhaseMolFraction": 1.0 - case.co2_fraction - case.water_fraction,
                "CarrierGas": True,
                "MolecularWeight": 0.0280134,
            },
            {
                "Name": "CO2",
                "GasPhaseMolFraction": case.co2_fraction,
                "MassTransferCoefficient": co2_ldf,
                "AxialDispersionCoefficient": axial_dispersion,
                "MolecularWeight": 0.0440095,
                # The column input currently requires every non-carrier to
                # declare a physisorption model. A zero-Henry placeholder
                # leaves all uptake in the chemisorption state.
                "PhysisorptionSites": [{"Type": "Henry", "Parameters": [0.0]}],
                "ChemisorptionSites": [
                    general_chemisorption_site(
                        # A strongly saturated Langmuir term supplies the
                        # approximately constant qm/2 free-site ceiling used
                        # by the paper's dry two-amine pathway. The kinetic
                        # reverse term sets the actual feed equilibrium.
                        {"Type": "Langmuir", "Parameters": [co2_capacity_loading, 1.0e3]},
                        maximum_loading=co2_capacity_loading,
                        heat=co2_heat,
                        adsorption_rate_coefficient=co2_adsorption_rate_coefficient,
                        capacity_order=2,
                        desorption_rate_coefficient=co2_desorption_rate_coefficient,
                    )
                ],
            },
            {
                "Name": "H2O",
                "GasPhaseMolFraction": case.water_fraction,
                "MassTransferCoefficient": min(case.water_ldf, WATER_LDF_CAP),
                "AxialDispersionCoefficient": axial_dispersion,
                "MolecularWeight": 0.01801528,
                "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [0.0, 1.0]}],
                "ChemisorptionSites": water_sites,
            },
        ],
        "PaperMapping": {
            "model": case.model,
            "flow_ml_min": case.flow_ml_min,
            "water_vol_percent": 100.0 * case.water_fraction,
            "nominal_water_percent": case.nominal_humidity.replace("p", "."),
            "co2_ppm": 1.0e6 * case.co2_fraction,
            "co2_toth_projection_relative_rms": 0.0
            if case.model == "toth_ldf"
            else projected_co2_isotherm(case.model, case.temperature, case.water_fraction * pressure_bar)[3],
            "gab_feed_loading_mol_kg": gab_loading(case.temperature, case.water_fraction * pressure_bar),
            "empty_column_feed_loading_mol_kg": empty_column_water_loading(case.water_fraction * pressure_bar),
            "paper_water_ldf_per_s": case.water_ldf,
            "ruptura_water_ldf_per_s": min(case.water_ldf, WATER_LDF_CAP),
            "co2_capacity_order": 2,
            "water_capacity_order": 1,
            "co2_adsorption_rate_coefficient": co2_adsorption_rate_coefficient,
            "co2_desorption_rate_coefficient": co2_desorption_rate_coefficient,
            "co2_feed_loading_mol_kg": co2_feed_loading,
            "co2_two_site_capacity_mol_kg": co2_capacity_loading,
            "water_rate_coefficient_per_s": water_rate_coefficient,
            "kinetic_mapping": "CO2 reversible General chemisorption with qm/2 capacity and squared free-site term; H2O FirstOrder chemisorption",
        },
    }


def write_case_json(path: Path, value: dict) -> None:
    # PaperMapping documents provenance but is not a Ruptura input key.
    mapping = value.pop("PaperMapping")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")
    (path.parent / "paper_mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")


def generate_breakthrough_cases() -> None:
    for case in CASES:
        case_dir = ROOT / case.directory_name
        config = case_config(case)
        write_case_json(case_dir / "simulation.json", config)
    print(f"Wrote {len(CASES)} case directories under {ROOT}")


# --------------------------------------------------------------------------
# Equilibrium reproductions (Figs. 2-4, S3, S5) - no Ruptura binary required
# --------------------------------------------------------------------------


EQUILIBRIUM_OUTPUT = ROOT / "equilibrium"
EQUILIBRIUM_COLORS = {298.15: "#169c31", 308.15: "#d1a700", 323.15: "#d8322b"}


def write_csv(name: str, header: list[str], rows: list[list[object]]) -> None:
    EQUILIBRIUM_OUTPUT.mkdir(parents=True, exist_ok=True)
    with (EQUILIBRIUM_OUTPUT / name).open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def svg_path(points: list[tuple[float, float]], x0: float, y0: float, width: float, height: float,
             xmax: float, ymax: float) -> str:
    values = []
    for x_value, y_value in points:
        x = x0 + width * x_value / xmax
        y = y0 + height * (1.0 - y_value / ymax)
        values.append(("M" if not values else "L") + f"{x:.2f},{y:.2f}")
    return " ".join(values)


def equilibrium_panel(x0: float, y0: float, width: float, height: float, title: str, x_label: str, y_label: str) -> str:
    return f'''<g><rect x="{x0}" y="{y0}" width="{width}" height="{height}" fill="white" stroke="#444"/>
<text x="{x0}" y="{y0-14}" class="panel-title">{title}</text>
<text x="{x0+width/2}" y="{y0+height+42}" text-anchor="middle" class="label">{x_label}</text>
<text x="{x0-52}" y="{y0+height/2}" transform="rotate(-90 {x0-52} {y0+height/2})" text-anchor="middle" class="label">{y_label}</text></g>'''


def generate_equilibrium_reproductions() -> None:
    temperatures = (298.15, 308.15, 323.15)
    humidities = (0.00165, 0.0100, 0.0195)
    co2_pressures = [index * 0.0042 / 100.0 for index in range(101)]

    gab_rows: list[list[object]] = []
    gab_curves: dict[float, list[tuple[float, float]]] = {}
    for temperature in temperatures:
        saturation = water_saturation_pressure_bar(temperature)
        points = []
        for index in range(1, 81):
            relative_pressure = index / 100.0
            pressure = relative_pressure * saturation
            loading = gab_loading(temperature, pressure)
            gab_rows.append([temperature, relative_pressure, pressure, loading])
            points.append((relative_pressure, loading))
        gab_curves[temperature] = points
    write_csv("figure2_gab.csv", ["temperature_K", "relative_pressure", "p_h2o_bar", "q_h2o_mol_kg"], gab_rows)

    model_rows: list[list[object]] = []
    model_curves: dict[tuple[str, float, float], list[tuple[float, float]]] = {}
    for model_name in ("five", "seven_bounded", "seven_relaxed"):
        model = MODELS[model_name]
        for temperature in temperatures:
            for humidity in humidities:
                points = []
                for pressure in co2_pressures:
                    q1, q2 = coadsorption_loading(model, temperature, pressure, humidity)
                    total = q1 + q2
                    model_rows.append([model_name, temperature, humidity, pressure, q1, q2, total])
                    points.append((pressure, total))
                model_curves[(model_name, temperature, humidity)] = points
    write_csv(
        "figures3_4_s3_coadsorption.csv",
        ["model", "temperature_K", "p_h2o_bar", "p_co2_bar", "q1_mol_kg", "q2_mol_kg", "q_total_mol_kg"],
        model_rows,
    )

    toth_rows: list[list[object]] = []
    toth_curves: dict[float, list[tuple[float, float]]] = {}
    for temperature in temperatures:
        qmax, affinity, exponent = toth_supplement_parameters(temperature)
        points = []
        for pressure in co2_pressures:
            loading = toth_loading(qmax, affinity, exponent, pressure)
            toth_rows.append([temperature, pressure, loading, qmax, affinity, exponent])
            points.append((pressure, loading))
        toth_curves[temperature] = points
    write_csv(
        "figure_s5_toth.csv",
        ["temperature_K", "p_co2_bar", "q_co2_mol_kg", "qmax_mol_kg", "b_per_bar", "toth_exponent"],
        toth_rows,
    )

    width, height = 1220, 920
    panel_width, panel_height = 470, 300
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#f7f8fa"/><style>text{{font-family:Arial,sans-serif;fill:#222}}.title{{font-size:24px;font-weight:700}}.panel-title{{font-size:18px;font-weight:700}}.label{{font-size:14px}}.legend{{font-size:12px}}</style>
<text x="610" y="36" text-anchor="middle" class="title">Elfving and Sainio (2021) equilibrium-model reproduction</text>''']
    specs = [
        (80, 95, "Figure 2 - H2O GAB", "Relative pressure, p/psat", "q_H2O [mol/kg]"),
        (675, 95, "Figure 3 - 5-parameter CO2", "p_CO2 [bar]", "q_CO2 [mol/kg]"),
        (80, 545, "Figure S3 - bounded 7-parameter CO2", "p_CO2 [bar]", "q_CO2 [mol/kg]"),
        (675, 545, "Figure S5 - Toth at 2 vol-% H2O", "p_CO2 [bar]", "q_CO2 [mol/kg]"),
    ]
    for item in specs:
        parts.append(equilibrium_panel(item[0], item[1], panel_width, panel_height, item[2], item[3], item[4]))
    # Compact endpoint labels keep the four panels readable without a plotting
    # dependency while still making their numerical ranges explicit.
    parts.extend([
        '<text x="80" y="414" class="legend">0</text><text x="550" y="414" text-anchor="end" class="legend">0.8</text><text x="70" y="99" text-anchor="end" class="legend">14</text>',
        '<text x="675" y="414" class="legend">0</text><text x="1145" y="414" text-anchor="end" class="legend">0.0042</text><text x="665" y="99" text-anchor="end" class="legend">1.6</text>',
        '<text x="80" y="864" class="legend">0</text><text x="550" y="864" text-anchor="end" class="legend">0.0042</text><text x="70" y="549" text-anchor="end" class="legend">1.6</text>',
        '<text x="675" y="864" class="legend">0</text><text x="1145" y="864" text-anchor="end" class="legend">0.0042</text><text x="665" y="549" text-anchor="end" class="legend">1.6</text>',
    ])

    for temperature, points in gab_curves.items():
        parts.append(f'<path d="{svg_path(points,80,95,panel_width,panel_height,0.8,14.0)}" fill="none" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="3"/>')
    for temperature in temperatures:
        for humidity, dash in zip(humidities, ("", "8 5", "2 4")):
            points = model_curves[("five", temperature, humidity)]
            parts.append(f'<path d="{svg_path(points,675,95,panel_width,panel_height,0.0042,1.6)}" fill="none" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="2.4" stroke-dasharray="{dash}"/>')
            points = model_curves[("seven_bounded", temperature, humidity)]
            parts.append(f'<path d="{svg_path(points,80,545,panel_width,panel_height,0.0042,1.6)}" fill="none" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="2.4" stroke-dasharray="{dash}"/>')
    for temperature, points in toth_curves.items():
        parts.append(f'<path d="{svg_path(points,675,545,panel_width,panel_height,0.0042,1.6)}" fill="none" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="3"/>')

    for index, temperature in enumerate(temperatures):
        label = f"{temperature-273.15:.0f} C"
        y = 120 + 20 * index
        parts.append(f'<line x1="430" y1="{y}" x2="460" y2="{y}" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="3"/><text x="468" y="{y+4}" class="legend">{label}</text>')
    parts.append('<text x="675" y="500" class="legend">CO2 line style: solid 0.2%, dashed 1%, dotted 2% H2O</text>')
    parts.append('</svg>')
    (EQUILIBRIUM_OUTPUT / "equilibrium_reproduction.svg").write_text("\n".join(parts))

    relaxed = ['''<svg xmlns="http://www.w3.org/2000/svg" width="680" height="520" viewBox="0 0 680 520">
<rect width="100%" height="100%" fill="#f7f8fa"/><style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:22px;font-weight:700}.label{font-size:14px}.legend{font-size:12px}</style>
<text x="340" y="32" text-anchor="middle" class="title">Figure 4 - relaxed 7-parameter CO2 model</text>
<rect x="90" y="75" width="510" height="340" fill="white" stroke="#444"/>
<text x="345" y="460" text-anchor="middle" class="label">p_CO2 [bar]</text>
<text x="30" y="245" transform="rotate(-90 30 245)" text-anchor="middle" class="label">q_CO2 [mol/kg]</text>
<text x="90" y="432" class="legend">0</text><text x="600" y="432" text-anchor="end" class="legend">0.0042</text><text x="80" y="79" text-anchor="end" class="legend">1.6</text>''']
    for temperature in temperatures:
        for humidity, dash in zip(humidities, ("", "8 5", "2 4")):
            points = model_curves[("seven_relaxed", temperature, humidity)]
            relaxed.append(f'<path d="{svg_path(points,90,75,510,340,0.0042,1.6)}" fill="none" stroke="{EQUILIBRIUM_COLORS[temperature]}" stroke-width="2.4" stroke-dasharray="{dash}"/>')
    relaxed.append('<text x="90" y="490" class="legend">Solid 0.2%, dashed 1%, dotted 2% H2O; green 25 C, gold 35 C, red 50 C</text></svg>')
    (EQUILIBRIUM_OUTPUT / "figure4_relaxed_reproduction.svg").write_text("\n".join(relaxed))
    print(f"Wrote equilibrium reproductions to {EQUILIBRIUM_OUTPUT}")


def main() -> None:
    generate_breakthrough_cases()
    generate_equilibrium_reproductions()


if __name__ == "__main__":
    main()
