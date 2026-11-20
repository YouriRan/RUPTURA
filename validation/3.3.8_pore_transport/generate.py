#!/usr/bin/env python3
"""Generate the pore-transport validation: CO2 from air on PEI/fumed silica, four particle sizes.

Breakthrough data of Goeppert et al., ChemSusChem 7 (2014) 1386, with the dual-porosity model and fitted parameters of
Auton et al., Appl. Math. Model. 130 (2024) 827 (Tables 5 and 8; arXiv:2311.10161). Their kinetics in the pores,
dq/dt = k+ c_p (q_max - q)^2 - k- q^2, are the General rate equation with orders 1, 2, 2. Their single exchange
coefficient k_p is split over the film and pore resistances of RUPTURA, 1/k_p = 1/k_f + R/(5 D_p). See conversions.md
for every number and README.md for the studies:

- grid:       the column length that closes Auton's mass balance (L_fit), no axial dispersion, 50 to 800 cells;
- split:      the resistance 1/k_p as film only, pore only, or half each (L_fit, no dispersion, 100 cells);
- beta:       the smallest size, whose published k_p means no transport limit, at beta = 100 and 1000 against direct uptake;
- affinity:   the site affinity that turns the site isotherm into the capacity q_max, 1e3 against 1e6 1/Pa;
- wave:       a column of 70 L_ref (Auton's length scale), 14 profiles, for the convergence to their travelling wave
              (Eq. 62);
- experiment: the published column (0.1 m) and L_fit, both with the axial dispersion of Table 5, against the data.

    python generate.py            # write simulations/ and simulations/index.json
    python generate.py --clean    # delete simulations/ first, including any output

All quantities are SI: Pa, K, s, m, mol/kg.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SIMULATION_DIR = ROOT / "simulations"

GAS_CONSTANT = 8.31446261815324  # J/(mol K), src/utils.h

TEMPERATURE = 298.15  # K
PRESSURE = 101325.0  # Pa
FEED_MOLE_FRACTION = 405.0e-6  # 405 ppm
VOLUMETRIC_FLOW = 335.0e-6 / 60.0  # m3/s, 335 mL/min
SAMPLE_MASS = 2.72e-3  # kg
COLUMN_DIAMETER = 8.0e-3  # m
PUBLISHED_LENGTH = 0.1  # m
VOID_FRACTION = 0.48
VISCOSITY = 1.8e-5  # Pa s; FixedVelocity computes no pressure drop

AREA = math.pi * COLUMN_DIAMETER**2 / 4.0  # m2
VELOCITY = VOLUMETRIC_FLOW / (VOID_FRACTION * AREA)  # m/s, interstitial; Table 5: 0.2314
BULK_DENSITY = SAMPLE_MASS / (AREA * PUBLISHED_LENGTH)  # kg/m3; Table 5: 541.127
PARTICLE_DENSITY = BULK_DENSITY / (1.0 - VOID_FRACTION)  # kg/m3, solid mass per solid volume in RUPTURA
FEED_CONCENTRATION = FEED_MOLE_FRACTION * PRESSURE / (GAS_CONSTANT * TEMPERATURE)  # mol/m3; Table 5: 0.0166

# Tables 5 and 8 of Auton et al.: R in m, D_ax in m2/s, half time in h, q_max in mol/kg, k+ in m3 kg/(mol2 s),
# k- in kg/(mol s), k_p in m/s.
SIZES = {
    "R1p25e-4m": {"sieve": "<0.25 mm", "radius": 1.25e-4, "dispersion": 0.6e-4, "half_time_h": 10.730,
                  "maximum_loading": 1.662, "k_plus": 0.0348, "k_minus": 5.13e-7, "k_p": 2.33e3},
    "R1p875e-4m": {"sieve": "0.25-0.5 mm", "radius": 1.875e-4, "dispersion": 0.9e-4, "half_time_h": 12.188,
                   "maximum_loading": 1.725, "k_plus": 0.0166, "k_minus": 2.45e-7, "k_p": 0.0294},
    "R5p5e-4m": {"sieve": "0.5-1.7 mm", "radius": 5.5e-4, "dispersion": 2.6e-4, "half_time_h": 12.080,
                 "maximum_loading": 1.751, "k_plus": 0.0145, "k_minus": 2.14e-7, "k_p": 0.0045},
    "R8p5e-4m": {"sieve": ">1.7 mm", "radius": 8.5e-4, "dispersion": 4.1e-4, "half_time_h": 14.028,
                 "maximum_loading": 1.915, "k_plus": 0.00440, "k_minus": 6.49e-8, "k_p": 0.0237},
}
# beta = 5.9e5 of the smallest size (k_p = 2.33e3 m/s) makes film and pore exchange 1e8 1/s, which CVODE does not get
# through; Auton et al. read it as "no intra-particle limit", which is direct uptake from the bulk in RUPTURA.
DIRECT_SIZE = "R1p25e-4m"

# Share of the resistance 1/k_p that is the film: k_f = k_p / share, D_p = R k_p / (5 (1 - share)).
SPLITS = {"film-limited": 0.999, "equal-split": 0.5, "pore-limited": 0.001}
SITE_AFFINITY = 1.0e3  # 1/Pa, as in elfving/; the site isotherm only supplies the capacity q_max of the General law

SIMULATION_TIME = 180_000.0  # s, 50 h
WRITE_INTERVAL = 900.0  # s
WAVE_REDUCED_LENGTH = 70.0  # wave column L = 70 L_ref, run for 70 T (the time its front needs to cross it)
WAVE_GRID_POINTS = 200
WAVE_OUTPUTS = 14  # column profiles per wave run
RELATIVE_TOLERANCE = 1.0e-7
ABSOLUTE_TOLERANCE_CONCENTRATION = 1.0e-9  # mol/m3
ABSOLUTE_TOLERANCE_LOADING = 1.0e-8  # mol/kg

SHOWN_GRID_POINTS = 100  # RMS 4e-4 or less from 800 cells, t_1/2 within 13 s (grid study)
GRID_POINTS = (50, 100, 200, 400, 800)
DENSE_GRID_POINTS = 200


def derived(size: dict) -> dict:
    """The dimensionless groups of Auton et al. (Eqs. 19, 21, 24) and the column length that closes their mass balance."""
    k_plus, k_minus, k_p = size["k_plus"], size["k_minus"], size["k_p"]
    mu = 1.0 + math.sqrt(k_minus / (k_plus * FEED_CONCENTRATION))
    equilibrium_loading = size["maximum_loading"] / mu
    time_scale = 1.0 / (equilibrium_loading * k_plus * FEED_CONCENTRATION)  # s
    length_scale = VELOCITY * time_scale * VOID_FRACTION * FEED_CONCENTRATION / (BULK_DENSITY * equilibrium_loading)  # m
    damkohler = length_scale / (VELOCITY * time_scale)
    beta = k_p * (3.0 / size["radius"]) * time_scale * damkohler * (1.0 - VOID_FRACTION) / VOID_FRACTION

    # Eq. (62) as coded by the authors: t - t_1/2 = (T/beta) (G(C) - G(1/2)), G(C) = A ln C - B ln(1 - C) + E ln(mu^2 - C).
    # The mean of t over C from 0 to 1 is the stoichiometric time, delta after the half time.
    mu2 = mu * mu
    a = 1.0 + beta / mu2
    b = (beta + (mu - 1.0) ** 2) / (mu2 - 1.0)
    e = (beta + mu2 * (mu - 1.0) ** 2) / (mu2 * (mu2 - 1.0))
    mean_g = -a + b + e * (mu2 * math.log(mu2) - (mu2 - 1.0) * math.log(mu2 - 1.0) - 1.0)
    g_half = (a - b) * math.log(0.5) + e * math.log(mu2 - 0.5)
    delta = time_scale / beta * (mean_g - g_half)  # s

    half_time = size["half_time_h"] * 3600.0  # s
    return {
        "mu": mu, "equilibrium_loading_mol_per_kg": equilibrium_loading, "time_scale_s": time_scale,
        "length_scale_m": length_scale, "beta": beta, "half_time_s": half_time, "delta_s": delta,
        "fit_length_m": (half_time + delta) * length_scale / time_scale,
        "published_stoichiometric_time_s": PUBLISHED_LENGTH * time_scale / length_scale,
        "peclet_inverse": size["dispersion"] / (VELOCITY * length_scale),
    }


def simulation(size: dict, length: float, dispersion: float, grid_points: int, transport: str, site_affinity: float,
               simulation_time: float, write_interval: float) -> dict:
    """One run. transport is a key of SPLITS, "direct" (uptake from the bulk, Auton's beta -> infinity), or
    "beta<value>" (half-and-half split with k_p scaled to that beta)."""
    radius = size["radius"]
    if transport == "direct":
        film_coefficient, pore_diffusivity = 0.0, 0.0
    else:
        share = SPLITS.get(transport, 0.5)
        k_p = size["k_p"] if transport in SPLITS else size["k_p"] * float(transport[len("beta"):]) / size["beta"]
        film_coefficient = k_p / share
        pore_diffusivity = radius * k_p / (5.0 * (1.0 - share))
    # Up to DENSE_GRID_POINTS cells CVODE uses the dense solver, whose difference-quotient Jacobian holds the film and
    # pore reservoirs. On finer grids it uses SPGMR with the banded sorption Jacobian as preconditioner, which leaves
    # the reservoirs out; that converges for the runs without axial dispersion, but stalls with the Table 5 dispersion.
    linear_solver = ({"CVODELinearSolver": "Dense"} if grid_points <= DENSE_GRID_POINTS
                     else {"CVODELinearSolver": "SPGMR", "CVODEPreconditioner": "SorptionJacobian"})
    return {
        "SimulationType": "Breakthrough",
        "DisplayName": f"CO2 on PEI/fumed silica, R = {radius} m, L = {length:.5f} m, D_ax = {dispersion} m2/s, "
                       f"{transport}, {grid_points} cells",
        "BreakthroughIntegrator": "CVODE",
        # The automatic CVODE cap, four cells at the gas velocity (0.007 s at 400 cells), is meant for fronts that move
        # with the gas; this front moves 1e5 times slower. The step is left to the error control.
        "CVODEMaximumTimeStep": 0.0,
        **linear_solver,
        "CVODERelativeTolerance": RELATIVE_TOLERANCE,
        "CVODEAbsoluteToleranceConcentration": ABSOLUTE_TOLERANCE_CONCENTRATION,
        "CVODEAbsoluteToleranceLoading": ABSOLUTE_TOLERANCE_LOADING,
        "ConvectionScheme": "WENO3",
        "InletCondition": "Danckwerts",
        "BoundaryCondition": "FixedVelocity",
        "MixturePredictionMethod": "SPI",
        "Temperature": TEMPERATURE,
        "InletPressure": PRESSURE,
        "ColumnEntranceVelocity": VELOCITY,
        "ColumnLength": length,
        "ColumnVoidFraction": VOID_FRACTION,
        "ParticleDensity": PARTICLE_DENSITY,
        "ParticleDiameter": 2.0 * radius,
        "DynamicViscosity": VISCOSITY,
        "NumberOfGridPoints": grid_points,
        "SimulationTime": simulation_time,
        "WriteInterval": write_interval,
        "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": VOID_FRACTION, "ParticleDiameter": 2.0 * radius},
        "Components": [
            {"Name": "Air", "GasPhaseMolFraction": 1.0 - FEED_MOLE_FRACTION, "MolecularWeight": 0.02897,
             "CarrierGas": True},
            {"Name": "CO2", "GasPhaseMolFraction": FEED_MOLE_FRACTION, "MolecularWeight": 0.04401,
             "MassTransferCoefficient": 1.0, "AxialDispersionCoefficient": dispersion,
             "PhysisorptionSites": [{"Type": "Henry", "Parameters": [0.0]}],
             "ChemisorptionSites": [{
                 "Type": "Langmuir",
                 "Parameters": [size["maximum_loading"], site_affinity],
                 "RateEquation": "General",
                 "RateEquationParameters": {
                     "AdsorptionRateCoefficient": size["k_plus"],
                     "DesorptionRateCoefficient": size["k_minus"],
                     "PoreConcentrationOrder": 1,
                     "CapacityOrder": 2,
                     "DesorptionOrder": 2,
                     "FilmMassTransferCoefficient": film_coefficient,
                     "PoreDiffusivity": pore_diffusivity,
                     "UsePoreSurfaceTransport": transport != "direct",
                 },
                 "HeatOfAdsorption": 0.0,
                 "MaximumLoading": size["maximum_loading"],
             }]},
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true",
                        help="delete simulations/ first, discarding any Ruptura output already there")
    arguments = parser.parse_args()

    if arguments.clean and SIMULATION_DIR.exists():
        shutil.rmtree(SIMULATION_DIR)
    SIMULATION_DIR.mkdir(parents=True, exist_ok=True)

    sizes = {name: {**size, **derived(size)} for name, size in SIZES.items()}

    # Each run is (size, column, dispersion, grid points, transport, site affinity); a run can serve several studies.
    studies: dict[tuple, list[str]] = {}

    def add(study: str, *run) -> None:
        studies.setdefault(run, []).append(study)

    for name in SIZES:
        default = "direct" if name == DIRECT_SIZE else "equal-split"
        for grid_points in GRID_POINTS:
            add("grid", name, "fit", "none", grid_points, default, SITE_AFFINITY)
        for column in ("fit", "published"):
            add("experiment", name, column, "table5", SHOWN_GRID_POINTS, default, SITE_AFFINITY)
        add("wave", name, "wave", "none", WAVE_GRID_POINTS, default, SITE_AFFINITY)
        transports = ("direct", "beta1000", "beta100") if name == DIRECT_SIZE else (*SPLITS, "direct")
        for transport in transports:
            add("beta" if name == DIRECT_SIZE else "split", name, "fit", "none", SHOWN_GRID_POINTS, transport,
                SITE_AFFINITY)
    add("affinity", "R5p5e-4m", "fit", "none", SHOWN_GRID_POINTS, "equal-split", SITE_AFFINITY)
    add("affinity", "R5p5e-4m", "fit", "none", SHOWN_GRID_POINTS, "equal-split", 1.0e6)

    index = []
    for (name, column, dispersion_case, grid_points, transport, site_affinity), run_studies in studies.items():
        size = sizes[name]
        directory = SIMULATION_DIR / (f"{name}_L-{column}_Dax-{dispersion_case}_{transport}_N{grid_points}"
                                      + ("" if site_affinity == SITE_AFFINITY else f"_affinity{site_affinity:g}"))
        directory.mkdir(parents=True, exist_ok=True)
        wave_length = WAVE_REDUCED_LENGTH * size["length_scale_m"]
        length = {"fit": size["fit_length_m"], "published": PUBLISHED_LENGTH, "wave": wave_length}[column]
        dispersion = size["dispersion"] if dispersion_case == "table5" else 0.0
        if column == "wave":
            write_interval = round(WAVE_REDUCED_LENGTH * size["time_scale_s"] / WAVE_OUTPUTS, -2)
            data = simulation(size, length, dispersion, grid_points, transport, site_affinity,
                              WAVE_OUTPUTS * write_interval, write_interval)
        else:
            data = simulation(size, length, dispersion, grid_points, transport, site_affinity, SIMULATION_TIME,
                              WRITE_INTERVAL)
        (directory / "simulation.json").write_text(json.dumps(data, indent=4) + "\n")
        index.append({"directory": str(directory.relative_to(ROOT)), "studies": run_studies, "size": name,
                      "column": column, "dispersion": dispersion_case, "grid_points": grid_points,
                      "transport": transport, "site_affinity_per_Pa": site_affinity, "length_m": length})

    parameters = {
        "temperature_K": TEMPERATURE, "pressure_Pa": PRESSURE, "feed_mole_fraction": FEED_MOLE_FRACTION,
        "feed_concentration_mol_per_m3": FEED_CONCENTRATION, "velocity_m_per_s": VELOCITY,
        "void_fraction": VOID_FRACTION, "bulk_density_kg_per_m3": BULK_DENSITY,
        "particle_density_kg_per_m3": PARTICLE_DENSITY, "published_length_m": PUBLISHED_LENGTH,
        "wave_reduced_length": WAVE_REDUCED_LENGTH, "shown_grid_points": SHOWN_GRID_POINTS, "direct_size": DIRECT_SIZE,
        "default_site_affinity_per_Pa": SITE_AFFINITY,
    }
    (SIMULATION_DIR / "index.json").write_text(
        json.dumps({"parameters": parameters, "sizes": sizes, "runs": index}, indent=2) + "\n")
    print(f"wrote {len(index)} simulation directories to {SIMULATION_DIR}")
    for name, size in sizes.items():
        print(f"{name}: mu {size['mu']:.5f}, T {size['time_scale_s']:.1f} s, L_ref {size['length_scale_m']:.3e} m, "
              f"beta {size['beta']:.3g}, Pe^-1 {size['peclet_inverse']:.3f}, delta {size['delta_s'] / 3600.0:.3f} h, "
              f"t_st(0.1 m) {size['published_stoichiometric_time_s'] / 3600.0:.2f} h, L_fit {size['fit_length_m']:.4f} m")


if __name__ == "__main__":
    main()
