#!/usr/bin/env python3
"""Fit the binary nC7 / C6m2 adsorption model to the BEA GCMC data.

The GCMC results in ``examples/BEA-alkanes-C7/fitting`` are single-temperature
(552 K) absolute-loading isotherms plus, in column 20, the simulated heat of
desorption.  Two things are fitted here and one is read:

* ``q_sat`` -- one saturation capacity **per site family, shared by both
  components**.  A competitive multinomial site model cannot give two species
  different capacities on the same sites, so the constraint is imposed during
  the fit rather than patched in afterwards.
* ``b``     -- one affinity per component and site family.
* ``Delta H`` is *not* fitted.  A single-temperature isotherm carries no
  information about it; the value is taken from the Henry-regime average of the
  GCMC heat of desorption, ``Q_eq = -R <heat of desorption>``, which is the
  positive quantity Ruptura's ``HeatOfAdsorption`` expects.

The result is written to ``fitted_parameters.json``, which
``generate_particle_distribution.py`` and ``generate_simulations.py`` both read.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

GAS_CONSTANT = 8.31446261815324  # J mol^-1 K^-1
AVOGADRO = 6.02214076e23  # mol^-1

TEMPERATURE = 552.0  # K, the temperature of the GCMC data
NUMBER_OF_SITES = 2
TOTAL_SITE_COUNT = 200  # particle-number window shared by both components

DATA_FILES = {
    "nC7": "Results.dat-BEA-Repeat-552K-nC7",
    "C6m2": "Results.dat-BEA-Repeat-552K-2mC6",
}
FUGACITY_COLUMN = 2  # column 3: fugacity [Pa]
LOADING_COLUMN = 7  # column 8: absolute loading [mol/kg framework]
LOADING_ERROR_COLUMN = 8  # column 9
HEAT_COLUMN = 19  # column 20: heat of desorption [K]

# Henry regime for the heat average: loadings below this fraction of the
# largest simulated loading.
HENRY_LOADING_FRACTION = 0.05

# Starting point: the existing single-component BEA fit (fitting/BEA_fitted.json).
INITIAL_CAPACITIES = (1.2, 0.4)
INITIAL_AFFINITIES = {
    "nC7": (6.55857e-05, 8.90731e-07),
    "C6m2": (3.90895e-05, 9.64046e-08),
}


def load_gcmc(fitting_dir: Path) -> dict[str, dict[str, np.ndarray]]:
    data = {}
    for component, file_name in DATA_FILES.items():
        table = np.loadtxt(fitting_dir / file_name, comments="#")
        data[component] = {
            "fugacity": table[:, FUGACITY_COLUMN],
            "loading": table[:, LOADING_COLUMN],
            "loading_error": table[:, LOADING_ERROR_COLUMN],
            "heat_of_desorption": table[:, HEAT_COLUMN],
        }
    return data


def dual_site_langmuir(capacities, affinities, fugacity):
    return sum(
        capacity * affinity * fugacity / (1.0 + affinity * fugacity)
        for capacity, affinity in zip(capacities, affinities)
    )


def nelder_mead(function, start, step, iterations=40000, tolerance=1e-15):
    """Plain Nelder-Mead, so the example needs nothing beyond NumPy."""
    dimension = len(start)
    simplex = [np.array(start, dtype=float)]
    for axis in range(dimension):
        vertex = np.array(start, dtype=float)
        vertex[axis] += step
        simplex.append(vertex)
    simplex = np.array(simplex)
    values = np.array([function(vertex) for vertex in simplex])

    for _ in range(iterations):
        order = np.argsort(values)
        simplex, values = simplex[order], values[order]
        if abs(values[-1] - values[0]) <= tolerance * (abs(values[0]) + tolerance):
            break
        centroid = simplex[:-1].mean(axis=0)
        reflected = centroid + (centroid - simplex[-1])
        reflected_value = function(reflected)
        if reflected_value < values[0]:
            expanded = centroid + 2.0 * (centroid - simplex[-1])
            expanded_value = function(expanded)
            simplex[-1], values[-1] = (
                (expanded, expanded_value) if expanded_value < reflected_value else (reflected, reflected_value)
            )
        elif reflected_value < values[-2]:
            simplex[-1], values[-1] = reflected, reflected_value
        else:
            contracted = centroid + 0.5 * (simplex[-1] - centroid)
            contracted_value = function(contracted)
            if contracted_value < values[-1]:
                simplex[-1], values[-1] = contracted, contracted_value
            else:
                simplex[1:] = simplex[0] + 0.5 * (simplex[1:] - simplex[0])
                values[1:] = [function(vertex) for vertex in simplex[1:]]

    order = np.argsort(values)
    return simplex[order][0], float(values[order][0])


def unpack(parameters: np.ndarray):
    """Parameters are fitted as logarithms, which keeps them positive."""
    values = np.exp(parameters)
    capacities = tuple(values[:NUMBER_OF_SITES])
    affinities = {
        component: tuple(values[NUMBER_OF_SITES + index * NUMBER_OF_SITES :
                                NUMBER_OF_SITES + (index + 1) * NUMBER_OF_SITES])
        for index, component in enumerate(DATA_FILES)
    }
    return capacities, affinities


def fit(data) -> tuple[tuple[float, ...], dict[str, tuple[float, ...]], float]:
    def residual_sum(parameters: np.ndarray) -> float:
        capacities, affinities = unpack(parameters)
        total = 0.0
        for component, measured in data.items():
            residual = dual_site_langmuir(capacities, affinities[component], measured["fugacity"])
            residual = residual - measured["loading"]
            total += float(residual @ residual)
        return total

    start = np.log(np.array(list(INITIAL_CAPACITIES) + [b for component in DATA_FILES for b in INITIAL_AFFINITIES[component]]))
    best, value = nelder_mead(residual_sum, start, 0.3)
    capacities, affinities = unpack(best)

    # Order the site families from strong to weak affinity for the first component.
    order = sorted(range(NUMBER_OF_SITES), key=lambda site: -affinities[next(iter(DATA_FILES))][site])
    capacities = tuple(capacities[site] for site in order)
    affinities = {component: tuple(values[site] for site in order) for component, values in affinities.items()}
    return capacities, affinities, value


def henry_regime_heat(measured) -> tuple[float, int]:
    """Q_eq [J/mol] from the low-loading GCMC heat of desorption."""
    threshold = HENRY_LOADING_FRACTION * measured["loading"].max()
    selection = measured["loading"] <= threshold
    if selection.sum() < 2:
        selection = measured["loading"] <= np.sort(measured["loading"])[2]
    return float(-GAS_CONSTANT * measured["heat_of_desorption"][selection].mean()), int(selection.sum())


def main() -> None:
    here = Path(__file__).resolve().parent
    data = load_gcmc(here.parent / "fitting")
    capacities, affinities, residual = fit(data)

    total_capacity = sum(capacities)
    # Site counts must be integers, and the two components share one particle
    # window.  Split TOTAL_SITE_COUNT in proportion to the fitted capacities and
    # let the framework mass follow, so the largest macrostate is exactly
    # TOTAL_SITE_COUNT particles.
    site_counts = [int(round(TOTAL_SITE_COUNT * capacity / total_capacity)) for capacity in capacities]
    site_counts[-1] = TOTAL_SITE_COUNT - sum(site_counts[:-1])
    framework_mass = TOTAL_SITE_COUNT / (AVOGADRO * total_capacity)
    represented_capacities = [count / (AVOGADRO * framework_mass) for count in site_counts]

    heats = {}
    print(f"GCMC fit of {' and '.join(DATA_FILES)} on BEA at {TEMPERATURE:g} K")
    print(f"  weighted residual sum of squares: {residual:.6e} (mol/kg)^2\n")
    print(f"{'site':>5} {'q_sat fitted':>13} {'q_sat used':>11} {'sites':>6} "
          + " ".join(f"{'b_' + component:>14}" for component in DATA_FILES))
    for site in range(NUMBER_OF_SITES):
        print(f"{site + 1:>5} {capacities[site]:>13.5f} {represented_capacities[site]:>11.5f} {site_counts[site]:>6} "
              + " ".join(f"{affinities[component][site]:>14.5e}" for component in DATA_FILES))
    print(f"\nFramework mass: {framework_mass:.10e} kg "
          f"({AVOGADRO * framework_mass:.3f} particles per mol/kg)")

    print("\nHeat of adsorption from the GCMC heat of desorption (Henry regime):")
    for component, measured in data.items():
        heat, count = henry_regime_heat(measured)
        heats[component] = heat
        print(f"  {component:>5}: Q_eq = {heat / 1000.0:7.3f} kJ/mol   (mean of {count} low-loading points)")

    print("\nFit quality:")
    for component, measured in data.items():
        predicted = dual_site_langmuir(capacities, affinities[component], measured["fugacity"])
        deviation = np.abs(predicted - measured["loading"])
        print(f"  {component:>5}: max |q_fit - q_GCMC| = {deviation.max():.4f} mol/kg "
              f"({deviation.max() / measured['loading'].max():.2%} of the largest loading)")

    payload = {
        "source": "examples/BEA-alkanes-C7/fitting (GCMC, BEA, 552 K)",
        "components": list(DATA_FILES),
        "referenceTemperature": TEMPERATURE,
        "totalSiteCount": TOTAL_SITE_COUNT,
        "frameworkMass": framework_mass,
        "siteFamilies": [
            {
                "sites": site_counts[site],
                "saturationCapacityFitted": capacities[site],
                "saturationCapacity": represented_capacities[site],
                "affinities": {component: affinities[component][site] for component in DATA_FILES},
            }
            for site in range(NUMBER_OF_SITES)
        ],
        "heatOfAdsorption": heats,
        "residualSumOfSquares": residual,
    }
    target = here / "fitted_parameters.json"
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote {target}")


if __name__ == "__main__":
    main()
