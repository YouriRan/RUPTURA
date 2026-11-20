"""Minimal plug-flow/LDF model used by Liu et al. (2024).

The paper solves the gas-phase balance with the method of lines.  Gas residence
times in all cases here are seconds while adsorption takes minutes to hours, so
this implementation uses the corresponding quasi-steady plug-flow limit.  Each
time step is conservative: a cell cannot adsorb more material than enters it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1
STANDARD_TEMPERATURE = 273.15  # K
STANDARD_PRESSURE = 100_000.0  # Pa


def toth_loading(pressure: float | np.ndarray, n_s: float, affinity: float, exponent: float):
    """Toth loading in mol/kg."""

    bp = np.maximum(affinity * np.asarray(pressure), 0.0)
    return n_s * bp / np.power(1.0 + np.power(bp, exponent), 1.0 / exponent)


def gab_loading(relative_humidity: float | np.ndarray, c_m: float, c_g: float, k_ads: float):
    """GAB water loading in mol/kg."""

    phi = np.clip(np.asarray(relative_humidity), 0.0, 0.999999)
    numerator = c_m * c_g * k_ads * phi
    denominator = (1.0 - k_ads * phi) * (1.0 + (c_g - 1.0) * k_ads * phi)
    return numerator / denominator


def buck_saturation_pressure_pa(temperature_k: float | np.ndarray):
    """Buck equation as printed in equation 9, returned in pascals."""

    temperature_c = np.asarray(temperature_k) - 273.15
    return 611.21 * np.exp(
        (18.678 - temperature_c / 234.5) * temperature_c / (257.14 + temperature_c)
    )


def temperature_dependent_toth(
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
    """Equations 4-6: return ``n_s, b, t_T`` at temperature ``T``."""

    t0 = reference_temperature
    n_s = n_s0 * np.exp(chi * (1.0 - temperature / t0))
    affinity = affinity0 * np.exp(heat / (GAS_CONSTANT * t0) * (t0 / temperature - 1.0))
    exponent = exponent0 + alpha * (1.0 - t0 / temperature)
    return float(n_s), float(affinity), float(exponent)


def standard_flow_to_actual_velocity(
    standard_flow_m3_s: float,
    area_m2: float,
    temperature: float,
    pressure_pa: float = STANDARD_PRESSURE,
) -> float:
    """Convert a standard volumetric flow in m3/s to superficial m/s."""

    actual_flow_m3_s = standard_flow_m3_s * temperature / STANDARD_TEMPERATURE
    actual_flow_m3_s *= STANDARD_PRESSURE / pressure_pa
    return actual_flow_m3_s / area_m2


@dataclass(frozen=True)
class Bed:
    length_m: float
    superficial_velocity_m_s: float
    bulk_density_kg_m3: float
    void_fraction: float = 0.5
    cells: int = 80


def _advance_species(
    q: np.ndarray,
    c_feed: float,
    bed: Bed,
    rate_constant_s: float,
    dt: float,
    equilibrium: Callable[[float, int], float],
) -> tuple[np.ndarray, float, np.ndarray]:
    """Advance one adsorbing species with a conservative cell sweep."""

    dz = bed.length_m / bed.cells
    c_in = float(c_feed)
    updated = q.copy()
    local_concentration = np.empty_like(q)
    relaxation = 1.0 - np.exp(-rate_constant_s * dt)
    cell_supply_factor = bed.superficial_velocity_m_s * dt / (bed.bulk_density_kg_m3 * dz)

    for index in range(bed.cells):
        q_equilibrium = max(float(equilibrium(c_in, index)), 0.0)
        potential = max((q_equilibrium - q[index]) * relaxation, 0.0)
        supply_limited = max(c_in * cell_supply_factor, 0.0)
        uptake = min(potential, supply_limited)
        updated[index] += uptake
        c_out = c_in - uptake / cell_supply_factor
        c_in = max(c_out, 0.0)
        local_concentration[index] = c_in

    return updated, c_in, local_concentration


def simulate_dry(
    times_s: np.ndarray,
    *,
    bed: Bed,
    feed_concentration_mol_m3: float,
    rate_constant_s: float,
    equilibrium: Callable[[float], float],
    internal_step_s: float = 10.0,
    initial_loading: float = 0.0,
) -> dict[str, np.ndarray]:
    """Simulate one dry adsorption step and return outlet/loading histories."""

    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1 or len(times) == 0 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times_s must be a strictly increasing vector beginning at zero")

    q = np.full(bed.cells, initial_loading, dtype=float)
    outlet = np.zeros_like(times)
    average_loading = np.full_like(times, initial_loading)
    current_time = 0.0
    current_outlet = 0.0

    for output_index in range(1, len(times)):
        target = times[output_index]
        while current_time < target - 1.0e-12:
            dt = min(internal_step_s, target - current_time)
            q, current_outlet, _ = _advance_species(
                q,
                feed_concentration_mol_m3,
                bed,
                rate_constant_s,
                dt,
                lambda concentration, _index: equilibrium(concentration),
            )
            current_time += dt
        outlet[output_index] = current_outlet / feed_concentration_mol_m3
        average_loading[output_index] = np.mean(q)

    return {
        "time_s": times,
        "outlet_fraction": np.clip(outlet, 0.0, 1.0),
        "average_loading_mol_kg": average_loading,
        "final_loading_profile_mol_kg": q,
    }


def simulate_humid(
    times_s: np.ndarray,
    *,
    bed: Bed,
    co2_feed_concentration_mol_m3: float,
    water_feed_concentration_mol_m3: float,
    co2_rate_constant_s: float,
    water_rate_constant_s: float,
    water_equilibrium: Callable[[float], float],
    co2_equilibrium: Callable[[float, float, float], float],
    feed_relative_humidity: float,
    internal_step_s: float = 10.0,
    initial_co2_loading: float = 0.0,
    initial_water_loading: float = 0.0,
) -> dict[str, np.ndarray]:
    """Simulate coupled H2O/CO2 adsorption in the paper's minimal model."""

    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1 or len(times) == 0 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times_s must be a strictly increasing vector beginning at zero")

    q_co2 = np.full(bed.cells, initial_co2_loading, dtype=float)
    q_water = np.full(bed.cells, initial_water_loading, dtype=float)
    co2_outlet = np.zeros_like(times)
    water_outlet = np.zeros_like(times)
    co2_average = np.full_like(times, initial_co2_loading)
    water_average = np.full_like(times, initial_water_loading)
    current_time = 0.0
    current_co2_outlet = 0.0
    current_water_outlet = 0.0

    for output_index in range(1, len(times)):
        target = times[output_index]
        while current_time < target - 1.0e-12:
            dt = min(internal_step_s, target - current_time)
            q_water, current_water_outlet, local_water = _advance_species(
                q_water,
                water_feed_concentration_mol_m3,
                bed,
                water_rate_constant_s,
                dt,
                lambda concentration, _index: water_equilibrium(
                    feed_relative_humidity * concentration / water_feed_concentration_mol_m3
                ),
            )

            def local_co2_equilibrium(concentration: float, index: int) -> float:
                local_rh = feed_relative_humidity * local_water[index] / water_feed_concentration_mol_m3
                return co2_equilibrium(concentration, q_water[index], local_rh)

            q_co2, current_co2_outlet, _ = _advance_species(
                q_co2,
                co2_feed_concentration_mol_m3,
                bed,
                co2_rate_constant_s,
                dt,
                local_co2_equilibrium,
            )
            current_time += dt

        co2_outlet[output_index] = current_co2_outlet / co2_feed_concentration_mol_m3
        water_outlet[output_index] = current_water_outlet / water_feed_concentration_mol_m3
        co2_average[output_index] = np.mean(q_co2)
        water_average[output_index] = np.mean(q_water)

    return {
        "time_s": times,
        "co2_outlet_fraction": np.clip(co2_outlet, 0.0, 1.0),
        "water_outlet_fraction": np.clip(water_outlet, 0.0, 1.0),
        "co2_average_loading_mol_kg": co2_average,
        "water_average_loading_mol_kg": water_average,
        "final_co2_loading_profile_mol_kg": q_co2,
        "final_water_loading_profile_mol_kg": q_water,
    }


def crossing_time(time: np.ndarray, response: np.ndarray, level: float) -> float:
    """Linearly interpolated first upward crossing, or NaN if not crossed."""

    indices = np.flatnonzero(np.asarray(response) >= level)
    if len(indices) == 0:
        return float("nan")
    upper = int(indices[0])
    if upper == 0:
        return float(time[0])
    lower = upper - 1
    x0, x1 = response[lower], response[upper]
    if x1 == x0:
        return float(time[upper])
    fraction = (level - x0) / (x1 - x0)
    return float(time[lower] + fraction * (time[upper] - time[lower]))
