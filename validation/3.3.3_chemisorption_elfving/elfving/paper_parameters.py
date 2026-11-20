"""Parameters and equations from Elfving and Sainio (2021).

All pressures passed to the co-adsorption and GAB helpers are in bar.  The
resulting loadings are in mol/kg, numerically equal to mmol/g.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


GAS_CONSTANT = 8.31446261815324
REFERENCE_TEMPERATURE = 298.15
GAB_PARAMETERS = (2.58, 0.155, 0.871, 6_600.0, 0.0)


@dataclass(frozen=True)
class CoAdsorptionModel:
    name: str
    qm: float
    b01: float
    b02: float
    t1: float
    t2: float
    heat1: float
    heat2: float


MODELS = {
    "five": CoAdsorptionModel("5-parameter", 1.91, 3.199e3, 7.461e4, 2.0, 1.0, 98_000.0, 111_500.0),
    "seven_bounded": CoAdsorptionModel(
        "7-parameter (t1=t2=3)", 2.61, 4.200e2, 2.540e4, 3.0, 3.0, 85_000.0, 123_400.0
    ),
    "seven_relaxed": CoAdsorptionModel(
        "7-parameter (relaxed)", 6.98, 9.224e-6, 7.372e-7, 10.56, 14.42, 81_700.0, 126_700.0
    ),
}


def temperature_affinity(reference_affinity: float, heat: float, temperature: float) -> float:
    """Paper Eqs. 12-13, with ``heat`` storing the positive value of -Delta H."""
    exponent = heat / (GAS_CONSTANT * REFERENCE_TEMPERATURE) * (
        REFERENCE_TEMPERATURE / temperature - 1.0
    )
    return reference_affinity * math.exp(exponent)


def coadsorption_loading(
    model: CoAdsorptionModel, temperature: float, p_co2_bar: float, p_h2o_bar: float
) -> tuple[float, float]:
    """Return the dry- and humid-pathway equilibrium CO2 loadings.

    Setting Eqs. 14-15 to zero leaves one monotone scalar balance for the
    remaining free-amine concentration.  Bisection avoids integrating an ODE
    merely to reach the same equilibrium.
    """
    if p_co2_bar <= 0.0:
        return 0.0, 0.0
    b1 = temperature_affinity(model.b01, model.heat1, temperature)
    b2 = temperature_affinity(model.b02, model.heat2, temperature)

    def balance(free_amine: float) -> float:
        q1 = b1 * p_co2_bar * free_amine**model.t1
        q2 = b2 * p_co2_bar * p_h2o_bar * free_amine**model.t2
        return free_amine + 2.0 * q1 + q2 - model.qm

    low, high = 0.0, model.qm
    for _ in range(160):
        midpoint = 0.5 * (low + high)
        if balance(midpoint) > 0.0:
            high = midpoint
        else:
            low = midpoint
    free_amine = 0.5 * (low + high)
    return (
        b1 * p_co2_bar * free_amine**model.t1,
        b2 * p_co2_bar * p_h2o_bar * free_amine**model.t2,
    )


def water_saturation_pressure_bar(temperature: float) -> float:
    """Kelvin form of the Antoine correlation for liquid water."""
    pressure_mmhg = 10.0 ** (8.07131 - 1730.63 / (temperature - 39.724))
    return pressure_mmhg * 133.32236842105263 / 100_000.0


def gab_loading(temperature: float, p_h2o_bar: float) -> float:
    """Paper Eqs. 5-7 and Table 2."""
    if p_h2o_bar <= 0.0:
        return 0.0
    qm_mono, c0, k0, delta_h_c, delta_h_k = GAB_PARAMETERS
    c_value = c0 * math.exp(delta_h_c / (GAS_CONSTANT * temperature))
    k_value = k0 * math.exp(delta_h_k / (GAS_CONSTANT * temperature))
    relative_pressure = p_h2o_bar / water_saturation_pressure_bar(temperature)
    numerator = qm_mono * c_value * k_value * relative_pressure
    denominator = (1.0 - k_value * relative_pressure) * (
        1.0 + k_value * relative_pressure * (c_value - 1.0)
    )
    return numerator / denominator


def empty_column_water_loading(p_h2o_bar: float) -> float:
    """Supplementary Fig. S2 Langmuir correction, normalized to 0.5 g resin."""
    bp = 67.1217 * max(0.0, p_h2o_bar)
    return 1.1216 * bp / (1.0 + bp)


def toth_supplement_parameters(temperature: float) -> tuple[float, float, float]:
    """Supplementary Eqs. S17-S19 and Table S1 (qmax, b in 1/bar, t)."""
    qmax = 1.84 * math.exp(3.52 * (1.0 - temperature / REFERENCE_TEMPERATURE))
    affinity = temperature_affinity(9.459e4, 112_800.0, temperature)
    exponent = 0.340 + 0.395 * (1.0 - REFERENCE_TEMPERATURE / temperature)
    return qmax, affinity, exponent


def toth_loading(qmax: float, affinity: float, exponent: float, pressure: float) -> float:
    bp = max(0.0, affinity * pressure)
    return qmax * bp / (1.0 + bp**exponent) ** (1.0 / exponent)


def fit_toth(pressures_bar: list[float], loadings: list[float]) -> tuple[float, float, float, float]:
    """Fit a three-parameter Toth projection without a SciPy dependency.

    A deterministic multi-start pattern search is sufficient for these smooth
    paper curves.  The fourth return value is the relative RMS error.
    """
    positive = [(p, q) for p, q in zip(pressures_bar, loadings) if p > 0.0 and q > 0.0]
    if not positive:
        return 1.0, 1.0, 1.0, 0.0

    scale = max(q for _, q in positive)

    def objective(log_parameters: tuple[float, float, float]) -> float:
        qmax, affinity, exponent = (math.exp(value) for value in log_parameters)
        # Exponents below 0.5 make the native Toth derivative singular enough
        # near a clean-bed pressure of zero to overwhelm CVODE in long columns.
        # A 0.5 floor regularizes the dynamic projection while retaining the
        # paper curves to within the declared fitting tolerance.
        if exponent < 0.5 or exponent > 8.0 or qmax < scale:
            return 1.0e30
        return sum(
            ((toth_loading(qmax, affinity, exponent, pressure) - target) / scale) ** 2
            for pressure, target in positive
        ) / len(positive)

    best = (1.0e30, (math.log(scale * 1.2), math.log(1.0e4), math.log(0.7)))
    for q_factor in (1.01, 1.2, 1.8, 3.0, 6.0):
        for affinity in (1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6):
            for exponent in (0.5, 0.7, 1.0, 1.7):
                point = (math.log(scale * q_factor), math.log(affinity), math.log(exponent))
                value = objective(point)
                if value < best[0]:
                    best = value, point

    value, point = best
    steps = [0.8, 1.2, 0.6]
    for _ in range(100):
        improved = False
        for dimension in range(3):
            for direction in (-1.0, 1.0):
                candidate = list(point)
                candidate[dimension] += direction * steps[dimension]
                candidate_tuple = tuple(candidate)
                candidate_value = objective(candidate_tuple)
                if candidate_value < value:
                    value, point = candidate_value, candidate_tuple
                    improved = True
        if not improved:
            steps = [step * 0.62 for step in steps]
            if max(steps) < 1.0e-7:
                break

    qmax, affinity, exponent = (math.exp(item) for item in point)
    return qmax, affinity, exponent, math.sqrt(value)


def effective_co2_heat(model: CoAdsorptionModel, q1: float, q2: float) -> float:
    total = q1 + q2
    return model.heat1 if total <= 0.0 else (model.heat1 * q1 + model.heat2 * q2) / total
