#!/usr/bin/env python3
"""Regenerate the two-component reference MPD used by this example."""

from math import exp
from pathlib import Path
from sys import float_info


N_MAX = 120
N_0 = 35.0
SIGMA_N = 3.0
SIGMA_D = 12.0


def unnormalized_probability(n_co2: int, n_n2: int) -> float:
    total = n_co2 + n_n2
    difference = n_co2 - n_n2
    preferred_co2_fraction = 0.70 - 0.60 / (1.0 + exp(-(total - 40.0) / 4.0))
    preferred_difference = (2.0 * preferred_co2_fraction - 1.0) * total
    exponent = -0.5 * ((total - N_0) / SIGMA_N) ** 2
    exponent -= 0.5 * ((difference - preferred_difference) / SIGMA_D) ** 2
    return exp(exponent)


def main() -> None:
    # C ordering: the last component (N2) varies fastest.
    probabilities = [
        unnormalized_probability(n_co2, n_n2)
        for n_co2 in range(N_MAX + 1)
        for n_n2 in range(N_MAX + 1)
    ]
    normalization = sum(probabilities)
    output = Path(__file__).with_name("particle_distribution.data")
    with output.open("w", encoding="utf-8") as stream:
        stream.write("# Pi_ref(n_CO2, n_N2); C order with n_N2 varying fastest\n")
        for probability in probabilities:
            normalized = probability / normalization
            # Some C++ stream implementations report an input error for
            # subnormal decimal values. These weights are below numerical
            # significance for the normalized distribution, so store zero.
            stored = normalized if normalized >= float_info.min else 0.0
            stream.write(f"{stored:.17g}\n")


if __name__ == "__main__":
    main()
