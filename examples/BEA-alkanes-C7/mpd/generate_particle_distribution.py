#!/usr/bin/env python3
"""Build the binary nC7/C6m2 macrostate particle distribution for BEA at 552 K.

The dual-site Langmuir isotherm of a component is a statement about independent
adsorption sites.  If a site family is shared by every component and has a
common saturation capacity, the occupancy of that family is multinomial.  This
script turns the fit produced by ``fit_isotherms.py`` into exactly that object,
convolves the two site families, and writes the dense C-ordered file consumed by
Ruptura's ``MPDSettings``.

Site family ``s`` with ``N_s`` sites and reference activities
``a_i = b_i,s f_ref`` has

    Pi_s(n) = N_s! / [(N_s - sum_i n_i)! prod_i n_i!] * prod_i a_i^n_i
              / (1 + sum_i a_i)^N_s ,

whose mean occupancy is ``<n_i> = N_s a_i / (1 + sum_j a_j)``.  Dividing by
``N_A m_fw`` reproduces the competitive Langmuir expression
``q_i = q_sat,s b_i,s f_i / (1 + sum_j b_j,s f_j)`` term by term.
"""

from __future__ import annotations

import json
from pathlib import Path
from sys import float_info

import numpy as np

AVOGADRO = 6.02214076e23  # mol^-1, exact SI definition

REFERENCE_FUGACITY = 1.0e5  # Pa, common component fugacity
N_MIN = 0
DELTA_N = 2

_PARAMETERS = json.loads((Path(__file__).with_name("fitted_parameters.json")).read_text(encoding="utf-8"))

COMPONENTS = tuple(_PARAMETERS["components"])
REFERENCE_TEMPERATURE = _PARAMETERS["referenceTemperature"]
REFERENCE_FRAMEWORK_MASS = _PARAMETERS["frameworkMass"]
HEAT_OF_ADSORPTION = _PARAMETERS["heatOfAdsorption"]
N_MAX = _PARAMETERS["totalSiteCount"]

#: ``(number of sites, (affinity per component, in COMPONENTS order))`` per family.
SITE_FAMILIES = tuple(
    (family["sites"], tuple(family["affinities"][component] for component in COMPONENTS))
    for family in _PARAMETERS["siteFamilies"]
)
SATURATION_CAPACITIES = tuple(family["saturationCapacity"] for family in _PARAMETERS["siteFamilies"])


def log_factorial(n: int) -> np.ndarray:
    """Return ``log(k!)`` for ``k = 0 .. n``."""
    return np.concatenate(([0.0], np.cumsum(np.log(np.arange(1, n + 1)))))


def site_family_distribution(number_of_sites: int, affinities: tuple[float, ...]) -> np.ndarray:
    """Normalized multinomial occupancy of one site family at the reference state."""
    activities = np.array([affinity * REFERENCE_FUGACITY for affinity in affinities])
    log_gamma = log_factorial(number_of_sites)

    counts = np.arange(number_of_sites + 1)
    first = counts[:, None]
    second = counts[None, :]
    occupied = first + second
    allowed = occupied <= number_of_sites

    log_weight = np.where(
        allowed,
        log_gamma[number_of_sites]
        - log_gamma[np.where(allowed, number_of_sites - occupied, 0)]
        - log_gamma[first]
        - log_gamma[second],
        -np.inf,
    )
    for axis, activity in enumerate(activities):
        index = first if axis == 0 else second
        if activity > 0.0:
            log_weight = log_weight + index * np.log(activity)
        else:
            log_weight = np.where(index > 0, -np.inf, log_weight)

    weight = np.exp(log_weight - log_weight.max())
    return weight / weight.sum()


def convolve_families() -> np.ndarray:
    """Convolve the site-resolved multinomials into total component counts."""
    distribution: np.ndarray | None = None
    for number_of_sites, affinities in SITE_FAMILIES:
        family = site_family_distribution(number_of_sites, affinities)
        if distribution is None:
            distribution = family
            continue
        width = distribution.shape[0] - 1
        extent = width + family.shape[0]
        convolved = np.zeros((extent, extent))
        for first, second in np.argwhere(family > 0.0):
            convolved[first : first + width + 1, second : second + width + 1] += (
                family[first, second] * distribution
            )
        distribution = convolved
    assert distribution is not None
    return distribution / distribution.sum()


def reference_distribution() -> tuple[np.ndarray, np.ndarray]:
    """Return the mesh Ruptura addresses and the distribution restricted to it."""
    mesh = np.arange(N_MIN, N_MAX + 1, DELTA_N)
    sampled = convolve_families()[np.ix_(mesh, mesh)]
    return mesh, sampled / sampled.sum()


def main() -> None:
    total_sites = sum(number_of_sites for number_of_sites, _ in SITE_FAMILIES)
    if total_sites != N_MAX:
        raise RuntimeError(f"site families hold {total_sites} sites but the window is [0, {N_MAX}]")

    # Ruptura addresses macrostates on the [N_MIN, N_MAX] mesh with spacing
    # DELTA_N.  Restricting the multinomial to that mesh and renormalizing keeps
    # the reweighted mean particle numbers exact wherever the distribution is
    # wide compared with DELTA_N; see the notebook for the residual error.
    mesh, sampled = reference_distribution()

    output = Path(__file__).with_name("particle_distribution.data")
    with output.open("w", encoding="utf-8") as stream:
        stream.write(f"# Pi_ref(n_{COMPONENTS[0]}, n_{COMPONENTS[1]}) for BEA at {REFERENCE_TEMPERATURE:g} K\n")
        stream.write(
            f"# Product of two multinomial site families "
            f"({', '.join(str(sites) for sites, _ in SITE_FAMILIES)} sites), "
            f"collapsed to total counts\n"
        )
        stream.write(
            f"# C order with n_{COMPONENTS[1]} varying fastest; "
            f"n = {N_MIN}, {N_MIN + DELTA_N}, ..., {N_MAX}\n"
        )
        for row in sampled:
            for probability in row:
                # Subnormal decimals upset some C++ stream implementations and
                # are far below numerical significance here.
                stream.write(f"{probability if probability >= float_info.min else 0.0:.17g}\n")

    counts = mesh.astype(float)
    mean_first = float(sampled.sum(axis=1) @ counts)
    mean_second = float(sampled.sum(axis=0) @ counts)
    scale = AVOGADRO * REFERENCE_FRAMEWORK_MASS

    print(f"Wrote {output} ({sampled.size} macrostates, {int((sampled > 0.0).sum())} nonzero)")
    print(f"Probability sum: {sampled.sum():.17g}")
    print(
        f"Reference mean particle numbers: "
        f"{COMPONENTS[0]} = {mean_first:.4f} ({mean_first / scale:.5f} mol/kg), "
        f"{COMPONENTS[1]} = {mean_second:.4f} ({mean_second / scale:.5f} mol/kg)"
    )


if __name__ == "__main__":
    main()
