#!/usr/bin/env python3
"""Build the BEA alkane MPD from independent multinomial site families."""

from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import factorial
from pathlib import Path


COMPONENTS = ("nC7", "C6m2", "C6m3", "C5m2m2", "C5m2m3", "C5m2m4", "C5m3m3")
REFERENCE_FUGACITY = 1.0e5  # Pa, common component fugacity

# The affinities are taken from examples/BEA-alkanes-C7/iast. A common
# saturation capacity is required within each competitive site family. Three
# primary sites and one secondary site represent q_sat = 1.2 and 0.4 mol/kg,
# respectively, for the framework mass used in simulation.json.
SITE_FAMILIES = (
    (
        3,
        (6.55857e-5, 3.90895e-5, 2.61740e-5, 6.99287e-6, 1.41992e-5, 2.11233e-5, 4.39800e-6),
    ),
    (1, (8.90731e-7, 9.64046e-8, 7.99897e-8, 0.0, 0.0, 0.0, 0.0)),
)


Macrostate = tuple[int, ...]


def multinomial_site_distribution(number_of_sites: int, affinities: tuple[float, ...]) -> dict[Macrostate, float]:
    """Return the normalized occupancy distribution for one site family."""
    activities = tuple(affinity * REFERENCE_FUGACITY for affinity in affinities)
    site_partition = 1.0 + sum(activities)
    distribution: dict[Macrostate, float] = {}

    for counts in product(range(number_of_sites + 1), repeat=len(COMPONENTS)):
        occupied = sum(counts)
        if occupied > number_of_sites:
            continue

        coefficient = factorial(number_of_sites) / factorial(number_of_sites - occupied)
        weight = coefficient / site_partition**number_of_sites
        for count, activity in zip(counts, activities, strict=True):
            weight *= activity**count / factorial(count)
        distribution[counts] = weight

    return distribution


def collapsed_distribution() -> dict[Macrostate, float]:
    """Convolve site-resolved multinomials into total component counts."""
    zero = (0,) * len(COMPONENTS)
    totals: dict[Macrostate, float] = {zero: 1.0}

    for number_of_sites, affinities in SITE_FAMILIES:
        family = multinomial_site_distribution(number_of_sites, affinities)
        convolved: defaultdict[Macrostate, float] = defaultdict(float)
        for previous, previous_probability in totals.items():
            for occupancy, occupancy_probability in family.items():
                macrostate = tuple(a + b for a, b in zip(previous, occupancy, strict=True))
                convolved[macrostate] += previous_probability * occupancy_probability
        totals = dict(convolved)

    normalization = sum(totals.values())
    return {macrostate: probability / normalization for macrostate, probability in totals.items()}


def main() -> None:
    distribution = collapsed_distribution()
    maximum_count = sum(number_of_sites for number_of_sites, _ in SITE_FAMILIES)
    output = Path(__file__).with_name("mpd") / "particle_distribution.data"

    with output.open("w", encoding="utf-8") as stream:
        stream.write("# Pi_ref(n_nC7,n_C6m2,n_C6m3,n_C5m2m2,n_C5m2m3,n_C5m2m4,n_C5m3m3)\n")
        stream.write("# Product of site-family multinomials, collapsed to total counts; last component varies fastest.\n")
        for macrostate in product(range(maximum_count + 1), repeat=len(COMPONENTS)):
            stream.write(f"{distribution.get(macrostate, 0.0):.17g}\n")

    print(f"Wrote {output} ({(maximum_count + 1) ** len(COMPONENTS)} dense macrostates)")
    print(f"Nonzero macrostates: {len(distribution)}; probability sum: {sum(distribution.values()):.17g}")


if __name__ == "__main__":
    main()
