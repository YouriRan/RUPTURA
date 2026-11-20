#!/usr/bin/env python3
"""Set the particle density of the ternary runs T1-T24 so that the adsorbent per unit fluid volume matches Marx et al.

RUPTURA writes the mass balance per unit fluid volume with the interstitial velocity v and the adsorbent term
rho_p (1 - eps_b)/eps_b dq/dt.  Marx et al. (2014), Table 3, write it with the total porosity eps_t, the superficial
velocity u = Q/A, and the bulk density rho_b = 502 kg/m3 (Table 4):

    eps_t dc/dt + d(u c)/dz = -rho_b dq/dt.

marx2 sets v = u/eps_t, so both balances agree when

    rho_p (1 - eps_b)/eps_b = rho_b/eps_t,    rho_p = rho_b eps_b / ((1 - eps_b) eps_t),

with eps_t = eps_b + (1 - eps_b) eps_p (see conversions.md).  The inputs used rho_b itself (502 kg/m3), i.e.
10% more adsorbent per unit fluid volume than the reference model.

    python particle_density.py            # print the table and write ParticleDensity into reference_settings_hassan_chung/T*/simulation.json
    python particle_density.py --dry-run  # print the table only

Only the key ParticleDensity is changed.  All quantities are SI.

These are the frozen reference settings of Hassan and Chung; the Marx validation in ternary/ uses the total-porosity
mapping of generate.py (conversions.md, section 6) instead.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

BULK_DENSITY = 502.0  # kg/m3, Marx et al. (2014), Table 4
# eps_p from eps_t = 0.742 and eps_b = 0.403 (Table 4), as in the conversion of the flow rate to ColumnEntranceVelocity
PARTICLE_POROSITY = 0.568


def particle_density(simulation: dict) -> float:
    """Particle density [kg/m3] that reproduces rho_b/eps_t per unit fluid volume."""
    bed_porosity = simulation["ColumnVoidFraction"]
    total_porosity = bed_porosity + (1.0 - bed_porosity) * PARTICLE_POROSITY
    return BULK_DENSITY * bed_porosity / ((1.0 - bed_porosity) * total_porosity)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="print the table without changing the inputs")
    arguments = parser.parse_args()

    print(f"{'run':4s} {'eps_b':>6s} {'rho_p/(kg/m3)':>14s} {'previous':>9s}")
    for index in range(1, 25):
        path = ROOT / "reference_settings_hassan_chung" / f"T{index}" / "simulation.json"
        simulation = json.loads(path.read_text())
        density = particle_density(simulation)
        previous = simulation["ParticleDensity"]
        print(f"T{index:<3d} {simulation['ColumnVoidFraction']:6.3f} {density:14.1f} {previous:9.1f}")
        if not arguments.dry_run:
            simulation["ParticleDensity"] = round(density, 1)
            path.write_text(json.dumps(simulation, indent=4))


if __name__ == "__main__":
    main()
