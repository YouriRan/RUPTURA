"""Stitch the pure-component windows of one temperature into one-component distributions read by RUPTURA.

Usage: python build_pure_distribution.py [300]  (temperature in K, default 300)

CO2 comes from the slice without propane, propane from the edge. These states are part of the binary distribution of
build_distribution.py but do not need the slices with propane, so the pure isotherms can be checked on their own.
Writes particle_distribution_T<T>K_<component>.data and .json in the format of build_distribution.py.
"""

import json
import sys
from pathlib import Path

import numpy as np

from build_distribution import BOLTZMANN, FRAMEWORK_MASS, HERE, REFERENCE_PRESSURE, TMMC, stitch

WINDOWS = {"CO2": "slice_C3H8_000_CO2_N*", "C3H8": "edge_C3H8_N*"}

if __name__ == "__main__":
    temperature = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    label = f"T{temperature:.0f}K"
    for name, pattern in WINDOWS.items():
        lnpi, _, energy, _, intramolecular = stitch(list((TMMC / label).glob(pattern)))
        molecules = np.arange(len(lnpi))
        fugacities = np.logspace(2, 9, 701)
        ranges = [np.ptp(lnpi + molecules * np.log(f / REFERENCE_PRESSURE)) for f in fugacities]
        reference_fugacity = float(fugacities[np.argmin(ranges)])
        if min(ranges) > 700.0:
            raise RuntimeError(f"{name}: ln Pi spans {min(ranges):.0f} at every reference fugacity; doubles underflow")
        shifted = lnpi + molecules * np.log(reference_fugacity / REFERENCE_PRESSURE)
        shifted -= shifted.max()
        shifted -= np.log(np.exp(shifted).sum())
        # <H> only enters away from the reference temperature, but there a state with <H> = 0 gains exp(-<H>/k_B
        # (1/T - 1/T_ref)), up to e^200, over its neighbours. States without production samples therefore get <H> by
        # linear interpolation in N between sampled states and, beyond them, by a linear fit to the nearest five.
        sampled = np.isfinite(energy)
        without_energy = ~sampled
        interaction = (energy - intramolecular)[sampled]
        filled = np.interp(molecules, molecules[sampled], interaction)
        below, above = molecules < molecules[sampled][0], molecules > molecules[sampled][-1]
        filled[below] = np.polyval(np.polyfit(molecules[sampled][:5], interaction[:5], 1), molecules[below])
        filled[above] = np.polyval(np.polyfit(molecules[sampled][-5:], interaction[-5:], 1), molecules[above])
        hamiltonian = BOLTZMANN * filled

        file_name = f"particle_distribution_{label}_{name}.data"
        with open(HERE / file_name, "w") as data:
            data.write(f"# pure {name} in MOR (2x2x4 unit cells), TMMC windows {pattern} (../tmmc), T = {temperature} K\n")
            data.write(f"# reference fugacity {reference_fugacity:.10g} Pa; framework mass {FRAMEWORK_MASS:.10g} kg\n")
            data.write(f"# rows: n_{name} = 0..{molecules[-1]}; Pi_ref, <H> = <U - U_intra> [J]\n")
            for p, h in zip(np.exp(shifted), hamiltonian):
                data.write(f"{p:.17g} {h:.17g}\n")
        settings = {
            "FileName": file_name,
            "ReferenceTemperature": temperature,
            "ReferenceFugacity": reference_fugacity,
            "ReferenceFrameworkMass": FRAMEWORK_MASS,
            "ComponentBounds": [{"Component": name, "NMin": 0, "NMax": int(molecules[-1]), "DeltaN": 1}],
        }
        (HERE / f"particle_distribution_{label}_{name}.json").write_text(json.dumps(settings, indent=2))
        print(f"{label} {name}: {len(lnpi)} states, reference fugacity {reference_fugacity:.4g} Pa, "
              f"ln Pi range {min(ranges):.1f}, <H> filled on {without_energy.sum()} states without production energies")
