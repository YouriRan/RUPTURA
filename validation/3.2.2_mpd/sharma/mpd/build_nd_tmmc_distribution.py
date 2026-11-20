"""Binary 300 K distribution from the N-d TMMC ln Pi(n, m) computed on another machine.

Usage: python build_nd_tmmc_distribution.py

Reads ../lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt (18 overlapping tiles of one TMMC over both counts, pooled by
C-weighted least squares; see its header). It holds 9865 states of the grid n + m <= N_TOTAL_MAX, m <= N_PROPANE_MAX
(n CO2, m propane molecules); all other states get Pi_ref = 0. Per its header, a reweighted condition is trustworthy
while the states on the edge of that domain (frontier = 1) carry less than about 1e-3 of the probability: every
composition up to a total fugacity of about 1e6 Pa, and 94 % of them up to 1e7 Pa.

The reference fugacity of the file, 1e5 Pa for both components, is kept (ln Pi spans 433 there; doubles underflow
below e^-708). <H> is the reweighting energy of the file header, U - 3/2 k_B T per propane molecule (the ideal-gas
mean intramolecular energy of propane), instead of the sampled U - U_intra of build_distribution.py; the two differ by
up to 3.5 % at the highest propane loadings. <H> only enters away from 300 K. The ideal-gas term is stored as
'intramolecular' in lnpi_T300K.npz, so that energy - intramolecular is <H>/k_B there too.

Writes particle_distribution_T300K.data and .json and lnpi_T300K.npz in the format of build_distribution.py (which
stops at missing states for the slice windows, README.md); runs/generate.py and ../analysis.ipynb read these files.
"""

import json

import numpy as np

from build_distribution import BOLTZMANN, FRAMEWORK_MASS, HERE, N_PROPANE_MAX, N_TOTAL_MAX

SOURCE = HERE.parent / "lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt"
TEMPERATURE = 300.0
REFERENCE_FUGACITY = 1.0e5  # f_CO2 = f_C3H8 of the file [Pa]

if __name__ == "__main__":
    # columns: N_CO2 N_C3H8 lnpi lnpi_error frontier U_mean[K] U_variance[K^2] U_intra_mean[K] visits
    table = np.loadtxt(SOURCE, comments="#", ndmin=2)
    n, m = table[:, 0].astype(int), table[:, 1].astype(int)
    shape = (N_TOTAL_MAX + 1, N_PROPANE_MAX + 1)
    lnpi = np.full(shape, -np.inf)
    energy, energy_squared = np.full(shape, np.nan), np.full(shape, np.nan)
    lnpi[n, m], energy[n, m], energy_squared[n, m] = table[:, 2], table[:, 5], table[:, 6]
    populated = np.isfinite(lnpi)
    intramolecular = np.where(populated, 1.5 * TEMPERATURE * np.arange(shape[1])[None, :], np.nan)
    label = f"T{TEMPERATURE:.0f}K"
    np.savez(HERE / f"lnpi_{label}.npz", lnpi=lnpi, energy=energy, energy_squared=energy_squared,
             intramolecular=intramolecular, temperature=TEMPERATURE, reference_pressure=REFERENCE_FUGACITY,
             framework_mass=FRAMEWORK_MASS)

    if np.ptp(lnpi[populated]) > 700.0:
        raise RuntimeError(f"ln Pi spans {np.ptp(lnpi[populated]):.0f} at {REFERENCE_FUGACITY} Pa; doubles underflow")
    shifted = lnpi - lnpi[populated].max()
    shifted -= np.log(np.exp(shifted[populated]).sum())
    probability = np.where(populated, np.exp(shifted), 0.0)
    # states without production samples (visits 0, energies nan) get <H> = 0
    without_energy = populated & ~np.isfinite(energy)
    hamiltonian = np.where(populated & ~without_energy, BOLTZMANN * (energy - intramolecular), 0.0)

    with open(HERE / f"particle_distribution_{label}.data", "w") as data:
        data.write(f"# CO2/propane in MOR (2x2x4 unit cells), N-d TMMC ({SOURCE.name}), T = {TEMPERATURE} K\n")
        data.write(f"# reference fugacity {REFERENCE_FUGACITY:.10g} Pa for both components; "
                   f"framework mass {FRAMEWORK_MASS:.10g} kg\n")
        data.write(f"# rows: n_CO2 = 0..{N_TOTAL_MAX}, n_C3H8 = 0..{N_PROPANE_MAX} fastest; Pi_ref, "
                   f"<H> = <U> - 3/2 k_B T n_C3H8 [J]\n")
        for p, h in zip(probability.ravel(), hamiltonian.ravel()):
            data.write(f"{p:.17g} {h:.17g}\n")
    settings = {
        "FileName": f"particle_distribution_{label}.data",
        "ReferenceTemperature": TEMPERATURE,
        "ReferenceFugacity": REFERENCE_FUGACITY,
        "ReferenceFrameworkMass": FRAMEWORK_MASS,
        "ComponentBounds": [{"Component": "CO2", "NMin": 0, "NMax": N_TOTAL_MAX, "DeltaN": 1},
                            {"Component": "C3H8", "NMin": 0, "NMax": N_PROPANE_MAX, "DeltaN": 1}],
    }
    (HERE / f"particle_distribution_{label}.json").write_text(json.dumps(settings, indent=2))
    print(f"{label}: {populated.sum()} states, reference fugacity {REFERENCE_FUGACITY:.4g} Pa, "
          f"ln Pi range {np.ptp(lnpi[populated]):.1f}")
    if without_energy.any():
        print(f"warning: {without_energy.sum()} states have no production energy samples (<H> = 0 written); the "
              f"probability they carry at the reference fugacity is {probability[without_energy].sum():.2e}")
