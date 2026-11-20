"""Stitch the TMMC windows of one temperature into the MPD reference distribution read by RUPTURA.

Usage: python build_distribution.py [300]  (temperature in K, default 300)

With n CO2 and m propane molecules (../tmmc/generate.py),

  ln Pi(n, m) = ln Pi_C3H8(m) + ln Pi_CO2|m(n) - ln Pi_CO2|m(0),

at the fugacity REFERENCE_PRESSURE of the TMMC runs for both components. The windows of one TMMC are joined by the
mean ln Pi difference over their overlapping states. Energies are averaged over the production samples of all
windows that contain a state; the mean Hamiltonian <H> of RUPTURA is the intermolecular and framework part,
U - U_intra, because the intramolecular energy of flexible propane is the same as in the ideal gas.

Writes particle_distribution_T<T>K.data (Pi_ref and <H> in J; C order, propane count fastest; states with
n + m > N_TOTAL_MAX have Pi_ref = 0) and particle_distribution_T<T>K.json (the MPDSettings for RUPTURA). The
reference fugacity is chosen such that no state underflows a double. Also writes lnpi_T<T>K.npz for the analysis.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
TMMC = HERE.parent / "tmmc"
sys.path.insert(0, str(TMMC))
from generate import N_PROPANE_MAX, N_TOTAL_MAX, REFERENCE_PRESSURE  # noqa: E402

BOLTZMANN = 1.380649e-23
AVOGADRO = 6.02214076e23
# 2x2x4 unit cells of Si48O96, with the RASPA3 pseudo-atom masses
FRAMEWORK_MASS = 16 * (48 * 28.0855 + 96 * 15.9994) * 1.0e-3 / AVOGADRO


def stitch(directories):
    """ln Pi, sample count, <U>, <dU^2> and <U_intra> over 0..top from the overlapping windows of one TMMC."""
    lnpi, samples, energy, energy_squared, intramolecular = {}, {}, {}, {}, {}
    for directory in sorted(directories):
        lines = (directory / "tmmc" / "tmmc_statistics.txt").read_text().splitlines()
        header = next(i for i, line in enumerate(lines) if line.startswith("N CM[-1]"))
        statistics = np.array([[float(x) for x in line.split()] for line in lines[header + 1:] if line.strip()])
        window = dict(zip(statistics[:, 0].astype(int), statistics[:, lines[header].split().index("lnpi")]))
        if not all(np.isfinite(list(window.values()))):
            raise RuntimeError(f"{directory}: ln Pi is not measured on every state")
        overlap = [n for n in window if n in lnpi]
        if lnpi and not overlap:
            raise RuntimeError(f"{directory}: no overlap with the previous window")
        offset = np.mean([lnpi[n] - window[n] for n in overlap]) if overlap else 0.0
        for n, value in window.items():
            lnpi[n] = 0.5 * (lnpi[n] + value + offset) if n in lnpi else value + offset

        energies = np.genfromtxt(directory / "tmmc" / "tmmc_energies.txt", comments="#", skip_header=4, ndmin=2)
        for n, count, mean, variance, _, mean_intramolecular in energies:
            if count == 0:
                continue
            n = int(n)
            total = samples.get(n, 0.0) + count
            # pooled mean and variance of the windows that share the state
            old_mean = energy.get(n, mean)
            new_mean = (samples.get(n, 0.0) * old_mean + count * mean) / total
            energy_squared[n] = (samples.get(n, 0.0) * (energy_squared.get(n, 0.0) + (old_mean - new_mean) ** 2)
                                 + count * (variance + (mean - new_mean) ** 2)) / total
            intramolecular[n] = (samples.get(n, 0.0) * intramolecular.get(n, 0.0) + count * mean_intramolecular) / total
            energy[n], samples[n] = new_mean, total
    states = np.arange(max(lnpi) + 1)
    if sorted(lnpi) != list(states):
        raise RuntimeError(f"{directories[0].parent}: missing states in {directories[0].name.rsplit('_N', 1)[0]}")
    # a state that no window visited in production has ln Pi (equilibration counts) but no energies: nan
    return (np.array([lnpi[n] for n in states]), np.array([samples.get(n, 0.0) for n in states]),
            np.array([energy.get(n, np.nan) for n in states]), np.array([energy_squared.get(n, np.nan) for n in states]),
            np.array([intramolecular.get(n, np.nan) for n in states]))


def distribution(temperature):
    """ln Pi[n, m], <U>[n, m], <dU^2>[n, m] and <U_intra>[n, m] in K at REFERENCE_PRESSURE; -inf/nan outside."""
    runs = TMMC / f"T{temperature:.0f}K"
    shape = (N_TOTAL_MAX + 1, N_PROPANE_MAX + 1)
    lnpi = np.full(shape, -np.inf)
    energy, energy_squared, intramolecular = np.full(shape, np.nan), np.full(shape, np.nan), np.full(shape, np.nan)
    edge = stitch(list(runs.glob("edge_C3H8_N*")))[0]
    for m in range(N_PROPANE_MAX + 1):
        slice_lnpi, _, slice_energy, slice_energy_squared, slice_intramolecular = stitch(
            list(runs.glob(f"slice_C3H8_{m:03d}_CO2_N*")))
        top = N_TOTAL_MAX - m
        lnpi[:top + 1, m] = edge[m] + slice_lnpi - slice_lnpi[0]
        energy[:top + 1, m] = slice_energy
        energy_squared[:top + 1, m] = slice_energy_squared
        intramolecular[:top + 1, m] = slice_intramolecular
    return lnpi, energy, energy_squared, intramolecular


if __name__ == "__main__":
    temperature = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    lnpi, energy, energy_squared, intramolecular = distribution(temperature)
    label = f"T{temperature:.0f}K"
    np.savez(HERE / f"lnpi_{label}.npz", lnpi=lnpi, energy=energy, energy_squared=energy_squared,
             intramolecular=intramolecular, temperature=temperature, reference_pressure=REFERENCE_PRESSURE,
             framework_mass=FRAMEWORK_MASS)

    # the reference fugacity with the smallest ln Pi range keeps every populated state above the double underflow
    populated = np.isfinite(lnpi)
    molecules = np.add.outer(np.arange(lnpi.shape[0]), np.arange(lnpi.shape[1]))
    fugacities = np.logspace(2, 9, 701)
    ranges = [np.ptp(lnpi[populated] + molecules[populated] * np.log(f / REFERENCE_PRESSURE)) for f in fugacities]
    reference_fugacity = float(fugacities[np.argmin(ranges)])
    if min(ranges) > 700.0:
        raise RuntimeError(f"ln Pi spans {min(ranges):.0f} at every reference fugacity; doubles underflow")
    shifted = lnpi + molecules * np.log(reference_fugacity / REFERENCE_PRESSURE)
    shifted -= shifted[populated].max()
    shifted -= np.log(np.exp(shifted[populated]).sum())
    probability = np.where(populated, np.exp(shifted), 0.0)
    # <H> only enters away from the reference temperature; states without production samples get 0
    without_energy = populated & ~np.isfinite(energy)
    hamiltonian = np.where(populated & ~without_energy, BOLTZMANN * (energy - intramolecular), 0.0)

    with open(HERE / f"particle_distribution_{label}.data", "w") as data:
        data.write(f"# CO2/propane in MOR (2x2x4 unit cells), TMMC slices with an edge (../tmmc), T = {temperature} K\n")
        data.write(f"# reference fugacity {reference_fugacity:.10g} Pa for both components; "
                   f"framework mass {FRAMEWORK_MASS:.10g} kg\n")
        data.write(f"# rows: n_CO2 = 0..{N_TOTAL_MAX}, n_C3H8 = 0..{N_PROPANE_MAX} fastest; Pi_ref, <H> = <U - U_intra> [J]\n")
        for p, h in zip(probability.ravel(), hamiltonian.ravel()):
            data.write(f"{p:.17g} {h:.17g}\n")
    settings = {
        "FileName": f"particle_distribution_{label}.data",
        "ReferenceTemperature": temperature,
        "ReferenceFugacity": reference_fugacity,
        "ReferenceFrameworkMass": FRAMEWORK_MASS,
        "ComponentBounds": [{"Component": "CO2", "NMin": 0, "NMax": N_TOTAL_MAX, "DeltaN": 1},
                            {"Component": "C3H8", "NMin": 0, "NMax": N_PROPANE_MAX, "DeltaN": 1}],
    }
    (HERE / f"particle_distribution_{label}.json").write_text(json.dumps(settings, indent=2))
    print(f"{label}: {populated.sum()} states, reference fugacity {reference_fugacity:.4g} Pa, "
          f"ln Pi range {min(ranges):.1f}")
    if without_energy.any():
        print(f"warning: {without_energy.sum()} states have no production energy samples (<H> = 0 written); the "
              f"probability they carry at the reference fugacity is {probability[without_energy].sum():.2e}")
