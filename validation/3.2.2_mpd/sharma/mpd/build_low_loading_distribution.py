"""Binary distribution from the well-sampled low-loading states of the 300 K windows.

Usage: python build_low_loading_distribution.py

Most slices with propane have not converged (README.md), so build_distribution.py stops. This keeps only the first
window of the edge and of the slices with m = 0..M_MAX propane molecules. Each slice is used from n = 0 up to its first
link with zero recorded acceptance or with a state visited fewer than MINIMUM_VISITS times; all other states get
Pi_ref = 0. Valid only at the reference temperature (<H> is written as 0) and only where the distribution stays away
from these cuts: with 300 to 10000 minimum visits and M_MAX from 15 to 42, the equimolar loadings change by less than
1 % up to 400 Pa and by 15 % at 500 Pa.

Writes particle_distribution_T300K_low_loading.data and .json in the format of build_distribution.py.
"""

import json

import numpy as np

from build_distribution import FRAMEWORK_MASS, HERE, REFERENCE_PRESSURE, TMMC

M_MAX = 42
N_MAX = 42
MINIMUM_VISITS = 1000

if __name__ == "__main__":
    statistics = {}
    for name in ["edge_C3H8_N000-042"] + [f"slice_C3H8_{m:03d}_CO2_N000-042" for m in range(M_MAX + 1)]:
        lines = (TMMC / "T300K" / name / "tmmc" / "tmmc_statistics.txt").read_text().splitlines()
        header = next(i for i, line in enumerate(lines) if line.startswith("N CM[-1]"))
        statistics[name] = np.array([[float(x) for x in line.split()] for line in lines[header + 1:] if line.strip()])

    edge = statistics["edge_C3H8_N000-042"][:, 5]
    lnpi = np.full((N_MAX + 1, M_MAX + 1), -np.inf)
    for m in range(M_MAX + 1):
        # columns: 1 CM[-1], 3 CM[+1], 5 lnpi, 8 histogram
        window = statistics[f"slice_C3H8_{m:03d}_CO2_N000-042"]
        defects = np.flatnonzero((window[:-1, 3] <= 1.0e-6) | (window[1:, 1] <= 1.0e-6)
                                 | (window[:-1, 8] < MINIMUM_VISITS) | (window[1:, 8] < MINIMUM_VISITS))
        top = defects[0] if len(defects) else N_MAX
        lnpi[:top + 1, m] = edge[m] - edge[0] + window[:top + 1, 5] - window[0, 5]

    populated = np.isfinite(lnpi)
    molecules = np.add.outer(np.arange(N_MAX + 1), np.arange(M_MAX + 1))
    fugacities = np.logspace(2, 9, 701)
    ranges = [np.ptp(lnpi[populated] + molecules[populated] * np.log(f / REFERENCE_PRESSURE)) for f in fugacities]
    reference_fugacity = float(fugacities[np.argmin(ranges)])
    shifted = lnpi + molecules * np.log(reference_fugacity / REFERENCE_PRESSURE)
    shifted -= shifted[populated].max()
    shifted -= np.log(np.exp(shifted[populated]).sum())
    probability = np.where(populated, np.exp(shifted), 0.0)

    file_name = "particle_distribution_T300K_low_loading.data"
    with open(HERE / file_name, "w") as data:
        data.write("# CO2/propane in MOR (2x2x4 unit cells), well-sampled low-loading states only, T = 300 K\n")
        data.write(f"# reference fugacity {reference_fugacity:.10g} Pa for both components; "
                   f"framework mass {FRAMEWORK_MASS:.10g} kg\n")
        data.write(f"# rows: n_CO2 = 0..{N_MAX}, n_C3H8 = 0..{M_MAX} fastest; Pi_ref, <H> = 0 (T_ref only)\n")
        for p in probability.ravel():
            data.write(f"{p:.17g} 0\n")
    settings = {
        "FileName": file_name,
        "ReferenceTemperature": 300.0,
        "ReferenceFugacity": reference_fugacity,
        "ReferenceFrameworkMass": FRAMEWORK_MASS,
        "ComponentBounds": [{"Component": "CO2", "NMin": 0, "NMax": N_MAX, "DeltaN": 1},
                            {"Component": "C3H8", "NMin": 0, "NMax": M_MAX, "DeltaN": 1}],
    }
    (HERE / "particle_distribution_T300K_low_loading.json").write_text(json.dumps(settings, indent=2))
    print(f"{populated.sum()} states kept, reference fugacity {reference_fugacity:.4g} Pa, ln Pi range {min(ranges):.1f}")
