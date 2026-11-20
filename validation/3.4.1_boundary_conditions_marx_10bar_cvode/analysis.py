#!/usr/bin/env python3
"""How far each CVODE run got, and how its total gas concentration compares with the ROS34PW2 run up to there.

    python analysis.py

The total gas concentration is the sum of the component concentrations at every node; the difference is the largest
over the nodes and the output times both runs share, relative to the feed concentration P_in / (R T_in).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
GAS_CONSTANT = 8.31446261815324  # J mol-1 K-1, the value in src/utils.h


def total_concentration(directory: Path) -> dict[float, np.ndarray]:
    """Output time [s] -> total gas concentration at every node [mol/m3]."""
    total = {}
    for path in sorted(directory.glob("component_*_*.data")):
        data = np.loadtxt(path)
        for time in np.unique(data[:, 1]):
            rows = data[data[:, 1] == time]
            key = round(60.0 * time, 6)
            total[key] = total.get(key, 0.0) + rows[:, 3]
    return total


def main() -> None:
    for entry in json.loads((ROOT / "simulations" / "index.json").read_text()):
        directory = ROOT / entry["directory"]
        reference = ROOT.parent / entry["reference"]
        log = (directory / "run.log").read_text()
        times = [float(value) for value in re.findall(r"Output \d+, time: ([0-9.]+)", log)]
        reached = times[-1] if times else 0.0
        simulation = json.loads((directory / "simulation.json").read_text())
        feed = simulation["InletPressure"] / (GAS_CONSTANT * simulation["InfluxTemperature"])
        cvode, ros = total_concentration(directory), total_concentration(reference)
        shared = sorted(set(cvode) & set(ros))
        deviation = max((np.abs(cvode[t] - ros[t]).max() / feed for t in shared), default=float("nan"))
        print(f"{entry['directory']:70s} reached {reached:7.1f} s; max |dc_tot| / c_feed {deviation:.1e} over "
              f"{len(shared)} outputs")


if __name__ == "__main__":
    main()
