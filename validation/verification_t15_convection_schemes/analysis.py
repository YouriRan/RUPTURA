#!/usr/bin/env python3
"""Print the CO2 front width t90 - t10 at z = L and the ROS34PW2 step counts of every scheme.

    python analysis.py

t10 and t90: the first output times at which y_CO2(L) reaches 10 and 90 % of its feed value, interpolated linearly.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent / "simulations"
FEED_CO2 = 0.4


def crossing(time: np.ndarray, value: np.ndarray, level: float) -> float:
    i = np.nonzero((value[:-1] < level) & (value[1:] >= level))[0][0]
    return time[i] + (level - value[i]) * (time[i + 1] - time[i]) / (value[i + 1] - value[i])


def main() -> None:
    for scheme in ("Upwind", "TVD", "WENO3", "WENO5"):
        line = []
        for integrator in ("ROS34PW2", "CVODE"):
            directory = ROOT / scheme / integrator
            data = np.loadtxt(directory / "component_1_CO2.data")
            outlet = data[np.isclose(data[:, 2], data[:, 2].max())]
            time, value = outlet[:, 1] * 60.0, outlet[:, 5] / FEED_CO2
            width = crossing(time, value, 0.9) - crossing(time, value, 0.1)
            steps = re.search(r"^\s+steps:\s+(\d+)", (directory / "run.log").read_text(), re.M)
            line.append(f"{integrator} {width:7.2f} s" + (f" ({steps.group(1)} steps)" if steps else ""))
        print(f"{scheme:7s} " + "; ".join(line))


if __name__ == "__main__":
    main()
