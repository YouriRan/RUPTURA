#!/usr/bin/env python3
"""Print the numbers of the T15 checks: breakthrough times, outlet deviations, temperature peaks, pressure drops,
step counts, and the times at which the three integrators end the phase.

    python analysis.py

Breakthrough time t50: the first output time at which y_j(L) reaches half its feed value, interpolated linearly.
Pressure drop: between the first and the last interior node at the end of the phase.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent / "simulations"
COMPONENTS = {"CO2": (1, 0.4), "N2": (2, 0.3), "H2": (3, 0.3)}


def outlet(directory: Path, component: str) -> tuple[np.ndarray, np.ndarray]:
    index, _ = COMPONENTS[component]
    data = np.loadtxt(directory / f"component_{index}_{component}.data")
    rows = data[np.isclose(data[:, 2], data[:, 2].max())]
    return rows[:, 1] * 60.0, rows[:, 5]


def t50(time: np.ndarray, value: np.ndarray, half: float) -> float:
    i = np.nonzero((value[:-1] < half) & (value[1:] >= half))[0][0]
    return time[i] + (half - value[i]) * (time[i + 1] - time[i]) / (value[i + 1] - value[i])


def statistics(directory: Path) -> str:
    log = (directory / "run.log").read_text()
    steps = re.search(r"^\s+steps:\s+(\d+)", log, re.M)
    total = re.search(r"Total time:\s+([0-9.]+)", log)
    return f"steps {steps.group(1) if steps else '-'}, {total.group(1) if total else '-'} s"


def main() -> None:
    for thermal in ("isothermal", "energy_balance"):
        print(f"== free gas pressure against the continuity closure, {thermal}")
        runs = {name: ROOT / "free_pressure" / name / thermal for name in ("momentum_balance", "continuity_closure")}
        deviation = 0.0
        for name, directory in runs.items():
            column = np.loadtxt(directory / "column.data")
            last = column[np.isclose(column[:, 1], column[:, 1].max())]
            # Between the first and the last interior node: the end rows hold the port (free gas pressure) or the
            # boundary value (closure), which are not the same quantity.
            drop = last[1, 4] - last[-2, 4]
            times = ", ".join(f"{component} {t50(*outlet(directory, component), 0.5 * feed):.2f}"
                              for component, (_, feed) in COMPONENTS.items())
            print(f"  {name:20s} t50 {times} s; peak T_g {column[:, 5].max():.2f} K, T_s {column[:, 7].max():.2f} K; "
                  f"pressure drop {drop:.2f} Pa; {statistics(directory)}")
        for component in COMPONENTS:
            (time_a, y_a), (time_b, y_b) = (outlet(directory, component) for directory in runs.values())
            n = min(len(time_a), len(time_b))
            deviation = max(deviation, float(np.abs(y_a[:n] - y_b[:n]).max()))
        print(f"  largest outlet mole-fraction difference {deviation:.1e}")
    klu = ROOT / "free_pressure" / "momentum_balance" / "isothermal_cvode_klu"
    print(f"== CVODE KLU, free gas pressure, isothermal: {statistics(klu)}; last line: "
          f"{(klu / 'run.log').read_text().strip().splitlines()[-1][:160]}")

    print("== phase end at y_CO2(L) = 0.2")
    ends = {}
    for directory in sorted((ROOT / "phase_end").iterdir()):
        phases = np.atleast_2d(np.loadtxt(directory / "phases.data"))
        ends[directory.name] = phases[0, 1]
        print(f"  {directory.name:12s} {ends[directory.name]:.4f} s; {statistics(directory)}")
    print(f"  spread {1e3 * (max(ends.values()) - min(ends.values())):.2f} ms")


if __name__ == "__main__":
    main()
