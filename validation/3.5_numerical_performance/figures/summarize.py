#!/usr/bin/env python3
"""Summarize repeated timings and matched CVODE/RK3 mole-fraction errors."""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path


FIGURE_ROOT = Path(__file__).resolve().parent
SIMULATION_ROOT = FIGURE_ROOT / "simulations"


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def final_mole_fractions(directory: Path) -> dict[tuple[int, float], float]:
    values = {}
    for component_path in sorted(directory.glob("component_*.data")):
        component_index = int(component_path.name.split("_", 2)[1])
        observations = []
        with component_path.open(encoding="utf-8") as stream:
            for line in stream:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split()
                observations.append((float(fields[1]), float(fields[2]), float(fields[5])))
        if not observations:
            continue
        final_time = max(time for time, _, _ in observations)
        for time, position, mole_fraction in observations:
            if math.isclose(time, final_time, rel_tol=0.0, abs_tol=1e-12):
                values[(component_index, position)] = mole_fraction
    return values


def main() -> None:
    manifest = read_rows(FIGURE_ROOT / "manifest.csv")
    results = read_rows(FIGURE_ROOT / "results.csv")
    manifest_by_directory = {row["directory"]: row for row in manifest}

    elapsed = defaultdict(list)
    for row in results:
        if int(row["return_code"]) == 0:
            elapsed[row["directory"]].append(float(row["elapsed_seconds"]))

    timing_rows = []
    for directory, repetitions in elapsed.items():
        timing_rows.append(
            {
                **manifest_by_directory[directory],
                "successful_repetitions": len(repetitions),
                "median_wall_seconds": statistics.median(repetitions),
                "min_wall_seconds": min(repetitions),
                "max_wall_seconds": max(repetitions),
            }
        )
    timing_rows.sort(key=lambda row: row["directory"])
    timing_path = FIGURE_ROOT / "timing_summary.csv"
    with timing_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=timing_rows[0].keys())
        writer.writeheader()
        writer.writerows(timing_rows)

    paired = defaultdict(dict)
    for row in timing_rows:
        if row["suite"] != "cvode-rk3":
            continue
        key = (row["case"], row["grid_points"], row["ncomp"], row["time_step"])
        paired[key][row["integrator"]] = row

    error_rows = []
    for key, pair in sorted(paired.items()):
        if set(pair) != {"rk3", "cvode"}:
            continue
        rk3_values = final_mole_fractions(SIMULATION_ROOT / pair["rk3"]["directory"])
        cvode_values = final_mole_fractions(SIMULATION_ROOT / pair["cvode"]["directory"])
        common = sorted(set(rk3_values) & set(cvode_values))
        differences = [abs(cvode_values[index] - rk3_values[index]) for index in common]
        error_rows.append(
            {
                "case": key[0],
                "grid_points": key[1],
                "ncomp": key[2],
                "time_step": key[3],
                "compared_mole_fractions": len(differences),
                "max_abs_mole_fraction_error": max(differences) if differences else "",
                "rms_mole_fraction_error": (
                    math.sqrt(sum(value * value for value in differences) / len(differences)) if differences else ""
                ),
            }
        )

    error_path = FIGURE_ROOT / "errors.csv"
    with error_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=error_rows[0].keys())
        writer.writeheader()
        writer.writerows(error_rows)
    print(f"Wrote {len(timing_rows)} timing rows and {len(error_rows)} paired error rows")


if __name__ == "__main__":
    main()
