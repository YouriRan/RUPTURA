#!/usr/bin/env python3
"""Compare Ruptura's BEA multinomial-MPD and Langmuir output files."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MPD_DIRECTORY = ROOT / "mpd"
LANGMUIR_DIRECTORY = ROOT / "langmuir"


def component_names() -> list[str]:
    with (MPD_DIRECTORY / "simulation.json").open(encoding="utf-8") as stream:
        simulation = json.load(stream)
    return [component["Name"] for component in simulation["Components"] if not component.get("CarrierGas", False)]


def read_mixture_loading(directory: Path, component_index: int, name: str) -> tuple[list[float], list[float]]:
    pressures: list[float] = []
    loadings: list[float] = []
    path = directory / f"component_{component_index}_{name}.data"
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            columns = line.split()
            pressures.append(float(columns[0]))
            loadings.append(float(columns[2]))
    return pressures, loadings


def main() -> None:
    names = component_names()
    rows: list[list[float]] = []
    series: dict[str, tuple[list[float], list[float], list[float]]] = {}
    maximum_absolute_error = 0.0
    maximum_scaled_error = 0.0

    for component_index, name in enumerate(names, start=1):
        pressures, mpd = read_mixture_loading(MPD_DIRECTORY, component_index, name)
        langmuir_pressures, langmuir = read_mixture_loading(LANGMUIR_DIRECTORY, component_index, name)
        if pressures != langmuir_pressures:
            raise RuntimeError(f"Pressure grids differ for {name}")

        absolute_errors = [abs(left - right) for left, right in zip(mpd, langmuir, strict=True)]
        scale = max(max(map(abs, langmuir)), 1.0)
        component_absolute_error = max(absolute_errors)
        component_scaled_error = component_absolute_error / scale
        maximum_absolute_error = max(maximum_absolute_error, component_absolute_error)
        maximum_scaled_error = max(maximum_scaled_error, component_scaled_error)
        series[name] = (pressures, mpd, langmuir)
        print(f"{name:8s} max |MPD-Langmuir| = {component_absolute_error:.3e} mol/kg")

    header = ["fugacity_Pa"]
    for name in names:
        header.extend((f"{name}_MPD_mol_kg", f"{name}_Langmuir_mol_kg", f"{name}_absolute_error"))
    pressures = series[names[0]][0]
    for row_index, pressure in enumerate(pressures):
        row = [pressure]
        for name in names:
            mpd = series[name][1][row_index]
            langmuir = series[name][2][row_index]
            row.extend((mpd, langmuir, abs(mpd - langmuir)))
        rows.append(row)

    with (ROOT / "comparison.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipped comparison.png")
    else:
        figure, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
        for name in names:
            pressure, mpd, langmuir = series[name]
            line = axis.plot(pressure, mpd, linewidth=2, label=f"{name} MPD")[0]
            axis.plot(
                pressure,
                langmuir,
                linestyle="none",
                marker="o",
                markersize=3,
                markevery=4,
                color=line.get_color(),
                label=f"{name} Langmuir",
            )
        axis.set_xscale("log")
        axis.set_xlabel("Total gas fugacity, f [Pa]")
        axis.set_ylabel("Mixture component loading, qᵢ [mol kg⁻¹]")
        axis.set_title("BEA C7 alkanes: multinomial MPD vs. competitive Langmuir")
        axis.grid(alpha=0.2)
        axis.legend(ncol=2, fontsize=8)
        figure.savefig(ROOT / "comparison.png", dpi=180)
        plt.close(figure)

    print(f"Overall maximum absolute error: {maximum_absolute_error:.3e} mol/kg")
    print(f"Overall maximum scaled error:   {maximum_scaled_error:.3e}")
    if maximum_scaled_error > 1.0e-10:
        raise RuntimeError("MPD and competitive Langmuir results do not agree within tolerance")


if __name__ == "__main__":
    main()
