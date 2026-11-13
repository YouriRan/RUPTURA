#!/usr/bin/env python3
"""Headless figures from the original 12/14-column breakthrough text files."""
from __future__ import annotations

import json
from pathlib import Path


def plot_legacy_results(input_path, output_dir=None, formats=("png", "pdf"), dpi=160):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    input_path = Path(input_path).resolve()
    document = json.loads(input_path.read_text())
    if document.get("SimulationType", "Breakthrough") != "Breakthrough":
        raise ValueError("Legacy static plotting supports Breakthrough inputs")
    source = input_path.parent
    output = Path(output_dir) if output_dir is not None else Path(
        document.get("Output", {}).get("Directory", "results"))
    if not output.is_absolute():
        output = source / output
    plots = output / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    def read(path, columns):
        values = np.loadtxt(path, comments="#", ndmin=2)
        if values.shape[1] < columns or not len(values):
            raise ValueError(f"{path} requires at least {columns} columns and one recorded row")
        if not np.isfinite(values).all():
            raise ValueError(f"{path} contains nonfinite recorded values")
        return values

    def outlet(values):
        return values[np.isclose(values[:, 2], np.max(values[:, 2]), rtol=1e-12, atol=1e-12)]

    def final(values):
        selected = values[values[:, 1] == np.max(values[:, 1])]
        return selected[np.argsort(selected[:, 2], kind="stable")]

    column = read(source / "column.data", 12)
    column_outlet, column_final = outlet(column), final(column)
    components = []
    for index, component in enumerate(document["Components"]):
        name = component.get("Name", f"component {index}")
        values = read(source / f"component_{index}_{name}.data", 14)
        components.append((component, values, outlet(values), final(values)))
    if not components:
        raise ValueError("The legacy input defines no components to plot")
    liquid = document.get("FluidPhase", "Gas") == "Liquid"
    saved = []
    context = {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
               "figure.constrained_layout.use": True, "savefig.dpi": dpi}

    def save(figure, name):
        figure.suptitle(document.get("DisplayName", "Legacy breakthrough"))
        figure.text(.5, .006, "Legacy breakthrough transient; figures use the recorded text results.",
                    ha="center", va="bottom", fontsize=7)
        figure.get_layout_engine().set(rect=(0, .07, 1, .93))
        for extension in formats:
            target = plots / f"{name}.{extension}"
            figure.savefig(target)
            saved.append(str(target))
        plt.close(figure)

    with plt.rc_context(context):
        figure, axis = plt.subplots(figsize=(8, 4))
        for component, _, end, _ in components:
            axis.plot(end[:, 1], end[:, 3 if liquid else 5], label=component["Name"])
        axis.set(xlabel="Time [min]",
                 ylabel="Solute concentration [mol/m³ liquid]" if liquid else "Outlet mole fraction [-]")
        axis.legend()
        save(figure, "legacy_breakthrough")

        fed = [(component, values, end, profile) for component, values, end, profile in components
               if float(component.get("LiquidPhaseConcentration" if liquid else "GasPhaseMolFraction", 0)) > 0]
        if fed:
            figure, axis = plt.subplots(figsize=(8, 4))
            for component, _, end, _ in fed:
                axis.plot(end[:, 1], end[:, 12], label=component["Name"])
            axis.set(xlabel="Time [min]", ylabel="Outlet c/cfeed [-]" if liquid else "Outlet y/yfeed [-]")
            axis.legend()
            save(figure, "legacy_normalized_breakthrough")

        figure, axes = plt.subplots(1, 2, figsize=(11, 4))
        for component, _, _, profile in components:
            line, = axes[0].plot(profile[:, 2], profile[:, 3 if liquid else 5], label=component["Name"])
            axes[1].plot(profile[:, 2], profile[:, 6] + profile[:, 8],
                         color=line.get_color(), label=component["Name"])
        axes[0].set(xlabel="Axial position [m]",
                    ylabel="Solute concentration [mol/m³ liquid]" if liquid else "Mole fraction [-]")
        axes[1].set(xlabel="Axial position [m]", ylabel="Total loading [mol/kg adsorbent]")
        for axis in axes:
            axis.legend()
        save(figure, "legacy_final_profiles")

        figure, axis = plt.subplots(figsize=(8, 4))
        for index, label in [(5, "Liquid" if liquid else "Gas"), (7, "Solid"), (9, "Wall")]:
            axis.plot(column_outlet[:, 1], column_outlet[:, index], label=label)
        axis.set(xlabel="Time [min]", ylabel="Outlet temperature [K]")
        axis.legend()
        save(figure, "legacy_temperature_history")

        figure, axis = plt.subplots(figsize=(8, 4))
        for index, label in [(5, "Liquid" if liquid else "Gas"), (7, "Solid"), (9, "Wall")]:
            axis.plot(column_final[:, 2], column_final[:, index], label=label)
        axis.set(xlabel="Axial position [m]", ylabel="Final temperature [K]")
        axis.legend()
        save(figure, "legacy_temperature_profiles")

        figure, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].plot(column_outlet[:, 1], column_outlet[:, 4] / 1e5)
        axes[0].set(xlabel="Time [min]", ylabel="Outlet pressure [bar]")
        axes[1].plot(column_outlet[:, 1], column_outlet[:, 3])
        axes[1].set(xlabel="Time [min]", ylabel="Outlet interstitial velocity [m/s]")
        save(figure, "legacy_pressure_velocity")

    manifest = {
        "backend": "Matplotlib", "legacy": True, "formats": list(formats),
        "files": [Path(path).name for path in saved],
        "source_input": input_path.name,
        "source_outputs": ["column.data"] + [
            f"component_{index}_{item['Name']}.data" for index, item in enumerate(document["Components"])],
        "final_time_min": float(np.max(column[:, 1])),
        "notes": [
            "Original text output schema: time is in minutes; axial position is in metres.",
            "Mole fractions, normalized driving forces and loadings are read from the recorded component files.",
            "The numerical solver and physical balances are not executed by this plotting routine.",
        ],
    }
    (plots / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return saved

