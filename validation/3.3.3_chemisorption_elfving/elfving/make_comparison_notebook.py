#!/usr/bin/env python3
"""Build the Figure 4-8 Matplotlib overlay notebook."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "elfving_numbered_examples.ipynb"


def source(text: str) -> list[str]:
    return textwrap.dedent(text).lstrip("\n").splitlines(keepends=True)


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source(text)}


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source(text),
    }


cells = [
    markdown(
        r"""
        # Elfving and Sainio (2021): Ruptura overlays for Figures 4-8

        This notebook is organized in the same order as Figures 4-8 of:

        > J. Elfving and T. Sainio, *Kinetic approach to modelling CO$_2$ adsorption from humid air using amine-functionalized resin: Equilibrium isotherms and column dynamics*, Chemical Engineering Science 246 (2021) 116885. https://doi.org/10.1016/j.ces.2021.116885

        Figures 5-8 focus on the requested CO$_2$ breakthrough overlays. Experimental controls are recovered from the red pixels in the published panels, calibrated to the printed axes, and stored as CSV files under `experimental_digitized/`. Markers are the digitized paper controls and lines are RUPTURA. They are graphical reconstructions rather than the authors' original numerical samples.

        The figures label the humidity-generator setpoints as 0.2, 1, and 2 vol-% H$_2$O. The paper also reports that the actual measured averages were 0.16, 0.95, and 1.90 vol-%, respectively. RUPTURA uses the measured values. Plot legends show only the nominal values; the measured values remain documented in the surrounding text and tables.

        Figure 4 is an equilibrium figure rather than a breakthrough figure. Its paper equation and RUPTURA-compatible Toth projection are compared on a logarithmic CO$_2$-pressure scale.

        The dynamic cases use RUPTURA's native chemisorption interface without modifying the solver. CO$_2$ uses the reversible `General` rate $k_a c_{\mathrm{CO_2}}(q_m/2-q)^2-k_dq$, approximating consumption of two free amine sites. The forward and reverse coefficients reproduce the paper model's clean-bed feed rate and feed loading for each experiment. H$_2$O uses the `FirstOrder` rate $k(q_\mathrm{eq}-q)$ and therefore depends linearly on the free surface capacity. It retains the native GAB equilibrium model, the reported LDF coefficient, and the paper's empty-column Langmuir correction.
        """
    ),
    code(
        r"""
        from pathlib import Path
        import html
        import json
        import math
        import subprocess
        import sys

        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        import numpy as np
        from IPython.display import HTML, display
        from matplotlib.lines import Line2D

        ROOT = Path.cwd()
        if (ROOT / "elfving").is_dir():
            ROOT = ROOT / "elfving"
        elif not (ROOT / "generate.py").is_file():
            raise RuntimeError("Open this notebook from the repository root or elfving/")

        sys.path.insert(0, str(ROOT))
        from paper_parameters import MODELS, coadsorption_loading, fit_toth, toth_loading
        from generate import CASES

        RUN_DIRECTORIES = {case.run: case.directory_name for case in CASES}

        if not (ROOT / RUN_DIRECTORIES[11] / "component_1_CO2.data").exists():
            subprocess.run([sys.executable, str(ROOT / "generate.py")], cwd=ROOT, check=True)
            subprocess.run([sys.executable, str(ROOT / "run.py")], cwd=ROOT, check=True)

        plt.style.use("default")
        plt.rcParams.update({
            "figure.dpi": 120,
            "axes.grid": False,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.frameon": True,
            "legend.fontsize": 10,
            "legend.title_fontsize": 10,
        })
        """
    ),
    code(
        r"""
        def read_blocks(path):
            blocks, rows = [], []
            for raw in Path(path).read_text().splitlines():
                line = raw.strip()
                if not line:
                    if rows:
                        blocks.append(rows)
                        rows = []
                elif not line.startswith("#"):
                    rows.append([float(value) for value in line.split()])
            if rows:
                blocks.append(rows)
            return blocks


        def crossing_time(x, y, level):
            for index, value in enumerate(y):
                if value >= level:
                    if index == 0:
                        return x[0]
                    previous = y[index - 1]
                    fraction = 0.0 if value == previous else (level - previous) / (value - previous)
                    return x[index - 1] + fraction * (x[index] - x[index - 1])
            return math.nan


        def load_case(run):
            directory = ROOT / RUN_DIRECTORIES[run]
            simulation = json.loads((directory / "simulation.json").read_text())
            mapping = json.loads((directory / "paper_mapping.json").read_text())
            co2_blocks = read_blocks(directory / "component_1_CO2.data")
            h2o_blocks = read_blocks(directory / "component_2_H2O.data")
            column_blocks = read_blocks(directory / "column.data")
            outlet = [block[-1] for block in co2_blocks]
            h2o_outlet = [block[-1] for block in h2o_blocks]
            time = np.array([60.0 * row[1] for row in outlet])
            normalized = np.array([row[12] for row in outlet])
            h2o_normalized = np.array([row[12] for row in h2o_outlet])
            column_midpoint = 0.5 * simulation["ColumnLength"]
            midpoint_temperature = np.array([
                min(block, key=lambda row: abs(row[2] - column_midpoint))[7]
                for block in column_blocks
            ])
            average_co2_loading = np.array([
                np.mean([row[6] + row[8] for row in block])
                for block in co2_blocks
            ])
            return {
                "run": run,
                "simulation": simulation,
                "mapping": mapping,
                "time": time,
                "normalized": normalized,
                "h2o_normalized": h2o_normalized,
                "midpoint_temperature": midpoint_temperature,
                "average_co2_loading": average_co2_loading,
                "feed_ppm": 1e6 * simulation["Components"][1]["GasPhaseMolFraction"],
                "t01": crossing_time(time, normalized, 0.01),
                "t50": crossing_time(time, normalized, 0.50),
            }


        CASES = {run: load_case(run) for run in range(1, 14)}
        """
    ),
    code(
        r"""
        def html_table(headers, rows):
            head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
            body = "".join(
                "<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>"
                for row in rows
            )
            return HTML(f'''<style>
            .elfving-table{{border-collapse:collapse;font-family:Arial;font-size:13px}}
            .elfving-table th{{background:#25364a;color:white}}
            .elfving-table td,.elfving-table th{{padding:6px 9px;border:1px solid #c7cbd1;text-align:right}}
            .elfving-table td:nth-child(2),.elfving-table th:nth-child(2){{text-align:left}}
            </style><table class="elfving-table"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>''')


        def load_control(filename):
            data = np.genfromtxt(ROOT / "experimental_digitized" / filename, delimiter=",", names=True)
            names = data.dtype.names
            result = {}
            if "time_s" in names:
                result["time_s"] = np.asarray(data["time_s"])
            else:
                result["time_s"] = 60.0 * np.asarray(data["time_min"])
            for name in names:
                if name.startswith("co2_mole_fraction_"):
                    suffix = name.removeprefix("co2_mole_fraction_")
                    result[f"co2_ppm_{suffix}"] = 1e6 * np.asarray(data[name])
                elif name.startswith("co2_ppm_"):
                    result[name] = np.asarray(data[name])
            return result


        NOMINAL_LABELS = {"0p2": "0.2%", "1": "1%", "2": "2%", "experimental": "2%"}
        SHADE_LEVELS = {"0p2": "light", "1": "normal", "2": "dark", "experimental": "dark"}


        def shade_level_color(base, level):
            rgb = np.array(mcolors.to_rgb(base))
            if level == "light":
                rgb = rgb + 0.42 * (1.0 - rgb)
            elif level == "dark":
                rgb = 0.68 * rgb
            return tuple(rgb)


        def shaded_color(base, column):
            return shade_level_color(base, SHADE_LEVELS[column])


        def run_line_style(run_index):
            return ("-", "--", ":")[run_index % 3]


        def breakthrough_overlay(control_file, specifications, extent, legend_location="best"):
            control = load_control(control_file)
            fig, ax = plt.subplots(figsize=(8, 6))
            elfving_handles = []
            ruptura_handles = []
            comparison_rows = []
            for column, condition, runs, marker in specifications:
                co2_color = shaded_color("C0", column)
                h2o_color = shaded_color("C1", column)
                nominal = NOMINAL_LABELS[column]
                marker_stride = max(1, len(control["time_s"]) // 55)
                experimental = ax.scatter(
                    control["time_s"][::marker_stride],
                    (control[f"co2_ppm_{column}"] / CASES[runs[0]]["feed_ppm"])[::marker_stride],
                    edgecolors=co2_color, marker=marker, s=30, linewidths=1.1,
                    facecolors="none",
                    label=f"Elfving: {nominal}", zorder=3,
                )
                elfving_handles.append(experimental)
                for run_index, run in enumerate(runs):
                    case = CASES[run]
                    mark_every = max(1, len(case["time"]) // 30)
                    line, = ax.plot(
                        case["time"], case["normalized"],
                        color=co2_color, linewidth=2.0, linestyle=run_line_style(run_index),
                        marker=marker, markersize=4, markevery=mark_every,
                        label=f"RUPTURA {nominal}", zorder=2,
                    )
                    ruptura_handles.append(line)
                    ax.plot(
                        case["time"], case["h2o_normalized"],
                        color=h2o_color, linewidth=2.0, linestyle=run_line_style(run_index),
                        marker=marker, markersize=4, markevery=mark_every,
                        label="_nolegend_", zorder=1,
                    )
                    experimental_normalized = control[f"co2_ppm_{column}"] / case["feed_ppm"]
                    comparison_rows.append([
                        condition,
                        run,
                        f'{crossing_time(control["time_s"], experimental_normalized, 0.10):.1f}',
                        f'{crossing_time(case["time"], case["normalized"], 0.10):.1f}',
                        f'{crossing_time(control["time_s"], experimental_normalized, 0.50):.1f}',
                        f'{case["t50"]:.1f}',
                    ])
            ax.set(
                xlabel="Time, $t$ / s",
                ylabel="Normalized concentration, $c_i/c_{i,0}$ / -",
                xlim=(extent[0], extent[1]),
                ylim=(0, 1.05),
            )
            handles = elfving_handles + ruptura_handles
            source_legend = ax.legend(
                handles=handles, loc=legend_location, ncol=3 if len(handles) <= 3 else 2
            )
            ax.add_artist(source_legend)
            ax.legend(
                handles=[
                    Line2D([], [], color="C0", linewidth=2, label="CO$_2$"),
                    Line2D([], [], color="C1", linewidth=2, label="H$_2$O"),
                ],
                loc="upper left",
            )
            fig.tight_layout()
            plt.show()
            display(html_table(
                ["Condition", "Run", "Pixel t10 [s]", "RUPTURA t10 [s]", "Pixel t50 [s]", "RUPTURA t50 [s]"],
                comparison_rows,
            ))


        def midpoint_temperature_plot(specifications, time_limit):
            fig, ax = plt.subplots(figsize=(8, 6))
            for column, _condition, runs, marker in specifications:
                color = shaded_color("C2", column)
                nominal = NOMINAL_LABELS[column]
                for run_index, run in enumerate(runs):
                    case = CASES[run]
                    ax.plot(
                        case["time"], case["midpoint_temperature"],
                        color=color, linewidth=2.0, linestyle=run_line_style(run_index),
                        marker=marker, markersize=4,
                        markevery=max(1, len(case["time"]) // 30),
                        label=f"RUPTURA {nominal}",
                    )
            ax.set(
                xlabel="Time, $t$ / s",
                ylabel="Solid temperature, $T_{\\mathrm{s}}$ / K",
                xlim=(0, time_limit),
            )
            ax.legend(loc="best")
            fig.tight_layout()
            plt.show()


        def average_co2_loading_plot(specifications, time_limit):
            fig, ax = plt.subplots(figsize=(8, 6))
            for column, _condition, runs, marker in specifications:
                color = shaded_color("C0", column)
                nominal = NOMINAL_LABELS[column]
                for run_index, run in enumerate(runs):
                    case = CASES[run]
                    ax.plot(
                        case["time"], case["average_co2_loading"],
                        color=color, linewidth=2.0, linestyle=run_line_style(run_index),
                        marker=marker, markersize=4,
                        markevery=max(1, len(case["time"]) // 30),
                        label=f"RUPTURA {nominal}",
                    )
            ax.set(
                xlabel="Time, $t$ / s",
                ylabel="Average CO$_2$ loading, $\\overline{q}_{\\mathrm{CO_2}}$ / mol kg$^{-1}$",
                xlim=(0, time_limit),
            )
            ax.set_ylim(bottom=0)
            ax.legend(loc="best")
            fig.tight_layout()
            plt.show()
        """
    ),
    markdown(
        r"""
        ## Numbered-example map

        The main-paper overlays use runs 1-8 and 10-11. Runs 9 and 12-13 are retained in the directory as supplementary comparisons, but do not correspond to a distinct main-paper breakthrough panel.
        """
    ),
    code(
        r"""
        figure_runs = {
            "Figure 5a": [1, 2, 3],
            "Figure 6a": [4, 5],
            "Figure 7a": [6, 7, 8],
            "Figure 8a": [10, 11],
        }
        rows = []
        for figure, runs in figure_runs.items():
            for run in runs:
                case = CASES[run]
                rows.append([
                    run,
                    figure,
                    case["mapping"]["model"],
                    f'{case["simulation"]["Temperature"] - 273.15:.0f}',
                    f'{case["mapping"]["water_vol_percent"]:.2f}',
                    f'{case["mapping"]["flow_ml_min"]:.0f}',
                    f'{case["t01"]:.1f}',
                    f'{case["t50"]:.1f}',
                ])
        display(html_table(
            ["Run", "Paper panel", "Model", "T [C]", "H₂O [vol-%]", "Flow [ml/min]", "t01 [s]", "t50 [s]"],
            rows,
        ))
        """
    ),
    markdown(
        r"""
        ## Figure 4 - 7-parameter equilibrium isotherms

        Figure 4 contains equilibrium isotherms, not column breakthrough. Solid lines evaluate the paper's relaxed 7-parameter equations (Table 3: `t1=10.56`, `t2=14.42`). Dashed lines are Toth projections compatible with RUPTURA. All three temperatures share one panel: C0/C1/C2 identify the 0.2/1/2 vol-% H$_2$O curves, and light/normal/dark shades of each color identify 25/35/50 °C. CO$_2$ pressure is logarithmic as requested.
        """
    ),
    code(
        r"""
        pressure_bar = np.geomspace(1e-6, 4.2e-3, 220)
        temperatures = [298.15, 308.15, 323.15]
        temperature_shades = {298.15: "light", 308.15: "normal", 323.15: "dark"}
        humidities = [0.0016, 0.0095, 0.0190]
        nominal_humidities = ["0.2%", "1%", "2%"]

        fig, ax = plt.subplots(figsize=(8, 6))
        condition_handles = []
        for index, (humidity, nominal) in enumerate(zip(humidities, nominal_humidities)):
            base = f"C{index}"
            for temperature in temperatures:
                color = shade_level_color(base, temperature_shades[temperature])
                paper = np.array([
                    sum(coadsorption_loading(MODELS["seven_relaxed"], temperature, pressure, humidity * 1.01325))
                    for pressure in pressure_bar
                ])
                qmax, affinity, exponent, rms = fit_toth(pressure_bar.tolist(), paper.tolist())
                ruptura = np.array([toth_loading(qmax, affinity, exponent, pressure) for pressure in pressure_bar])
                ax.plot(pressure_bar, paper, color=color, linewidth=2.0)
                ax.plot(pressure_bar, ruptura, color=color, linewidth=2.0, linestyle="--")
                condition_handles.append(Line2D(
                    [], [], color=color, linewidth=2,
                    label=f"{nominal} H$_2$O, {temperature - 273.15:.0f} °C",
                ))
        ax.set_xscale("log")
        ax.set_xlim(1e-6, 4.2e-3)
        ax.set_ylim(0, 1.6)
        ax.set_xlabel("CO$_2$ partial pressure, $p_{\\mathrm{CO_2}}$ / bar")
        ax.set_ylabel("CO$_2$ loading, $q_{\\mathrm{CO_2}}$ / mol kg$^{-1}$")
        condition_legend = ax.legend(
            handles=condition_handles, loc="upper left", ncol=3, fontsize=8,
            title="C0/C1/C2 = 0.2/1/2% H$_2$O; light/normal/dark = 25/35/50 °C",
        )
        ax.add_artist(condition_legend)
        ax.legend(
            handles=[
                Line2D([], [], color="0.2", linewidth=2, label="Elfving (relaxed 7-parameter)"),
                Line2D([], [], color="0.2", linewidth=2, linestyle="--", label="RUPTURA (Toth projection)"),
            ],
            loc="lower right",
        )
        fig.tight_layout()
        plt.show()
        """
    ),
    markdown(
        r"""
        ## Figure 5 - CO$_2$ breakthrough at 25-26 °C and 500 ml/min

        Runs 1-3 use the measured humidities corresponding to the paper's three nominal setpoints. The temperature diagnostic uses the solid temperature at the grid point nearest the bed midpoint. The loading diagnostic is the arithmetic mean of total adsorbed CO$_2$ over every output grid point.

        ### Figure 5a-b reference - normalized CO$_2$ and H$_2$O breakthrough
        """
    ),
    code(
        r"""
        figure_5_specs = [
            ("0p2", "0.2% nominal / 0.16% measured", [1], "o"),
            ("1", "1% nominal / 0.95% measured", [2], "s"),
            ("2", "2% nominal / 1.90% measured", [3], "D"),
        ]
        breakthrough_overlay(
            "figure_5a.csv",
            figure_5_specs,
            (0, 7200, 0, 400),
            legend_location="lower right",
        )
        """
    ),
    markdown(
        r"""
        ### Figure 5c reference - midpoint solid temperature
        """
    ),
    code(
        r"""
        midpoint_temperature_plot(figure_5_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ### Figure 5d reference - grid-averaged CO$_2$ loading
        """
    ),
    code(
        r"""
        average_co2_loading_plot(figure_5_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ## Figure 6 - CO$_2$ breakthrough at 35 °C and 500 ml/min

        Runs 4 and 5 represent the measured humidities corresponding to the nominal 1 and 2 vol-% setpoints. Temperature is sampled at the bed midpoint and CO$_2$ loading is averaged over all grid points.

        ### Figure 6a-b reference - normalized CO$_2$ and H$_2$O breakthrough
        """
    ),
    code(
        r"""
        figure_6_specs = [
            ("1", "1% nominal / 0.95% measured", [4], "s"),
            ("2", "2% nominal / 1.90% measured", [5], "D"),
        ]
        breakthrough_overlay(
            "figure_6a.csv",
            figure_6_specs,
            (0, 7200, 0, 400),
            legend_location="lower right",
        )
        """
    ),
    markdown(
        r"""
        ### Figure 6c reference - midpoint solid temperature
        """
    ),
    code(
        r"""
        midpoint_temperature_plot(figure_6_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ### Figure 6d reference - grid-averaged CO$_2$ loading
        """
    ),
    code(
        r"""
        average_co2_loading_plot(figure_6_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ## Figure 7 - CO$_2$ breakthrough at 25-26 °C and 500 ml/min

        Runs 6-8 use the bounded 7-parameter projection with `t1=t2=3` and the measured humidity values. Figure 7 in the paper does not include H$_2$O or temperature panels, so those two RUPTURA diagnostics are included here as requested.

        ### Figure 7a reference - normalized CO$_2$ and H$_2$O breakthrough
        """
    ),
    code(
        r"""
        figure_7_specs = [
            ("0p2", "0.2% nominal / 0.16% measured", [6], "o"),
            ("1", "1% nominal / 0.95% measured", [7], "s"),
            ("2", "2% nominal / 1.90% measured", [8], "D"),
        ]
        breakthrough_overlay(
            "figure_7a.csv",
            figure_7_specs,
            (0, 7200, 0, 400),
            legend_location="lower right",
        )
        """
    ),
    markdown(
        r"""
        ### Figure 7 reference - midpoint solid temperature
        """
    ),
    code(
        r"""
        midpoint_temperature_plot(figure_7_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ### Figure 7b reference - grid-averaged CO$_2$ loading
        """
    ),
    code(
        r"""
        average_co2_loading_plot(figure_7_specs, time_limit=7200)
        """
    ),
    markdown(
        r"""
        ## Figure 8 - CO$_2$ breakthrough at 25 °C and 100 ml/min

        Runs 10-11 compare the 5- and bounded 7-parameter RUPTURA projections against the same pixel-derived experiment. The breakthrough plot retains the paper inset's 0-27000 s window. The temperature and loading plots use the full available RUPTURA output horizon of 57600 s.

        ### Figure 8a-b reference - normalized CO$_2$ and H$_2$O breakthrough
        """
    ),
    code(
        r"""
        figure_8_specs = [
            ("experimental", "2% nominal / about 2% measured", [10, 11], "D"),
        ]
        breakthrough_overlay(
            "figure_8a.csv",
            figure_8_specs,
            (0, 27000, 0, 400),
            legend_location="lower right",
        )
        """
    ),
    markdown(
        r"""
        ### Figure 8c reference - midpoint solid temperature
        """
    ),
    code(
        r"""
        midpoint_temperature_plot(figure_8_specs, time_limit=57600)
        """
    ),
    markdown(
        r"""
        ### Figure 8d reference - grid-averaged CO$_2$ loading
        """
    ),
    code(
        r"""
        average_co2_loading_plot(figure_8_specs, time_limit=57600)
        """
    ),
    markdown(
        r"""
        ## Reading the overlays

        - Markers are pixel-digitized Elfving CO$_2$ controls; lines are RUPTURA CO$_2$ and H$_2$O results.
        - C0 identifies CO$_2$ and C1 identifies H$_2$O. Within each color family, 0.2%, 1%, and 2% humidity use light, normal, and dark shades with circle, square, and diamond markers.
        - The temperature curves use solid temperature at the grid point nearest the bed midpoint and remain in kelvin.
        - The adsorption curves average total CO$_2$ loading over all grid points at each output time.
        - Figures 5-7 make the differences in breakthrough timing and front shape directly visible against reconstructed paper data.
        - Figure 8 uses the paper's inset because the full 162000 s axis compresses the breakthrough transition into a narrow region.
        - The tables compare 10% and 50% breakthrough times. The 1% level is intentionally omitted because it lies too close to the rendered baseline for reliable pixel recovery.
        - Treat the recovered values as digitization-based controls, not original experimental precision.
        """
    ),
]

for index, cell in enumerate(cells):
    cell["id"] = f"elfving-{index:02d}"


notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

OUTPUT.write_text(json.dumps(notebook, indent=1) + "\n")
print(f"Wrote {OUTPUT}")
