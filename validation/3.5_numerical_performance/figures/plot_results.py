"""Plot helpers for the figure benchmark results."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


FIGURE_ROOT = Path(__file__).resolve().parent
CASE_LABELS = {
    "bea-alkanes-C7": "BEA-alkanes-C7",
    "marx-403030_298_1e6_2e-1": r"marx-403030_298_1e6_2e-1",
}
DEFAULT_NCOMP = {"bea-alkanes-C7": 7, "marx-403030_298_1e6_2e-1": 3}
NCOMP_SWEEP = {"bea-alkanes-C7": (2, 4, 6), "marx-403030_298_1e6_2e-1": (1, 2, 3)}
GRID_SWEEP = (10, 20, 100)
METHODS = ("IAST", "SIAST", "SEI", "SPI")


def load_figure_data(root: Path = FIGURE_ROOT) -> tuple[pd.DataFrame, pd.DataFrame]:
    timing = pd.read_csv(root / "timing_summary.csv")
    errors = pd.read_csv(root / "errors.csv")
    for frame in (timing, errors):
        for column in ("grid_points", "ncomp", "time_step"):
            if column in frame:
                frame[column] = pd.to_numeric(frame[column])
    for column in ("median_wall_seconds", "min_wall_seconds", "max_wall_seconds"):
        timing[column] = pd.to_numeric(timing[column])
    for column in ("max_abs_mole_fraction_error", "rms_mole_fraction_error"):
        errors[column] = pd.to_numeric(errors[column])
    return timing, errors


def settings_for_case(case: str) -> list[dict]:
    default_ncomp = DEFAULT_NCOMP[case]
    settings = [
        {
            "sweep": "grid",
            "value": grid,
            "grid_points": grid,
            "ncomp": default_ncomp,
            "color": f"C{index}",
            "marker": "o",
            "label": rf"$N_{{grid}}={grid}$, default $n_{{comp}}$",
        }
        for index, grid in enumerate(GRID_SWEEP)
    ]
    settings.extend(
        {
            "sweep": "ncomp",
            "value": ncomp,
            "grid_points": 100,
            "ncomp": ncomp,
            "color": f"C{index}",
            "marker": "s",
            "label": rf"$n_{{comp}}={ncomp}$, $N_{{grid}}=100$",
        }
        for index, ncomp in enumerate(NCOMP_SWEEP[case])
    )
    return settings


def setting_subset(frame: pd.DataFrame, case: str, setting: dict) -> pd.DataFrame:
    return frame[
        frame["case"].eq(case)
        & frame["grid_points"].eq(setting["grid_points"])
        & frame["ncomp"].eq(setting["ncomp"])
    ].sort_values("time_step")


def add_setting_legend(axis: plt.Axes, case: str) -> None:
    handles = [
        Line2D(
            [0],
            [0],
            color=setting["color"],
            marker=setting["marker"],
            linestyle="none",
            markersize=5,
            label=setting["label"],
        )
        for setting in settings_for_case(case)
    ]
    handles.extend(
        [
            Line2D([0], [0], color="0.2", linestyle="-", label="RK3"),
            Line2D([0], [0], color="0.2", linestyle="--", label="CVODE"),
        ]
    )
    axis.legend(handles=handles, fontsize=7.2, ncol=2, loc="best", frameon=False)


def make_cvode_rk3_figure(
    timing: pd.DataFrame, errors: pd.DataFrame
) -> tuple[plt.Figure, np.ndarray]:
    comparison = timing[timing["suite"].eq("cvode-rk3")].copy()
    pivot = comparison.pivot_table(
        index=["case", "grid_points", "ncomp", "time_step"],
        columns="integrator",
        values="median_wall_seconds",
        aggfunc="last",
    ).reset_index()
    pivot["cvode_over_rk3"] = pivot["cvode"] / pivot["rk3"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    cases = ("bea-alkanes-C7", "marx-403030_298_1e6_2e-1")
    for axis, case, panel in zip(axes[0], cases, ("(a)", "(b)")):
        case_data = comparison[comparison["case"].eq(case)]
        for setting in settings_for_case(case):
            values = setting_subset(case_data, case, setting)
            for integrator, linestyle in (("rk3", "-"), ("cvode", "--")):
                series = values[values["integrator"].eq(integrator)]
                axis.loglog(
                    series["time_step"],
                    series["median_wall_seconds"],
                    color=setting["color"],
                    marker=setting["marker"],
                    linestyle=linestyle,
                    linewidth=1.35,
                    markersize=4.2,
                    alpha=0.92,
                )
        axis.set_title(f"{panel} {CASE_LABELS[case]}")
        axis.set_xlabel("Time step, Δt [s]")
        axis.set_ylabel("Median wall time [s]")
        axis.grid(True, which="both", alpha=0.25)
        add_setting_legend(axis, case)

    speedup_axis = axes[1, 0]
    for case, linestyle in ((cases[0], "-"), (cases[1], "--")):
        for setting in settings_for_case(case):
            values = setting_subset(pivot, case, setting)
            speedup_axis.loglog(
                values["time_step"],
                values["cvode_over_rk3"],
                color=setting["color"],
                marker=setting["marker"],
                linestyle=linestyle,
                linewidth=1.25,
                markersize=4.2,
                alpha=0.9,
            )
    speedup_axis.axhline(1.0, color="0.2", linestyle=":", linewidth=1)
    speedup_axis.set_title("(c) Integrator wall-time ratio")
    speedup_axis.set_xlabel("Time step, Δt [s]")
    speedup_axis.set_ylabel("CVODE wall time / RK3 wall time")
    speedup_axis.grid(True, which="both", alpha=0.25)
    speedup_axis.legend(
        handles=[
            Line2D([0], [0], color="0.2", linestyle="-", label="BEA"),
            Line2D([0], [0], color="0.2", linestyle="--", label="Marx"),
            Line2D([0], [0], color="0.2", marker="o", linestyle="none", label="grid sweep"),
            Line2D([0], [0], color="0.2", marker="s", linestyle="none", label="component sweep"),
        ],
        fontsize=8,
        frameon=False,
    )

    error_axis = axes[1, 1]
    for case, filled in ((cases[0], True), (cases[1], False)):
        for setting in settings_for_case(case):
            values = setting_subset(errors, case, setting)
            error_axis.scatter(
                values["time_step"],
                values["max_abs_mole_fraction_error"],
                marker=setting["marker"],
                s=28,
                facecolors=setting["color"] if filled else "none",
                edgecolors=setting["color"],
                linewidths=1.1,
                alpha=0.92,
            )
    error_axis.set_xscale("log")
    error_axis.set_yscale("log")
    error_axis.set_title("(d) CVODE error against RK3")
    error_axis.set_xlabel("Time step, Δt [s]")
    error_axis.set_ylabel(r"max $|y_{CVODE}-y_{RK3}|$ over all mole fractions")
    error_axis.grid(True, which="both", alpha=0.25)
    error_axis.legend(
        handles=[
            Line2D([0], [0], color="0.2", marker="o", markerfacecolor="0.2", linestyle="none", label="BEA"),
            Line2D([0], [0], color="0.2", marker="o", markerfacecolor="none", linestyle="none", label="Marx"),
            Line2D([0], [0], color="0.2", marker="o", linestyle="none", label="grid sweep"),
            Line2D([0], [0], color="0.2", marker="s", linestyle="none", label="component sweep"),
        ],
        fontsize=8,
        frameon=False,
    )
    fig.suptitle("CVODE versus RK3: scaling, speedup, and mole-fraction error", fontsize=16)
    return fig, axes


def make_method_bar_figure(timing: pd.DataFrame) -> tuple[plt.Figure, np.ndarray]:
    methods = timing[timing["suite"].eq("methods")].copy()
    methods["energy_balance"] = methods["energy_balance"].astype(str).str.lower().eq("true")
    methods["chemisorption_transport"] = (
        methods["chemisorption_transport"].astype(str).str.lower().eq("true")
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True, constrained_layout=True)
    x = np.arange(len(METHODS))
    width = 0.36
    for row, transport in enumerate((False, True)):
        for column, energy in enumerate((False, True)):
            axis = axes[row, column]
            subset = methods[
                methods["energy_balance"].eq(energy)
                & methods["chemisorption_transport"].eq(transport)
            ]
            for offset, (integrator, color, hatch) in enumerate(
                (("rk3", "C0", ""), ("cvode", "C1", "//"))
            ):
                values = (
                    subset[subset["integrator"].eq(integrator)]
                    .set_index("method")
                    .reindex(METHODS)["median_wall_seconds"]
                )
                bars = axis.bar(
                    x + (offset - 0.5) * width,
                    values,
                    width,
                    color=color,
                    edgecolor="0.2",
                    linewidth=0.6,
                    hatch=hatch,
                    label=integrator.upper(),
                )
                axis.bar_label(bars, fmt="%.2g", padding=2, fontsize=7.5)
            axis.set_yscale("log")
            axis.set_xticks(x, METHODS)
            axis.set_title(
                f"Energy {'on' if energy else 'off'} · "
                f"chemisorption + transport {'on' if transport else 'off'}"
            )
            axis.grid(True, axis="y", which="both", alpha=0.25)
            if column == 0:
                axis.set_ylabel("Median wall time [s]")
            if row == 1:
                axis.set_xlabel("Mixture prediction method")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.97))
    fig.suptitle("Marx method benchmark: eight RK3/CVODE, energy, and transport settings per method", fontsize=15)
    return fig, axes
