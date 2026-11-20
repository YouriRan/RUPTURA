#!/usr/bin/env python3
"""Digitize all experimental Marx et al. breakthrough panels from the PDFs.

The paper PDFs contain vector marker paths.  This script extracts marker centres
directly from those paths, calibrates them against the printed axes, converts
temperature from degrees Celsius to kelvin, and interpolates every trace onto
exactly 50 points.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "reference_digitized"
SUPPLEMENT = ROOT / "10450_2013_9593_MOESM1_ESM.pdf"
PAPER = ROOT / "s10450-013-9593-5-1.pdf"

N_POINTS = 50


def centre(obj):
    return 0.5 * (obj["x0"] + obj["x1"]), 0.5 * (obj["top"] + obj["bottom"])


def resample(points, x_transform, y_transform):
    """Sort, merge coincident abscissae, and interpolate to exactly 50 points."""
    xy = np.asarray(sorted(points), dtype=float)
    if len(xy) < 2:
        raise ValueError(f"Need at least two marker points, found {len(xy)}")
    # Sub-point line-width offsets can map the first marker a few seconds below
    # zero; physical experiment time is nonnegative.
    x = np.maximum(np.asarray(x_transform(xy[:, 0]), dtype=float), 0.0)
    y = np.asarray(y_transform(xy[:, 1]), dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    unique_x, inverse = np.unique(x, return_inverse=True)
    if len(unique_x) != len(x):
        totals = np.bincount(inverse, weights=y)
        counts = np.bincount(inverse)
        x, y = unique_x, totals / counts
    target = np.linspace(x[0], x[-1], N_POINTS)
    return np.column_stack((target, np.interp(target, x, y)))


def write_trace(run_id, name, values, unit):
    run_dir = OUT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    header = f"time_s,{name}_{unit}"
    np.savetxt(run_dir / f"exp-{name}.csv", values, delimiter=",", header=header, comments="# ", fmt="%.10g")


def numeric_words_below_axis(page, x0, x1, y_bottom):
    pairs = []
    for word in page.extract_words():
        text = word["text"].strip()
        xc = 0.5 * (word["x0"] + word["x1"])
        if re.fullmatch(r"\d+(?:\.\d+)?", text) and x0 - 4 <= xc <= x1 + 4 and y_bottom + 1 <= word["top"] <= y_bottom + 18:
            pairs.append((xc, float(text)))
    if len(pairs) < 2:
        raise ValueError(f"Could not calibrate time axis at x={x0}, y={y_bottom}: {pairs}")
    pairs = sorted(pairs)
    return np.polyfit([p[0] for p in pairs], [p[1] for p in pairs], 1)


def temperature_calibration(page, x0, x1, y_top, y_bottom):
    horizontal = sorted(
        {
            round(obj["top"], 4)
            for obj in page.lines
            if abs(obj["bottom"] - obj["top"]) < 0.15
            and 0.85 * (x1 - x0) <= obj["x1"] - obj["x0"] <= 1.05 * (x1 - x0)
            and y_top - 1 <= obj["top"] <= y_bottom + 1
        }
    )
    pairs = []
    for y in horizontal:
        chars = []
        for char in page.chars:
            xc, yc = centre(char)
            if x0 - 17 <= xc <= x0 - 1 and abs(yc - y) <= 3.5 and char["text"] in "-0123456789.":
                chars.append(char)
        text = "".join(c["text"] for c in sorted(chars, key=lambda c: c["x0"]))
        if re.fullmatch(r"-?\d+(?:\.\d+)?", text):
            pairs.append((y, float(text)))
    # Duplicate path segments can generate duplicate y coordinates.
    pairs = sorted(dict(pairs).items())
    if len(pairs) < 2:
        raise ValueError(f"Could not calibrate temperature axis at x={x0}: {pairs}")
    return np.polyfit([p[0] for p in pairs], [p[1] for p in pairs], 1)


def in_box(obj, box):
    x, y = centre(obj)
    x0, x1, y0, y1 = box
    return x0 <= x <= x1 and y0 <= y <= y1


def digitize_binary():
    columns = [(149.20, 238.55), (262.64, 351.95), (376.04, 465.38)]
    rows = [((194.31, 277.17), (279.16, 362.02)), ((404.57, 487.43), (489.41, 572.28))]
    species_colours = {"CO2": (0.0, 0.0, 1.0), "N2": (0.0, 0.0, 0.0)}
    temperature_colour = (0.92549, 0.466667, 0.00784314)

    import pdfplumber  # only the extraction needs it; digitize_model.py imports the panel constants

    with pdfplumber.open(SUPPLEMENT) as pdf:
        for page_offset, page_index in enumerate(range(4, 10)):
            page = pdf.pages[page_index]
            for panel_index in range(6):
                row, col = divmod(panel_index, 3)
                run_id = f"B{page_offset * 6 + panel_index + 1}"
                x0, x1 = columns[col]
                (bt_top, bt_bottom), (t_top, t_bottom) = rows[row]
                x_fit = numeric_words_below_axis(page, x0, x1, t_bottom)
                t_fit = temperature_calibration(page, x0, x1, t_top, t_bottom)
                x_transform = lambda x, fit=x_fit: np.polyval(fit, x)
                fraction_transform = lambda y, top=bt_top, bottom=bt_bottom: np.clip((bottom - y) / (bottom - top), 0.0, 1.0)
                kelvin_transform = lambda y, fit=t_fit: np.polyval(fit, y) + 273.15

                bt_box = (x0, x1, bt_top, bt_bottom)
                temp_box = (x0, x1, t_top, t_bottom)
                raw_counts = {}
                for species, colour in species_colours.items():
                    objects = [
                        obj
                        for obj in page.curves
                        if obj.get("stroking_color") == colour
                        and obj.get("fill") is False
                        and 1.0 <= obj["x1"] - obj["x0"] <= 2.2
                        and obj["bottom"] - obj["top"] <= 2.2
                        and in_box(obj, bt_box)
                    ]
                    raw_counts[species] = len(objects)
                    write_trace(run_id, species, resample([centre(o) for o in objects], x_transform, fraction_transform), "mole_fraction")

                temp_objects = [
                    obj
                    for obj in page.rects
                    if obj.get("stroking_color") == temperature_colour
                    and obj.get("fill") is False
                    and 1.0 <= obj["x1"] - obj["x0"] <= 2.2
                    and 1.0 <= obj["bottom"] - obj["top"] <= 2.2
                    and in_box(obj, temp_box)
                ]
                raw_counts["Temp"] = len(temp_objects)
                write_trace(run_id, "Temp", resample([centre(o) for o in temp_objects], x_transform, kelvin_transform), "K")
                (OUT / run_id / "metadata.json").write_text(
                    json.dumps(
                        {
                            "source_pdf": SUPPLEMENT.name,
                            "pdf_page": page_index + 1,
                            "source_figure": f"S{page_offset + 2}",
                            "source_panel": panel_index + 1,
                            "points_per_trace": N_POINTS,
                            "raw_marker_counts": raw_counts,
                            "temperature_sensor": "T5 (outlet-side)",
                        },
                        indent=2,
                    )
                    + "\n"
                )


MAIN_COLUMNS_TOP = [(111.69, 223.61), (253.51, 365.40), (395.30, 507.23)]
MAIN_COLUMNS_BOTTOM = [(114.09, 226.02), (255.92, 367.81), (397.71, 509.64)]
MAIN_ROWS = [((168.88, 272.81), (274.97, 378.79)), ((433.27, 537.20), (539.36, 643.18))]

# Panel order is top-left to top-right, then bottom-left to bottom-right.
MAIN_PAGE_SETTINGS = {
    12: {
        "figure": 8,
        "runs": ["T1", "T2", "T6", "T5", "T4", "T3"],
        "x_step": [500, 500, 200, 100, 100, 200],
        "temp_bottom_tick": [24, 25, 20, 25, 25, 20],
        "temp_step": [2, 5, 10, 5, 5, 5],
    },
    13: {
        "figure": 9,
        "runs": ["T13", "T14", "T18", "T17", "T16", "T15"],
        "x_step": [500, 200, 100, 100, 100, 200],
        "temp_bottom_tick": [22, 20, 20, 30, 20, 20],
        "temp_step": [2, 10, 10, 10, 10, 10],
    },
    14: {
        "figure": 10,
        "runs": ["T7", "T8", "T12", "T11", "T10", "T9"],
        "x_step": [500, 200, 100, 100, 100, 200],
        "temp_bottom_tick": [45, 45, 50, 40, 40, 40],
        "temp_step": [5, 5, 10, 5, 5, 5],
    },
    15: {
        "figure": 11,
        "runs": ["T19", "T20", "T24", "T23", "T22", "T21"],
        "x_step": [500, 200, 100, 100, 100, 100],
        "temp_bottom_tick": [40, 40, 40, 40, 40, 40],
        "temp_step": [5, 5, 10, 10, 10, 10],
    },
}


def frequent_grid_positions(page, box, orientation):
    x0, x1, y0, y1 = box
    if orientation == "vertical":
        positions = [
            round(centre(obj)[0], 2)
            for obj in page.rects
            if obj.get("non_stroking_color") == (0.07451, 0.078431, 0.07451)
            and in_box(obj, box)
            and obj["x1"] - obj["x0"] < 0.2
            and 0.2 < obj["bottom"] - obj["top"] < 0.4
        ]
    else:
        positions = [
            round(centre(obj)[1], 2)
            for obj in page.rects
            if obj.get("non_stroking_color") == (0.07451, 0.078431, 0.07451)
            and in_box(obj, box)
            and 0.2 < obj["x1"] - obj["x0"] < 0.4
            and obj["bottom"] - obj["top"] < 0.2
        ]
    values, counts = np.unique(positions, return_counts=True)
    return np.sort(values[counts > 45])


def digitize_ternary():
    species_colours = {
        "H2": (0.87451, 0.19216, 0.13725),
        "N2": (0.07451, 0.078431, 0.07451),
        "CO2": (0.18431, 0.26667, 0.63137),
    }
    temperature_colour = (0.49804, 0.20392, 0.011765)

    import pdfplumber

    with pdfplumber.open(PAPER) as pdf:
        for page_index, settings in MAIN_PAGE_SETTINGS.items():
            page = pdf.pages[page_index]
            for panel_index, run_id in enumerate(settings["runs"]):
                row, col = divmod(panel_index, 3)
                x0, x1 = (MAIN_COLUMNS_TOP if row == 0 else MAIN_COLUMNS_BOTTOM)[col]
                (bt_top, bt_bottom), (t_top, t_bottom) = MAIN_ROWS[row]
                # Markers lying exactly on y=0 or y=1 are clipped slightly outside
                # the nominal axes path in this publisher-generated PDF.
                bt_box = (x0 - 1.25, x1 + 1.25, bt_top - 1.25, bt_bottom + 1.25)
                temp_box = (x0 - 0.25, x1 + 0.25, t_top - 0.25, t_bottom + 0.25)

                x_grid = frequent_grid_positions(page, bt_box, "vertical")
                if len(x_grid) < 2:
                    raise ValueError(f"Could not detect x grid for {run_id}")
                spacing = float(np.median(np.diff(x_grid)))
                x_origin = x_grid[0] if abs(x_grid[0] - x0) < 5 else x0
                x_scale = settings["x_step"][panel_index] / spacing
                x_transform = lambda x, origin=x_origin, scale=x_scale: (x - origin) * scale

                y_grid = frequent_grid_positions(page, temp_box, "horizontal")
                if len(y_grid) < 2:
                    raise ValueError(f"Could not detect temperature grid for {run_id}")
                y_spacing = float(np.median(np.diff(y_grid)))
                bottom_grid = float(y_grid[-1])
                bottom_value = settings["temp_bottom_tick"][panel_index]
                y_step = settings["temp_step"][panel_index]
                kelvin_transform = lambda y, yg=bottom_grid, base=bottom_value, step=y_step, spacing=y_spacing: base + (yg - y) * step / spacing + 273.15
                fraction_transform = lambda y, top=bt_top, bottom=bt_bottom: np.clip((bottom - y) / (bottom - top), 0.0, 1.0)

                raw_counts = {}
                for species, colour in species_colours.items():
                    objects = [
                        obj
                        for obj in page.curves
                        if obj.get("non_stroking_color") == colour
                        # Every plotted marker contains a short, coloured centre
                        # stroke. Unlike the outer path, it survives clipping at
                        # the axes boundaries, so it recovers baseline points too.
                        and 0.7 <= obj["x1"] - obj["x0"] <= 1.1
                        and obj["bottom"] - obj["top"] < 0.2
                        and in_box(obj, bt_box)
                    ]
                    raw_counts[species] = len(objects)
                    write_trace(run_id, species, resample([centre(o) for o in objects], x_transform, fraction_transform), "mole_fraction")

                temp_objects = [
                    obj
                    for obj in page.curves
                    if obj.get("non_stroking_color") == temperature_colour
                    and 1.8 <= obj["x1"] - obj["x0"] <= 2.7
                    and 1.8 <= obj["bottom"] - obj["top"] <= 2.7
                    and in_box(obj, temp_box)
                ]
                raw_counts["Temp"] = len(temp_objects)
                write_trace(run_id, "Temp", resample([centre(o) for o in temp_objects], x_transform, kelvin_transform), "K")
                (OUT / run_id / "metadata.json").write_text(
                    json.dumps(
                        {
                            "source_pdf": PAPER.name,
                            "pdf_page": page_index + 1,
                            "source_figure": settings["figure"],
                            "source_panel": panel_index + 1,
                            "points_per_trace": N_POINTS,
                            "raw_marker_counts": raw_counts,
                            "temperature_sensor": "T5 (outlet-side)",
                        },
                        indent=2,
                    )
                    + "\n"
                )


def validate():
    expected = [f"B{i}" for i in range(1, 37)] + [f"T{i}" for i in range(1, 25)]
    for run_id in expected:
        species = ["CO2", "N2"] if run_id.startswith("B") else ["H2", "N2", "CO2"]
        for name in [*species, "Temp"]:
            path = OUT / run_id / f"exp-{name}.csv"
            data = np.loadtxt(path, delimiter=",")
            if data.shape != (N_POINTS, 2) or not np.isfinite(data).all():
                raise ValueError(f"Invalid digitized data: {path} {data.shape}")
            if np.any(np.diff(data[:, 0]) <= 0):
                raise ValueError(f"Non-increasing time axis: {path}")


def write_readme():
    (OUT / "README.md").write_text(
        """# Digitized experimental reference data

These CSV files were extracted from the vector marker paths in the Marx et al.
paper and supporting-information PDFs. Each trace contains exactly 50 points.
The breakthrough values are mole fractions and temperatures are converted from
degrees Celsius in the source figures to kelvin. `Temp` is the outlet-side T5
sensor, matching the Ruptura outlet comparison. Every run directory includes
provenance metadata with the source PDF, page, figure, panel, and raw marker
count. Run `digitize_paper.py` to reproduce the extraction.
"""
    )


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    digitize_binary()
    digitize_ternary()
    validate()
    write_readme()
    print(f"Digitized 60 experiments to {OUT} with {N_POINTS} points per trace")
