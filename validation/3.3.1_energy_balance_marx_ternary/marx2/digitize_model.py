#!/usr/bin/env python3
"""Digitize the model curves of Marx et al. (2014), Figs. 8-11, for the ternary runs T1-T24.

The PDF draws every model line as a thin filled outline. Within each panel of digitize_paper.py, this script takes
the outlines in the colours of H2, N2, CO2, and of the outlet-side thermocouple (110 cm, whose markers
digitize_paper.py took), and reduces each to its centre line: the midpoint of the outline in bins of 0.2 pt along
the time axis. The axes are calibrated on the dashed grid lines, as for the markers.

    python digitize_model.py   # writes reference_digitized/T*/model-{H2,N2,CO2,Temp}.csv

Uses PyMuPDF. All quantities are SI (temperatures in K).
"""

from __future__ import annotations

import numpy as np
import pymupdf

from digitize_paper import MAIN_COLUMNS_BOTTOM, MAIN_COLUMNS_TOP, MAIN_PAGE_SETTINGS, MAIN_ROWS, OUT, PAPER

GRID = (0.07451, 0.078431, 0.07451)  # also the colour of N2
SPECIES = {"H2": (0.87451, 0.19216, 0.13725), "N2": GRID, "CO2": (0.18431, 0.26667, 0.63137)}
TEMPERATURE = (0.49804, 0.20392, 0.011765)  # thermocouple at 110 cm
BIN = 0.2  # pt


def same_colour(a, b) -> bool:
    return a is not None and len(a) == 3 and max(abs(x - y) for x, y in zip(a, b)) < 2.0e-3


def outline_points(drawing) -> np.ndarray:
    points = []
    for item in drawing["items"]:
        if item[0] == "re":
            r = item[1]
            points += [(r.x0, r.y0), (r.x1, r.y1)]
        elif item[0] == "qu":
            points += [(p.x, p.y) for p in (item[1].ul, item[1].ur, item[1].ll, item[1].lr)]
        else:  # "l" and "c": the points of the segment
            points += [(p.x, p.y) for p in item[1:] if isinstance(p, pymupdf.Point)]
    return np.asarray(points, dtype=float)


def centre_line(drawings, box) -> np.ndarray:
    """Midpoint of the outlines in each bin along x, ordered in x."""
    points = np.vstack([outline_points(d) for d in drawings])
    points = points[(points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)]
    bins = np.floor((points[:, 0] - box.x0) / BIN).astype(int)
    line = [(box.x0 + (b + 0.5) * BIN, 0.5 * (points[bins == b, 1].min() + points[bins == b, 1].max()))
            for b in np.unique(bins)]
    return np.asarray(line)


def grid_positions(drawings, box, vertical: bool) -> np.ndarray:
    """Centres of the dashed grid lines (paths of many dashes) in a box."""
    lines = [d["rect"] for d in drawings
             if same_colour(d.get("fill"), GRID) and len(d["items"]) > 40 and box.intersects(d["rect"])
             and ((d["rect"].width < 0.5) if vertical else (d["rect"].height < 0.5))]
    return np.sort([(r.x0 + r.x1) / 2 if vertical else (r.y0 + r.y1) / 2 for r in lines])


def model_paths(drawings, box, colour):
    """Wide filled outlines of one colour: the model lines (markers are 2-3 pt wide, grid lines are dashed).

    The panel titles are drawn in the colour of N2 just above the axes; requiring the outline to lie within the
    panel vertically leaves them out.
    """
    return [d for d in drawings
            if same_colour(d.get("fill"), colour) and box.intersects(d["rect"]) and d["rect"].y0 >= box.y0
            and d["rect"].y1 <= box.y1 and d["rect"].width > 30.0 and d["rect"].height > 0.6
            and len(d["items"]) > 30]


def main() -> None:
    document = pymupdf.open(PAPER)
    for page_index, settings in MAIN_PAGE_SETTINGS.items():
        drawings = document[page_index].get_drawings()
        for panel, run in enumerate(settings["runs"]):
            row, column = divmod(panel, 3)
            x0, x1 = (MAIN_COLUMNS_TOP if row == 0 else MAIN_COLUMNS_BOTTOM)[column]
            (bt_top, bt_bottom), (t_top, t_bottom) = MAIN_ROWS[row]
            bt_box = pymupdf.Rect(x0 - 1.25, bt_top - 1.25, x1 + 1.25, bt_bottom + 1.25)
            temp_box = pymupdf.Rect(x0 - 0.25, t_top - 0.25, x1 + 0.25, t_bottom + 0.25)

            x_grid = grid_positions(drawings, bt_box, vertical=True)
            x_origin = x_grid[0] if abs(x_grid[0] - x0) < 5 else x0
            x_scale = settings["x_step"][panel] / float(np.median(np.diff(x_grid)))
            y_grid = grid_positions(drawings, temp_box, vertical=False)
            y_scale = settings["temp_step"][panel] / float(np.median(np.diff(y_grid)))

            traces = {name: model_paths(drawings, bt_box, colour) for name, colour in SPECIES.items()}
            traces["Temp"] = model_paths(drawings, temp_box, TEMPERATURE)
            for name, paths in traces.items():
                if not paths:
                    raise ValueError(f"{run}: no model line for {name}")
                line = centre_line(paths, bt_box if name != "Temp" else temp_box)
                time = (line[:, 0] - x_origin) * x_scale
                if name == "Temp":
                    value = settings["temp_bottom_tick"][panel] + (y_grid[-1] - line[:, 1]) * y_scale + 273.15
                    header = "time_s,Temp_K"
                else:
                    value = np.clip((bt_bottom - line[:, 1]) / (bt_bottom - bt_top), 0.0, 1.0)
                    header = f"time_s,{name}_mole_fraction"
                keep = time >= 0.0
                np.savetxt(OUT / run / f"model-{name}.csv", np.column_stack((time[keep], value[keep])), delimiter=",",
                           header=header, comments="# ", fmt="%.10g")
            print(f"{run}: " + ", ".join(f"{name} {sum(len(p['items']) for p in paths)} segments"
                                          for name, paths in traces.items()))


if __name__ == "__main__":
    main()
