#!/usr/bin/env python3
"""Vector digitisation of Figures 4-8 of Stadler, Bender and Pfeifer (2022).

Figures 4-8 of the article PDF (fceng-04-1000064.pdf) are vector graphics, so
nothing is read from pixels.  Each marker is its own filled path, and each curve
is a single stroked polyline (a dashed curve carries a dash pattern but still
lists every vertex).  A marker is located at the centre of its bounding box and a
curve by its vertices.  Series are told apart by colour and, within a colour, by
whether the path is a marker or a solid or dashed line.

Every panel frame coincides with its axis limits (checked against the tick marks),
so the frame rectangle and the axis limits are the whole calibration.  Markers in
a panel's legend box are discarded.

Writes one CSV per series into ``reference_digitized/``, in SI units (pressures in
Pa, loadings in mol/kg, times in s, fractions [-]):

    figure4_<experiment|model|equilibrium>_<T>K.csv  steam_to_gas, co_conversion
    figure5_<experiment|freundlich_fit>.csv           pressure_Pa, loading_mol_per_kg
    figure6_cycle<N>_<CO|CO2>.csv                     time_s, dry_fraction
    figure7_cycle<N>_<CO|CO2>.csv                     time_s, dry_fraction
    figure8<a|b>_<experiment|simulation>_<CO|CO2>.csv time_s, dry_fraction

Requires PyMuPDF (``pip install pymupdf``).
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pymupdf

ROOT = Path(__file__).resolve().parent
PAPER = ROOT / "fceng-04-1000064.pdf"
OUTPUT = ROOT / "reference_digitized"

MMHG = 133.322387415          # Pa
MMOL_PER_G = 1.0              # mol/kg

# ------------------------------------------------------------ colours in the PDF
BLUE = (0.189, 0.392, 0.677)
GREEN = (0.0, 0.590, 0.506)
ORANGE = (0.954, 0.610, 0.0)
LIGHT_BLUE = (0.534, 0.637, 0.807)
LIGHT_ORANGE = (0.977, 0.766, 0.446)

# ------------------------------------------------------------------- panels
# page is zero-based; frame = (x0, y0, x1, y1) in PDF points; x and y are the axis
# limits at the frame edges; legend = (x0, y0, x1, y1) or None.
FIGURE_4 = dict(page=8, frame=(91.52, 474.54, 279.36, 595.33), x=(0.0, 4.5), y=(0.0, 1.0),
                legend=(135.8, 573.2, 272.5, 588.4))
FIGURE_5 = dict(page=8, frame=(347.13, 470.51, 529.95, 586.43), x=(0.0, 800.0), y=(0.0, 0.3),
                legend=(372.3, 476.6, 504.8, 491.9))

# Figures 6 and 7: five panels each, 0-200 s by 0-0.08 dry volume fraction.
FIGURE_6_CYCLES = (1, 2, 3, 4, 5)
FIGURE_6_FRAMES = [(154.9, 533.7, 197.2, 661.3), (214.3, 533.7, 256.7, 661.3),
                   (273.8, 533.7, 315.6, 661.3), (333.2, 533.7, 375.0, 661.3),
                   (392.1, 533.7, 434.5, 661.3)]
FIGURE_7_CYCLES = (1, 3, 13, 43, 63)
FIGURE_7_FRAMES = [(153.8, 99.2, 196.3, 227.2), (213.5, 99.2, 255.9, 227.2),
                   (273.1, 99.2, 315.0, 227.2), (332.8, 99.2, 374.7, 227.2),
                   (391.9, 99.2, 434.3, 227.2)]
FIGURE_8 = {
    "a": dict(page=10, frame=(166.5, 371.7, 283.7, 449.4), x=(0.0, 200.0), y=(0.0, 0.08),
              legend=(210.6, 421.4, 279.3, 444.9)),
    "b": dict(page=10, frame=(331.7, 371.7, 448.9, 449.4), x=(0.0, 200.0), y=(0.0, 0.08),
              legend=(375.8, 410.1, 444.5, 434.1)),
}

# Figure 4 legend: blue circles 225 C, green squares 250 C, orange triangles 275 C.
FIGURE_4_SERIES = {498: BLUE, 523: GREEN, 548: ORANGE}


def same_colour(rgb, reference, tolerance=0.01):
    return rgb is not None and max(abs(a - b) for a, b in zip(rgb, reference)) <= tolerance


def inside(x, y, box, margin=0.0):
    x0, y0, x1, y1 = box
    return x0 - margin <= x <= x1 + margin and y0 - margin <= y <= y1 + margin


def to_data(points, panel):
    """PDF points (x right, y down) to axis values."""
    x0, y0, x1, y1 = panel["frame"]
    (xmin, xmax), (ymin, ymax) = panel["x"], panel["y"]
    points = np.asarray(points, dtype=float).reshape(-1, 2)
    x = xmin + (points[:, 0] - x0) / (x1 - x0) * (xmax - xmin)
    y = ymin + (y1 - points[:, 1]) / (y1 - y0) * (ymax - ymin)
    return x, y


def markers(drawings, colour, panel, stroked=False):
    """Bounding-box centres of the filled (or stroked) markers drawn in `colour`.

    Figure 8 draws its open triangles as one stroked path of three segments per
    triangle; those are split back into triangles here.
    """
    centres = []
    for drawing in drawings:
        if stroked:
            if drawing["type"] != "s" or not same_colour(drawing["color"], colour):
                continue
            segments = [item for item in drawing["items"] if item[0] == "l"]
            boxes = []
            for i in range(0, len(segments) - 2, 3):
                xs = [point.x for item in segments[i:i + 3] for point in item[1:]]
                ys = [point.y for item in segments[i:i + 3] for point in item[1:]]
                boxes.append((min(xs), min(ys), max(xs), max(ys)))
        else:
            if drawing["type"] != "f" or not same_colour(drawing["fill"], colour):
                continue
            boxes = [tuple(drawing["rect"])]
        for x0, y0, x1, y1 in boxes:
            centre = (0.5 * (x0 + x1), 0.5 * (y0 + y1))
            if not inside(*centre, panel["frame"], margin=2.0):
                continue
            if panel.get("legend") and inside(*centre, panel["legend"]):
                continue
            centres.append(centre)
    centres.sort()
    return to_data(centres, panel)


def polyline(drawings, colour, panel, dashed=False):
    """Vertices of the stroked line in `colour` that lies in the panel frame."""
    for drawing in drawings:
        if drawing["type"] != "s" or not same_colour(drawing["color"], colour):
            continue
        if (drawing["dashes"] not in (None, "[] 0")) != dashed:
            continue
        segments = [item for item in drawing["items"] if item[0] == "l"]
        if len(segments) < 10:  # legend sample lines
            continue
        rect = drawing["rect"]
        if not inside(0.5 * (rect.x0 + rect.x1), 0.5 * (rect.y0 + rect.y1), panel["frame"]):
            continue
        # A path can hold several subpaths (Figure 7, cycle 3 draws the start of its curve last); split it
        # where a segment does not start at the previous end, and order the pieces along the x axis.
        pieces = []
        for item in segments:
            start, end = (item[1].x, item[1].y), (item[2].x, item[2].y)
            if not pieces or abs(pieces[-1][-1][0] - start[0]) + abs(pieces[-1][-1][1] - start[1]) > 1.0e-3:
                pieces.append([start])
            pieces[-1].append(end)
        pieces.sort(key=lambda piece: piece[0][0])
        return to_data([vertex for piece in pieces for vertex in piece], panel)
    raise LookupError(f"no {'dashed' if dashed else 'solid'} line of colour {colour} in {panel['frame']}")


def write(name, header, x, y):
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT / f"{name}.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for row in zip(x, y):
            writer.writerow([f"{value:.6g}" for value in row])
    print(f"  {name:40s} {len(x):4d} points")


def main():
    document = pymupdf.open(PAPER)

    # Figure 4: steady-state CO conversion over the steam-to-gas ratio.
    drawings = document[FIGURE_4["page"]].get_drawings()
    for temperature, colour in FIGURE_4_SERIES.items():
        header = ("steam_to_gas", "co_conversion")
        write(f"figure4_experiment_{temperature}K", header, *markers(drawings, colour, FIGURE_4))
        write(f"figure4_model_{temperature}K", header, *polyline(drawings, colour, FIGURE_4))
        write(f"figure4_equilibrium_{temperature}K", header,
              *polyline(drawings, colour, FIGURE_4, dashed=True))

    # Figure 5: CO2 isotherm at 523 K; mmHg and mmol/g converted to Pa and mol/kg.
    drawings = document[FIGURE_5["page"]].get_drawings()
    header = ("pressure_Pa", "loading_mol_per_kg")
    pressure, loading = markers(drawings, BLUE, FIGURE_5)
    write("figure5_experiment", header, pressure * MMHG, loading * MMOL_PER_G)
    pressure, loading = polyline(drawings, BLUE, FIGURE_5)
    write("figure5_freundlich_fit", header, pressure * MMHG, loading * MMOL_PER_G)

    # Figures 6 and 7: measured outlet curves, one panel per cycle.
    header = ("time_s", "dry_fraction")
    for figure, page, cycles, frames in ((6, 9, FIGURE_6_CYCLES, FIGURE_6_FRAMES),
                                         (7, 10, FIGURE_7_CYCLES, FIGURE_7_FRAMES)):
        drawings = document[page].get_drawings()
        for cycle, frame in zip(cycles, frames):
            panel = dict(frame=frame, x=(0.0, 200.0), y=(0.0, 0.08))
            write(f"figure{figure}_cycle{cycle}_CO2", header, *polyline(drawings, BLUE, panel))
            write(f"figure{figure}_cycle{cycle}_CO", header, *polyline(drawings, ORANGE, panel))

    # Figure 8: third cycle, measured (open markers) against the authors' model (lines).
    for key, panel in FIGURE_8.items():
        drawings = document[panel["page"]].get_drawings()
        write(f"figure8{key}_experiment_CO2", header, *markers(drawings, LIGHT_BLUE, panel))
        write(f"figure8{key}_experiment_CO", header,
              *markers(drawings, LIGHT_ORANGE, panel, stroked=True))
        write(f"figure8{key}_simulation_CO2", header, *polyline(drawings, BLUE, panel))
        write(f"figure8{key}_simulation_CO", header, *polyline(drawings, ORANGE, panel))


if __name__ == "__main__":
    main()
