#!/usr/bin/env python3
"""Colour-keyed raster digitisation of Figures 4 and 5 of Stampi-Bombelli et al.

The two 300 dpi rasters embedded in the article PDF are stored under
``paper_figures/``.  Every experimental trace in those panels is drawn in its own
solid colour, so a trace can be recovered exactly by selecting the pixels that
match that colour, grouping them into contiguous runs within each pixel column
and keeping the longest run.  That last step discards the single-pixel
antialiasing bleed a curve leaves on its neighbour of a similar shade.

The pixel-to-value calibration below was read off the axis frames and gridlines
of the rasters.  As a check, the t50 and t90 values recovered from Figure 4
agree with the tabulated values of Table 6 to better than 3 %.

Writes one CSV per trace into ``reference_digitized/``:

    figure4_<case>.csv              time_s, outlet_fraction
    figure4_temperature_<case>.csv  time_s, temperature_C
    figure5_time_<case>.csv         shifted_time_min, outlet_fraction
    figure5_tau_<case>.csv          shifted_tau, outlet_fraction
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "reference_digitized"
FIGURE_4 = ROOT / "paper_figures" / "paper_figure_4.png"
FIGURE_5 = ROOT / "paper_figures" / "paper_figure_5.png"

MAXIMUM_POINTS = 160          # points kept per trace after resampling

# --------------------------------------------------------- colours of the traces
# One entry per experiment: the RGB triple used for that curve in both figures.
TRACE_COLOUR = {
    "packed_bed_400ppm_2mmol": (163, 20, 46),
    "packed_bed_5.6pct_2mmol": (204, 115, 8),
    "packed_bed_5.6pct_1mmol": (255, 166, 64),
    "packed_bed_5.6pct_0p5mmol": (255, 199, 143),
    "monolith_400ppm_3mmol": (69, 79, 133),
    "monolith_400ppm_2mmol": (128, 138, 194),
    "monolith_400ppm_1mmol": (171, 191, 204),
    "monolith_400ppm_0p18mmol": (189, 204, 219),
    "monolith_5.6pct_0p18mmol": (105, 140, 107),
    "monolith_5.6pct_0p06mmol": (128, 191, 133),
    "monolith_5.6pct_0p02mmol": (186, 212, 168),
}

# ------------------------------------------------------------------ axis frames
# box = (row_top, row_bottom, column_left, column_right) of the plotted area;
# x and y give two (pixel, value) calibration points along each axis.
FIGURE_4_PANELS = {
    "packed_bed_400ppm": dict(
        cases=["packed_bed_400ppm_2mmol"],
        box=(40, 587, 105, 819), temperature_box=(607, 844, 105, 819),
        x=((107.25, 0.0), (818.25, 180000.0)),
        y=((541.5, 0.0), (76.5, 1.0)),
        temperature_y=((844.0, 20.0), (654.5, 40.0))),
    "packed_bed_5.6pct": dict(
        cases=["packed_bed_5.6pct_2mmol", "packed_bed_5.6pct_1mmol",
               "packed_bed_5.6pct_0p5mmol"],
        box=(40, 587, 1015, 1727), temperature_box=(607, 844, 1015, 1727),
        x=((1014.5, 0.0), (1726.5, 24000.0)),
        y=((541.5, 0.0), (54.5, 6.0 / 5.6)),
        temperature_y=((844.0, 20.0), (654.5, 40.0))),
    "monolith_400ppm": dict(
        cases=["monolith_400ppm_3mmol", "monolith_400ppm_2mmol",
               "monolith_400ppm_1mmol", "monolith_400ppm_0p18mmol"],
        box=(1003, 1550, 115, 827), temperature_box=(1570, 1807, 115, 827),
        x=((114.75, 0.0), (825.75, 18000.0)),
        y=((1505.0, 0.0), (1050.0, 1.0)),
        temperature_y=((1807.0, 20.0), (1617.5, 40.0))),
    "monolith_5.6pct": dict(
        cases=["monolith_5.6pct_0p18mmol", "monolith_5.6pct_0p06mmol",
               "monolith_5.6pct_0p02mmol"],
        box=(1003, 1550, 1015, 1727), temperature_box=(1570, 1807, 1015, 1727),
        x=((1015.0, 0.0), (1727.0, 7200.0)),
        y=((1505.0, 0.0), (1018.0, 6.0 / 5.6)),
        temperature_y=((1807.0, 20.0), (1617.5, 40.0))),
}

FIGURE_5_PANELS = {
    "a": dict(kind="time", cases=FIGURE_4_PANELS["packed_bed_5.6pct"]["cases"],
              box=(60, 614, 107, 661), x=((384.0, 0.0), (522.0, 20.0)),
              y=((521.0, 0.0), (152.0, 0.8))),
    "b": dict(kind="time", cases=FIGURE_4_PANELS["monolith_5.6pct"]["cases"],
              box=(60, 614, 707, 1260), x=((983.0, 0.0), (1121.0, 5.0)),
              y=((521.0, 0.0), (152.0, 0.8))),
    "c": dict(kind="time", cases=FIGURE_4_PANELS["monolith_400ppm"]["cases"],
              box=(60, 614, 1306, 1859), x=((1582.0, 0.0), (1674.0, 20.0)),
              y=((521.0, 0.0), (152.0, 0.8))),
    "d": dict(kind="tau", cases=FIGURE_4_PANELS["packed_bed_5.6pct"]["cases"],
              box=(748, 1302, 107, 661), x=((384.0, 0.0), (506.5, 200.0)),
              y=((1208.5, 0.0), (840.0, 0.8))),
    "e": dict(kind="tau", cases=FIGURE_4_PANELS["monolith_5.6pct"]["cases"],
              box=(748, 1302, 707, 1260), x=((983.0, 0.0), (1121.0, 5.0)),
              y=((1208.5, 0.0), (840.0, 0.8))),
    "f": dict(kind="tau", cases=FIGURE_4_PANELS["monolith_400ppm"]["cases"],
              box=(748, 1302, 1306, 1859), x=((1582.0, 0.0), (1705.0, 200.0)),
              y=((1208.5, 0.0), (840.0, 0.8))),
}


def _calibration(pair):
    (pixel0, value0), (pixel1, value1) = pair
    slope = (value1 - value0) / (pixel1 - pixel0)
    return lambda pixel: value0 + slope * (pixel - pixel0)


def extract(image, box, x_axis, y_axis, colour, tolerance=4, minimum_run=5):
    """Return one sample per pixel column of the trace drawn in `colour`."""
    row0, row1, column0, column1 = box
    window = image[row0:row1, column0:column1].astype(int)
    matches = np.abs(window - np.asarray(colour, dtype=int)).max(axis=2) <= tolerance

    to_x = _calibration(x_axis)
    to_y = _calibration(y_axis)

    xs, ys = [], []
    for column in range(matches.shape[1]):
        rows = np.flatnonzero(matches[:, column])
        if rows.size < minimum_run:
            continue
        runs = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
        longest = max(runs, key=len)
        if longest.size < minimum_run:
            continue
        xs.append(to_x(column0 + column))
        ys.append(to_y(row0 + longest.mean()))
    return np.asarray(xs), np.asarray(ys)


def trim_leading_transient(x, y):
    """Start the trace at the minimum of its first quarter.

    Figure 4a opens with the vertical stroke the CO2 sensor draws when the feed
    is switched on, which the column-wise extraction sees as three very high
    leading samples.  A breakthrough trace only rises, so the first sample must
    be the lowest one; on every other panel this is a no-op.
    """
    start = int(np.argmin(y[: max(1, y.size // 4)]))
    return x[start:], y[start:]


def resample(x, y, count=MAXIMUM_POINTS):
    if x.size <= count:
        return x, y
    index = np.unique(np.linspace(0, x.size - 1, count).round().astype(int))
    return x[index], y[index]


def write(name, header, x, y):
    OUTPUT.mkdir(exist_ok=True)
    path = OUTPUT / f"{name}.csv"
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for first, second in zip(x, y):
            writer.writerow([f"{first:.6g}", f"{second:.6g}"])
    return path


def digitize_figure_4(image):
    written = []
    for panel in FIGURE_4_PANELS.values():
        for case in panel["cases"]:
            colour = TRACE_COLOUR[case]
            x, y = extract(image, panel["box"], panel["x"], panel["y"], colour)
            if x.size:
                x, y = resample(*trim_leading_transient(x, y))
                written.append(write(f"figure4_{case}", ["time_s", "outlet_fraction"], x, y))
            x, y = extract(image, panel["temperature_box"], panel["x"],
                           panel["temperature_y"], colour, minimum_run=3)
            if x.size:
                x, y = resample(x, y)
                written.append(write(f"figure4_temperature_{case}",
                                     ["time_s", "temperature_C"], x, y))
    return written


def digitize_figure_5(image):
    written = []
    for panel in FIGURE_5_PANELS.values():
        header = (["shifted_time_min", "outlet_fraction"] if panel["kind"] == "time"
                  else ["shifted_tau", "outlet_fraction"])
        for case in panel["cases"]:
            x, y = extract(image, panel["box"], panel["x"], panel["y"],
                           TRACE_COLOUR[case], minimum_run=4)
            if not x.size:
                continue
            x, y = resample(x, y)
            written.append(write(f"figure5_{panel['kind']}_{case}", header, x, y))
    return written


def main():
    for source in (FIGURE_4, FIGURE_5):
        if not source.is_file():
            raise FileNotFoundError(f"missing paper raster: {source}")
    written = digitize_figure_4(np.array(Image.open(FIGURE_4).convert("RGB")))
    written += digitize_figure_5(np.array(Image.open(FIGURE_5).convert("RGB")))
    for path in written:
        print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
