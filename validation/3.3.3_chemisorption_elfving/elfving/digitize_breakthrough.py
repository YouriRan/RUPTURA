#!/usr/bin/env python3
"""Recover experimental CO2, H2O, and temperature traces from the red pixels in
paper Figures 5, 6, and 8. Figure 7 has no paper H2O or temperature panels."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parent
IMAGE_DIR = ROOT / "paper_figures"
OUTPUT_DIR = ROOT / "experimental_digitized"


@dataclass(frozen=True)
class Panel:
    figure: str
    image: str
    output: str
    time_max: float
    value_column: str
    series: tuple[str, ...]
    sample_step: float
    value_max: float = 400.0
    value_min: float = 0.0
    baseline_end: float = 0.0
    zero_baseline: bool = True
    mask_boxes: tuple[tuple[int, int, int, int], ...] = ()


# Legend/inset regions are given as (x0, y0, x1, y1) pixel boxes within each
# cropped panel image, found by masking everything but red ink and checking
# what is left (see git history for the exploratory crops).
PANELS = (
    Panel("5a", "figure_5a.png", "figure_5a.csv", 120.0, "co2_mole_fraction",
          ("0p2", "1", "2"), 1.0, baseline_end=8.0, mask_boxes=((190, 120, 389, 318),)),
    Panel("6a", "figure_6a.png", "figure_6a.csv", 120.0, "co2_mole_fraction",
          ("1", "2"), 1.0, baseline_end=4.0, mask_boxes=((195, 100, 385, 316),)),
    Panel("7a", "figure_7a.png", "figure_7a.csv", 120.0, "co2_mole_fraction",
          ("0p2", "1", "2"), 1.0, baseline_end=8.0, mask_boxes=((180, 120, 376, 308),)),
    Panel("8a", "figure_8a_zoom.png", "figure_8a.csv", 450.0, "co2_mole_fraction",
          ("experimental",), 2.0, baseline_end=160.0),
    Panel("5b", "figure_5b.png", "figure_5b.csv", 120.0, "h2o_vol_percent",
          ("0p2", "1", "2"), 1.0, value_max=2.0, baseline_end=8.0,
          mask_boxes=((455, 85, 900, 420), (410, 470, 1044, 780))),
    Panel("5c", "figure_5c.png", "figure_5c.csv", 120.0, "temperature_C",
          ("0p2", "1", "2"), 1.0, value_min=24.0, value_max=36.0, zero_baseline=False,
          mask_boxes=((225, 50, 915, 265), (470, 270, 940, 660))),
    Panel("6b", "figure_6b.png", "figure_6b.csv", 120.0, "h2o_vol_percent",
          ("1", "2"), 1.0, value_max=2.0, baseline_end=4.0,
          mask_boxes=((430, 90, 1035, 395), (430, 510, 1035, 700))),
    Panel("6c", "figure_6c.png", "figure_6c.csv", 120.0, "temperature_C",
          ("1", "2"), 1.0, value_min=34.0, value_max=44.0, zero_baseline=False,
          mask_boxes=((550, 10, 1036, 180), (450, 280, 940, 630))),
    Panel("8b", "figure_8b.png", "figure_8b.csv", 2700.0, "h2o_vol_percent",
          ("experimental",), 5.0, value_max=2.0, baseline_end=2.0,
          mask_boxes=((425, 95, 1059, 295), (200, 290, 810, 755))),
    Panel("8c", "figure_8c.png", "figure_8c.csv", 2700.0, "temperature_C",
          ("experimental",), 5.0, value_min=24.0, value_max=30.0, zero_baseline=False,
          mask_boxes=((435, 5, 1059, 225), (395, 220, 965, 655))),
)


def red_mask(image: np.ndarray) -> np.ndarray:
    """Select the red experimental ink while rejecting black model curves."""
    rgb = image.astype(np.int16)
    red, green, blue = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    return (red > 155) & (red - green > 55) & (red - blue > 55)


def one_dimensional_kmeans(values: np.ndarray, clusters: int) -> np.ndarray:
    """Return ordered vertical centers for the requested number of traces."""
    if clusters == 1:
        return np.array([np.median(values)])

    centers = np.quantile(values, np.linspace(0.08, 0.92, clusters))
    for _ in range(30):
        assignment = np.argmin(np.abs(values[:, None] - centers[None, :]), axis=1)
        updated = np.array(
            [np.median(values[assignment == index]) if np.any(assignment == index) else centers[index]
             for index in range(clusters)]
        )
        if np.max(np.abs(updated - centers)) < 0.02:
            break
        centers = updated
    return np.sort(centers)


def median_filter(values: np.ndarray, radius: int = 2) -> np.ndarray:
    result = np.empty_like(values)
    for index in range(len(values)):
        result[index] = np.median(values[max(0, index - radius): index + radius + 1])
    return result


def digitize(panel: Panel) -> tuple[np.ndarray, np.ndarray]:
    image = np.asarray(Image.open(IMAGE_DIR / panel.image).convert("RGB"))
    mask = red_mask(image)
    for x0, y0, x1, y1 in panel.mask_boxes:
        mask[y0:y1, x0:x1] = False

    height, width = mask.shape
    times = np.arange(0.0, panel.time_max + 0.5 * panel.sample_step, panel.sample_step)
    rows = np.full((len(times), len(panel.series)), np.nan)

    for time_index, time in enumerate(times):
        x_center = time / panel.time_max * (width - 1)
        for radius in (2, 4, 7, 11):
            low = max(0, int(round(x_center)) - radius)
            high = min(width, int(round(x_center)) + radius + 1)
            y_values = np.flatnonzero(mask[:, low:high]) // (high - low)
            if len(y_values) >= max(4, 2 * len(panel.series)):
                rows[time_index] = one_dimensional_kmeans(y_values.astype(float), len(panel.series))
                break

    for series_index in range(rows.shape[1]):
        valid = np.isfinite(rows[:, series_index])
        if valid.sum() < 2:
            raise RuntimeError(f"Not enough red pixels recovered for Figure {panel.figure}")
        rows[:, series_index] = np.interp(times, times[valid], rows[valid, series_index])
        rows[:, series_index] = median_filter(rows[:, series_index])

    # Smaller image y means a higher plotted value. The panel crops are bounded
    # by the published axes, so this affine transform is the calibration.
    span = panel.value_max - panel.value_min
    values = panel.value_min + span * (1.0 - rows / (height - 1))
    if panel.zero_baseline:
        # Remove the sub-pixel vertical offset of the rendered zero line using
        # the pre-breakthrough baseline. This preserves low-value information
        # needed for a 1% breakthrough estimate without treating line
        # thickness as a real outlet reading. Not meaningful for temperature,
        # which has no zero baseline to snap to.
        baseline = np.median(values[times <= panel.baseline_end], axis=0)
        values = np.clip(values - baseline, 0.0, span)
        values[values < 0.0025 * span] = 0.0
        values[times <= panel.baseline_end] = 0.0
    return times, values


def write_csv(panel: Panel, times: np.ndarray, values: np.ndarray) -> None:
    output = OUTPUT_DIR / panel.output
    scale = 1.0e-6 if panel.value_column == "co2_mole_fraction" else 1.0
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", *(f"{panel.value_column}_{name}" for name in panel.series)])
        for time, row in zip(times, values):
            writer.writerow([f"{60.0 * time:.3f}", *(f"{scale * value:.9f}" for value in row)])


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    for panel in PANELS:
        times, concentrations = digitize(panel)
        write_csv(panel, times, concentrations)
        print(f"Wrote {OUTPUT_DIR / panel.output} ({len(times)} calibrated samples)")


if __name__ == "__main__":
    main()
