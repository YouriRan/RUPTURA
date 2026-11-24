#!/usr/bin/env python3
"""Compare the Danckwerts runs of ``generate.py`` with the Dirichlet runs, and the Dirichlet runs with the stored
results of the source studies.

For every component that enters the initially empty column, the outlet mole fraction normalized by the inlet one,
y_out/y_in, gives the breakthrough times t10, t50, and t90 (first crossings, linear interpolation). For every case,
the table also lists the largest difference between the outlet mole fractions of two runs, relative to the largest
inlet mole fraction of that component, over the common output times and all components; it also covers the cyclic
runs, whose first-cycle breakthrough times say little.

    python compare.py             # write comparison.csv and print a summary per study
    python compare.py -j 4        # four worker processes (default 4)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parent.parent
LEVELS = (0.1, 0.5, 0.9)


def read_outlets(directory: Path) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Component name -> (time in s, inlet mole fraction, outlet mole fraction) at every output."""
    simulation = {key.lower(): value for key, value in json.loads((directory / "simulation.json").read_text()).items()}
    feed = {}
    for component in simulation["components"]:
        component = {key.lower(): value for key, value in component.items()}
        feed[component["name"]] = component.get("gasphasemolfraction", 0.0)
    total = sum(feed.values())  # Ruptura normalizes the feed (the Marx inputs give percentages)
    feed = {name: value / total for name, value in feed.items()}
    outlets = {}
    for path in sorted(directory.glob("component_*_*.data")):
        name = path.stem.split("_", 2)[2]
        frame = pd.read_csv(path, comment="#", sep=r"\s+", header=None, usecols=[1, 2, 5], names=["t", "z", "y"])
        z = frame["z"].to_numpy()
        restarts = np.flatnonzero(z[1:] < z[:-1])
        # A WriteOutletOnly file has one row per output, all at z = L.
        # shortcut: a reverse-flow phase would write node 1 instead; detect by the key if such a source is added.
        nodes = int(restarts[0]) + 1 if restarts.size else (1 if np.ptp(z) == 0.0 else z.size)
        outputs = z.size // nodes
        t = frame["t"].to_numpy()[: outputs * nodes].reshape(outputs, nodes)[:, 0] * 60.0
        y = frame["y"].to_numpy()[: outputs * nodes].reshape(outputs, nodes)
        # Without the z = 0 row the inlet is the feed, which node 0 holds in a full file.
        y_in = np.full(outputs, feed[name]) if nodes == 1 else y[:, 0]
        outlets[name] = (t, y_in, y[:, -1])
    return outlets


def crossing(t: np.ndarray, ratio: np.ndarray, level: float) -> float:
    above = np.flatnonzero(ratio >= level)
    if not above.size:
        return math.nan
    i = int(above[0])
    if i == 0:
        return float(t[0])
    return float(t[i - 1] + (level - ratio[i - 1]) * (t[i] - t[i - 1]) / (ratio[i] - ratio[i - 1]))


def breakthrough_times(outlets) -> dict[str, tuple[float, float, float]]:
    """t10, t50, t90 of every component that is fed but not initially present at the outlet."""
    times = {}
    for name, (t, y_in, y_out) in outlets.items():
        if y_in.max() <= 0.0:
            continue
        ratio = np.divide(y_out, y_in, out=np.zeros_like(y_out), where=y_in > 0.0)
        if ratio[0] >= 0.5:
            continue
        times[name] = tuple(crossing(t, ratio, level) for level in LEVELS)
    return times


def largest_outlet_difference(first, second) -> float:
    """Largest |y_out difference| / max(y_in) over common output times and components; nan without overlap."""
    largest = math.nan
    for name in first.keys() & second.keys():
        t_a, y_in, y_a = first[name]
        t_b, _, y_b = second[name]
        common, i_a, i_b = np.intersect1d(np.round(t_a, 6), np.round(t_b, 6), return_indices=True)
        if not common.size or y_in.max() <= 0.0:
            continue
        difference = float(np.max(np.abs(y_a[i_a] - y_b[i_b])) / y_in.max())
        largest = difference if math.isnan(largest) else max(largest, difference)
    return largest


def finished(directory: Path) -> str:
    for log in ("run.log", "ruptura.log"):
        if (directory / log).is_file():
            return "ok" if "Final time:" in (directory / log).read_text(errors="ignore") else "stopped"
    return "missing"


def compare_case(entry: dict) -> list[dict]:
    dirichlet_dir = ROOT / entry["dirichlet"]
    danckwerts_dir = ROOT / entry["danckwerts"]
    stored_dir = REPOSITORY / entry["source"]
    outlets = {
        "stored": read_outlets(stored_dir),
        "dirichlet": read_outlets(dirichlet_dir),
        "danckwerts": read_outlets(danckwerts_dir),
    }
    times = {run: breakthrough_times(data) for run, data in outlets.items()}
    case = {
        "study": entry["study"],
        "case": entry["case"],
        "status_dirichlet": finished(dirichlet_dir),
        "status_danckwerts": finished(danckwerts_dir),
        "max_dy_danckwerts_vs_dirichlet": largest_outlet_difference(outlets["danckwerts"], outlets["dirichlet"]),
        "max_dy_dirichlet_vs_stored": largest_outlet_difference(outlets["dirichlet"], outlets["stored"]),
    }
    rows = []
    for name in sorted(times["dirichlet"].keys() | times["danckwerts"].keys()):
        row = dict(case, component=name)
        for run in ("stored", "dirichlet", "danckwerts"):
            for level, value in zip(LEVELS, times[run].get(name, (math.nan,) * 3)):
                row[f"t{round(100 * level)}_{run}"] = value
        rows.append(row)
    return rows or [dict(case, component="")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-j", "--workers", type=int, default=4)
    arguments = parser.parse_args()

    pairs = {}
    for entry in json.loads((ROOT / "simulations" / "index.json").read_text()):
        pair = pairs.setdefault((entry["study"], entry["case"]), {"study": entry["study"], "case": entry["case"]})
        pair[entry["inlet"].lower()] = entry["directory"]
        pair["source"] = entry["source"]

    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        rows = [row for case_rows in executor.map(compare_case, pairs.values()) for row in case_rows]

    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (ROOT / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    frame = pd.DataFrame(rows)
    frame["dt50_rel"] = (frame["t50_danckwerts"] - frame["t50_dirichlet"]) / frame["t50_dirichlet"]
    spread_dirichlet = frame["t90_dirichlet"] - frame["t10_dirichlet"]
    frame["dspread_rel"] = (frame["t90_danckwerts"] - frame["t10_danckwerts"] - spread_dirichlet) / spread_dirichlet
    frame["dt50_rel_stored"] = (frame["t50_dirichlet"] - frame["t50_stored"]) / frame["t50_stored"]
    pd.set_option("display.width", 200)
    for study, group in frame.groupby("study", sort=False):
        cases = group.drop_duplicates("case")
        print(f"\n{study}: {len(cases)} cases; Dirichlet {dict(cases['status_dirichlet'].value_counts())}, "
              f"Danckwerts {dict(cases['status_danckwerts'].value_counts())}")
        print(f"  Danckwerts vs Dirichlet: max |dt50|/t50 {group['dt50_rel'].abs().max():.2e}, "
              f"max |d(t90-t10)|/(t90-t10) {group['dspread_rel'].abs().max():.2e}, "
              f"max |dy_out|/y_in {cases['max_dy_danckwerts_vs_dirichlet'].max():.2e}")
        print(f"  Dirichlet vs stored:     max |dt50|/t50 {group['dt50_rel_stored'].abs().max():.2e}, "
              f"max |dy_out|/y_in {cases['max_dy_dirichlet_vs_stored'].max():.2e}")


if __name__ == "__main__":
    main()
