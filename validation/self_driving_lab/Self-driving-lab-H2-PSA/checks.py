#!/usr/bin/env python3
"""Fidelity and ablation checks of the two campaigns, quoted in the paper (Results, self-driving laboratory).

  * one column: the two ends of the Pareto front and the front setup with the largest backfill, on 160 cells with CVODE
    at a relative tolerance of 1e-6; and that setup without its backfill, at the fidelity of the campaign;
  * two beds: every setup of the Pareto front on 160 cells with ROS34PW2 at a relative tolerance of 1e-6; and the most
    productive front setup without its pressure equalization, at the fidelity of the campaign.

Each check is an experiment directory below checks/, written by the robot and run by `python -m ruptura.selfdriving`,
8 at a time, as in the campaigns.  Run after the campaigns:

    python checks.py          # write and run checks/, then print the comparison with the campaign
    python checks.py --report # print the comparison only
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

import pandas as pd

from ruptura.selfdriving import PressureSwingCycle, TwoBedPressureSwingCycle, pareto_front

ROOT = Path(__file__).resolve().parent
EXAMPLES = ROOT.parents[2] / "examples"
CHECKS = ROOT / "checks"


def front(campaign: Path, robot) -> pd.DataFrame:
    """The Pareto front of the experiments at the purity target, in productivity and recovery (as analysis.ipynb)."""
    data = pd.read_csv(campaign / "experiments.csv")
    feasible = data[(data["status"] != "failed") & (data["purity"] >= robot.purity_target)]
    return feasible[pareto_front(feasible[["productivity", "recovery"]].to_numpy())].sort_values("productivity")


def setup_of(campaign: Path, name: str) -> dict:
    return json.loads(next(campaign.glob(f"round_*/{name}/setup.json")).read_text())


def plan() -> list[tuple[str, object, dict, str]]:
    """(check directory, robot, setup, name of the campaign experiment it compares with)."""
    checks = []
    one = PressureSwingCycle(EXAMPLES / "PSA-H2-activated-carbon" / "simulation.json", product="H2", purity_target=0.999)
    fine_one = replace(one, grid_points=160, integrator="CVODE", relative_tolerance=1.0e-6)
    pareto = front(ROOT / "campaign", one)
    largest_backfill = pareto.sort_values("backfill").iloc[-1]
    for row in {pareto.iloc[0]["name"], pareto.iloc[-1]["name"], largest_backfill["name"]}:
        checks.append((f"one_column_160_cells_cvode_1e-6/{row}", fine_one, setup_of(ROOT / "campaign", row), row))
    setup = setup_of(ROOT / "campaign", largest_backfill["name"])
    checks.append((f"one_column_without_backfill/{largest_backfill['name']}", one, {**setup, "backfill": 0.0},
                   largest_backfill["name"]))

    two = TwoBedPressureSwingCycle(EXAMPLES / "PSA-H2-activated-carbon-two-beds" / "simulation.json", product="H2",
                                   purity_target=0.999)
    fine_two = replace(two, grid_points=160, relative_tolerance=1.0e-6)
    pareto = front(ROOT / "campaign_two_beds", two)
    for row in pareto["name"]:
        checks.append((f"two_beds_160_cells_ros34pw2_1e-6/{row}", fine_two, setup_of(ROOT / "campaign_two_beds", row), row))
    best = pareto.iloc[-1]["name"]
    setup = setup_of(ROOT / "campaign_two_beds", best)
    checks.append((f"two_beds_without_equalization/{best}", two, {**setup, "equalization": 0.0}, best))
    return checks


def report(checks) -> None:
    campaigns = {"one_column": ROOT / "campaign", "two_beds": ROOT / "campaign_two_beds"}
    rows = []
    for directory, _, _, name in checks:
        campaign = campaigns["two_beds" if directory.startswith("two_beds") else "one_column"]
        reference = json.loads(next(campaign.glob(f"round_*/{name}/kpis.json")).read_text())
        check = json.loads((CHECKS / directory / "kpis.json").read_text())
        rows.append({"check": directory, "status": check.get("status"),
                     "N2 / campaign": check.get("product_fraction_N2", float("nan")) / reference["product_fraction_N2"],
                     "N2": check.get("product_fraction_N2"), "N2 campaign": reference["product_fraction_N2"],
                     "recovery": check.get("recovery"), "recovery campaign": reference["recovery"],
                     "productivity / campaign": check.get("productivity", float("nan")) / reference["productivity"]})
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.max_colwidth", 80):
        print(pd.DataFrame(rows).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--report", action="store_true", help="print the comparison of existing checks only")
    arguments = parser.parse_args()
    checks = plan()
    if not arguments.report:
        for directory, robot, setup, _ in checks:
            robot.write(CHECKS / directory, setup)
        run = lambda directory: subprocess.run([sys.executable, "-m", "ruptura.selfdriving", str(CHECKS / directory)],
                                               stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        with ThreadPoolExecutor(max_workers=8) as pool:
            list(pool.map(run, [directory for directory, *_ in checks]))
    report(checks)


if __name__ == "__main__":
    main()
