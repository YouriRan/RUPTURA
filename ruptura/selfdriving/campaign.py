"""The closed loop of the self-driving laboratory: propose a batch, run it, read the results, and propose the next batch.

Round 0 is the baseline of the cycle and space-filling (Sobol) setups; every later round is a batch from
optimizer.propose_batch. Each experiment gets its own directory, <directory>/round_<r>/<name>/, and runs in its own
Python process: a simulation keeps the interpreter lock while it runs, and a failed integration stops only its own
experiment. <directory>/experiments.csv is the lab notebook, one row per experiment with its setup and its results; a
campaign continues from it after an interruption.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from .optimizer import propose_batch


class Campaign:
    """A campaign of the self-driving laboratory with a robot such as psa.PressureSwingCycle.

    It maximizes the productivity and the recovery of the cycles whose purity reaches the robot's purity_target:
    rounds batches of batch_size after a start of start_size setups, workers experiments at a time.
    """

    def __init__(self, robot, directory, rounds=15, batch_size=8, start_size=16, workers=8, seed=2026):
        self.robot = robot
        self.directory = Path(directory)
        self.rounds = rounds
        self.batch_size = batch_size
        self.start_size = start_size
        self.workers = workers
        self.seed = seed
        self.log = self.directory / "experiments.csv"

    def results(self):
        """The rows of the lab notebook, as dicts of strings."""
        return list(csv.DictReader(self.log.open(encoding="utf-8"))) if self.log.exists() else []

    def run(self, verbose=True):
        """Runs the rounds that the lab notebook does not have yet, and returns its rows (see results)."""
        rows = self.results()
        rng = np.random.default_rng(self.seed + len(rows))
        finished_rounds = {int(row["round"]) for row in rows}
        target = self.robot.purity_target
        for round_index in range(self.rounds + 1):
            if round_index in finished_rounds:
                continue
            if round_index == 0:
                start = qmc.Sobol(len(self.robot.variables), seed=rng).random(self.start_size)[1:]
                setups = [self.robot.baseline] + self.robot.from_unit(start)
            else:
                setups = self.robot.from_unit(self._propose(rows, rng))

            directories = []
            for setup in setups:
                name = self.robot.name(len(rows) + len(directories), setup)
                directories.append(self.directory / f"round_{round_index:02d}" / name)
                self.robot.write(directories[-1], setup)
            self._run(directories)

            columns = self.robot.result_columns()
            for directory, setup in zip(directories, setups):
                kpis = {"status": "failed"}
                if (directory / "kpis.json").exists():
                    kpis = json.loads((directory / "kpis.json").read_text(encoding="utf-8"))
                rows.append({"index": len(rows), "round": round_index, "name": directory.name, **setup,
                             **{key: kpis.get(key, "") for key in columns}})
            with self.log.open("w", newline="", encoding="utf-8") as file:
                writer = csv.DictWriter(file, fieldnames=["index", "round", "name", *self.robot.variables, *columns])
                writer.writeheader()
                writer.writerows(rows)

            if verbose:
                feasible = [row for row in rows if row["status"] != "failed" and float(row["purity"]) >= target]
                best = max(feasible, key=lambda row: float(row["productivity"]), default=None)
                print(f"round {round_index:2d}: {len(rows)} experiments, {len(feasible)} at purity >= {target}"
                      + (f", most productive {1e3 * float(best['productivity']):.3f} mmol/(kg s) at recovery "
                         f"{float(best['recovery']):.3f} ({best['name']})" if best else ""), flush=True)
        return self.results()

    def write_inputs(self):
        """Writes the directory of every experiment in the lab notebook again, from its setup."""
        for row in self.results():
            setup = {name: float(row[name]) for name in self.robot.variables}
            self.robot.write(self.directory / f"round_{int(row['round']):02d}" / row["name"], setup)

    def rerun(self):
        """Runs every experiment directory of the campaign again, workers at a time, for instance after write_inputs."""
        self._run(sorted(path.parent for path in self.directory.glob("round_*/*/simulation.json")))

    def _propose(self, rows, rng):
        """The next batch on the unit cube, from the setups and results of the finished experiments."""
        target = -np.log10(1.0 - self.robot.purity_target)
        # A failed run counts as a purity of 0, so that the batch stays away from where runs fail. The constraint is the
        # number of nines of the purity, up to 4: what lies far above the target only distorts the model.
        purity = np.clip([0.0 if row["status"] == "failed" else float(row["purity"]) for row in rows], 0.0, 1.0)
        nines = -np.log10(np.maximum(1.0 - purity, 1.0e-4))
        # A cycle whose purge and backfill take more product than it makes counts as making none: keeping its negative
        # values sent the batches further to the bounds of the cube. A run below the target keeps (nines / target)^4 of
        # its objectives: the cycles that are far off do not pull the batch to them, while those close to the boundary,
        # where the Pareto front lies, still count.
        objectives = np.array([[max(float(row["productivity"]), 0.0), max(float(row["recovery"]), 0.0)]
                               if row["status"] != "failed" else [0.0, 0.0] for row in rows])
        objectives *= np.minimum(1.0, nines / target)[:, None] ** 4
        return propose_batch(self.robot.to_unit(rows), objectives, nines, target, self.batch_size, rng)

    def _run(self, directories):
        """Runs experiments in their own processes, workers at a time; each writes run.log and kpis.json."""

        def run(directory):
            with open(directory / "run.log", "w", encoding="utf-8") as log:
                subprocess.run([sys.executable, "-m", "ruptura.selfdriving", str(directory)], stdout=log,
                               stderr=subprocess.STDOUT, env={**os.environ, "OMP_NUM_THREADS": "1"}, check=False)

        with ThreadPoolExecutor(self.workers) as pool:
            list(pool.map(run, directories))
