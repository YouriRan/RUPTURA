"""The robot and the analytics of the self-driving laboratory: one pressure-swing cycle, run to its cyclic steady state.

A setup sets the temperature T of the feed and the surroundings, the adsorption and purge pressures P_H and P_L, the
interstitial feed velocity v_F, the adsorption time t_ads (that of the purge as well), the purge-to-feed ratio G, the
column length L, and the backfill share b. The cycle has five steps:

1. backfill: product enters at z = L, the feed end closed, while the pressure rises linearly from P_L to
   P_L + b (P_H - P_L) in b t_r;
2. feed pressurization: the feed of the base input enters at z = 0, the product end closed, while the pressure rises
   linearly to P_H in (1 - b) t_r;
3. adsorption: the feed enters at z = 0 at v_F and P_H for t_ads, and the product leaves at z = L;
4. blowdown: gas leaves at z = 0, the product end closed, while the pressure falls linearly to P_L in t_r;
5. purge: product enters at z = L at P_L for t_ads with v_P = G v_F P_H / P_L, so that G is the product purged per
   mole of feed.

A backfill or feed pressurization shorter than 0.1 s is left out; with b = 0 the cycle is the Skarstrom cycle. The
ramp time t_r = max(10 s, t_ref (L/L_0)^2 (P_H - P_L) / (9 P_L)) keeps the pressure drop of the gas that a ramp moves,
which grows as L^2 (dP/dt) / P_L, small against P_L: the closures follow the pressure of the open end only. t_ref = 40 s
holds it below 1 % for the 1.2 m column of Marx T15 between 10 and 1 bar. The column keeps the sorbent mass of the base
input, D = D_0 (L_0 / L)^(1/2), and its wall thickness.

The analytics read the phase records, the gas and adsorbed amounts n_j per unit cross-section at the start and end of
each phase; their change dn_j is what entered minus what left. With the feed N_F,j = eps v_F P_H / (R T) t_ads y_F,j,
the product of the adsorption step is N_P,j = max(N_F,j - dn_j(adsorption), 0), the product used is
G sum_j N_F,j + dn_product(backfill), through whose closed feed end nothing leaves, and

    purity       = N_P,product / sum_j N_P,j
    recovery     = (N_P,product - used) / (N_F,product + dn_product(feed pressurization))
    productivity = (N_P,product - used) / ((1 - eps) rho_p L t_cycle)      [mol/(kg s)]

A run goes on in blocks of cycles, each from the column state at the end of the previous one, until in the last cycle
the amount of every component in the column changes by less than the tolerance times the feed of a cycle, the recovery
and the productivity by less than the tolerance relative to their values, and the purity by less than the tolerance
times its margin to the target.
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field, fields
from pathlib import Path

import numpy as np

GAS_CONSTANT = 8.31446261815324  # J/(mol K)

# The setup the optimizer may choose, name: (lower bound, upper bound, scale of the search).
VARIABLES = {
    "temperature": (288.0, 328.0, "linear"),  # K, of the feed and the column surroundings
    "high_pressure": (4.0e5, 1.0e6, "linear"),  # Pa
    "low_pressure": (1.0e5, 2.0e5, "linear"),  # Pa
    "feed_velocity": (0.01, 0.06, "log"),  # m/s, interstitial
    "adsorption_time": (20.0, 200.0, "log"),  # s, also that of the purge
    "purge_to_feed": (0.02, 0.6, "log"),  # [-], purge over feed in mol during the adsorption time
    "length": (0.6, 3.6, "log"),  # m, at the sorbent mass of the base column, which sets the diameter
    "backfill": (0.0, 1.0, "linear"),  # [-], share of the pressure rise made with product from z = L
}
# A hand-tuned Skarstrom cycle for the T15 column: 10/1 bar, 0.03 m/s, 30 s steps, a purge of 0.08 m/s.
BASELINE = {
    "temperature": 298.0,
    "high_pressure": 1.0e6,
    "low_pressure": 1.0e5,
    "feed_velocity": 0.03,
    "adsorption_time": 30.0,
    "purge_to_feed": 0.08 * 1.0e5 / (0.03 * 1.0e6),
    "length": 1.2,
    "backfill": 0.0,
}


@dataclass
class Robot:
    """What the process of one experiment needs besides its input: the product, the target, and when to stop."""

    product: str = "H2"
    purity_target: float = 0.999
    cycles_per_block: int = 10
    maximum_blocks: int = 10
    cyclic_steady_state_tolerance: float = 2.0e-3


@dataclass
class PressureSwingCycle:
    """A five-step pressure-swing cycle on the column of a base input (see the module description).

    base_input is a RUPTURA swing-adsorption input whose feed is the gas to separate; product names its component that
    the cycle purifies, and purity_target the mole fraction of the product that a cycle must reach. variables and
    baseline override the search space and the start cycle of a campaign; grid_points, integrator, and
    relative_tolerance set the fidelity of the robot.
    """

    base_input: Path | str
    product: str = "H2"
    purity_target: float = 0.999
    variables: dict = field(default_factory=lambda: dict(VARIABLES))
    baseline: dict = field(default_factory=lambda: dict(BASELINE))
    grid_points: int = 80
    integrator: str = "ROS34PW2"
    relative_tolerance: float = 1.0e-5
    cycles_per_block: int = 10
    maximum_blocks: int = 10
    cyclic_steady_state_tolerance: float = 2.0e-3
    reference_ramp_time: float = 40.0  # s, for the base length and P_H / P_L = 10
    shortest_ramp_time: float = 10.0  # s
    shortest_phase_time: float = 0.1  # s

    def __post_init__(self):
        self.base_input = Path(self.base_input)
        self.base = json.loads(self.base_input.read_text(encoding="utf-8"))
        self.logarithmic = np.array([bounds[2] == "log" for bounds in self.variables.values()])
        lower = np.array([bounds[0] for bounds in self.variables.values()])
        upper = np.array([bounds[1] for bounds in self.variables.values()])
        self.low = np.where(self.logarithmic, np.log(np.where(self.logarithmic, lower, 1.0)), lower)
        self.high = np.where(self.logarithmic, np.log(upper), upper)

    def to_unit(self, setups):
        """Points on the unit cube of the optimizer from setups (dicts or rows with the variables), the logarithmic
        variables on a log scale."""
        values = np.array([[float(setup[name]) for name in self.variables] for setup in setups])
        values[:, self.logarithmic] = np.log(values[:, self.logarithmic])
        return (values - self.low) / (self.high - self.low)

    def from_unit(self, points):
        """Setups (dicts in SI units) from points on the unit cube."""
        values = self.low + np.atleast_2d(points) * (self.high - self.low)
        values[:, self.logarithmic] = np.exp(values[:, self.logarithmic])
        return [dict(zip(self.variables, map(float, row))) for row in values]

    def result_columns(self):
        """The columns of the lab notebook that the KPIs of an experiment fill."""
        others = [component["Name"] for component in self.base["Components"]
                  if component["Name"] != self.product and component["GasPhaseMolFraction"] > 0.0]
        return ["status", "purity", *[f"product_fraction_{name}" for name in others], "recovery", "productivity",
                "backfill_to_feed", "cycle_time", "gas_deviation", "cycle_amount_change", "cyclic_steady_state_change",
                "cycles", "wall_time"]

    @staticmethod
    def name(index, setup):
        """An explicit directory name: the experiment number and its setup."""
        return (
            f"exp_{index:03d}_T{setup['temperature']:.0f}K_PH{setup['high_pressure'] / 1e5:.2f}bar"
            f"_PL{setup['low_pressure'] / 1e5:.2f}bar_vf{setup['feed_velocity']:.3f}ms"
            f"_tads{setup['adsorption_time']:.0f}s"
            f"_PF{setup['purge_to_feed']:.3f}_L{setup['length']:.2f}m_BF{setup['backfill']:.2f}"
        )

    def simulation_input(self, setup):
        """The RUPTURA input of a setup: the base input with what the setup and the fidelity change."""
        data = json.loads(json.dumps(self.base))
        swing = (setup["high_pressure"] - setup["low_pressure"]) / setup["low_pressure"]
        ramp_time = max(self.shortest_ramp_time,
                        self.reference_ramp_time * (setup["length"] / data["ColumnLength"]) ** 2 * swing / 9.0)
        # The diameter of a column with the sorbent mass of the base column, and the wall thickness of the base column.
        geometry = data["Geometry"]
        diameter = geometry["InternalDiameter"] * math.sqrt(data["ColumnLength"] / setup["length"])
        wall = geometry["OuterDiameter"] - geometry["InternalDiameter"]
        geometry.update(InternalDiameter=diameter, OuterDiameter=diameter + wall)
        for key in [key for key in data if key.lower().startswith(("cvode", "ros34pw2"))]:
            del data[key]
        data.update(
            Temperature=setup["temperature"],
            InfluxTemperature=setup["temperature"],
            InletPressure=setup["low_pressure"],
            ColumnEntranceVelocity=setup["feed_velocity"],
            ColumnLength=setup["length"],
            NumberOfGridPoints=self.grid_points,
            BreakthroughIntegrator=self.integrator,
            NumberOfCycles=self.cycles_per_block,
            WriteInterval=1.0e7,
        )
        data[f"{self.integrator}RelativeTolerance"] = self.relative_tolerance
        # The backfill raises the pressure from P_L by the share b of P_H - P_L with product that enters at z = L, the
        # feed end closed; the feed then raises it to P_H from z = 0, the product end closed. Both keep the ramp's rate.
        backfill_time = setup["backfill"] * ramp_time
        swing_pressure = setup["high_pressure"] - setup["low_pressure"]
        intermediate_pressure = setup["low_pressure"] + setup["backfill"] * swing_pressure
        product_gas = {self.product: 1.0}
        phases = [
            {"Name": "Backfill", "ClosedEnd": "FeedEnd", "InletPressure": intermediate_pressure,
             "Temperature": setup["temperature"], "GasPhaseMolFractions": product_gas, "Duration": backfill_time,
             "PressureRamp": {"Shape": "Linear", "TimeConstant": backfill_time}},
            {"Name": "Pressurization", "ClosedEnd": "ProductEnd", "InletPressure": setup["high_pressure"],
             "Temperature": setup["temperature"], "Duration": ramp_time - backfill_time,
             "PressureRamp": {"Shape": "Linear", "TimeConstant": ramp_time - backfill_time}},
            {"Name": "Adsorption", "InletPressure": setup["high_pressure"], "Temperature": setup["temperature"],
             "ColumnEntranceVelocity": setup["feed_velocity"], "Duration": setup["adsorption_time"]},
            {"Name": "Blowdown", "ClosedEnd": "ProductEnd", "InletPressure": setup["low_pressure"],
             "Temperature": setup["temperature"], "Duration": ramp_time,
             "PressureRamp": {"Shape": "Linear", "TimeConstant": ramp_time}},
            {"Name": "Purge", "InletPressure": setup["low_pressure"], "Temperature": setup["temperature"],
             "ColumnEntranceVelocity": setup["purge_to_feed"] * setup["feed_velocity"] * setup["high_pressure"]
             / setup["low_pressure"],
             "FlowDirection": "Reverse", "GasPhaseMolFractions": product_gas, "Duration": setup["adsorption_time"]},
        ]
        data["SwingAdsorptionPhases"] = [phase for phase in phases if phase["Duration"] >= self.shortest_phase_time]
        return data

    def write(self, directory, setup):
        """Writes one experiment: its setup, its input, and what its process needs to run it (robot.json)."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "setup.json").write_text(json.dumps(setup, indent=2) + "\n", encoding="utf-8")
        (directory / "simulation.json").write_text(json.dumps(self.simulation_input(setup), indent=2) + "\n",
                                                   encoding="utf-8")
        robot = {item.name: getattr(self, item.name) for item in fields(Robot)}
        (directory / "robot.json").write_text(json.dumps(robot, indent=2) + "\n", encoding="utf-8")


def key_performance(records, setup, data, robot):
    """The KPIs of the last cycle of a run, and how far it is from the cyclic steady state.

    The records hold the amounts per unit cross-section, in mol/m^2 (see the module description for the equations).
    """
    names = [component["Name"] for component in data["Components"]]
    product = names.index(robot.product)
    fractions = np.array([component["GasPhaseMolFraction"] for component in data["Components"]])
    void_fraction = data["ColumnVoidFraction"]
    sorbent_mass = (1.0 - void_fraction) * data["ParticleDensity"] * setup["length"]  # kg/m^2
    molar_density = setup["high_pressure"] / (GAS_CONSTANT * setup["temperature"])  # mol/m^3
    phases_per_cycle = len(data["SwingAdsorptionPhases"])
    amount = lambda record, end: (np.array(record.gas_at_end if end else record.gas_at_start)
                                  + (record.adsorbed_at_end if end else record.adsorbed_at_start))
    change = lambda record: amount(record, True) - amount(record, False)

    def performance(cycle):
        phase = {record.name.split(" (cycle ")[0]: record for record in cycle}  # "Adsorption (cycle 3)"
        fed = void_fraction * setup["feed_velocity"] * molar_density * setup["adsorption_time"] * fractions
        purged = setup["purge_to_feed"] * fed.sum()
        backfilled = change(phase["Backfill"])[product] if "Backfill" in phase else 0.0
        pressurized = change(phase["Pressurization"])[product] if "Pressurization" in phase else 0.0
        # A component that the column takes up beyond its feed (gas then flows back in at z = L) is not in the product;
        # without a product the purity is 0.
        out = np.maximum(fed - change(phase["Adsorption"]), 0.0)
        net = out[product] - purged - backfilled
        total = out.sum() if out.sum() > 0.0 else np.inf
        cycle_time = cycle[-1].end_time - cycle[0].start_time
        return {
            "purity": float(out[product] / total),
            **{f"product_fraction_{name}": float(out[k] / total) for k, name in enumerate(names)
               if k != product and fractions[k] > 0.0},
            "recovery": float(net / (fed[product] + pressurized)),
            "productivity": float(net / (sorbent_mass * cycle_time)),  # mol/(kg s)
            "backfill_to_feed": float(backfilled / fed.sum()),
            "cycle_time": float(cycle_time),
            "gas_deviation": float(max(record.gas_deviation for record in cycle)),
            "cycle_amount_change": float(np.max(np.abs(amount(cycle[-1], True) - amount(cycle[0], False))) / fed.sum()),
        }

    last = performance(records[-phases_per_cycle:])
    previous = performance(records[-2 * phases_per_cycle:-phases_per_cycle])
    relative = lambda now, before: abs(now - before) / max(abs(now), 1.0e-9)
    last["cyclic_steady_state_change"] = max(
        last["cycle_amount_change"],
        abs(last["purity"] - previous["purity"]) / (1.0 - robot.purity_target),
        relative(last["recovery"], previous["recovery"]),
        relative(last["productivity"], previous["productivity"]),
    )
    return last


def run_experiment(directory):
    """Runs the experiment of a directory through the Python interface and writes its KPIs to kpis.json.

    The run goes on in blocks of cycles, each from the column state at the end of the previous one, until the cyclic
    steady state or the last block. A failed integration is a result of the experiment.
    """
    import ruptura

    directory = Path(directory).resolve()
    start = time.perf_counter()
    setup = json.loads((directory / "setup.json").read_text(encoding="utf-8"))
    data = json.loads((directory / "simulation.json").read_text(encoding="utf-8"))
    robot = Robot(**json.loads((directory / "robot.json").read_text(encoding="utf-8")))
    kpis = {"status": "failed"}
    working_directory = Path.cwd()
    try:
        os.chdir(directory)  # RUPTURA writes its output next to the input
        input_name = "simulation.json"
        for block in range(1, robot.maximum_blocks + 1):
            simulation = ruptura.load_simulation(input_name)
            simulation.run()
            kpis = key_performance(simulation.phase_records, setup, data, robot)
            kpis["cycles"] = block * robot.cycles_per_block
            kpis["status"] = "converged"
            if kpis["cyclic_steady_state_change"] < robot.cyclic_steady_state_tolerance:
                break
            kpis["status"] = "not_converged"
            Path("continued.json").write_text(json.dumps({**data, "ReadColumnFile": "column.json"}), encoding="utf-8")
            input_name = "continued.json"
    except Exception as error:
        kpis = {"status": "failed", "error": str(error)}
    finally:
        os.chdir(working_directory)
    kpis["wall_time"] = time.perf_counter() - start
    (directory / "kpis.json").write_text(json.dumps(kpis, indent=2) + "\n", encoding="utf-8")
    return kpis
