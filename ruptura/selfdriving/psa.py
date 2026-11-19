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

The analytics read the stream amounts of the phase records: what entered and what left at each end of the column in
each phase, per component and unit cross-section. With N_P,j what left at z = L in the adsorption step, N_F,j what
entered at z = 0 in it and in the feed pressurization, and the product used that which entered at z = L in the purge
and the backfill,

    purity       = N_P,product / sum_j N_P,j
    recovery     = (N_P,product - used) / N_F,product
    productivity = (N_P,product - used) / ((1 - eps) rho_p L t_cycle)      [mol/(kg s)]

A run ends at its cyclic steady state (CyclicSteadyState of the input): the first cycle whose state, amounts, and
streams are within the tolerance of those of the cycle before, or after maximum_cycles.

TwoBedPressureSwingCycle runs the cycle on two such columns, A and B, half a cycle apart, on a base input with
BoundaryCondition MomentumBalance. Valves between their product ends and a product tank take the place of the product
that the one-column cycle takes from outside for purge and backfill, and a pressure equalization, with the share e of
P_H - P_L, takes the place of the backfill. Each half cycle has three phases, here for A adsorbing:

1. pressurization of A and blowdown of B through their feed ends, linearly to P_H and P_L in t_r;
2. adsorption in A, the feed at v_F and P_H, while B is purged at P_L: the product end of A feeds the tank (held at
   P_H, a check valve) and the product end of B (a check valve) through valves that pass the molar feed flow
   F = eps v_F P_H A / (R T) at 2 % above P_H and G F between P_H and P_L;
3. equalization: a valve from the product end of A to that of B until the mean pressure of B is P_L + e (P_H - P_L),
   at most t_r; it evens out the gas of a column, eps A L P / (R T), with a time constant of a tenth of t_r. An
   equalization below 1 % is left out; with e = 0 the cycle is the Skarstrom cycle on two beds.

The product is what reaches the tank, the net amount through the product ends of both columns, N_P,j = sum over the
phases and columns of (left - entered) at z = L; nothing is used from outside, and the productivity is per kg of the
sorbent of both columns.
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field, fields
from pathlib import Path

import numpy as np

from ..webapp.config import canonicalize

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
# Two beds: the pressure equalization, as a share of P_H - P_L, in place of the backfill; the baseline is the same
# Skarstrom cycle on two beds.
TWO_BED_VARIABLES = {**{name: bounds for name, bounds in VARIABLES.items() if name != "backfill"},
                     "equalization": (0.0, 0.5, "linear")}
TWO_BED_BASELINE = {**{name: value for name, value in BASELINE.items() if name != "backfill"}, "equalization": 0.0}


@dataclass
class Robot:
    """What the process of one experiment needs besides its input: the product and its target."""

    product: str = "H2"
    purity_target: float = 0.999


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
    maximum_cycles: int = 100
    cyclic_steady_state_tolerance: float = 2.0e-3
    reference_ramp_time: float = 40.0  # s, for the base length and P_H / P_L = 10
    shortest_ramp_time: float = 10.0  # s
    shortest_phase_time: float = 0.1  # s

    def __post_init__(self):
        self.base_input = Path(self.base_input)
        self.base = canonicalize(json.loads(self.base_input.read_text(encoding="utf-8")))  # keys in any case
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
            NumberOfCycles=self.maximum_cycles,
            CyclicSteadyState={"Tolerance": self.cyclic_steady_state_tolerance},
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


@dataclass
class TwoBedPressureSwingCycle(PressureSwingCycle):
    """The cycle on two columns with valves, a product tank, and pressure equalization (see the module description).

    base_input is a swing-adsorption input with BoundaryCondition MomentumBalance and Columns, such as
    examples/PSA-H2-activated-carbon-two-beds; its phases, columns, and buffers are replaced.
    """

    variables: dict = field(default_factory=lambda: dict(TWO_BED_VARIABLES))
    baseline: dict = field(default_factory=lambda: dict(TWO_BED_BASELINE))

    def result_columns(self):
        """The columns of the lab notebook that the KPIs of an experiment fill: those of one column without backfill."""
        return [name for name in super().result_columns() if name != "backfill_to_feed"]

    @staticmethod
    def name(index, setup):
        """An explicit directory name: the experiment number and its setup."""
        return (
            f"exp_{index:03d}_T{setup['temperature']:.0f}K_PH{setup['high_pressure'] / 1e5:.2f}bar"
            f"_PL{setup['low_pressure'] / 1e5:.2f}bar_vf{setup['feed_velocity']:.3f}ms"
            f"_tads{setup['adsorption_time']:.0f}s"
            f"_PF{setup['purge_to_feed']:.3f}_L{setup['length']:.2f}m_EQ{setup['equalization']:.2f}"
        )

    def simulation_input(self, setup):
        """The RUPTURA input of a setup: that of the one-column cycle without backfill, with the columns, the tank, and
        the six phases of the module description."""
        data = super().simulation_input({**setup, "backfill": 0.0})
        ramp_time = data["SwingAdsorptionPhases"][0]["Duration"]  # the feed pressurization, without backfill
        temperature, high, low = setup["temperature"], setup["high_pressure"], setup["low_pressure"]
        del data["InletPressure"]
        data["OutletPressure"] = low  # the columns start at rest at P_L
        area = math.pi * data["Geometry"]["InternalDiameter"] ** 2 / 4.0
        feed = data["ColumnVoidFraction"] * setup["feed_velocity"] * high / (GAS_CONSTANT * temperature) * area  # mol/s
        holdup = data["ColumnVoidFraction"] * area * setup["length"] / (GAS_CONSTANT * temperature)  # mol/Pa
        data["Columns"] = ["A", "B"]
        data["Buffers"] = {"Tank": {"Volume": area * setup["length"], "Pressure": high, "MaximumPressure": high,
                                    "Temperature": temperature, "GasPhaseMolFractions": {self.product: 1.0}}}
        ramp = {"Shape": "Linear", "TimeConstant": ramp_time}
        phases = []
        for bed, other in (("A", "B"), ("B", "A")):
            phases += [
                {"Name": f"Pressurization {bed}, blowdown {other}", "Temperature": temperature, "Duration": ramp_time,
                 "Columns": {bed: {"FeedEnd": {"Mode": "Pressure", "Pressure": high, "PressureRamp": ramp}},
                             other: {"FeedEnd": {"Mode": "Pressure", "Pressure": low, "PressureRamp": ramp}}}},
                {"Name": f"Adsorption {bed}, purge {other}", "Temperature": temperature,
                 "Duration": setup["adsorption_time"],
                 "Columns": {bed: {"FeedEnd": {"Mode": "Velocity", "Velocity": setup["feed_velocity"],
                                               "Pressure": high}},
                             other: {"FeedEnd": {"Mode": "Pressure", "Pressure": low}}},
                 "Connections": [
                     {"Name": "Product", "From": f"{bed}.ProductEnd", "To": "Tank", "Conductance": feed / (0.02 * high),
                      "CheckValve": True},
                     {"Name": "Purge", "From": f"{bed}.ProductEnd", "To": f"{other}.ProductEnd",
                      "Conductance": setup["purge_to_feed"] * feed / (high - low), "CheckValve": True}]},
                {"Name": f"Equalization {bed} to {other}", "Temperature": temperature, "Duration": ramp_time,
                 "Connections": [{"Name": "Equalization", "From": f"{bed}.ProductEnd", "To": f"{other}.ProductEnd",
                                  "Conductance": 5.0 * holdup / ramp_time, "CheckValve": True}],
                 "StopCriterion": {"Column": other, "Metric": "AveragePressure", "Operator": ">=",
                                   "Value": low + setup["equalization"] * (high - low)}},
            ]
        data["SwingAdsorptionPhases"] = [phase for phase in phases if
                                         not phase["Name"].startswith("Equalization") or setup["equalization"] >= 0.01]
        return data


def key_performance(records, setup, data, robot):
    """The KPIs of the last cycle of a run (see the module description for the equations).

    The stream amounts of a record are what entered at z = 0, left at z = 0, entered at z = L, and left at z = L in its
    phase, each per component in mol/m^2 of column cross-section.
    """
    names = [component["Name"] for component in data["Components"]]
    product = names.index(robot.product)
    beds = max(1, len(data.get("Columns", [])))  # the records of a phase are those of every column in turn
    sorbent_mass = beds * (1.0 - data["ColumnVoidFraction"]) * data["ParticleDensity"] * setup["length"]  # kg/m^2
    cycle = records[-beds * len(data["SwingAdsorptionPhases"]):]
    streams = [np.reshape(record.stream_amounts, (4, len(names))) for record in cycle]
    fed = sum(stream[0] for stream in streams)  # the feed enters at z = 0 only
    kpis = {}
    if beds > 1:
        # The product ends of the columns feed each other and the tank: what reaches the tank is their net outflow.
        out = sum(stream[3] - stream[2] for stream in streams)
        net = out[product]
    else:
        # The product leaves in the adsorption phase; purge and backfill take product from outside.
        nothing = np.zeros((4, len(names)))
        phase = {record.name.split(" (cycle ")[0]: stream for record, stream in zip(cycle, streams)}  # "Adsorption (cycle 3)"
        out = phase["Adsorption"][3]
        backfilled = phase.get("Backfill", nothing)[2][product]
        net = out[product] - phase.get("Purge", nothing)[2][product] - backfilled
        kpis["backfill_to_feed"] = float(backfilled / phase["Adsorption"][0].sum())
    total = out.sum() if out.sum() > 0.0 else np.inf  # without a product the purity is 0
    cycle_time = cycle[-1].end_time - cycle[0].start_time
    return {
        "purity": float(out[product] / total),
        **{f"product_fraction_{name}": float(out[k] / total) for k, name in enumerate(names)
           if k != product and data["Components"][k]["GasPhaseMolFraction"] > 0.0},
        "recovery": float(net / fed[product]),
        "productivity": float(net / (sorbent_mass * cycle_time)),  # mol/(kg s)
        **kpis,
        "cycle_time": float(cycle_time),
        "gas_deviation": float(max(record.gas_deviation for record in cycle)),
    }


def run_experiment(directory):
    """Runs the experiment of a directory through the Python interface and writes its KPIs to kpis.json.

    cycles.data of the run holds, for its last cycle, the changes from the cycle before, which the cyclic steady state
    is measured by. A failed integration is a result of the experiment.
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
        simulation = ruptura.load_simulation("simulation.json")
        simulation.run()
        kpis = key_performance(simulation.phase_records, setup, data, robot)
        last = np.loadtxt("cycles.data", ndmin=2)[-1]  # cycle, time, change of state, amounts, streams, mass balance
        kpis["cycle_amount_change"] = float(last[3])
        kpis["cyclic_steady_state_change"] = float(last[2:5].max())
        kpis["cycles"] = int(last[0])
        converged = kpis["cyclic_steady_state_change"] <= data["CyclicSteadyState"]["Tolerance"]
        kpis["status"] = "converged" if converged else "not_converged"
    except Exception as error:
        kpis = {"status": "failed", "error": str(error)}
    finally:
        os.chdir(working_directory)
    kpis["wall_time"] = time.perf_counter() - start
    (directory / "kpis.json").write_text(json.dumps(kpis, indent=2) + "\n", encoding="utf-8")
    return kpis
