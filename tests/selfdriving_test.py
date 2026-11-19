"""The self-driving laboratory: brain, robot input, analytics, and a small campaign through the Python interface."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("scipy")

from ruptura.selfdriving import (Campaign, GaussianProcess, PressureSwingCycle, TwoBedPressureSwingCycle, hypervolume,
                                 key_performance, pareto_front, propose_batch)
from ruptura.selfdriving.psa import GAS_CONSTANT, Robot

BASE_INPUT = Path(__file__).resolve().parents[1] / "examples" / "PSA-H2-activated-carbon" / "simulation.json"
TWO_BED_INPUT = Path(__file__).resolve().parents[1] / "examples" / "PSA-H2-activated-carbon-two-beds" / "simulation.json"


def test_pareto_front_and_hypervolume_of_a_known_set():
    points = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 1.0], [1.5, 1.5], [2.0, 1.0]])
    assert pareto_front(points).tolist() == [True, True, True, False, False]
    # The staircase of the front above (0, 0): 3 x 1 + 2 x (2 - 1) + 1 x (3 - 2); above (1, 1) only (2, 2) counts.
    assert hypervolume(points, np.zeros(2)) == pytest.approx(6.0)
    assert hypervolume(points, np.array([1.0, 1.0])) == pytest.approx(1.0)


def test_gaussian_process_interpolates_and_believes():
    rng = np.random.default_rng(0)
    inputs = rng.random((30, 2)) * np.array([0.5, 1.0])  # nothing beyond x = 0.5
    targets = np.sin(3.0 * inputs[:, 0]) + inputs[:, 1] ** 2
    model = GaussianProcess(inputs, targets)
    mean, deviation = model.predict(inputs[:5])
    assert mean == pytest.approx(targets[:5], abs=2.0e-2)
    # The Kriging believer: the model takes its own mean at a new point as measured, which keeps the mean and removes
    # the uncertainty there, so that the next proposal of a batch goes elsewhere.
    point = np.array([[0.9, 0.5]])
    before_mean, before_deviation = model.predict(point)
    model.believe(point[0], before_mean[0])
    after_mean, after_deviation = model.predict(point)
    assert after_mean[0] == pytest.approx(before_mean[0], abs=1.0e-3)
    assert after_deviation[0] < 0.2 * before_deviation[0]


def test_propose_batch_stays_in_the_cube_and_spreads_the_batch():
    rng = np.random.default_rng(1)
    inputs = rng.random((20, 3))
    objectives = np.column_stack([inputs[:, 0], 1.0 - inputs[:, 0] + 0.1 * inputs[:, 1]])
    constraint = 3.0 * inputs[:, 2]
    batch = propose_batch(inputs, objectives, constraint, 1.5, 4, rng)
    assert batch.shape == (4, 3)
    assert np.all((batch >= 0.0) & (batch <= 1.0))
    assert np.min(np.linalg.norm(batch[:, None, :] - batch[None, :, :], axis=-1) + np.eye(4)) > 1.0e-3


def test_cycle_input_has_the_five_steps_their_ramps_and_the_purge_flow():
    cycle = PressureSwingCycle(BASE_INPUT, product="H2")
    phases = cycle.simulation_input(cycle.baseline)["SwingAdsorptionPhases"]
    # Without backfill the cycle is the Skarstrom cycle, with 40 s ramps on the 1.2 m column between 10 and 1 bar.
    assert [phase["Name"] for phase in phases] == ["Pressurization", "Adsorption", "Blowdown", "Purge"]
    assert phases[0]["Duration"] == pytest.approx(40.0)
    assert phases[3]["ColumnEntranceVelocity"] == pytest.approx(0.08)
    assert phases[3]["GasPhaseMolFractions"] == {"H2": 1.0}

    setup = {**cycle.baseline, "backfill": 0.25, "length": 2.4}
    data = cycle.simulation_input(setup)
    backfill, pressurization = data["SwingAdsorptionPhases"][:2]
    # Four times the ramp time for twice the length; the diameter keeps the sorbent mass.
    assert backfill["Duration"] + pressurization["Duration"] == pytest.approx(160.0)
    assert backfill["Duration"] == pytest.approx(40.0)
    assert backfill["ClosedEnd"] == "FeedEnd" and pressurization["ClosedEnd"] == "ProductEnd"
    assert backfill["InletPressure"] == pytest.approx(1.0e5 + 0.25 * 9.0e5)
    assert data["Geometry"]["InternalDiameter"] == pytest.approx(0.025 / np.sqrt(2.0))
    assert data["BreakthroughIntegrator"] == "ROS34PW2" and data["ROS34PW2RelativeTolerance"] == 1.0e-5


def test_unit_cube_round_trip():
    cycle = PressureSwingCycle(BASE_INPUT)
    unit = cycle.to_unit([cycle.baseline])
    assert np.all((unit >= 0.0) & (unit <= 1.0))
    assert cycle.from_unit(unit)[0] == pytest.approx(cycle.baseline)


def test_key_performance_from_the_stream_amounts():
    # One component of each kind: helium (no feed), a strong adsorbate A, and the product H2, fed 0.5 / 0.5.
    data = {"Components": [{"Name": "Helium", "GasPhaseMolFraction": 0.0}, {"Name": "A", "GasPhaseMolFraction": 0.5},
                           {"Name": "H2", "GasPhaseMolFraction": 0.5}],
            "ColumnVoidFraction": 0.4, "ParticleDensity": 500.0,
            "SwingAdsorptionPhases": [{"Name": name} for name in ("Backfill", "Pressurization", "Adsorption", "Purge")]}

    def record(name, start, **streams):  # entered or left at z = 0 or z = L, per component
        amounts = [streams.get(key, [0.0, 0.0, 0.0]) for key in ("entered_0", "left_0", "entered_L", "left_L")]
        return SimpleNamespace(name=name + " (cycle 2)", start_time=start, end_time=start + 10.0,
                               stream_amounts=np.ravel(amounts).tolist(), gas_deviation=1.0e-3)

    cycle = [record("Backfill", 0.0, entered_L=[0.0, 0.0, 10.0]),
             record("Pressurization", 10.0, entered_0=[0.0, 25.0, 25.0]),
             record("Adsorption", 20.0, entered_0=[0.0, 200.0, 200.0], left_L=[0.0, 1.0, 180.0]),
             record("Purge", 30.0, entered_L=[0.0, 0.0, 40.0], left_0=[0.0, 224.0, 75.0])]
    kpis = key_performance(cycle + cycle, {"length": 1.0}, data, Robot(product="H2", purity_target=0.99))
    net = 180.0 - 40.0 - 10.0
    assert kpis["purity"] == pytest.approx(180.0 / 181.0)
    assert kpis["product_fraction_A"] == pytest.approx(1.0 / 181.0)
    assert kpis["recovery"] == pytest.approx(net / (200.0 + 25.0))
    assert kpis["productivity"] == pytest.approx(net / (0.6 * 500.0 * 40.0))
    assert kpis["backfill_to_feed"] == pytest.approx(10.0 / 400.0)
    assert kpis["cycle_time"] == pytest.approx(40.0)


def test_small_campaign_runs_resumes_and_rewrites_its_inputs(tmp_path):
    pytest.importorskip("ruptura._ruptura")
    cycle = PressureSwingCycle(BASE_INPUT, grid_points=12, maximum_cycles=2)
    campaign = Campaign(cycle, tmp_path / "campaign", rounds=1, batch_size=2, start_size=4, workers=4)
    rows = campaign.run(verbose=False)
    assert [int(row["round"]) for row in rows] == [0, 0, 0, 0, 1, 1]
    assert all(row["status"] in ("converged", "not_converged") for row in rows)
    assert 0.99 < float(rows[0]["purity"]) <= 1.0  # the baseline cycle
    # A finished campaign runs nothing more, and its inputs can be written again.
    assert campaign.run(verbose=False) == rows
    first = tmp_path / "campaign" / "round_00" / rows[0]["name"] / "simulation.json"
    before = first.read_text(encoding="utf-8")
    first.unlink()
    campaign.write_inputs()
    assert first.read_text(encoding="utf-8") == before
    assert json.loads((first.parent / "robot.json").read_text(encoding="utf-8"))["product"] == "H2"


def test_two_bed_input_sizes_its_valves_from_the_feed():
    cycle = TwoBedPressureSwingCycle(TWO_BED_INPUT)
    data = cycle.simulation_input({**cycle.baseline, "equalization": 0.3})
    phases = data["SwingAdsorptionPhases"]
    assert [phase["Name"] for phase in phases] == [
        "Pressurization A, blowdown B", "Adsorption A, purge B", "Equalization A to B",
        "Pressurization B, blowdown A", "Adsorption B, purge A", "Equalization B to A"]
    assert data["Columns"] == ["A", "B"] and data["Buffers"]["Tank"]["MaximumPressure"] == 1.0e6
    assert "InletPressure" not in data and data["OutletPressure"] == 1.0e5
    assert phases[0]["Duration"] == pytest.approx(40.0) and phases[1]["Duration"] == pytest.approx(30.0)
    # The product valve passes the feed flow at 2 % over the tank, the purge valve G times it between P_H and P_L.
    feed = 0.403 * 0.03 * 1.0e6 / (GAS_CONSTANT * 298.0) * np.pi * 0.025 ** 2 / 4.0
    product, purge = phases[1]["Connections"]
    assert product["Conductance"] * 0.02e6 == pytest.approx(feed) and product["To"] == "Tank"
    assert purge["Conductance"] * 9.0e5 == pytest.approx(cycle.baseline["purge_to_feed"] * feed)
    assert purge["From"] == "A.ProductEnd" and purge["To"] == "B.ProductEnd"
    assert phases[2]["StopCriterion"]["Column"] == "B"
    assert phases[2]["StopCriterion"]["Value"] == pytest.approx(1.0e5 + 0.3 * 9.0e5)
    assert phases[4]["Columns"]["B"]["FeedEnd"] == {"Mode": "Velocity", "Velocity": 0.03, "Pressure": 1.0e6}
    # Without equalization it is the Skarstrom cycle on two beds.
    assert len(cycle.simulation_input(cycle.baseline)["SwingAdsorptionPhases"]) == 4
    assert "backfill_to_feed" not in cycle.result_columns() and "equalization" in cycle.variables
    assert cycle.name(3, {**cycle.baseline, "equalization": 0.3}).endswith("_EQ0.30")


def test_key_performance_of_two_beds_counts_what_reaches_the_tank():
    data = {"Columns": ["A", "B"],
            "Components": [{"Name": "Helium", "GasPhaseMolFraction": 0.0}, {"Name": "A", "GasPhaseMolFraction": 0.5},
                           {"Name": "H2", "GasPhaseMolFraction": 0.5}],
            "ColumnVoidFraction": 0.4, "ParticleDensity": 500.0, "SwingAdsorptionPhases": [{}, {}, {}, {}]}
    zero = [0.0, 0.0, 0.0]

    def record(phase, entered_0=zero, left_0=zero, entered_L=zero, left_L=zero):
        return SimpleNamespace(name=f"Phase {phase} (cycle 2)", start_time=10.0 * phase, end_time=10.0 * phase + 10.0,
                               stream_amounts=np.ravel([entered_0, left_0, entered_L, left_L]).tolist(),
                               gas_deviation=0.0)

    # Pressurization of one bed and blowdown of the other, then adsorption in it while it purges the other with 10 % of
    # what leaves it; the second half the other way round.
    cycle = []
    for half in (0, 1):
        pressurized = record(2 * half, entered_0=[0.0, 10.0, 10.0])
        blown_down = record(2 * half, left_0=[0.0, 8.0, 2.0])
        adsorbing = record(2 * half + 1, entered_0=[0.0, 200.0, 200.0], left_L=[0.0, 2.0, 180.0])
        purged = record(2 * half + 1, entered_L=[0.0, 0.2, 18.0], left_0=[0.0, 150.0, 40.0])
        cycle += [pressurized, blown_down, adsorbing, purged] if half == 0 else [blown_down, pressurized, purged,
                                                                                adsorbing]
    kpis = key_performance(cycle + cycle, {"length": 1.0}, data, Robot(product="H2", purity_target=0.99))
    product = 2.0 * (180.0 - 18.0)
    assert kpis["purity"] == pytest.approx(product / (product + 2.0 * 1.8))
    assert kpis["recovery"] == pytest.approx(product / 420.0)
    assert kpis["productivity"] == pytest.approx(product / (2.0 * 0.6 * 500.0 * 40.0))
    assert "backfill_to_feed" not in kpis


def test_small_two_bed_campaign_runs(tmp_path):
    pytest.importorskip("ruptura._ruptura")
    cycle = TwoBedPressureSwingCycle(TWO_BED_INPUT, grid_points=12, maximum_cycles=2)
    rows = Campaign(cycle, tmp_path / "campaign", rounds=1, batch_size=2, start_size=4, workers=4).run(verbose=False)
    assert [int(row["round"]) for row in rows] == [0, 0, 0, 0, 1, 1]
    assert all(row["status"] in ("converged", "not_converged") for row in rows)
    assert 0.9 < float(rows[0]["purity"]) <= 1.0 and float(rows[0]["recovery"]) > 0.0  # the baseline cycle
    assert (tmp_path / "campaign" / "round_00" / rows[0]["name"] / "A" / "phases.data").exists()
