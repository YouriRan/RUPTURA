"""The self-driving laboratory: brain, robot input, analytics, and a small campaign through the Python interface."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("scipy")

from ruptura.selfdriving import (Campaign, GaussianProcess, PressureSwingCycle, hypervolume, key_performance,
                                 pareto_front, propose_batch)
from ruptura.selfdriving.psa import Robot

BASE_INPUT = Path(__file__).resolve().parents[1] / "examples" / "PSA-H2-activated-carbon" / "simulation.json"


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


def test_key_performance_from_the_phase_records():
    # One component of each kind: helium (no feed), a strong adsorbate A, and the product H2, fed 0.5 / 0.5.
    data = {"Components": [{"Name": "Helium", "GasPhaseMolFraction": 0.0}, {"Name": "A", "GasPhaseMolFraction": 0.5},
                           {"Name": "H2", "GasPhaseMolFraction": 0.5}],
            "ColumnVoidFraction": 0.4, "ParticleDensity": 500.0,
            "SwingAdsorptionPhases": [{"Name": name} for name in ("Backfill", "Pressurization", "Adsorption")]}
    setup = {"length": 1.0, "high_pressure": 8.31446261815324e5, "temperature": 100.0, "feed_velocity": 0.1,
             "adsorption_time": 10.0, "purge_to_feed": 0.1}
    fed = 0.4 * 0.1 * 1000.0 * 10.0 * 0.5  # 200 mol/m^2 of each

    def record(name, change, start=0.0):
        return SimpleNamespace(name=name, start_time=start, end_time=start + 10.0, gas_at_start=[0.0, 0.0, 0.0],
                               adsorbed_at_start=[0.0, 0.0, 0.0], gas_at_end=change, adsorbed_at_end=[0.0, 0.0, 0.0],
                               gas_deviation=1.0e-3)

    cycle = [record("Backfill", [0.0, 0.0, 10.0]), record("Pressurization", [0.0, 25.0, 25.0], 10.0),
             record("Adsorption", [0.0, 199.0, 20.0], 20.0)]
    kpis = key_performance(cycle + cycle, setup, data, Robot(product="H2", purity_target=0.99))
    net = 180.0 - 0.1 * 2.0 * fed - 10.0
    assert kpis["purity"] == pytest.approx(180.0 / 181.0)
    assert kpis["product_fraction_A"] == pytest.approx(1.0 / 181.0)
    assert kpis["recovery"] == pytest.approx(net / (fed + 25.0))
    assert kpis["productivity"] == pytest.approx(net / (0.6 * 500.0 * 30.0))
    assert kpis["backfill_to_feed"] == pytest.approx(10.0 / (2.0 * fed))
    assert kpis["cyclic_steady_state_change"] == pytest.approx(kpis["cycle_amount_change"])


def test_small_campaign_runs_resumes_and_rewrites_its_inputs(tmp_path):
    pytest.importorskip("ruptura._ruptura")
    cycle = PressureSwingCycle(BASE_INPUT, grid_points=12, cycles_per_block=2, maximum_blocks=1)
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
