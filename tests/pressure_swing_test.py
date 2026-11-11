"""Phase records of a PSA cycle through the Python interface."""

import json
from pathlib import Path

import numpy as np
import pytest

import ruptura

pytest.importorskip("ruptura._ruptura")

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "PSA-H2-activated-carbon" / "simulation.json"


def test_phase_records_hold_the_cycle_mass_balance(tmp_path, monkeypatch):
    data = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    data.update(NumberOfGridPoints=12, NumberOfCycles=2, WriteInterval=1000.0)
    (tmp_path / "simulation.json").write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    simulation = ruptura.load_simulation("simulation.json")
    simulation.run()

    records = simulation.phase_records
    assert [record.name for record in records[:2]] == ["Pressurization (cycle 1)", "Adsorption (cycle 1)"]
    assert len(records) == 8
    assert records[0].pressure_at_zero == pytest.approx(1.0e6, rel=5.0e-2)
    assert records[2].pressure_at_zero == pytest.approx(1.0e5, rel=5.0e-2)

    # Only the feed enters a pressurization.
    pressurization = records[4]
    entered = (
        np.array(pressurization.gas_at_end)
        + pressurization.adsorbed_at_end
        - pressurization.gas_at_start
        - pressurization.adsorbed_at_start
    )
    assert entered[1:] / entered.sum() == pytest.approx([0.4, 0.3, 0.3], abs=1.0e-3)

    # phases.data holds the same numbers, one row per phase.
    table = np.loadtxt("phases.data")
    assert table.shape == (8, 5 + 4 * 4)
    assert table[4, 1] == pytest.approx(pressurization.end_time)
    assert table[4, 4] == pytest.approx(pressurization.gas_deviation)
    assert 0.0 < pressurization.gas_deviation < 1.0e-2
    assert table[4, 5:9] == pytest.approx(
        [
            pressurization.gas_at_start[0],
            pressurization.adsorbed_at_start[0],
            pressurization.gas_at_end[0],
            pressurization.adsorbed_at_end[0],
        ]
    )
