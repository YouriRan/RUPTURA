import sys
import types

import numpy as np
import pytest

try:
    import plotly.graph_objects  # noqa: F401
except ModuleNotFoundError:
    plotly = types.ModuleType("plotly")
    graph_objects = types.ModuleType("plotly.graph_objects")
    graph_objects.Figure = object
    graph_objects.Scatter = object
    graph_objects.Frame = object
    graph_objects.Layout = object
    plotly.graph_objects = graph_objects
    sys.modules["plotly"] = plotly
    sys.modules["plotly.graph_objects"] = graph_objects

try:
    import pandas  # noqa: F401
except ModuleNotFoundError:
    pandas = types.ModuleType("pandas")
    sys.modules["pandas"] = pandas

from ruptura.plot_breakthrough import (
    BreakthroughPlotly,
    COMPONENT_METRICS,
    _column_explorer_y_options,
)


def test_breakthrough_y_units_match_component_writer_columns():
    plotter = BreakthroughPlotly(
        displayName="test",
        externalTemperature=298.15,
        externalPressure=101325.0,
        components=[],
    )
    data = np.zeros((2, 14), dtype=float)
    data[:, COMPONENT_METRICS["C"].col_0based] = [12.0, 24.0]
    data[:, COMPONENT_METRICS["Y"].col_0based] = [0.12, 0.24]

    np.testing.assert_allclose(plotter._y_values(data, "concentration", None), [12.0, 24.0])
    np.testing.assert_allclose(plotter._y_values(data, "molefraction", None), [0.12, 0.24])


def test_component_metric_indices_follow_component_data_header():
    assert COMPONENT_METRICS["C"].col_1based == 4
    assert COMPONENT_METRICS["Cdot"].col_1based == 5
    assert COMPONENT_METRICS["Y"].col_1based == 6
    assert COMPONENT_METRICS["Q"].col_1based == 7
    assert COMPONENT_METRICS["Dqdt"].col_1based == 8
    assert COMPONENT_METRICS["Qchem"].col_1based == 9
    assert COMPONENT_METRICS["Dqchemdt"].col_1based == 10
    assert COMPONENT_METRICS["P"].col_1based == 11
    assert COMPONENT_METRICS["Qeq"].col_1based == 12
    assert COMPONENT_METRICS["Pnorm"].col_1based == 13
    assert COMPONENT_METRICS["Qchemeq"].col_1based == 14


def test_explorer_y_options_use_full_column_header_names():
    options = _column_explorer_y_options()

    assert ("Physisorption, q_phy_i [mol/kg]", "Q") in options
    assert ("Total pressure, p_t [Pa]", "Pt") in options

    labels = [label for label, _value in options]
    assert labels.count("Concentration time derivative, dc_i/dt [mol/m^3/s]") == 1


def test_breakthrough_reader_rejects_old_component_schema(tmp_path):
    component_file = tmp_path / "component_0_test.data"
    component_file.write_text("0 0 0 1 2 3 4 5 6 7\n", encoding="utf-8")
    plotter = BreakthroughPlotly(
        displayName="test",
        externalTemperature=298.15,
        externalPressure=101325.0,
        components=[],
    )

    with pytest.raises(ValueError, match="expected the current 14-column schema"):
        plotter._read_component_blocks(component_file)
