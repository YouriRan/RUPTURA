import sys
import types

import numpy as np

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
    _canonicalize_component_blocks,
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


def test_legacy_multibed_component_blocks_are_upgraded_to_canonical_schema():
    legacy = np.asarray(
        [
            [0.0, 0.0, 0.0, 2.0, 0.1, 1.2, 0.2, 25.0, 1.3, 0.5],
            [0.0, 0.0, 1.0, 4.0, 0.2, 1.4, 0.3, 50.0, 1.5, 1.0],
        ]
    )
    column = np.zeros((2, 12), dtype=float)
    column[:, :3] = legacy[:, :3]
    column[:, 4] = 100.0

    [canonical] = _canonicalize_component_blocks([legacy], [column])

    assert canonical.shape == (2, 14)
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Y"].col_0based], [0.25, 0.5])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Q"].col_0based], legacy[:, 5])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Dqdt"].col_0based], legacy[:, 6])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["P"].col_0based], legacy[:, 7])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Qeq"].col_0based], legacy[:, 8])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Pnorm"].col_0based], legacy[:, 9])
    np.testing.assert_allclose(canonical[:, COMPONENT_METRICS["Qchem"].col_0based], 0.0)
