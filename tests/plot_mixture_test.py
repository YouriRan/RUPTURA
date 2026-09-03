import pytest

from ruptura.plot_mixture import MixturePredictionPlotly
from ruptura.utils import ComponentInfo


@pytest.fixture
def plotter(tmp_path):
    (tmp_path / "component_0_test.data").write_text(
        "1 2 3 0.5 0.25\n10 4 5 0.5 0.75\n",
        encoding="utf-8",
    )
    return MixturePredictionPlotly(
        displayName="test",
        temperature=298.15,
        components=[ComponentInfo(index=0, name="test", initialGasMoleFraction=0.5)],
        data_dir=tmp_path,
    )


@pytest.mark.parametrize(
    "plot_method",
    ["pure_components", "mixture_loading", "mixture_adsorbed_molfractions"],
)
def test_mixture_figures_can_hide_markers(plotter, plot_method):
    figure = getattr(plotter, plot_method)(show_markers=False)

    assert figure.data[0].mode == "lines"
    assert figure.data[0].marker.to_plotly_json() == {}


def test_mixture_figures_show_markers_by_default(plotter):
    figure = plotter.mixture_loading()

    assert figure.data[0].mode == "lines+markers"
    assert figure.data[0].marker.size == 8


@pytest.fixture
def chemisorption_plotter(tmp_path):
    """A MixturePrediction output carrying the chemisorption columns 8-11."""
    (tmp_path / "component_0_test.data").write_text(
        "# column 8: pure component chemisorption isotherm value\n"
        "1 2 3 0.5 0.25 1 0.9 0.5 0.4 2.5 3.4\n"
        "10 4 5 0.5 0.75 10 1.8 1.5 1.4 5.5 6.4\n",
        encoding="utf-8",
    )
    return MixturePredictionPlotly(
        displayName="test",
        temperature=298.15,
        components=[ComponentInfo(index=0, name="test", initialGasMoleFraction=0.5)],
        data_dir=tmp_path,
    )


def test_physisorption_only_output_has_no_chemisorption(plotter):
    assert plotter.has_chemisorption() is False


def test_chemisorption_output_is_detected(chemisorption_plotter):
    assert chemisorption_plotter.has_chemisorption() is True


@pytest.mark.parametrize(
    "plot_method, expected",
    [
        ("pure_components_chemisorption", [0.5, 1.5]),
        ("mixture_chemisorption_loading", [0.4, 1.4]),
        ("pure_components_total_loading", [2.5, 5.5]),
        ("mixture_total_loading", [3.4, 6.4]),
    ],
)
def test_chemisorption_figures_read_the_right_column(chemisorption_plotter, plot_method, expected):
    figure = getattr(chemisorption_plotter, plot_method)()

    assert list(figure.data[0].y) == expected


@pytest.mark.parametrize(
    "plot_method",
    [
        "pure_components_chemisorption",
        "mixture_chemisorption_loading",
        "pure_components_total_loading",
        "mixture_total_loading",
    ],
)
def test_chemisorption_figures_reject_physisorption_only_output(plotter, plot_method):
    with pytest.raises(ValueError, match="no chemisorption columns"):
        getattr(plotter, plot_method)()
