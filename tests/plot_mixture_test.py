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
