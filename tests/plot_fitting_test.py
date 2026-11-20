import numpy as np

from ruptura.plot_fitting import FittingPlotly, JsonIsothermSite, JsonMultiSiteIsotherm


def test_start_isotherm_with_zero_parameters_evaluates_to_nan():
    # A fit from scratch starts at all-zero parameters; the start curve must not stop the figure.
    toth = JsonMultiSiteIsotherm((JsonIsothermSite("Toth", (0.0, 0.0, 0.0)),))
    assert np.isnan(FittingPlotly._evaluate_isotherm(toth, np.array([1.0e3, 1.0e4]))).all()
