from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .base_plotly import BasePlotly
from .utils import ComponentInfo


class FittingPlotly(BasePlotly):
    """
    Plotly interface for the C++ fitting workflow.

    Mirrors the intent of:
      - Fitting::createPlotScripts(...)
      - Fitting::createPlotScript()

    For each component it plots:
      - the starting isotherm
      - the optimized fit
      - the raw data used in fitting
    """

    def __init__(
        self,
        displayName: str,
        components: List[ComponentInfo],
        start_isotherms: Sequence[object],
        fit_isotherms: Sequence[object],
        filenames: Sequence[Union[str, Path]],
        data_dir: Union[str, Path] = ".",
        pressureScale: str = "log",
        columnPressure: int = 1,
        columnLoading: int = 2,
        columnError: Optional[int] = None,
        carrierGasComponent: Optional[int] = None,
    ) -> None:
        super().__init__(
            displayName=displayName,
            components=components,
            data_dir=data_dir,
            carrierGasComponent=carrierGasComponent,
        )

        if len(start_isotherms) != self.Ncomp:
            raise ValueError("start_isotherms length must match number of components")
        if len(fit_isotherms) != self.Ncomp:
            raise ValueError("fit_isotherms length must match number of components")
        if len(filenames) != self.Ncomp:
            raise ValueError("filenames length must match number of components")

        self.start_isotherms = list(start_isotherms)
        self.fit_isotherms = list(fit_isotherms)
        self.filenames = [self.data_dir / Path(fileName) for fileName in filenames]
        self.pressureScale = str(pressureScale).lower()

        # Match C++ InputReader semantics: incoming columns are 1-based
        self.columnPressure = int(columnPressure) - 1
        self.columnLoading = int(columnLoading) - 1
        self.columnError = None if columnError is None else int(columnError) - 1

    def _plot_title(self, mode: str = "Fitted isotherm") -> str:
        return f"{self.displayName} — {mode}"

    def _read_fit_data(self, fileName: Union[str, Path]) -> pd.DataFrame:
        df = pd.read_csv(
            fileName,
            sep=r"\s+",
            comment="#",
            header=None,
            engine="python",
        )

        required_cols = max(self.columnPressure, self.columnLoading)
        if self.columnError is not None:
            required_cols = max(required_cols, self.columnError)

        if df.shape[1] <= required_cols:
            raise ValueError(f"{fileName} does not contain the requested fitting columns")

        out = pd.DataFrame(
            {
                "pressure": pd.to_numeric(df.iloc[:, self.columnPressure], errors="coerce"),
                "loading": pd.to_numeric(df.iloc[:, self.columnLoading], errors="coerce"),
            }
        )

        if self.columnError is not None:
            out["error"] = pd.to_numeric(df.iloc[:, self.columnError], errors="coerce")

        out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["pressure", "loading"])
        out = out.sort_values("pressure", kind="mergesort").reset_index(drop=True)

        if out.empty:
            raise ValueError(f"No valid pressure/loading points found in {fileName}")

        return out

    @staticmethod
    def _evaluate_isotherm(isotherm: object, pressure: np.ndarray) -> np.ndarray:
        if hasattr(isotherm, "value") and callable(isotherm.value):
            return np.asarray([isotherm.value(float(p)) for p in pressure], dtype=float)

        if callable(isotherm):
            return np.asarray([isotherm(float(p)) for p in pressure], dtype=float)

        raise TypeError("isotherm must either be callable or expose a value(pressure) method")

    def _pressure_grid(self, pressure_data: np.ndarray, npoints: int = 400) -> np.ndarray:
        pmin = float(np.nanmin(pressure_data))
        pmax = float(np.nanmax(pressure_data))

        if pmin <= 0.0 and self.pressureScale == "log":
            raise ValueError("Log-scale fitting plots require strictly positive pressures")

        if np.isclose(pmin, pmax):
            return np.asarray([pmin], dtype=float)

        if self.pressureScale == "log":
            return np.geomspace(pmin, pmax, npoints)

        return np.linspace(pmin, pmax, npoints)

    def fit_figure(self, comp_index: int) -> go.Figure:
        if comp_index < 0 or comp_index >= self.Ncomp:
            raise IndexError(f"comp_index {comp_index} out of range [0, {self.Ncomp - 1}]")

        comp = self.components[comp_index]
        df = self._read_fit_data(self.filenames[comp_index])

        pressure_data = df["pressure"].to_numpy(dtype=float)
        loading_data = df["loading"].to_numpy(dtype=float)
        pressure_grid = self._pressure_grid(pressure_data)

        start_loading = self._evaluate_isotherm(self.start_isotherms[comp_index], pressure_grid)
        fit_loading = self._evaluate_isotherm(self.fit_isotherms[comp_index], pressure_grid)

        ymax = 0.0
        for values in (loading_data, start_loading, fit_loading):
            if np.isfinite(values).any():
                ymax = max(ymax, float(np.nanmax(values)))

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=pressure_grid,
                y=start_loading,
                mode="lines",
                name="start f(x)",
                line=dict(width=2.4, dash="dash"),
                hovertemplate=("start f(x)<br>" + "P=%{x:.4g} Pa<br>" + "q=%{y:.6g} mol/kg<extra></extra>"),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=pressure_grid,
                y=fit_loading,
                mode="lines",
                name="fit f(x)",
                line=dict(width=3.0),
                hovertemplate=("fit f(x)<br>" + "P=%{x:.4g} Pa<br>" + "q=%{y:.6g} mol/kg<extra></extra>"),
            )
        )

        scatter_kwargs = {}
        if "error" in df.columns:
            error_values = df["error"].to_numpy(dtype=float)
            if np.isfinite(error_values).any():
                scatter_kwargs["error_y"] = dict(
                    type="data",
                    array=error_values,
                    visible=True,
                )

        fig.add_trace(
            go.Scatter(
                x=pressure_data,
                y=loading_data,
                mode="markers",
                name="raw data",
                marker=dict(size=7),
                hovertemplate=("raw data<br>" + "P=%{x:.4g} Pa<br>" + "q=%{y:.6g} mol/kg<extra></extra>"),
                **scatter_kwargs,
            )
        )

        y_upper = 1.08 * ymax if ymax > 0.0 else 1.0

        fig = self._publication_layout(
            fig,
            xaxis_title="Pressure, p [Pa]",
            yaxis_title="Absolute loading, q [mol/kg]",
            title_text=f"{self._plot_title('Fitted isotherm')} — {comp.name}",
            y_range=[0.0, y_upper],
            width=900,
            height=540,
            margin_left=105,
            margin_right=35,
            margin_top=92,
            margin_bottom=88,
            show_legend=True,
        )

        fig.update_layout(
            legend=dict(
                orientation="v",
                yanchor="bottom",
                y=0.02,
                xanchor="right",
                x=0.98,
                title=dict(text=comp.name),
            )
        )

        if self.pressureScale == "log":
            fig.update_xaxes(
                type="log",
                exponentformat="power",
                minor=dict(ticks="", showgrid=False),
                minorloglabels="none",
            )

        return fig

    def all_fit_figures(self) -> dict[str, go.Figure]:
        return {
            comp.name: self.fit_figure(comp.index) for comp in sorted(self.components, key=self._component_sort_key)
        }
