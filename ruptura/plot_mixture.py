from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import plotly.graph_objects as go

from .base_plotly import BasePlotly
from .utils import ComponentInfo, load_simulation_metadata


class MixturePredictionPlotly(BasePlotly):
    """
    Publication-oriented Plotly interface for MixturePrediction output files.
    """

    def __init__(
        self,
        displayName: str,
        temperature: float,
        components: List[ComponentInfo],
        data_dir: Union[str, Path] = ".",
        pressureScale: str = "log",
        pressureStart: Optional[float] = None,
        pressureEnd: Optional[float] = None,
        carrierGasComponent: Optional[int] = None,
    ) -> None:
        super().__init__(
            displayName=displayName,
            components=components,
            data_dir=data_dir,
            carrierGasComponent=carrierGasComponent,
        )
        self.temperature = temperature
        self.pressureScale = pressureScale.lower()
        self.pressureStart = pressureStart
        self.pressureEnd = pressureEnd

    @classmethod
    def from_simulation_json(
        cls,
        simulation_json: Union[str, Path],
        data_dir: Union[str, Path] = ".",
        carrierGasComponent: Optional[int] = None,
    ) -> "MixturePredictionPlotly":
        meta = load_simulation_metadata(simulation_json)
        return cls(
            displayName=str(meta["displayName"]),
            temperature=float(meta["temperature"]),
            components=meta["components"],
            data_dir=data_dir,
            pressureScale=str(meta["pressureScale"]),
            pressureStart=meta["pressureStart"],
            pressureEnd=meta["pressureEnd"],
            carrierGasComponent=carrierGasComponent,
        )

    def _read_component_data(self, fileName: Union[str, Path]):
        return super()._read_component_data(fileName, min_columns=5)

    def _plot_title(self, mode: str = "Mixture prediction") -> str:
        return f"{self.displayName} — {mode}  T={self.temperature:g} K"

    def _publication_layout(
        self,
        fig: go.Figure,
        xaxis_title: str,
        yaxis_title: str,
        title_text: str,
        y_range: Optional[List[float]] = None,
        show_legend: bool = True,
    ) -> go.Figure:
        fig = super()._publication_layout(
            fig,
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title,
            title_text=title_text,
            y_range=y_range,
            width=640,
            height=480,
            margin_left=105,
            margin_right=35,
            margin_top=60,
            margin_bottom=88,
            show_legend=show_legend,
        )

        if self.pressureScale == "log":
            if self.pressureStart is not None and self.pressureEnd is not None:
                fig.update_xaxes(
                    type="log",
                    tickmode="linear",
                    tick0=np.floor(np.log10(self.pressureStart)) if self.pressureStart is not None else 0,
                    dtick=1,
                    exponentformat="power",
                    minor=dict(ticks="", showgrid=False),
                    minorloglabels="none",
                )
        else:
            if self.pressureStart is not None and self.pressureEnd is not None:
                fig.update_xaxes(range=[self.pressureStart, self.pressureEnd])

        return fig

    def _build_figure(
        self,
        y_col: str,
        yaxis_title: str,
        include_yi: bool,
        include_carrier_gas: bool,
        title_mode: str,
        force_zero_ymin: bool = True,
    ) -> go.Figure:
        fig = go.Figure()
        ymax = -np.inf
        ymin = np.inf

        for comp in sorted(self.components, key=self._component_sort_key):
            if (not include_carrier_gas) and comp.isCarrierGas:
                continue

            fileName = self._component_file_name(comp.index, comp.name)
            df = self._read_component_data(fileName)

            x = df["col_1"].to_numpy(dtype=float)
            y = df[y_col].to_numpy(dtype=float)

            if np.isfinite(y).any():
                ymax = max(ymax, float(np.nanmax(y)))
                ymin = min(ymin, float(np.nanmin(y)))

            line_width = 3.6 if comp.isCarrierGas else 3.2
            label = self._component_label(comp, include_yi=include_yi)

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=label,
                    line=dict(width=line_width),
                    hovertemplate=(
                        f"{label}<br>"
                        + "P=%{x:.4g} Pa<br>"
                        + "y=%{y:.6g}<extra></extra>"
                    ),
                )
            )

        if not np.isfinite(ymax):
            ymax = 1.0
            ymin = 0.0

        if force_zero_ymin:
            y_lower = 0.0
            y_upper = 1.08 * ymax if ymax > 0.0 else 1.0
        else:
            if ymin >= 0.0:
                y_lower = 0.0
                y_upper = 1.08 * ymax if ymax > 0.0 else 1.0
            else:
                y_lower = 1.08 * ymin
                y_upper = 1.08 * ymax if ymax != 0.0 else 1.0

        fig = self._publication_layout(
            fig,
            xaxis_title="Total bulk fluid phase fugacity, f [Pa]",
            yaxis_title=yaxis_title,
            title_text=self._plot_title(title_mode),
            y_range=[y_lower, y_upper],
        )
        return fig

    def pure_components(self, include_carrier_gas: bool = True) -> go.Figure:
        return self._build_figure(
            y_col="col_2",
            yaxis_title="Absolute loading, q_i",
            include_yi=False,
            include_carrier_gas=include_carrier_gas,
            title_mode="Pure-component isotherms",
            force_zero_ymin=True,
        )

    def mixture_loading(self, include_carrier_gas: bool = True) -> go.Figure:
        return self._build_figure(
            y_col="col_3",
            yaxis_title="Absolute loading, q_i",
            include_yi=True,
            include_carrier_gas=include_carrier_gas,
            title_mode="Mixture prediction",
            force_zero_ymin=True,
        )

    def mixture_adsorbed_molfractions(self, include_carrier_gas: bool = True) -> go.Figure:
        fig = self._build_figure(
            y_col="col_5",
            yaxis_title="Adsorbed mol-fraction, x_i [-]",
            include_yi=True,
            include_carrier_gas=include_carrier_gas,
            title_mode="Mixture adsorbed mol-fractions",
            force_zero_ymin=True,
        )
        fig.update_yaxes(range=[0.0, 1.02])
        return fig