from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from .utils import ComponentInfo


class BasePlotly:
    """
    Shared helpers for Plotly-based simulation visualizers.
    """

    def __init__(
        self,
        displayName: str,
        components: List[ComponentInfo],
        data_dir: Union[str, Path] = ".",
        carrierGasComponent: Optional[int] = None,
    ) -> None:
        self.displayName = displayName
        self.components = components
        self.Ncomp = len(components)
        self.data_dir = Path(data_dir)
        self.carrierGasComponent = carrierGasComponent

        if self.carrierGasComponent is not None:
            for comp in self.components:
                comp.isCarrierGas = (comp.index == self.carrierGasComponent)

    def _component_file_name(self, i: int, name: str) -> Path:
        return self.data_dir / f"component_{i}_{name}.data"

    def _component_sort_key(self, comp: ComponentInfo):
        return (0 if comp.isCarrierGas else 1, comp.index)

    def _component_label(self, comp: ComponentInfo, include_yi: bool = True) -> str:
        if include_yi and comp.Yi0 is not None:
            return f"{comp.name} (y_i={comp.Yi0:g})"
        return comp.name

    def _read_component_data(
        self,
        fileName: Union[str, Path],
        min_columns: int,
    ) -> pd.DataFrame:
        df = pd.read_csv(
            fileName,
            sep=r"\s+",
            comment="#",
            header=None,
            engine="python",
        )
        if df.shape[1] < min_columns:
            raise ValueError(f"{fileName} must contain at least {min_columns} columns")
        df.columns = [f"col_{i + 1}" for i in range(df.shape[1])]
        return df

    def _publication_layout(
            self,
            fig: go.Figure,
            xaxis_title: str,
            yaxis_title: str,
            title_text: str,
            x_range: Optional[List[float]] = None,
            y_range: Optional[List[float]] = None,
            width: int = 640,
            height: int = 480,
            margin_left: int = 95,
            margin_right: int = 30,
            margin_top: int = 40,
            margin_bottom: int = 110,
            show_legend: bool = True,
        ) -> go.Figure:
            fig.update_layout(
                title=dict(
                    text=title_text,
                    x=0.5,
                    xanchor="center",
                    y=0.96,
                    yanchor="top",
                    font=dict(size=20),
                ),
                width=width,
                height=height,
                template="simple_white",
                paper_bgcolor="white",
                plot_bgcolor="white",
                margin=dict(
                    l=margin_left,
                    r=margin_right,
                    t=margin_top,
                    b=margin_bottom,
                ),
                font=dict(size=16),
                showlegend=show_legend,
                legend=dict(
                    orientation="h",
                    yanchor="top",
                    y=-0.3,
                    xanchor="center",
                    x=0.5,
                    bgcolor="rgba(255,255,255,0.95)",
                    bordercolor="black",
                    borderwidth=0.8,
                    font=dict(size=13),
                    title=dict(text="Components"),
                ),
            )

            fig.update_xaxes(
                title=xaxis_title,
                showline=True,
                linewidth=1.4,
                linecolor="black",
                mirror=True,
                ticks="outside",
                tickwidth=1.2,
                ticklen=6,
                showgrid=False,
                zeroline=False,
            )
            fig.update_yaxes(
                title=yaxis_title,
                showline=True,
                linewidth=1.4,
                linecolor="black",
                mirror=True,
                ticks="outside",
                tickwidth=1.2,
                ticklen=6,
                showgrid=False,
                zeroline=False,
            )

            if x_range is not None:
                fig.update_xaxes(range=x_range)
            if y_range is not None:
                fig.update_yaxes(range=y_range)

            return fig

    @staticmethod
    def _safe_minmax(values: np.ndarray) -> Optional[tuple[float, float]]:
        if not np.isfinite(values).any():
            return None
        return float(np.nanmin(values)), float(np.nanmax(values))