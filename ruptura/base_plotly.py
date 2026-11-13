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
        showMarkers: bool = False,
    ) -> None:
        self.displayName = displayName
        self.components = components
        self.numberOfComponents = len(components)
        self.data_dir = Path(data_dir)
        self.carrierGasComponent = carrierGasComponent
        self.showMarkers = showMarkers

        if self.carrierGasComponent is not None:
            for comp in self.components:
                comp.isCarrierGas = comp.index == self.carrierGasComponent

    def _component_file_name(self, i: int, name: str) -> Path:
        return self.data_dir / f"component_{i}_{name}.data"

    def _component_sort_key(self, comp: ComponentInfo):
        return (0 if comp.isCarrierGas else 1, comp.index)

    def _component_label(self, comp: ComponentInfo, include_yi: bool = True) -> str:
        if include_yi and comp.initialGasMoleFraction is not None:
            return f"{comp.name} (y_i={comp.initialGasMoleFraction:g})"
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

    def _trace_mode(self, with_lines: bool = True) -> str:
        if with_lines:
            return "lines+markers" if self.showMarkers else "lines"
        return "markers" if self.showMarkers else "lines"

    def _trace_kwargs(self, with_lines: bool = True) -> dict:
        kwargs = {"mode": self._trace_mode(with_lines=with_lines)}
        if self.showMarkers:
            kwargs["marker"] = {"size": 7}
        return kwargs

    def _publication_layout(
        self,
        fig: go.Figure,
        xaxis_title: str,
        yaxis_title: str,
        title_text: str,
        x_range: Optional[List[float]] = None,
        y_range: Optional[List[float]] = None,
        width: int = 800,
        height: int = 600,
        margin_left: int = 95,
        margin_right: int = 30,
        margin_top: int = 50,
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
                y=-0.2,
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

    def save_figure(
        self,
        fig: go.Figure,
        fileName: Union[str, Path],
        format: Optional[str] = None,
        width: Optional[int] = None,
        height: Optional[int] = None,
        scale: float = 1.0,
    ) -> Path:
        """
        Save a Plotly figure to disk.

        Static image export requires kaleido:
            pip install kaleido
        """
        path = Path(fileName)
        fmt = (format or path.suffix.lstrip(".") or "pdf").lower()

        if path.suffix == "":
            path = path.with_suffix(f".{fmt}")

        path.parent.mkdir(parents=True, exist_ok=True)

        fig.write_image(
            str(path),
            format=fmt,
            width=width,
            height=height,
            scale=scale,
        )
        return path

    @staticmethod
    def _safe_minmax(values: np.ndarray) -> Optional[tuple[float, float]]:
        if not np.isfinite(values).any():
            return None
        return float(np.nanmin(values)), float(np.nanmax(values))