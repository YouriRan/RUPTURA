from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import plotly.graph_objects as go

try:
    import ipywidgets as widgets
    from IPython.display import display
except Exception:  # pragma: no cover
    widgets = None
    display = None

from .base_plotly import BasePlotly
from .utils import (
    ComponentInfo,
    HeaderComponentInfo,
    load_simulation_metadata,
    parse_header_metadata,
    read_blocks,
)


class BreakthroughPlotly(BasePlotly):
    """
    Notebook-native Plotly interface for breakthrough and column-profile data.

    This class:
      - reads the existing simulation output files
      - constructs Plotly figures in memory
      - does not write HTML or other output files
    """

    def __init__(
        self,
        displayName: str,
        externalTemperature: float,
        externalPressure: float,
        components: List[ComponentInfo],
        data_dir: Union[str, Path] = ".",
        carrierGasComponent: Optional[int] = None,
    ) -> None:
        super().__init__(
            displayName=displayName,
            components=components,
            data_dir=data_dir,
            carrierGasComponent=carrierGasComponent,
        )
        self.externalTemperature = externalTemperature
        self.externalPressure = externalPressure

    @classmethod
    def from_simulation_json(
        cls,
        simulation_json: Union[str, Path],
        data_dir: Union[str, Path] = ".",
        carrierGasComponent: Optional[int] = None,
    ) -> "BreakthroughPlotly":
        meta = load_simulation_metadata(simulation_json)
        return cls(
            displayName=str(meta["displayName"]),
            externalTemperature=float(meta["temperature"]),
            externalPressure=float(meta["pressureEnd"]) if meta["pressureEnd"] is not None else 0.0,
            components=meta["components"],
            data_dir=data_dir,
            carrierGasComponent=carrierGasComponent,
        )

    def _read_component_data(self, fileName: Union[str, Path]):
        return super()._read_component_data(fileName, min_columns=3)

    def _read_column_data(self, fileName: Union[str, Path] = "column.data") -> List[np.ndarray]:
        return read_blocks(self.data_dir / fileName)

    def _plot_title(self, mode: str = "Breakthrough") -> str:
        return (
            f"{self.displayName} — {mode}  "
            f"T={self.externalTemperature:g} K, "
            f"p_t={self.externalPressure * 1e-3:g} kPa"
        )

    def _metric_info(self, metric: str) -> Dict[str, object]:
        metric_map: Dict[str, Dict[str, object]] = {
            "V": {"ylabel": "Interstitial velocity, v [m/s]", "mode": "single", "base_col_1based": 2},
            "Pt": {"ylabel": "Total pressure, p_t [Pa]", "mode": "single", "base_col_1based": 3},
            "Q": {"ylabel": "Concentration, c_i [mol/kg]", "mode": "component", "base_col_1based": 4},
            "Qeq": {"ylabel": "Concentration, c_i [mol/kg]", "mode": "component", "base_col_1based": 5},
            "P": {"ylabel": "Partial pressure, p_i [Pa]", "mode": "component", "base_col_1based": 6},
            "Pnorm": {"ylabel": "Partial pressure, p_i [-]", "mode": "component", "base_col_1based": 7},
            "Dpdt": {"ylabel": "Pressure derivative, dp_i/dt [Pa/s]", "mode": "component", "base_col_1based": 8},
            "Dqdt": {"ylabel": "Loading derivative, dq_i/dt [mol/kg/s]", "mode": "component", "base_col_1based": 9},
        }

        if metric not in metric_map:
            raise ValueError(f"Unknown metric '{metric}'")

        return metric_map[metric]

    def _column_index_0based(self, metric: str, comp_index: Optional[int] = None) -> int:
        info = self._metric_info(metric)
        base_col_1based = int(info["base_col_1based"])
        mode = str(info["mode"])

        if mode == "single":
            return base_col_1based - 1

        if comp_index is None:
            raise ValueError(f"metric '{metric}' requires comp_index")

        return (base_col_1based - 1) + 6 * comp_index

    def _breakthrough_figure(
        self,
        x_col: str,
        xaxis_title: str,
        include_carrier_gas: bool = True,
    ) -> go.Figure:
        fig = go.Figure()
        ymax = 0.0
        xmin = None
        xmax = None

        for comp in sorted(self.components, key=self._component_sort_key):
            if (not include_carrier_gas) and comp.isCarrierGas:
                continue

            fileName = self._component_file_name(comp.index, comp.name)
            df = self._read_component_data(fileName)

            x = df[x_col].to_numpy(dtype=float)
            y = df["col_3"].to_numpy(dtype=float)

            if np.isfinite(x).any():
                x_local_min = float(np.nanmin(x))
                x_local_max = float(np.nanmax(x))
                xmin = x_local_min if xmin is None else min(xmin, x_local_min)
                xmax = x_local_max if xmax is None else max(xmax, x_local_max)

            if np.isfinite(y).any():
                ymax = max(ymax, float(np.nanmax(y)))

            line_width = 3.6 if comp.isCarrierGas else 3.2
            label = self._component_label(comp, include_yi=True)

            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode="lines",
                    name=label,
                    line=dict(width=line_width),
                    hovertemplate=(f"{label}<br>" + "x=%{x:.3f}<br>" + "y=%{y:.3f}<extra></extra>"),
                )
            )

        if xmax is None:
            xmin = 0.0
            xmax = 1.0

        y_upper = 1.05 if ymax <= 1.02 else 1.08 * ymax
        fig = self._publication_layout(
            fig,
            xaxis_title=xaxis_title,
            yaxis_title="Concentration exit gas, c<sub>i</sub>/c<sub>i,0</sub> [-]",
            title_text=self._plot_title("Breakthrough"),
            x_range=[xmin, xmax],
            y_range=[0.0, y_upper],
        )
        return fig

    def breakthrough_dimensionless(self, include_carrier_gas: bool = True) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_1",
            xaxis_title="Dimensionless time, τ = tv/L [-]",
            include_carrier_gas=include_carrier_gas,
        )

    def breakthrough_time(self, include_carrier_gas: bool = True) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_2",
            xaxis_title="Time, t [min]",
            include_carrier_gas=include_carrier_gas,
        )

    def column_snapshot(self, metric: str, frame_index: int = 0) -> go.Figure:
        blocks = self._read_column_data("column.data")
        info = self._metric_info(metric)
        ylabel = str(info["ylabel"])
        mode = str(info["mode"])

        if frame_index < 0 or frame_index >= len(blocks):
            raise IndexError(f"frame_index {frame_index} out of range [0, {len(blocks) - 1}]")

        block = blocks[frame_index]
        fig = go.Figure()

        if mode == "single":
            col = self._column_index_0based(metric)
            fig.add_trace(
                go.Scatter(
                    x=block[:, 0],
                    y=block[:, col],
                    mode="lines+markers",
                    name=metric,
                )
            )
        else:
            for comp_index, comp in enumerate(self.components):
                col = self._column_index_0based(metric, comp_index)
                fig.add_trace(
                    go.Scatter(
                        x=block[:, 0],
                        y=block[:, col],
                        mode="lines+markers",
                        name=self._component_label(comp, include_yi=True),
                    )
                )

        fig.update_layout(
            title=f"{self._plot_title()} — {metric} — frame {frame_index}",
            xaxis_title="Adsorber position [m]",
            yaxis_title=ylabel,
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.18,
                xanchor="center",
                x=0.5,
            ),
            template="plotly_white",
            margin=dict(b=110),
        )
        return fig

    def column_animation(self, metric: str, every: int = 1) -> go.Figure:
        blocks = self._read_column_data("column.data")
        info = self._metric_info(metric)
        ylabel = str(info["ylabel"])
        mode = str(info["mode"])

        sampled_indices = list(range(0, len(blocks), max(1, every)))
        sampled_blocks = [blocks[i] for i in sampled_indices]

        x_min = min(block[:, 0].min() for block in sampled_blocks)
        x_max = max(block[:, 0].max() for block in sampled_blocks)

        if mode == "single":
            col = self._column_index_0based(metric)
            y_min = min(block[:, col].min() for block in sampled_blocks)
            y_max = max(block[:, col].max() for block in sampled_blocks)
        else:
            y_arrays = []
            for comp_index in range(self.Ncomp):
                col = self._column_index_0based(metric, comp_index)
                for block in sampled_blocks:
                    y_arrays.append(block[:, col])
            y_all = np.concatenate(y_arrays)
            y_min = float(np.min(y_all))
            y_max = float(np.max(y_all))

        if y_min >= 0.0:
            y_min = 0.0
            y_max = 1.1 * y_max if y_max > 0.0 else 1.0
        else:
            y_min = 1.1 * y_min
            y_max = 1.1 * y_max

        first_block = sampled_blocks[0]
        fig = go.Figure()

        if mode == "single":
            col = self._column_index_0based(metric)
            fig.add_trace(
                go.Scatter(
                    x=first_block[:, 0],
                    y=first_block[:, col],
                    mode="lines+markers",
                    name=metric,
                )
            )
        else:
            for comp_index, comp in enumerate(self.components):
                col = self._column_index_0based(metric, comp_index)
                fig.add_trace(
                    go.Scatter(
                        x=first_block[:, 0],
                        y=first_block[:, col],
                        mode="lines+markers",
                        name=self._component_label(comp, include_yi=True),
                    )
                )

        frames = []
        for idx, block in zip(sampled_indices, sampled_blocks):
            if mode == "single":
                col = self._column_index_0based(metric)
                data = [
                    go.Scatter(
                        x=block[:, 0],
                        y=block[:, col],
                        mode="lines+markers",
                        name=metric,
                    )
                ]
            else:
                data = []
                for comp_index, comp in enumerate(self.components):
                    col = self._column_index_0based(metric, comp_index)
                    data.append(
                        go.Scatter(
                            x=block[:, 0],
                            y=block[:, col],
                            mode="lines+markers",
                            name=self._component_label(comp, include_yi=True),
                        )
                    )

            frames.append(
                go.Frame(
                    name=str(idx),
                    data=data,
                    layout=go.Layout(title_text=f"{self._plot_title()} — {metric} — frame {idx}"),
                )
            )

        fig.frames = frames

        fig.update_layout(
            title=f"{self._plot_title()} — {metric} — frame {sampled_indices[0]}",
            xaxis_title="Adsorber position [m]",
            yaxis_title=ylabel,
            xaxis=dict(range=[x_min, x_max]),
            yaxis=dict(range=[y_min, y_max]),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.18,
                xanchor="center",
                x=0.5,
            ),
            template="plotly_white",
            margin=dict(b=110),
            updatemenus=[
                {
                    "type": "buttons",
                    "showactive": True,
                    "buttons": [
                        {
                            "label": "Play",
                            "method": "animate",
                            "args": [
                                None,
                                {
                                    "frame": {"duration": 80, "redraw": True},
                                    "transition": {"duration": 0},
                                    "fromcurrent": True,
                                },
                            ],
                        },
                        {
                            "label": "Pause",
                            "method": "animate",
                            "args": [
                                [None],
                                {
                                    "frame": {"duration": 0, "redraw": False},
                                    "transition": {"duration": 0},
                                    "mode": "immediate",
                                },
                            ],
                        },
                    ],
                }
            ],
            sliders=[
                {
                    "steps": [
                        {
                            "label": str(idx),
                            "method": "animate",
                            "args": [
                                [str(idx)],
                                {
                                    "frame": {"duration": 0, "redraw": True},
                                    "transition": {"duration": 0},
                                    "mode": "immediate",
                                },
                            ],
                        }
                        for idx in sampled_indices
                    ]
                }
            ],
        )

        return fig

    def all_column_figures(self, every: int = 1) -> Dict[str, go.Figure]:
        return {
            metric: self.column_animation(metric, every=every)
            for metric in ["V", "Pt", "Q", "Qeq", "P", "Pnorm", "Dpdt", "Dqdt"]
        }

    def explorer(
        self,
        datafile: Union[str, Path] = "column.data",
        dt: Optional[float] = None,
        time0: float = 0.0,
        component_names: Optional[Dict[int, str]] = None,
        legend_title: str = "",
    ) -> "ColumnDataExplorer":
        if component_names is None:
            component_names = {comp.index: self._component_label(comp, include_yi=True) for comp in self.components}

        if not legend_title:
            legend_title = self._plot_title()

        dt_use = self.dt_guess(dt)
        return ColumnDataExplorer(
            self.data_dir / datafile,
            dt=dt_use,
            time0=time0,
            component_names=component_names,
            legend_title=legend_title,
        )

    def dt_guess(self, dt: Optional[float]) -> Optional[float]:
        return dt


class ColumnDataExplorer:
    """
    Interactive explorer for column.data-like files (plotly + ipywidgets).
    """

    def __init__(
        self,
        datafile: Union[str, Path],
        dt: Optional[float] = None,
        time0: float = 0.0,
        component_names: Optional[Dict[int, str]] = None,
        legend_title: str = "",
    ):
        if widgets is None:
            raise ImportError("ipywidgets is required for ColumnDataExplorer (pip install ipywidgets).")

        self.datafile = Path(datafile)
        self.meta = parse_header_metadata(self.datafile)
        self.blocks = read_blocks(self.datafile)

        if len(self.blocks) == 0:
            raise ValueError(f"No blocks found in {self.datafile}")

        if component_names:
            for ci, nm in component_names.items():
                if ci in self.meta["components"]:
                    self.meta["components"][ci].name = nm

        self.dt = dt
        self.time0 = time0
        self.legend_title = legend_title

        self.z_grid = self.blocks[0][:, self.meta["z_col"]].astype(float)
        self.zmin = float(np.nanmin(self.z_grid))
        self.zmax = float(np.nanmax(self.z_grid))

        prop_set = set()
        for comp in self.meta["components"].values():
            prop_set.update(comp.props.keys())
        self.available_props = sorted(prop_set)

        self._build_widgets()
        self._build_figure()
        self._render()

    def _build_widgets(self):
        default_prop = (
            "P" if "P" in self.available_props else (self.available_props[0] if self.available_props else "P")
        )
        self.w_prop = widgets.Dropdown(options=self.available_props or ["P"], value=default_prop, description="y:")

        self.w_xmode = widgets.Dropdown(
            options=[("grid (z)", "grid"), ("time", "time")],
            value="grid",
            description="x:",
        )

        self.w_block = widgets.IntSlider(value=0, min=0, max=len(self.blocks) - 1, step=1, description="block")

        step = (self.zmax - self.zmin) / 200 if (self.zmax > self.zmin) else 1.0
        if step <= 0 or not np.isfinite(step):
            step = 1.0
        self.w_z = widgets.FloatSlider(value=self.zmin, min=self.zmin, max=self.zmax, step=step, description="z@")

        comp_options = [
            (f"{ci}: {self.meta['components'][ci].name}", ci) for ci in sorted(self.meta["components"].keys())
        ]
        self.w_comps = widgets.SelectMultiple(
            options=comp_options,
            value=tuple(ci for _, ci in comp_options),
            description="components",
            rows=min(8, max(3, len(comp_options))),
        )

        self.controls_top = widgets.HBox([self.w_xmode, self.w_prop, self.w_block])
        self.controls_bottom = widgets.HBox([self.w_z])
        self.controls = widgets.VBox([self.controls_top, self.controls_bottom, self.w_comps])

        for w in [self.w_prop, self.w_xmode, self.w_block, self.w_z, self.w_comps]:
            w.observe(self._on_change, names="value")

        self._sync_visibility()

    def _build_figure(self):
        self.fig = go.FigureWidget()
        self.fig.update_layout(
            height=520,
            margin=dict(l=70, r=20, t=70, b=110),
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.18,
                xanchor="center",
                x=0.5,
            ),
            template="plotly_white",
        )

    def _sync_visibility(self):
        if self.w_xmode.value == "grid":
            self.w_block.layout.display = ""
            self.w_z.layout.display = "none"
        else:
            self.w_block.layout.display = "none"
            self.w_z.layout.display = ""

    def _on_change(self, _change):
        self._sync_visibility()
        self._render()

    def _time_axis(self) -> np.ndarray:
        n = len(self.blocks)
        idx = np.arange(n, dtype=float)
        if self.dt is None:
            return self.time0 + idx
        return self.time0 + idx * float(self.dt)

    def _nearest_z_index(self, z_value: float) -> int:
        return int(np.nanargmin(np.abs(self.z_grid - z_value)))

    def _selected_components(self) -> List[int]:
        return list(self.w_comps.value)

    def _replace_traces(self, traces: List[go.Scatter]):
        with self.fig.batch_update():
            self.fig.data = ()
            for tr in traces:
                self.fig.add_trace(tr)

    def _render_grid(self, prop: str, block_index: int, comp_ids: List[int]):
        block = self.blocks[block_index]
        z = block[:, self.meta["z_col"]].astype(float)

        ymax = 0.0
        ymin = np.inf
        traces: List[go.Scatter] = []

        for ci in comp_ids:
            comp = self.meta["components"].get(ci)
            if not comp or prop not in comp.props:
                continue
            col = comp.props[prop]
            if col >= block.shape[1]:
                continue

            y = block[:, col].astype(float)
            if np.isfinite(y).any():
                ymax = max(ymax, float(np.nanmax(y)))
                ymin = min(ymin, float(np.nanmin(y)))

            traces.append(
                go.Scatter(
                    x=z,
                    y=y,
                    mode="lines+markers",
                    name=comp.name,
                )
            )

        self._replace_traces(traces)

        title = f"{prop} vs z (block {block_index})"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(title="Adsorber position / [m]", range=[self.zmin, self.zmax]),
                yaxis=dict(title=prop),
            )
            if np.isfinite(ymin) and np.isfinite(ymax):
                if ymin >= 0.0 and ymax > 0.0:
                    self.fig.update_yaxes(range=[0.0, 1.1 * ymax])
                elif ymin < 0.0:
                    self.fig.update_yaxes(range=[1.1 * ymin, 1.1 * ymax if ymax != 0 else 1.0])
                else:
                    self.fig.update_yaxes(autorange=True)
            else:
                self.fig.update_yaxes(autorange=True)

    def _render_time(self, prop: str, z_value: float, comp_ids: List[int]):
        t = self._time_axis()
        zi = self._nearest_z_index(z_value)
        z_used = float(self.z_grid[zi])

        traces: List[go.Scatter] = []

        for ci in comp_ids:
            comp = self.meta["components"].get(ci)
            if not comp or prop not in comp.props:
                continue
            col = comp.props[prop]

            ys = []
            for b in self.blocks:
                if zi < b.shape[0] and col < b.shape[1]:
                    ys.append(float(b[zi, col]))
                else:
                    ys.append(np.nan)
            ys = np.asarray(ys, dtype=float)

            traces.append(
                go.Scatter(
                    x=t,
                    y=ys,
                    mode="lines+markers",
                    name=comp.name,
                )
            )

        self._replace_traces(traces)

        xlab = "time" if self.dt is not None else "timestep (block index)"
        title = f"{prop} vs {xlab} (z≈{z_used:g} m)"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(title=xlab),
                yaxis=dict(title=f"{prop} at z≈{z_used:g} m"),
            )
            self.fig.update_yaxes(autorange=True)

    def _render(self):
        prop = self.w_prop.value
        comp_ids = self._selected_components()

        if self.w_xmode.value == "grid":
            self._render_grid(prop=prop, block_index=int(self.w_block.value), comp_ids=comp_ids)
        else:
            self._render_time(prop=prop, z_value=float(self.w_z.value), comp_ids=comp_ids)

    def widget(self):
        return widgets.VBox([self.controls, self.fig])

    def display(self):
        if display is None:
            raise RuntimeError("IPython.display.display not available.")
        display(self.widget())
        return self
