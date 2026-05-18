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
from .utils import ComponentInfo, load_simulation_metadata, read_blocks


class BreakthroughPlotly(BasePlotly):
    """
    Notebook-native Plotly interface for breakthrough and column-profile data.

    File layout:
      - column.data contains only column-level properties
      - each component file contains only that component's properties

    Breakthrough curves are reconstructed by selecting one grid row from each
    stored block in the component_*.data files. By default, the last grid point
    is used, i.e. the outlet.
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
        return super()._read_component_data(fileName, min_columns=10)

    def _read_component_blocks(self, fileName: Union[str, Path]) -> List[np.ndarray]:
        return read_blocks(self.data_dir / fileName)

    def _read_column_data(self, fileName: Union[str, Path] = "column.data") -> List[np.ndarray]:
        return read_blocks(self.data_dir / fileName)

    def _plot_title(self, mode: str = "Breakthrough") -> str:
        return (
            f"{self.displayName} — {mode}  "
            f"T={self.externalTemperature:g} K, "
            f"p_t={self.externalPressure * 1e-3:g} kPa"
        )

    def _trace_mode(self, show_markers: bool) -> str:
        return "lines+markers" if show_markers else "lines"

    def _metric_info(self, metric: str) -> Dict[str, object]:
        metric_map: Dict[str, Dict[str, object]] = {
            # column.data
            "V": {"ylabel": "Interstitial velocity, v [m/s]", "source": "column", "col_1based": 4},
            "Pt": {"ylabel": "Total pressure, p_t [Pa]", "source": "column", "col_1based": 5},
            "Tg": {"ylabel": "Gas temperature, T_g [K]", "source": "column", "col_1based": 6},
            "dTgdt": {"ylabel": "Gas temperature derivative, dT_g/dt [K/s]", "source": "column", "col_1based": 7},
            "Ts": {"ylabel": "Solid temperature, T_s [K]", "source": "column", "col_1based": 8},
            "dTsdt": {"ylabel": "Solid temperature derivative, dT_s/dt [K/s]", "source": "column", "col_1based": 9},
            "Tw": {"ylabel": "Wall temperature, T_w [K]", "source": "column", "col_1based": 10},
            "dTwdt": {"ylabel": "Wall temperature derivative, dT_w/dt [K/s]", "source": "column", "col_1based": 11},
            "rho": {"ylabel": "Gas density, ρ_g [kg/m³]", "source": "column", "col_1based": 12},

            # component_*.data
            "C": {"ylabel": "Concentration, c_i [mol/m³]", "source": "component", "col_1based": 4},
            "Dcdt": {
                "ylabel": "Concentration derivative, dc_i/dt [mol/m³/s]",
                "source": "component",
                "col_1based": 5,
            },
            "Q": {"ylabel": "Adsorption, q_i [mol/kg]", "source": "component", "col_1based": 6},
            "Dqdt": {"ylabel": "Adsorption derivative, dq_i/dt [mol/kg/s]", "source": "component", "col_1based": 7},
            "P": {"ylabel": "Partial pressure, p_i [Pa]", "source": "component", "col_1based": 8},
            "Qeq": {"ylabel": "Equilibrium adsorption, q_i* [mol/kg]", "source": "component", "col_1based": 9},
            "Pnorm": {
                "ylabel": "Normalized partial pressure, p_i/p_{t}y_{i,0} [-]",
                "source": "component",
                "col_1based": 10,
            },
        }

        if metric not in metric_map:
            raise ValueError(f"Unknown metric '{metric}'")

        return metric_map[metric]

    def _column_index_0based(self, metric: str) -> int:
        info = self._metric_info(metric)
        return int(info["col_1based"]) - 1

    def _resolve_grid_index(self, block: np.ndarray, grid_index: int) -> int:
        ngrid_local = block.shape[0]
        idx = grid_index if grid_index >= 0 else ngrid_local + grid_index
        if idx < 0 or idx >= ngrid_local:
            raise IndexError(f"grid_index {grid_index} out of range for block with {ngrid_local} grid points")
        return idx

    def _breakthrough_component_series(
        self,
        fileName: Union[str, Path],
        grid_index: int = -1,
    ) -> Dict[str, np.ndarray]:
        blocks = self._read_component_blocks(fileName)
        if len(blocks) == 0:
            raise ValueError(f"No blocks found in {fileName}")

        rows = []
        for block in blocks:
            idx = self._resolve_grid_index(block, grid_index)
            rows.append(block[idx, :])

        data = np.asarray(rows, dtype=float)

        return {
            "col_1": data[:, 0],  # dimensionless time
            "col_2": data[:, 1],  # time [min]
            "col_3": data[:, 9],  # normalized partial pressure [-]
            "col_4": data[:, 7],  # partial pressure [Pa]
            "col_5": data[:, 3],  # concentration [mol/m^3]
            "col_6": data[:, 4],  # dc_i/dt [mol/m^3/s]
            "col_7": data[:, 5],  # q_i [mol/kg]
            "col_8": data[:, 6],  # dq_i/dt [mol/kg/s]
            "col_9": data[:, 8],  # q_i^* [mol/kg]
            "z": data[:, 2],      # selected z, usually constant across blocks
        }

    def _breakthrough_figure(
        self,
        x_col: str,
        xaxis_title: str,
        y_col: str = "col_3",
        yaxis_title: str = "Concentration exit gas, c<sub>i</sub>/c<sub>i,0</sub> [-]",
        mode: str = "Breakthrough",
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        fig = go.Figure()
        ymax = 0.0
        ymin = np.inf
        xmin = None
        xmax = None
        z_selected = None

        for comp in sorted(self.components, key=self._component_sort_key):
            if (not include_carrier_gas) and comp.isCarrierGas:
                continue

            fileName = self._component_file_name(comp.index, comp.name)
            series = self._breakthrough_component_series(fileName, grid_index=grid_index)

            x = np.asarray(series[x_col], dtype=float)
            y = np.asarray(series[y_col], dtype=float)

            if z_selected is None and "z" in series and len(series["z"]) > 0:
                z_selected = float(series["z"][0])

            if np.isfinite(x).any():
                x_local_min = float(np.nanmin(x))
                x_local_max = float(np.nanmax(x))
                xmin = x_local_min if xmin is None else min(xmin, x_local_min)
                xmax = x_local_max if xmax is None else max(xmax, x_local_max)

            if np.isfinite(y).any():
                ymax = max(ymax, float(np.nanmax(y)))
                ymin = min(ymin, float(np.nanmin(y)))

            line_width = 3.6 if comp.isCarrierGas else 3.2
            label = self._component_label(comp, include_yi=True)

            trace_kwargs = {
                "x": x,
                "y": y,
                "mode": self._trace_mode(show_markers),
                "name": label,
                "line": dict(width=line_width),
                "hovertemplate": (f"{label}<br>" + "x=%{x:.3f}<br>" + "y=%{y:.3f}<extra></extra>"),
            }
            if show_markers:
                trace_kwargs["marker"] = dict(symbol="circle", size=8, line=dict(width=1, color="Black"))

            fig.add_trace(go.Scatter(**trace_kwargs))

        if xmax is None:
            xmin = 0.0
            xmax = 1.0

        if not np.isfinite(ymin):
            ymin = 0.0
        if not np.isfinite(ymax):
            ymax = 1.0

        if ymin >= 0.0:
            y_lower = 0.0
            y_upper = 1.05 if ymax <= 1.02 and y_col == "col_3" else (1.08 * ymax if ymax > 0.0 else 1.0)
        else:
            y_lower = 1.08 * ymin
            y_upper = 1.08 * ymax if ymax > 0.0 else 1.0

        title_text = self._plot_title(mode)
        if z_selected is not None:
            title_text += f"  z={z_selected:g} m"

        fig = self._publication_layout(
            fig,
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title,
            title_text=title_text,
            x_range=[xmin, xmax],
            y_range=[y_lower, y_upper],
        )
        return fig

    def breakthrough_dimensionless(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_1",
            xaxis_title="Dimensionless time, τ = tv/L [-]",
            y_col="col_3",
            yaxis_title="Concentration exit gas, c<sub>i</sub>/c<sub>i,0</sub> [-]",
            mode="Breakthrough",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def breakthrough_time(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_2",
            xaxis_title="Time, t [min]",
            y_col="col_3",
            yaxis_title="Concentration exit gas, c<sub>i</sub>/c<sub>i,0</sub> [-]",
            mode="Breakthrough",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def breakthrough_concentration_dimensionless(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_1",
            xaxis_title="Dimensionless time, τ = tv/L [-]",
            y_col="col_5",
            yaxis_title="Concentration, c<sub>i</sub> [mol/m³]",
            mode="Breakthrough (concentration)",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def breakthrough_concentration_time(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_2",
            xaxis_title="Time, t [min]",
            y_col="col_5",
            yaxis_title="Concentration, c<sub>i</sub> [mol/m³]",
            mode="Breakthrough (concentration)",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def breakthrough_partial_pressure_dimensionless(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_1",
            xaxis_title="Dimensionless time, τ = tv/L [-]",
            y_col="col_4",
            yaxis_title="Partial pressure, p<sub>i</sub> [Pa]",
            mode="Breakthrough (partial pressure)",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def breakthrough_partial_pressure_time(
        self,
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        return self._breakthrough_figure(
            x_col="col_2",
            xaxis_title="Time, t [min]",
            y_col="col_4",
            yaxis_title="Partial pressure, p<sub>i</sub> [Pa]",
            mode="Breakthrough (partial pressure)",
            include_carrier_gas=include_carrier_gas,
            show_markers=show_markers,
            grid_index=grid_index,
        )

    def temperature_triplet_time(
        self,
        grid_index: int = -1,
        show_markers: bool = True,
    ) -> go.Figure:
        blocks = self._read_column_data("column.data")
        if len(blocks) == 0:
            raise ValueError("No blocks found in column.data")

        rows = []
        for block in blocks:
            idx = self._resolve_grid_index(block, grid_index)
            rows.append(block[idx, :])

        data = np.asarray(rows, dtype=float)

        t = data[:, 1]  # time [min]
        tg = data[:, 5]  # T_g [K]
        ts = data[:, 7]  # T_s [K]
        tw = data[:, 9]  # T_w [K]
        z_selected = float(data[0, 2]) if data.shape[0] > 0 else None

        fig = go.Figure()

        trace_specs = [
            ("Gas temperature", tg),
            ("Solid temperature", ts),
            ("Wall temperature", tw),
        ]

        ymin = np.inf
        ymax = -np.inf

        for label, y in trace_specs:
            if np.isfinite(y).any():
                ymin = min(ymin, float(np.nanmin(y)))
                ymax = max(ymax, float(np.nanmax(y)))

            trace_kwargs = {
                "x": t,
                "y": y,
                "mode": self._trace_mode(show_markers),
                "name": label,
                "line": dict(width=3.2),
                "hovertemplate": (f"{label}<br>" + "t=%{x:.3f} min<br>" + "T=%{y:.3f} K<extra></extra>"),
            }
            if show_markers:
                trace_kwargs["marker"] = dict(symbol="circle", size=8, line=dict(width=1, color="Black"))

            fig.add_trace(go.Scatter(**trace_kwargs))

        if not np.isfinite(ymin):
            ymin = 0.0
        if not np.isfinite(ymax):
            ymax = 1.0

        if ymin >= 0.0:
            y_lower = 0.0
            y_upper = 1.08 * ymax if ymax > 0.0 else 1.0
        else:
            y_lower = 1.08 * ymin
            y_upper = 1.08 * ymax if ymax > 0.0 else 1.0

        x_min = float(np.nanmin(t)) if np.isfinite(t).any() else 0.0
        x_max = float(np.nanmax(t)) if np.isfinite(t).any() else 1.0

        title_text = self._plot_title("Temperature history")
        if z_selected is not None:
            title_text += f"  z={z_selected:g} m"

        fig = self._publication_layout(
            fig,
            xaxis_title="Time, t [min]",
            yaxis_title="Temperature [K]",
            title_text=title_text,
            x_range=[x_min, x_max],
            y_range=[y_lower, y_upper],
        )
        return fig

    def column_snapshot(
        self,
        metric: str,
        frame_index: int = 0,
        show_markers: bool = True,
    ) -> go.Figure:
        info = self._metric_info(metric)
        ylabel = str(info["ylabel"])
        source = str(info["source"])
        fig = go.Figure()

        if source == "column":
            blocks = self._read_column_data("column.data")
            if frame_index < 0 or frame_index >= len(blocks):
                raise IndexError(f"frame_index {frame_index} out of range [0, {len(blocks) - 1}]")

            block = blocks[frame_index]
            z = block[:, 2]
            col = self._column_index_0based(metric)

            fig.add_trace(
                go.Scatter(
                    x=z,
                    y=block[:, col],
                    mode=self._trace_mode(show_markers),
                    name=metric,
                )
            )
        else:
            for comp in self.components:
                fileName = self._component_file_name(comp.index, comp.name)
                blocks = self._read_component_blocks(fileName)

                if frame_index < 0 or frame_index >= len(blocks):
                    raise IndexError(f"frame_index {frame_index} out of range for {fileName}")

                block = blocks[frame_index]
                z = block[:, 2]
                col = self._column_index_0based(metric)
                y = block[:, col]

                fig.add_trace(
                    go.Scatter(
                        x=z,
                        y=y,
                        mode=self._trace_mode(show_markers),
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

    def column_animation(
        self,
        metric: str,
        every: int = 1,
        show_markers: bool = True,
    ) -> go.Figure:
        info = self._metric_info(metric)
        ylabel = str(info["ylabel"])
        source = str(info["source"])

        sampled_indices: List[int]
        fig = go.Figure()

        if source == "column":
            blocks = self._read_column_data("column.data")
            sampled_indices = list(range(0, len(blocks), max(1, every)))
            sampled_blocks = [blocks[i] for i in sampled_indices]

            x_min = min(block[:, 2].min() for block in sampled_blocks)
            x_max = max(block[:, 2].max() for block in sampled_blocks)

            col = self._column_index_0based(metric)
            y_min = min(block[:, col].min() for block in sampled_blocks)
            y_max = max(block[:, col].max() for block in sampled_blocks)

            first_block = sampled_blocks[0]
            fig.add_trace(
                go.Scatter(
                    x=first_block[:, 2],
                    y=first_block[:, col],
                    mode=self._trace_mode(show_markers),
                    name=metric,
                )
            )

            frames = []
            for idx, block in zip(sampled_indices, sampled_blocks):
                frames.append(
                    go.Frame(
                        name=str(idx),
                        data=[
                            go.Scatter(
                                x=block[:, 2],
                                y=block[:, col],
                                mode=self._trace_mode(show_markers),
                                name=metric,
                            )
                        ],
                        layout=go.Layout(title_text=f"{self._plot_title()} — {metric} — frame {idx}"),
                    )
                )
        else:
            comp_blocks: Dict[int, List[np.ndarray]] = {}
            nframes = None

            for comp in self.components:
                fileName = self._component_file_name(comp.index, comp.name)
                blocks = self._read_component_blocks(fileName)
                comp_blocks[comp.index] = blocks
                if nframes is None:
                    nframes = len(blocks)
                elif len(blocks) != nframes:
                    raise ValueError("Component files do not have the same number of blocks")

            if nframes is None or nframes == 0:
                raise ValueError("No component blocks found")

            sampled_indices = list(range(0, nframes, max(1, every)))

            x_values = []
            y_values = []
            col = self._column_index_0based(metric)

            for comp in self.components:
                for idx in sampled_indices:
                    block = comp_blocks[comp.index][idx]
                    x_values.append(block[:, 2])
                    y_values.append(block[:, col])

            x_all = np.concatenate(x_values)
            y_all = np.concatenate(y_values)
            x_min = float(np.min(x_all))
            x_max = float(np.max(x_all))
            y_min = float(np.min(y_all))
            y_max = float(np.max(y_all))

            first_idx = sampled_indices[0]
            for comp in self.components:
                block = comp_blocks[comp.index][first_idx]
                fig.add_trace(
                    go.Scatter(
                        x=block[:, 2],
                        y=block[:, col],
                        mode=self._trace_mode(show_markers),
                        name=self._component_label(comp, include_yi=True),
                    )
                )

            frames = []
            for idx in sampled_indices:
                data = []
                for comp in self.components:
                    block = comp_blocks[comp.index][idx]
                    data.append(
                        go.Scatter(
                            x=block[:, 2],
                            y=block[:, col],
                            mode=self._trace_mode(show_markers),
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

        if y_min >= 0.0:
            y_min = 0.0
            y_max = 1.1 * y_max if y_max > 0.0 else 1.0
        else:
            y_min = 1.1 * y_min
            y_max = 1.1 * y_max

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

    def all_column_figures(self, every: int = 1, show_markers: bool = True) -> Dict[str, go.Figure]:
        return {
            metric: self.column_animation(metric, every=every, show_markers=show_markers)
            for metric in [
                "V",
                "Pt",
                "Tg",
                "dTgdt",
                "Ts",
                "dTsdt",
                "Tw",
                "dTwdt",
                "rho",
                "C",
                "Dcdt",
                "Q",
                "Dqdt",
                "P",
                "Qeq",
                "Pnorm",
            ]
        }

    def explorer(
        self,
        column_file: Union[str, Path] = "column.data",
        dt: Optional[float] = None,
        time0: float = 0.0,
        legend_title: str = "",
        show_markers: bool = True,
    ) -> "ColumnDataExplorer":
        if not legend_title:
            legend_title = self._plot_title()

        dt_use = self.dt_guess(dt)
        return ColumnDataExplorer(
            data_dir=self.data_dir,
            components=self.components,
            component_file_names={
                comp.index: self._component_file_name(comp.index, comp.name) for comp in self.components
            },
            column_file=column_file,
            dt=dt_use,
            time0=time0,
            legend_title=legend_title,
            carrierGasComponent=self.carrierGasComponent,
            show_markers=show_markers,
        )

    def dt_guess(self, dt: Optional[float]) -> Optional[float]:
        return dt


class ColumnDataExplorer:
    """
    Interactive explorer for the split file structure:
      - column.data for column-level variables
      - per-component files for component-level variables
    """

    def __init__(
        self,
        data_dir: Union[str, Path],
        components: List[ComponentInfo],
        component_file_names: Dict[int, Union[str, Path]],
        column_file: Union[str, Path] = "column.data",
        dt: Optional[float] = None,
        time0: float = 0.0,
        legend_title: str = "",
        carrierGasComponent: Optional[int] = None,
        show_markers: bool = True,
    ):
        if widgets is None:
            raise ImportError("ipywidgets is required for ColumnDataExplorer (pip install ipywidgets).")

        self.data_dir = Path(data_dir)
        self.components = list(components)
        self.component_file_names = {int(k): Path(v) for k, v in component_file_names.items()}
        self.column_file = self.data_dir / column_file

        self.column_blocks = read_blocks(self.column_file)
        if len(self.column_blocks) == 0:
            raise ValueError(f"No blocks found in {self.column_file}")

        self.component_blocks: Dict[int, List[np.ndarray]] = {}
        for comp in self.components:
            if comp.index not in self.component_file_names:
                raise KeyError(f"Missing component file name for component index {comp.index}")
            fileName = self.data_dir / self.component_file_names[comp.index]
            self.component_blocks[comp.index] = read_blocks(fileName)

        self.dt = dt
        self.time0 = time0
        self.legend_title = legend_title
        self.carrierGasComponent = carrierGasComponent
        self.show_markers = show_markers

        self.z_grid = self.column_blocks[0][:, 2].astype(float)
        self.zmin = float(np.nanmin(self.z_grid))
        self.zmax = float(np.nanmax(self.z_grid))

        self.column_props = {
            "V": {"ylabel": "Interstitial velocity, v [m/s]", "col_0based": 3, "source": "column"},
            "Pt": {"ylabel": "Total pressure, p_t [Pa]", "col_0based": 4, "source": "column"},
            "Tg": {"ylabel": "Gas temperature, T_g [K]", "col_0based": 5, "source": "column"},
            "dTgdt": {"ylabel": "Gas temperature derivative, dT_g/dt [K/s]", "col_0based": 6, "source": "column"},
            "Ts": {"ylabel": "Solid temperature, T_s [K]", "col_0based": 7, "source": "column"},
            "dTsdt": {"ylabel": "Solid temperature derivative, dT_s/dt [K/s]", "col_0based": 8, "source": "column"},
            "Tw": {"ylabel": "Wall temperature, T_w [K]", "col_0based": 9, "source": "column"},
            "dTwdt": {"ylabel": "Wall temperature derivative, dT_w/dt [K/s]", "col_0based": 10, "source": "column"},
            "rho": {"ylabel": "Gas density, ρ_g [kg/m³]", "col_0based": 11, "source": "column"},
        }

        self.component_props = {
            "C": {"ylabel": "Concentration, c_i [mol/m³]", "col_0based": 3, "source": "component"},
            "Dcdt": {
                "ylabel": "Concentration derivative, dc_i/dt [mol/m³/s]",
                "col_0based": 4,
                "source": "component",
            },
            "Q": {"ylabel": "Adsorption, q_i [mol/kg]", "col_0based": 5, "source": "component"},
            "Dqdt": {"ylabel": "Adsorption derivative, dq_i/dt [mol/kg/s]", "col_0based": 6, "source": "component"},
            "P": {"ylabel": "Partial pressure, p_i [Pa]", "col_0based": 7, "source": "component"},
            "Qeq": {"ylabel": "Equilibrium adsorption, q_i* [mol/kg]", "col_0based": 8, "source": "component"},
            "Pnorm": {
                "ylabel": "Normalized partial pressure, p_i/p_{t}y_{i,0} [-]",
                "col_0based": 9,
                "source": "component",
            },
        }

        self.available_props = sorted(list(self.column_props.keys()) + list(self.component_props.keys()))

        self._build_widgets()
        self._build_figure()
        self._render()

    def _trace_mode(self) -> str:
        return "lines+markers" if self.show_markers else "lines"

    def _build_widgets(self):
        default_prop = "C" if "C" in self.available_props else self.available_props[0]
        self.w_prop = widgets.Dropdown(options=self.available_props, value=default_prop, description="y:")

        self.w_xmode = widgets.Dropdown(
            options=[("grid (z)", "grid"), ("time", "time")],
            value="grid",
            description="x:",
        )

        self.w_block = widgets.IntSlider(value=0, min=0, max=len(self.column_blocks) - 1, step=1, description="block")

        step = (self.zmax - self.zmin) / 200 if (self.zmax > self.zmin) else 1.0
        if step <= 0 or not np.isfinite(step):
            step = 1.0
        self.w_z = widgets.FloatSlider(value=self.zmin, min=self.zmin, max=self.zmax, step=step, description="z@")

        comp_options = [(f"{comp.index}: {comp.name}", comp.index) for comp in self.components]
        self.w_comps = widgets.SelectMultiple(
            options=comp_options,
            value=tuple(ci for _, ci in comp_options),
            description="components",
            rows=min(8, max(3, len(comp_options))),
        )

        self.w_show_carrier = widgets.Checkbox(
            value=False,
            description="show carrier gas",
            indent=False,
        )

        self.w_show_markers = widgets.Checkbox(
            value=self.show_markers,
            description="show markers",
            indent=False,
        )

        self.w_save_path = widgets.Text(
            value="column_plot.pdf",
            description="save as",
            layout=widgets.Layout(width="420px"),
        )
        self.w_save_button = widgets.Button(
            description="Save PDF",
            tooltip="Save current figure to PDF",
            button_style="",
        )
        self.w_save_status = widgets.HTML(value="")

        self.w_save_button.on_click(self._save_current_pdf)

        self.controls_top = widgets.HBox([self.w_xmode, self.w_prop, self.w_block])
        self.controls_bottom = widgets.HBox([self.w_z, self.w_show_carrier, self.w_show_markers])
        self.controls_export = widgets.HBox([self.w_save_path, self.w_save_button, self.w_save_status])
        self.controls = widgets.VBox([self.controls_top, self.controls_bottom, self.w_comps, self.controls_export])

        for w in [self.w_prop, self.w_xmode, self.w_block, self.w_z, self.w_comps, self.w_show_carrier]:
            w.observe(self._on_change, names="value")
        self.w_show_markers.observe(self._on_markers_change, names="value")

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

        prop = str(self.w_prop.value)
        is_component_prop = prop in self.component_props
        self.w_comps.layout.display = "" if is_component_prop else "none"
        self.w_show_carrier.layout.display = "" if is_component_prop else "none"

    def _on_change(self, _change):
        self._sync_visibility()
        self._render()

    def _on_markers_change(self, _change):
        self.show_markers = bool(self.w_show_markers.value)
        self._render()

    def _time_axis(self) -> np.ndarray:
        n = len(self.column_blocks)
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

    def _default_pdf_name(self) -> str:
        prop = str(self.w_prop.value)
        if self.w_xmode.value == "grid":
            return f"{prop}_grid_block_{int(self.w_block.value)}.pdf"

        z_value = float(self.w_z.value)
        return f"{prop}_time_z_{z_value:.4g}.pdf"

    def _save_current_pdf(self, _button) -> None:
        out = self.w_save_path.value.strip() or self._default_pdf_name()
        path = Path(out)

        if path.suffix.lower() != ".pdf":
            path = path.with_suffix(".pdf")

        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            self.fig.write_image(str(path), format="pdf")
            self.w_save_status.value = f"<span style='color: green;'>Saved: {path}</span>"
        except Exception as exc:
            self.w_save_status.value = f"<span style='color: red;'>Save failed: {exc}</span>"

    def _include_component(self, ci: int) -> bool:
        if self.carrierGasComponent is not None and (not self.w_show_carrier.value) and ci == self.carrierGasComponent:
            return False
        return True

    def _render_grid(self, prop: str, block_index: int, comp_ids: List[int]):
        traces: List[go.Scatter] = []
        ymin = np.inf
        ymax = -np.inf

        if prop in self.column_props:
            info = self.column_props[prop]
            block = self.column_blocks[block_index]
            z = block[:, 2].astype(float)
            col = int(info["col_0based"])
            y = block[:, col].astype(float)

            if np.isfinite(y).any():
                ymin = min(ymin, float(np.nanmin(y)))
                ymax = max(ymax, float(np.nanmax(y)))

            traces.append(
                go.Scatter(
                    x=z,
                    y=y,
                    mode=self._trace_mode(),
                    name=prop,
                )
            )

            ylabel = str(info["ylabel"])
        else:
            info = self.component_props[prop]
            col = int(info["col_0based"])

            for ci in comp_ids:
                if not self._include_component(ci):
                    continue

                blocks = self.component_blocks[ci]
                block = blocks[block_index]
                z = block[:, 2].astype(float)
                y = block[:, col].astype(float)

                if np.isfinite(y).any():
                    ymin = min(ymin, float(np.nanmin(y)))
                    ymax = max(ymax, float(np.nanmax(y)))

                comp = next(comp for comp in self.components if comp.index == ci)
                traces.append(
                    go.Scatter(
                        x=z,
                        y=y,
                        mode=self._trace_mode(),
                        name=comp.name,
                    )
                )

            ylabel = str(info["ylabel"])

        self._replace_traces(traces)

        title = f"{prop} vs z (block {block_index})"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(
                    title="Adsorber position [m]",
                    range=[self.zmin, self.zmax],
                    autorange=False,
                ),
                yaxis=dict(title=ylabel),
            )

            if np.isfinite(ymin) and np.isfinite(ymax):
                if ymin >= 0.0 and ymax > 0.0:
                    self.fig.update_yaxes(range=[0.0, 1.1 * ymax], autorange=False)
                elif ymin < 0.0:
                    self.fig.update_yaxes(range=[1.1 * ymin, 1.1 * ymax if ymax != 0 else 1.0], autorange=False)
                else:
                    self.fig.update_yaxes(autorange=True)
            else:
                self.fig.update_yaxes(autorange=True)

    def _render_time(self, prop: str, z_value: float, comp_ids: List[int]):
        t = self._time_axis()
        zi = self._nearest_z_index(z_value)
        z_used = float(self.z_grid[zi])

        traces: List[go.Scatter] = []

        if prop in self.column_props:
            info = self.column_props[prop]
            col = int(info["col_0based"])

            ys = []
            for block in self.column_blocks:
                if zi < block.shape[0] and col < block.shape[1]:
                    ys.append(float(block[zi, col]))
                else:
                    ys.append(np.nan)
            ys = np.asarray(ys, dtype=float)

            traces.append(
                go.Scatter(
                    x=t,
                    y=ys,
                    mode=self._trace_mode(),
                    name=prop,
                )
            )
            ylabel = str(info["ylabel"])
        else:
            info = self.component_props[prop]
            col = int(info["col_0based"])

            for ci in comp_ids:
                if not self._include_component(ci):
                    continue

                ys = []
                for block in self.component_blocks[ci]:
                    if zi < block.shape[0] and col < block.shape[1]:
                        ys.append(float(block[zi, col]))
                    else:
                        ys.append(np.nan)
                ys = np.asarray(ys, dtype=float)

                comp = next(comp for comp in self.components if comp.index == ci)
                traces.append(
                    go.Scatter(
                        x=t,
                        y=ys,
                        mode=self._trace_mode(),
                        name=comp.name,
                    )
                )
            ylabel = str(info["ylabel"])

        self._replace_traces(traces)

        xlab = "time" if self.dt is not None else "timestep (block index)"
        title = f"{prop} vs {xlab} (z≈{z_used:g} m)"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(
                    title=xlab,
                    autorange=True,
                    range=None,
                ),
                yaxis=dict(title=ylabel),
            )
            self.fig.update_xaxes(autorange=True)
            self.fig.update_yaxes(autorange=True)

    def _render(self):
        prop = str(self.w_prop.value)
        comp_ids = self._selected_components()

        if self.w_xmode.value == "grid":
            self._render_grid(prop=prop, block_index=int(self.w_block.value), comp_ids=comp_ids)
        else:
            self._render_time(prop=prop, z_value=float(self.w_z.value), comp_ids=comp_ids)

        if not self.w_save_path.value.strip() or self.w_save_path.value == "column_plot.pdf":
            self.w_save_path.value = self._default_pdf_name()

    def widget(self):
        return widgets.VBox([self.controls, self.fig])

    def display(self):
        if display is None:
            raise RuntimeError("IPython.display.display not available.")
        display(self.widget())
        return self