from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

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


@dataclass(frozen=True)
class MetricSpec:
    ylabel: str
    source: str
    col_0based: int
    header_label: Optional[str] = None

    @property
    def col_1based(self) -> int:
        return self.col_0based + 1

    @property
    def selector_label(self) -> str:
        return self.header_label or self.ylabel


@dataclass(frozen=True)
class AxisSpec:
    title: str
    col_0based: int
    scale: float = 1.0


@dataclass(frozen=True)
class TemperatureUnitSpec:
    title: str
    hover_unit: str
    scale: float = 1.0
    offset: float = 0.0


@dataclass(frozen=True)
class BreakthroughYSpec:
    title: str
    mode: str
    component_col_0based: Optional[int]
    hover_label: str
    normalized: bool = False


COLUMN_METRICS: Dict[str, MetricSpec] = {
    "V": MetricSpec("Interstitial velocity, v [m/s]", "column", 3, "Interstitial gas velocity, v [m/s]"),
    "Pt": MetricSpec("Total pressure, p<sub>t</sub> [Pa]", "column", 4, "Total pressure, p_t [Pa]"),
    "Tg": MetricSpec("Gas temperature, T<sub>g</sub> [K]", "column", 5, "Gas temperature, T_g [K]"),
    "dTgdt": MetricSpec(
        "Gas temperature derivative, dT<sub>g</sub>/dt [K/s]",
        "column",
        6,
        "Gas temperature time derivative, dT_g/dt [K/s]",
    ),
    "Ts": MetricSpec("Solid temperature, T<sub>s</sub> [K]", "column", 7, "Solid temperature, T_s [K]"),
    "dTsdt": MetricSpec(
        "Solid temperature derivative, dT<sub>s</sub>/dt [K/s]",
        "column",
        8,
        "Solid temperature time derivative, dT_s/dt [K/s]",
    ),
    "Tw": MetricSpec("Wall temperature, T<sub>w</sub> [K]", "column", 9, "Wall temperature, T_w [K]"),
    "dTwdt": MetricSpec(
        "Wall temperature derivative, dT<sub>w</sub>/dt [K/s]",
        "column",
        10,
        "Wall temperature time derivative, dT_w/dt [K/s]",
    ),
    "rho": MetricSpec("Gas density, ρ<sub>g</sub> [kg/m³]", "column", 11, "Gas density, rho_g [kg/m^3]"),
}

COMPONENT_METRICS: Dict[str, MetricSpec] = {
    "C": MetricSpec("Concentration, c<sub>i</sub> [mol/m³]", "component", 3, "Concentration, c_i [mol/m^3]"),
    "Cdot": MetricSpec(
        "Concentration derivative, dc<sub>i</sub>/dt [mol/m³/s]",
        "component",
        4,
        "Concentration time derivative, dc_i/dt [mol/m^3/s]",
    ),
    "Dcdt": MetricSpec(
        "Concentration derivative, dc<sub>i</sub>/dt [mol/m³/s]",
        "component",
        4,
        "Concentration time derivative, dc_i/dt [mol/m^3/s]",
    ),
    "Y": MetricSpec("Mole fraction, y<sub>i</sub> [-]", "component", 5, "Mole fraction, y_i [-]"),
    "Q": MetricSpec("Physisorption, q<sub>phy,i</sub> [mol/kg]", "component", 6, "Physisorption, q_phy_i [mol/kg]"),
    "Dqdt": MetricSpec(
        "Physisorption derivative, dq<sub>phy,i</sub>/dt [mol/kg/s]",
        "component",
        7,
        "Physisorption time derivative, dq_phy_i/dt [mol/kg/s]",
    ),
    "Qchem": MetricSpec(
        "Chemisorption, q<sub>chem,i</sub> [mol/kg]",
        "component",
        8,
        "Chemisorption, q_chem_i [mol/kg]",
    ),
    "Dqchemdt": MetricSpec(
        "Chemisorption derivative, dq<sub>chem,i</sub>/dt [mol/kg/s]",
        "component",
        9,
        "Chemisorption time derivative, dq_chem_i/dt [mol/kg/s]",
    ),
    "P": MetricSpec("Partial pressure, p<sub>i</sub> [Pa]", "component", 10, "Partial pressure, p_i [Pa]"),
    "Qeq": MetricSpec(
        "Equilibrium physisorption, q<sub>phy,i</sub>* [mol/kg]",
        "component",
        11,
        "Equilibrium physisorption, q_phy_i^* [mol/kg]",
    ),
    "Pnorm": MetricSpec(
        "Normalized partial pressure, p<sub>i</sub>/(p<sub>t</sub>y<sub>i,0</sub>) [-]",
        "component",
        12,
        "Normalized partial pressure, p_i / (p_t y_i,0) [-]",
    ),
    "Qchemeq": MetricSpec(
        "Equilibrium chemisorption, q<sub>chem,i</sub>* [mol/kg]",
        "component",
        13,
        "Equilibrium chemisorption, q_chem_i^* [mol/kg]",
    ),
}

METRIC_SPECS: Dict[str, MetricSpec] = {**COLUMN_METRICS, **COMPONENT_METRICS}

BREAKTHROUGH_X_SPECS: Dict[str, AxisSpec] = {
    "dimensionless": AxisSpec("Dimensionless time, τ = tv/L [-]", 0, 1.0),
    "min": AxisSpec("Time, t [min]", 1, 1.0),
    "s": AxisSpec("Time, t [s]", 1, 60.0),
    "hr": AxisSpec("Time, t [hr]", 1, 1 / 60.0),
}

BREAKTHROUGH_Y_SPECS: Dict[str, BreakthroughYSpec] = {
    "concentration": BreakthroughYSpec(
        title="Concentration, c<sub>i</sub> [mol/m³]",
        mode="Breakthrough (concentration)",
        component_col_0based=COMPONENT_METRICS["C"].col_0based,
        hover_label="c<sub>i</sub>",
    ),
    "normalized_concentration": BreakthroughYSpec(
        title="Normalized concentration, c<sub>i</sub>/c<sub>i,0</sub> [-]",
        mode="Breakthrough (normalized concentration)",
        component_col_0based=COMPONENT_METRICS["Pnorm"].col_0based,
        hover_label="c<sub>i</sub>/c<sub>i,0</sub>",
        normalized=True,
    ),
    "molefraction": BreakthroughYSpec(
        title="Mole fraction, y<sub>i</sub> [-]",
        mode="Breakthrough (mole fraction)",
        component_col_0based=COMPONENT_METRICS["Y"].col_0based,
        hover_label="y<sub>i</sub>",
        normalized=True,
    ),
    "partial_pressure": BreakthroughYSpec(
        title="Partial pressure, p<sub>i</sub> [Pa]",
        mode="Breakthrough (partial pressure)",
        component_col_0based=COMPONENT_METRICS["P"].col_0based,
        hover_label="p<sub>i</sub>",
    ),
}

TEMPERATURE_Y_SPECS: Dict[str, TemperatureUnitSpec] = {
    "kelvin": TemperatureUnitSpec("Temperature [K]", "K"),
    "celsius": TemperatureUnitSpec("Temperature [°C]", "°C", offset=-273.15),
}

X_UNITS = set(BREAKTHROUGH_X_SPECS)
Y_UNITS = set(BREAKTHROUGH_Y_SPECS)
TEMPERATURE_Y_UNITS = set(TEMPERATURE_Y_SPECS)


def _canonicalize_component_blocks(
    blocks: List[np.ndarray],
    column_blocks: Optional[List[np.ndarray]] = None,
) -> List[np.ndarray]:
    """Upgrade legacy 10-column multibed output to the 14-column component schema."""

    canonical_blocks: List[np.ndarray] = []
    for block_index, block in enumerate(blocks):
        if block.ndim != 2 or block.shape[1] != 10:
            canonical_blocks.append(block)
            continue

        canonical = np.zeros((block.shape[0], 14), dtype=float)
        canonical[:, :5] = block[:, :5]
        canonical[:, 6] = block[:, 5]  # physisorption
        canonical[:, 7] = block[:, 6]  # physisorption derivative
        canonical[:, 10] = block[:, 7]  # partial pressure
        canonical[:, 11] = block[:, 8]  # equilibrium physisorption
        canonical[:, 12] = block[:, 9]  # normalized partial pressure
        canonical[:, 5] = np.nan

        if column_blocks is not None and block_index < len(column_blocks):
            column_block = column_blocks[block_index]
            matching_grid = (
                column_block.ndim == 2
                and column_block.shape[0] == block.shape[0]
                and column_block.shape[1] >= 5
                and np.allclose(column_block[:, :3], block[:, :3], equal_nan=True)
            )
            if matching_grid:
                total_pressure = column_block[:, COLUMN_METRICS["Pt"].col_0based]
                with np.errstate(divide="ignore", invalid="ignore"):
                    canonical[:, 5] = np.where(
                        total_pressure != 0.0,
                        canonical[:, COMPONENT_METRICS["P"].col_0based] / total_pressure,
                        np.nan,
                    )

        canonical_blocks.append(canonical)

    return canonical_blocks


def _canonical_key(value: str, allowed_values: set[str], option_name: str) -> str:
    key = str(value).strip().lower()
    if key in allowed_values:
        return key

    allowed = ", ".join(sorted(allowed_values))
    raise ValueError(f"Unknown {option_name} '{value}'. Allowed values: {allowed}")


def _props_from_specs(specs: Dict[str, MetricSpec]) -> Dict[str, Dict[str, object]]:
    return {
        name: {
            "ylabel": spec.ylabel,
            "header_label": spec.selector_label,
            "source": spec.source,
            "col_0based": spec.col_0based,
            "col_1based": spec.col_1based,
        }
        for name, spec in specs.items()
    }


def _selector_options_from_specs(specs: Dict[str, MetricSpec]) -> List[Tuple[str, str]]:
    options: List[Tuple[str, str]] = []
    seen_labels = set()
    for name, spec in specs.items():
        label = spec.selector_label
        if label in seen_labels:
            continue
        seen_labels.add(label)
        options.append((label, name))
    return options


def _column_explorer_y_options() -> List[Tuple[str, str]]:
    return _selector_options_from_specs(COLUMN_METRICS) + _selector_options_from_specs(COMPONENT_METRICS)


def _finite_min_max(values: List[np.ndarray]) -> Tuple[Optional[float], Optional[float]]:
    finite_chunks = []
    for value in values:
        arr = np.asarray(value, dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size:
            finite_chunks.append(finite)

    if not finite_chunks:
        return None, None

    all_values = np.concatenate(finite_chunks)
    return float(np.nanmin(all_values)), float(np.nanmax(all_values))


def _padded_y_range(ymin: Optional[float], ymax: Optional[float], normalized: bool = False) -> List[float]:
    if ymin is None or ymax is None:
        return [0.0, 1.0]

    if ymin >= 0.0:
        if normalized and ymax <= 1.02:
            return [0.0, 1.05]
        return [0.0, 1.08 * ymax if ymax > 0.0 else 1.0]

    span = ymax - ymin
    if span <= 0.0 or not np.isfinite(span):
        span = abs(ymax) if ymax != 0.0 else 1.0
    padding = 0.08 * span
    return [ymin - padding, ymax + padding]


class BreakthroughPlotly(BasePlotly):
    """
    Notebook-native Plotly interface for breakthrough and column-profile data.

    Expected file layout:
      - column.data contains column-level properties only.
      - component files contain component-level properties only.

    Breakthrough curves are reconstructed by selecting one grid row from each
    stored block in each component file. By default, the selected row is the
    outlet row, i.e. grid_index=-1.
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
        return super()._read_component_data(fileName, min_columns=11)

    def _read_component_blocks(self, fileName: Union[str, Path]) -> List[np.ndarray]:
        blocks = read_blocks(fileName)
        column_blocks = None
        if any(block.ndim == 2 and block.shape[1] == 10 for block in blocks):
            column_path = self.data_dir / "column.data"
            if column_path.exists():
                column_blocks = read_blocks(column_path)
        return _canonicalize_component_blocks(blocks, column_blocks)

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
        if metric not in METRIC_SPECS:
            allowed = ", ".join(sorted(METRIC_SPECS))
            raise ValueError(f"Unknown metric '{metric}'. Allowed metrics: {allowed}")

        spec = METRIC_SPECS[metric]
        return {
            "ylabel": spec.ylabel,
            "source": spec.source,
            "col_0based": spec.col_0based,
            "col_1based": spec.col_1based,
        }

    def _column_index_0based(self, metric: str) -> int:
        info = self._metric_info(metric)
        return int(info["col_0based"])

    def _resolve_grid_index(self, block: np.ndarray, grid_index: Optional[int]) -> int:
        ngrid_local = block.shape[0]
        if grid_index is None:
            return ngrid_local // 2

        idx = grid_index if grid_index >= 0 else ngrid_local + grid_index
        if idx < 0 or idx >= ngrid_local:
            raise IndexError(f"grid_index {grid_index} out of range for block with {ngrid_local} grid points")
        return idx

    def _rows_at_grid(
        self,
        blocks: List[np.ndarray],
        grid_index: Optional[int],
        fileName: Union[str, Path],
        min_columns: int,
    ) -> np.ndarray:
        if len(blocks) == 0:
            raise ValueError(f"No blocks found in {fileName}")

        blocks_with_columns = [
            (block_index, block)
            for block_index, block in enumerate(blocks)
            if block.shape[1] >= min_columns
        ]
        if len(blocks_with_columns) == 0:
            column_counts = ", ".join(str(block.shape[1]) for block in blocks)
            raise ValueError(
                f"No complete blocks found in {fileName}; expected at least {min_columns} columns "
                f"but found column counts: {column_counts}"
            )

        row_counts = {}
        for _, block in blocks_with_columns:
            row_counts[block.shape[0]] = row_counts.get(block.shape[0], 0) + 1
        expected_rows = max(row_counts, key=lambda count: (row_counts[count], count))

        rows = []
        for _, block in blocks_with_columns:
            if block.shape[0] != expected_rows:
                continue
            idx = self._resolve_grid_index(block, grid_index)
            rows.append(block[idx, :])

        if len(rows) == 0:
            raise ValueError(f"No complete blocks found in {fileName}")

        return np.asarray(rows, dtype=float)

    def _breakthrough_component_data(
        self,
        fileName: Union[str, Path],
        grid_index: Optional[int] = -1,
    ) -> np.ndarray:
        blocks = self._read_component_blocks(fileName)
        return self._rows_at_grid(blocks, grid_index=grid_index, fileName=fileName, min_columns=11)

    def _breakthrough_column_data(
        self,
        grid_index: Optional[int] = -1,
        fileName: Union[str, Path] = "column.data",
    ) -> np.ndarray:
        blocks = self._read_column_data(fileName)
        return self._rows_at_grid(blocks, grid_index=grid_index, fileName=fileName, min_columns=12)

    def _x_values(self, data: np.ndarray, x_units: str) -> np.ndarray:
        spec = BREAKTHROUGH_X_SPECS[x_units]
        return data[:, spec.col_0based].astype(float) * spec.scale

    def _y_values(
        self,
        data: np.ndarray,
        y_units: str,
        column_data: Optional[np.ndarray],
    ) -> np.ndarray:
        spec = BREAKTHROUGH_Y_SPECS[y_units]
        if spec.component_col_0based is not None:
            return data[:, spec.component_col_0based].astype(float)

        if y_units != "molefraction":
            raise ValueError(f"Unsupported derived breakthrough y-axis '{y_units}'")
        if column_data is None:
            raise ValueError("column.data is required to calculate mole fractions")
        if column_data.shape[0] != data.shape[0]:
            raise ValueError(
                "Cannot calculate mole fraction: component and column files have different block counts "
                f"({data.shape[0]} vs {column_data.shape[0]})"
            )

        p_i = data[:, COMPONENT_METRICS["P"].col_0based].astype(float)
        p_t = column_data[:, COLUMN_METRICS["Pt"].col_0based].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(p_t != 0.0, p_i / p_t, np.nan)

    def _temperature_values(self, values: np.ndarray, y_units: str) -> np.ndarray:
        spec = TEMPERATURE_Y_SPECS[y_units]
        return values.astype(float) * spec.scale + spec.offset

    def breakthrough(
        self,
        x_units: str = "dimensionless",
        y_units: str = "normalized_concentration",
        include_carrier_gas: bool = True,
        show_markers: bool = True,
        grid_index: int = -1,
    ) -> go.Figure:
        """
        Plot breakthrough curves with explicit axis selections.

        Parameters
        ----------
        x_units:
            One of: "dimensionless", "min", "s".
        y_units:
            One of: "concentration", "normalized_concentration",
            "molefraction", "partial_pressure".
        include_carrier_gas:
            If False, components marked as carrier gas are skipped.
        show_markers:
            If True, traces are rendered as lines plus markers.
        grid_index:
            Grid row to extract from each block. The default, -1, is the outlet.
        """
        x_key = _canonical_key(x_units, X_UNITS, "x_units")
        y_key = _canonical_key(y_units, Y_UNITS, "y_units")
        x_spec = BREAKTHROUGH_X_SPECS[x_key]
        y_spec = BREAKTHROUGH_Y_SPECS[y_key]

        column_data = None

        fig = go.Figure()
        x_values: List[np.ndarray] = []
        y_values: List[np.ndarray] = []
        z_selected: Optional[float] = None

        for comp in sorted(self.components, key=self._component_sort_key):
            if (not include_carrier_gas) and comp.isCarrierGas:
                continue

            fileName = self._component_file_name(comp.index, comp.name)
            data = self._breakthrough_component_data(fileName, grid_index=grid_index)

            x = self._x_values(data, x_key)
            y = self._y_values(data, y_key, column_data)
            x_values.append(x)
            y_values.append(y)

            if z_selected is None and data.shape[0] > 0:
                z_selected = float(data[0, 2])

            line_width = 3.6 if comp.isCarrierGas else 3.2
            label = self._component_label(comp, include_yi=True)
            trace_kwargs = {
                "x": x,
                "y": y,
                "mode": self._trace_mode(show_markers),
                "name": label,
                "line": dict(width=line_width),
                "hovertemplate": (
                    f"{label}<br>"
                    f"{x_spec.title}=%{{x:.4g}}<br>"
                    f"{y_spec.hover_label}=%{{y:.4g}}<extra></extra>"
                ),
            }
            if show_markers:
                trace_kwargs["marker"] = dict(symbol="circle", size=8, line=dict(width=1, color="Black"))

            fig.add_trace(go.Scatter(**trace_kwargs))

        xmin, xmax = _finite_min_max(x_values)
        ymin, ymax = _finite_min_max(y_values)
        x_range = [0.0, 1.0] if xmin is None or xmax is None else [xmin, xmax]
        y_range = _padded_y_range(ymin, ymax, normalized=y_spec.normalized)

        title_text = self._plot_title(y_spec.mode)
        if z_selected is not None:
            title_text += f"  z={z_selected:g} m"

        return self._publication_layout(
            fig,
            xaxis_title=x_spec.title,
            yaxis_title=y_spec.title,
            title_text=title_text,
            x_range=x_range,
            y_range=y_range,
        )

    def temperature_triplet_time(
        self,
        grid_index: Optional[int] = None,
        x_units: str = "min",
        y_units: str = "kelvin",
        show_markers: bool = True,
    ) -> go.Figure:
        """
        Plot gas, solid, and wall temperature histories at a selected grid row.

        Parameters
        ----------
        grid_index:
            Grid row to extract from each block. The default, None, selects
            the middle row of each block.
        x_units:
            One of: "dimensionless", "min", "s", "hr".
        y_units:
            One of: "kelvin", "celsius".
        show_markers:
            If True, traces are rendered as lines plus markers.
        """
        x_key = _canonical_key(x_units, X_UNITS, "x_units")
        y_key = _canonical_key(y_units, TEMPERATURE_Y_UNITS, "y_units")
        x_spec = BREAKTHROUGH_X_SPECS[x_key]
        y_spec = TEMPERATURE_Y_SPECS[y_key]

        data = self._breakthrough_column_data(grid_index=grid_index)

        x = self._x_values(data, x_key)
        tg = self._temperature_values(data[:, COLUMN_METRICS["Tg"].col_0based], y_key)
        ts = self._temperature_values(data[:, COLUMN_METRICS["Ts"].col_0based], y_key)
        tw = self._temperature_values(data[:, COLUMN_METRICS["Tw"].col_0based], y_key)
        z_selected = float(data[0, 2]) if data.shape[0] > 0 else None

        fig = go.Figure()
        y_values = [tg, ts, tw]
        for label, y in [
            ("Gas temperature", tg),
            ("Solid temperature", ts),
            ("Wall temperature", tw),
        ]:
            trace_kwargs = {
                "x": x,
                "y": y,
                "mode": self._trace_mode(show_markers),
                "name": label,
                "line": dict(width=3.2),
                "hovertemplate": (
                    f"{label}<br>"
                    f"{x_spec.title}=%{{x:.4g}}<br>"
                    f"T=%{{y:.4g}} {y_spec.hover_unit}<extra></extra>"
                ),
            }
            if show_markers:
                trace_kwargs["marker"] = dict(symbol="circle", size=8, line=dict(width=1, color="Black"))
            fig.add_trace(go.Scatter(**trace_kwargs))

        xmin, xmax = _finite_min_max([x])
        ymin, ymax = _finite_min_max(y_values)
        title_text = self._plot_title("Temperature history")
        if z_selected is not None:
            title_text += f"  z={z_selected:g} m"

        return self._publication_layout(
            fig,
            xaxis_title=x_spec.title,
            yaxis_title=y_spec.title,
            title_text=title_text,
            x_range=[0.0, 1.0] if xmin is None or xmax is None else [xmin, xmax],
            y_range=_padded_y_range(ymin, ymax),
        )

    def column_snapshot(
        self,
        metric: str,
        frame_index: int = 0,
        show_markers: bool = True,
    ) -> go.Figure:
        info = self._metric_info(metric)
        ylabel = str(info["ylabel"])
        source = str(info["source"])
        col = int(info["col_0based"])
        fig = go.Figure()
        y_values: List[np.ndarray] = []

        if source == "column":
            blocks = self._read_column_data("column.data")
            if frame_index < 0 or frame_index >= len(blocks):
                raise IndexError(f"frame_index {frame_index} out of range [0, {len(blocks) - 1}]")

            block = blocks[frame_index]
            z = block[:, 2]
            y = block[:, col]
            y_values.append(y)
            fig.add_trace(go.Scatter(x=z, y=y, mode=self._trace_mode(show_markers), name=metric))
        else:
            for comp in sorted(self.components, key=self._component_sort_key):
                fileName = self._component_file_name(comp.index, comp.name)
                blocks = self._read_component_blocks(fileName)
                if frame_index < 0 or frame_index >= len(blocks):
                    raise IndexError(f"frame_index {frame_index} out of range for {fileName}")

                block = blocks[frame_index]
                z = block[:, 2]
                y = block[:, col]
                y_values.append(y)
                fig.add_trace(
                    go.Scatter(
                        x=z,
                        y=y,
                        mode=self._trace_mode(show_markers),
                        name=self._component_label(comp, include_yi=True),
                    )
                )

        ymin, ymax = _finite_min_max(y_values)
        fig.update_layout(
            title=f"{self._plot_title()} — {metric} — frame {frame_index}",
            xaxis_title="Adsorber position [m]",
            yaxis_title=ylabel,
            yaxis=dict(range=_padded_y_range(ymin, ymax)),
            legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="center", x=0.5),
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
        col = int(info["col_0based"])
        every = max(1, int(every))
        fig = go.Figure()

        if source == "column":
            blocks = self._read_column_data("column.data")
            if len(blocks) == 0:
                raise ValueError("No blocks found in column.data")

            sampled_indices = list(range(0, len(blocks), every))
            sampled_blocks = [blocks[i] for i in sampled_indices]
            x_values = [block[:, 2] for block in sampled_blocks]
            y_values = [block[:, col] for block in sampled_blocks]
            x_min, x_max = _finite_min_max(x_values)
            y_min, y_max = _finite_min_max(y_values)

            first_block = sampled_blocks[0]
            fig.add_trace(
                go.Scatter(
                    x=first_block[:, 2],
                    y=first_block[:, col],
                    mode=self._trace_mode(show_markers),
                    name=metric,
                )
            )

            frames = [
                go.Frame(
                    name=str(idx),
                    data=[go.Scatter(x=block[:, 2], y=block[:, col], mode=self._trace_mode(show_markers), name=metric)],
                    layout=go.Layout(title_text=f"{self._plot_title()} — {metric} — frame {idx}"),
                )
                for idx, block in zip(sampled_indices, sampled_blocks)
            ]
        else:
            comp_blocks: Dict[int, List[np.ndarray]] = {}
            nframes: Optional[int] = None
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

            sampled_indices = list(range(0, nframes, every))
            x_values = []
            y_values = []
            for comp in self.components:
                for idx in sampled_indices:
                    block = comp_blocks[comp.index][idx]
                    x_values.append(block[:, 2])
                    y_values.append(block[:, col])

            x_min, x_max = _finite_min_max(x_values)
            y_min, y_max = _finite_min_max(y_values)

            first_idx = sampled_indices[0]
            for comp in sorted(self.components, key=self._component_sort_key):
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
                for comp in sorted(self.components, key=self._component_sort_key):
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

        fig.frames = frames
        fig.update_layout(
            title=f"{self._plot_title()} — {metric} — frame {sampled_indices[0]}",
            xaxis_title="Adsorber position [m]",
            yaxis_title=ylabel,
            xaxis=dict(range=[0.0, 1.0] if x_min is None or x_max is None else [x_min, x_max]),
            yaxis=dict(range=_padded_y_range(y_min, y_max)),
            legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="center", x=0.5),
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
        return {metric: self.column_animation(metric, every=every, show_markers=show_markers) for metric in METRIC_SPECS}

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

        return ColumnDataExplorer(
            data_dir=self.data_dir,
            components=self.components,
            component_file_names={comp.index: self._component_file_name(comp.index, comp.name) for comp in self.components},
            column_file=column_file,
            dt=dt,
            time0=time0,
            legend_title=legend_title,
            carrierGasComponent=self.carrierGasComponent,
            show_markers=show_markers,
        )

    def dt_guess(self, dt: Optional[float]) -> Optional[float]:
        return dt


class ColumnDataExplorer:
    """
    Interactive explorer for split column/component data files.
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
            self.component_blocks[comp.index] = _canonicalize_component_blocks(
                read_blocks(fileName), self.column_blocks
            )

        self.dt = dt
        self.time0 = time0
        self.legend_title = legend_title
        self.carrierGasComponent = carrierGasComponent
        self.show_markers = show_markers
        self._syncing_play = False

        self.z_grid = self.column_blocks[0][:, 2].astype(float)
        self.zmin = float(np.nanmin(self.z_grid))
        self.zmax = float(np.nanmax(self.z_grid))

        self.column_props = _props_from_specs(COLUMN_METRICS)
        self.component_props = _props_from_specs(COMPONENT_METRICS)
        self.property_options = _column_explorer_y_options()
        self.available_props = [value for _, value in self.property_options]

        self._build_widgets()
        self._build_figure()
        self._render()

    def _trace_mode(self) -> str:
        return "lines+markers" if self.show_markers else "lines"

    def _build_widgets(self):
        default_prop = "C" if "C" in self.available_props else self.available_props[0]
        self.w_prop = widgets.Dropdown(
            options=self.property_options,
            value=default_prop,
            description="y:",
            layout=widgets.Layout(width="520px"),
            style={"description_width": "24px"},
        )
        self.w_xmode = widgets.Dropdown(
            options=[("grid (z)", "grid"), ("time", "time")],
            value="grid",
            description="x:",
        )
        self.w_block = widgets.IntSlider(
            value=0,
            min=0,
            max=len(self.column_blocks) - 1,
            step=1,
            description="block",
            continuous_update=False,
            layout=widgets.Layout(width="380px"),
        )
        self.w_play = widgets.Play(
            value=0,
            min=0,
            max=len(self.column_blocks) - 1,
            step=1,
            interval=120,
            repeat=True,
            description="play",
            disabled=False,
            layout=widgets.Layout(width="150px"),
        )

        step = (self.zmax - self.zmin) / 200 if self.zmax > self.zmin else 1.0
        if step <= 0.0 or not np.isfinite(step):
            step = 1.0
        self.w_z = widgets.FloatSlider(value=self.zmin, min=self.zmin, max=self.zmax, step=step, description="z@")

        comp_options = [(f"{comp.index}: {comp.name}", comp.index) for comp in self.components]
        self.w_comps = widgets.SelectMultiple(
            options=comp_options,
            value=tuple(ci for _, ci in comp_options),
            description="components",
            rows=min(8, max(3, len(comp_options))),
        )
        self.w_show_carrier = widgets.Checkbox(value=False, description="show carrier gas", indent=False)
        self.w_show_markers = widgets.Checkbox(value=self.show_markers, description="show markers", indent=False)
        self.w_save_path = widgets.Text(value="column_plot.pdf", description="save as", layout=widgets.Layout(width="420px"))
        self.w_save_button = widgets.Button(description="Save PDF", tooltip="Save current figure to PDF", button_style="")
        self.w_save_status = widgets.HTML(value="")

        self.w_save_button.on_click(self._save_current_pdf)

        self.controls_top = widgets.HBox([self.w_xmode, self.w_prop])
        self.controls_x = widgets.HBox([self.w_block, self.w_play, self.w_z])
        self.controls_bottom = widgets.HBox([self.w_show_carrier, self.w_show_markers])
        self.controls_export = widgets.HBox([self.w_save_path, self.w_save_button, self.w_save_status])
        self.controls = widgets.VBox(
            [self.controls_top, self.controls_x, self.controls_bottom, self.w_comps, self.controls_export]
        )

        for w in [self.w_prop, self.w_xmode, self.w_z, self.w_comps, self.w_show_carrier]:
            w.observe(self._on_change, names="value")
        self.w_block.observe(self._on_block_change, names="value")
        self.w_play.observe(self._on_play_change, names="value")
        self.w_show_markers.observe(self._on_markers_change, names="value")
        self._sync_visibility()

    def _build_figure(self):
        self.fig = go.FigureWidget()
        self.fig.update_layout(
            height=520,
            margin=dict(l=70, r=20, t=70, b=110),
            legend=dict(orientation="h", yanchor="top", y=-0.18, xanchor="center", x=0.5),
            template="plotly_white",
        )

    def _sync_visibility(self):
        if self.w_xmode.value == "grid":
            self.w_block.layout.display = ""
            self.w_play.layout.display = ""
            self.w_play.disabled = False
            self.w_z.layout.display = "none"
        else:
            self.w_block.layout.display = "none"
            self.w_play.layout.display = "none"
            self.w_play.disabled = True
            self.w_z.layout.display = ""

        prop = str(self.w_prop.value)
        is_component_prop = prop in self.component_props
        self.w_comps.layout.display = "" if is_component_prop else "none"
        self.w_show_carrier.layout.display = "" if is_component_prop else "none"

    def _on_change(self, _change):
        self._sync_visibility()
        self._render()

    def _on_block_change(self, change):
        if not self._syncing_play:
            self._sync_play_to_block(int(change["new"]))
        self._on_change(change)

    def _on_play_change(self, change):
        if self._syncing_play or self.w_xmode.value != "grid":
            return
        block_index = int(change["new"])
        if self.w_block.value != block_index:
            self.w_block.value = block_index

    def _on_markers_change(self, _change):
        self.show_markers = bool(self.w_show_markers.value)
        self._render()

    def _sync_play_to_block(self, block_index: int) -> None:
        if self.w_play.value == block_index:
            return
        self._syncing_play = True
        try:
            self.w_play.value = block_index
        finally:
            self._syncing_play = False

    def _time_axis(self) -> np.ndarray:
        if self.dt is not None:
            idx = np.arange(len(self.column_blocks), dtype=float)
            return self.time0 + idx * float(self.dt)

        values = []
        for block in self.column_blocks:
            values.append(float(block[0, 1]) if block.shape[0] and block.shape[1] > 1 else np.nan)
        return np.asarray(values, dtype=float)

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

    def _component_by_index(self, ci: int) -> ComponentInfo:
        for comp in self.components:
            if comp.index == ci:
                return comp
        raise KeyError(f"Unknown component index {ci}")

    def _prop_header_label(self, prop: str) -> str:
        if prop in self.column_props:
            return str(self.column_props[prop]["header_label"])
        return str(self.component_props[prop]["header_label"])

    def _render_grid(self, prop: str, block_index: int, comp_ids: List[int]):
        traces: List[go.Scatter] = []
        y_values: List[np.ndarray] = []
        prop_label = self._prop_header_label(prop)

        if prop in self.column_props:
            info = self.column_props[prop]
            block = self.column_blocks[block_index]
            z = block[:, 2].astype(float)
            col = int(info["col_0based"])
            y = block[:, col].astype(float)
            y_values.append(y)
            traces.append(go.Scatter(x=z, y=y, mode=self._trace_mode(), name=prop_label))
            ylabel = str(info["ylabel"])
        else:
            info = self.component_props[prop]
            col = int(info["col_0based"])
            for ci in comp_ids:
                if not self._include_component(ci):
                    continue
                blocks = self.component_blocks[ci]
                if block_index < 0 or block_index >= len(blocks):
                    raise IndexError(f"block_index {block_index} out of range for component {ci}")

                block = blocks[block_index]
                z = block[:, 2].astype(float)
                y = block[:, col].astype(float)
                y_values.append(y)
                comp = self._component_by_index(ci)
                traces.append(go.Scatter(x=z, y=y, mode=self._trace_mode(), name=comp.name))
            ylabel = str(info["ylabel"])

        self._replace_traces(traces)
        title = f"{prop_label} vs z (block {block_index})"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        ymin, ymax = _finite_min_max(y_values)
        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(title="Adsorber position [m]", range=[self.zmin, self.zmax], autorange=False),
                yaxis=dict(title=ylabel, range=_padded_y_range(ymin, ymax)),
            )

    def _render_time(self, prop: str, z_value: float, comp_ids: List[int]):
        t = self._time_axis()
        zi = self._nearest_z_index(z_value)
        z_used = float(self.z_grid[zi])
        traces: List[go.Scatter] = []
        y_values: List[np.ndarray] = []
        prop_label = self._prop_header_label(prop)

        if prop in self.column_props:
            info = self.column_props[prop]
            col = int(info["col_0based"])
            ys = []
            for block in self.column_blocks:
                ys.append(float(block[zi, col]) if zi < block.shape[0] and col < block.shape[1] else np.nan)
            ys = np.asarray(ys, dtype=float)
            y_values.append(ys)
            traces.append(go.Scatter(x=t, y=ys, mode=self._trace_mode(), name=prop_label))
            ylabel = str(info["ylabel"])
        else:
            info = self.component_props[prop]
            col = int(info["col_0based"])
            for ci in comp_ids:
                if not self._include_component(ci):
                    continue
                blocks = self.component_blocks[ci]
                ys = []
                for block in blocks:
                    ys.append(float(block[zi, col]) if zi < block.shape[0] and col < block.shape[1] else np.nan)
                ys = np.asarray(ys, dtype=float)
                y_values.append(ys)
                comp = self._component_by_index(ci)
                x = t[: ys.shape[0]] if t.shape[0] >= ys.shape[0] else np.arange(ys.shape[0], dtype=float)
                traces.append(go.Scatter(x=x, y=ys, mode=self._trace_mode(), name=comp.name))
            ylabel = str(info["ylabel"])

        self._replace_traces(traces)
        xlab = "Time, t [min]" if self.dt is None else "time"
        title = f"{prop_label} vs {xlab} (z≈{z_used:g} m)"
        if self.legend_title:
            title = f"{self.legend_title}<br>{title}"

        ymin, ymax = _finite_min_max(y_values)
        with self.fig.batch_update():
            self.fig.update_layout(
                title=dict(text=title),
                xaxis=dict(title=xlab, autorange=True, range=None),
                yaxis=dict(title=ylabel, range=_padded_y_range(ymin, ymax)),
            )
            self.fig.update_xaxes(autorange=True)

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
