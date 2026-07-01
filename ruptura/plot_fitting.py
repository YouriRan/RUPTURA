from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, List, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pc

from .base_plotly import BasePlotly
from .utils import ComponentInfo, load_simulation_metadata


@dataclass(frozen=True)
class JsonIsothermSite:
    """Single isotherm site read from the fitted Components JSON."""

    type: str
    parameters: tuple[float, ...]

    @staticmethod
    def _type_key(typeName: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", typeName.lower())

    def _require_parameters(self, numberOfParameters: int) -> None:
        if len(self.parameters) < numberOfParameters:
            raise ValueError(
                f"Isotherm site '{self.type}' requires at least {numberOfParameters} parameters; "
                f"found {len(self.parameters)}"
            )

    def value(self, pressure: float) -> float:
        p = float(pressure)
        parameters = self.parameters
        type_key = self._type_key(self.type)

        if type_key == "langmuir":
            self._require_parameters(2)
            temp = parameters[1] * p
            return parameters[0] * temp / (1.0 + temp)

        if type_key == "antilangmuir":
            self._require_parameters(2)
            return parameters[0] * p / (1.0 - parameters[1] * p)

        if type_key == "bet":
            self._require_parameters(3)
            return (
                parameters[0]
                * parameters[1]
                * p
                / ((1.0 - parameters[2] * p) * (1.0 - parameters[2] + parameters[1] * p))
            )

        if type_key == "henry":
            self._require_parameters(1)
            return parameters[0] * p

        if type_key == "freundlich":
            self._require_parameters(2)
            return parameters[0] * math.pow(p, 1.0 / parameters[1])

        if type_key == "sips":
            self._require_parameters(3)
            temp = math.pow(parameters[1] * p, 1.0 / parameters[2])
            return parameters[0] * temp / (1.0 + temp)

        if type_key == "langmuirfreundlich":
            self._require_parameters(3)
            temp = parameters[1] * math.pow(p, parameters[2])
            return parameters[0] * temp / (1.0 + temp)

        if type_key == "redlichpeterson":
            self._require_parameters(3)
            return parameters[0] * p / (1.0 + parameters[1] * math.pow(p, parameters[2]))

        if type_key == "toth":
            self._require_parameters(3)
            temp = parameters[1] * p
            return parameters[0] * temp / math.pow(1.0 + math.pow(temp, parameters[2]), 1.0 / parameters[2])

        if type_key == "unilan":
            self._require_parameters(3)
            temp1 = 1.0 + parameters[1] * math.exp(parameters[2]) * p
            temp2 = 1.0 + parameters[1] * math.exp(-parameters[2]) * p
            return parameters[0] * (0.5 / parameters[2]) * math.log(temp1 / temp2)

        if type_key in {"obrienmyers", "obrianmyers"}:
            self._require_parameters(3)
            temp1 = parameters[1] * p
            temp2 = 1.0 + temp1
            return parameters[0] * (
                temp1 / temp2
                + parameters[2] * parameters[2] * temp1 * (1.0 - temp1) / (temp2 * temp2 * temp2)
            )

        if type_key == "quadratic":
            self._require_parameters(3)
            temp1 = parameters[1] * p
            temp2 = parameters[2] * p * p
            return parameters[0] * (temp1 + 2.0 * temp2) / (1.0 + temp1 + temp2)

        if type_key in {"temkin", "asymptotictemkin"}:
            self._require_parameters(3)
            temp = parameters[1] * p
            temp1 = temp / (1.0 + temp)
            return parameters[0] * (temp1 + parameters[2] * temp1 * temp1 * (temp1 - 1.0))

        if type_key == "bingelwalton":
            self._require_parameters(3)
            return parameters[0] * (1.0 - math.exp(-(parameters[1] + parameters[2]) * p)) / (
                1.0 + (parameters[2] / parameters[1]) * math.exp(-(parameters[1] + parameters[2]) * p)
            )

        raise ValueError(f"Unsupported isotherm type in fitted JSON: {self.type}")


@dataclass(frozen=True)
class JsonMultiSiteIsotherm:
    """Multi-site isotherm read from the fitted Components JSON."""

    sites: tuple[JsonIsothermSite, ...]
    isCarrierGas: bool = False

    @classmethod
    def from_component_json(cls, componentJson: Mapping[str, Any]) -> JsonMultiSiteIsotherm:
        isCarrierGas = bool(componentJson.get("CarrierGas", False))
        sitesJson = componentJson.get("IsothermSites", [])

        if sitesJson is None:
            sitesJson = []
        if not isinstance(sitesJson, list):
            raise ValueError("Component field 'IsothermSites' must be a list")

        sites: list[JsonIsothermSite] = []
        for siteJson in sitesJson:
            if not isinstance(siteJson, Mapping):
                raise ValueError("Each isotherm site in 'IsothermSites' must be an object")

            siteType = siteJson.get("Type")
            parametersJson = siteJson.get("Parameters")
            if not isinstance(siteType, str):
                raise ValueError("Each isotherm site requires a string 'Type'")
            if not isinstance(parametersJson, list):
                raise ValueError(f"Isotherm site '{siteType}' requires a list field 'Parameters'")

            sites.append(JsonIsothermSite(siteType, tuple(float(value) for value in parametersJson)))

        return cls(tuple(sites), isCarrierGas=isCarrierGas)

    def value(self, pressure: float) -> float:
        if self.isCarrierGas or not self.sites:
            raise ValueError("Carrier-gas components or components without isotherm sites cannot be evaluated")

        return float(sum(site.value(pressure) for site in self.sites))


class FittingPlotly(BasePlotly):
    """
    Plotly interface for the C++ fitting workflow.

    This mirrors the old gnuplot output generated by Fitting::createPlotScripts(...):
      - starting isotherm, labelled "start f(x)"
      - optimized isotherm, labelled "fit f(x)"
      - raw fitting data, labelled "raw data"

    The optimized isotherm is read from the fitted Components JSON emitted by the
    C++ fitting code after the Nelder-Mead step.
    """

    def __init__(
        self,
        displayName: str,
        components: List[ComponentInfo],
        start_isotherms: Sequence[object],
        fit_json: Union[str, Path],
        filenames: Sequence[Union[str, Path]],
        data_dir: Union[str, Path] = ".",
        pressureScale: str = "log",
        columnPressure: int = 1,
        columnLoading: int = 2,
        columnError: Optional[int] = None,
        carrierGasComponent: Optional[int] = None,
        showMarkers: bool = False,
    ) -> None:
        super().__init__(
            displayName=displayName,
            components=components,
            data_dir=data_dir,
            carrierGasComponent=carrierGasComponent,
            showMarkers=showMarkers,
        )

        if len(start_isotherms) != self.numberOfComponents:
            raise ValueError("start_isotherms length must match number of components")
        if len(filenames) != self.numberOfComponents:
            raise ValueError("filenames length must match number of components")

        self.start_isotherms = list(start_isotherms)
        self.filenames = [self.data_dir / Path(fileName) for fileName in filenames]
        self.fit_json = self.data_dir / Path(fit_json)
        self.fit_isotherms = self._read_fit_isotherms_json(self.fit_json)
        self.pressureScale = self._normalise_pressure_scale(pressureScale)

#Match C++ InputReader semantics : incoming columns are 1 - based.
        self.columnPressure = int(columnPressure) - 1
        self.columnLoading = int(columnLoading) - 1
        self.columnError = None if columnError is None else int(columnError) - 1


    @staticmethod
    def _normalise_pressure_scale(value: object) -> str:
        scale = str(value).strip().lower()
        if scale in {"0", "log", "logarithmic"}:
            return "log"
        if scale in {"1", "normal", "linear", "lin"}:
            return "normal"
        return scale

    @staticmethod
    def _component_filename_from_json(componentJson: Mapping[str, Any], fallback: str) -> str:
        for key in (
            "Filename",
            "FileName",
            "FittingFilename",
            "FittingFileName",
            "FittingDataFilename",
            "FittingDataFileName",
            "DataFilename",
            "DataFileName",
            "IsothermFilename",
            "IsothermFileName",
        ):
            value = componentJson.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()

        return fallback

    @classmethod
    def from_simulation_json(
        cls,
        simulation_json: Union[str, Path],
        fit_json: Optional[Union[str, Path]] = None,
        data_dir: Union[str, Path] = ".",
        carrierGasComponent: Optional[int] = None,
        columnPressure: Optional[int] = None,
        columnLoading: Optional[int] = None,
        columnError: Optional[int] = None,
        showMarkers: bool = False,
    ) -> "FittingPlotly":
        simulation_json_path = Path(simulation_json)
        with simulation_json_path.open("r", encoding="utf-8") as stream:
            simulationJson = json.load(stream)

        if not isinstance(simulationJson, Mapping):
            raise ValueError(f"{simulation_json_path} must contain a JSON object")

        meta = load_simulation_metadata(simulation_json_path)
        components = list(meta["components"])
        componentsJson = simulationJson.get("Components", [])
        if not isinstance(componentsJson, list):
            raise ValueError(f"{simulation_json_path} must contain a list field named 'Components'")
        if len(componentsJson) != len(components):
            raise ValueError(
                f"{simulation_json_path} contains {len(componentsJson)} component definitions, "
                f"but load_simulation_metadata returned {len(components)} components"
            )

        if carrierGasComponent is None:
            for componentIndex, componentJson in enumerate(componentsJson):
                if isinstance(componentJson, Mapping) and bool(componentJson.get("CarrierGas", False)):
                    carrierGasComponent = componentIndex
                    break

        if carrierGasComponent is not None:
            for comp in components:
                comp.isCarrierGas = comp.index == carrierGasComponent

        start_isotherms: list[JsonMultiSiteIsotherm] = []
        filenames: list[str] = []
        for componentIndex, componentJson in enumerate(componentsJson):
            if not isinstance(componentJson, Mapping):
                raise ValueError("Each item in 'Components' must be an object")

            start_isotherms.append(JsonMultiSiteIsotherm.from_component_json(componentJson))

            componentName = components[componentIndex].name
            fallbackFileName = f"component_{componentIndex}_{componentName}.data"
            filenames.append(cls._component_filename_from_json(componentJson, fallback=fallbackFileName))

        if fit_json is None:
            displayName = str(meta["displayName"])
            fit_json = f"{displayName}_fitted_components.json" if displayName else "fitted_components.json"

        pressureScale = cls._normalise_pressure_scale(meta.get("pressureScale", simulationJson.get("PressureScale", "log")))
        resolvedColumnPressure = int(
            columnPressure if columnPressure is not None else simulationJson.get("ColumnPressure", 1)
        )
        resolvedColumnLoading = int(
            columnLoading if columnLoading is not None else simulationJson.get("ColumnLoading", 2)
        )
        resolvedColumnError = columnError
        if resolvedColumnError is None and "ColumnError" in simulationJson:
            resolvedColumnError = int(simulationJson["ColumnError"])

        return cls(
            displayName=str(meta["displayName"]),
            components=components,
            start_isotherms=start_isotherms,
            fit_json=fit_json,
            filenames=filenames,
            data_dir=data_dir,
            pressureScale=pressureScale,
            columnPressure=resolvedColumnPressure,
            columnLoading=resolvedColumnLoading,
            columnError=resolvedColumnError,
            carrierGasComponent=carrierGasComponent,
            showMarkers=showMarkers,
        )

    def _plot_title(self, mode: str = "Fitted isotherm") -> str:
        return f"{self.displayName} — {mode}"

    def _read_fit_isotherms_json(self, fit_json: Union[str, Path]) -> list[JsonMultiSiteIsotherm]:
        fit_json_path = Path(fit_json)
        with fit_json_path.open("r", encoding="utf-8") as stream:
            payload = json.load(stream)

        if isinstance(payload, list):
            componentsJson = payload
        elif isinstance(payload, Mapping):
            componentsJson = payload.get("Components")
        else:
            raise ValueError(f"{fit_json_path} must contain either a JSON object or a component list")

        if not isinstance(componentsJson, list):
            raise ValueError(f"{fit_json_path} must contain a list field named 'Components'")

        componentJsonByName: dict[str, Mapping[str, Any]] = {}
        for componentJson in componentsJson:
            if not isinstance(componentJson, Mapping):
                raise ValueError("Each item in 'Components' must be an object")

            name = componentJson.get("Name")
            if not isinstance(name, str):
                raise ValueError("Each item in 'Components' requires a string field 'Name'")
            if name in componentJsonByName:
                raise ValueError(f"Duplicate component name in fitted JSON: {name}")

            componentJsonByName[name] = componentJson

        fit_isotherms: list[JsonMultiSiteIsotherm] = []
        for comp in self.components:
            componentJson = componentJsonByName.get(comp.name)
            if componentJson is None:
                raise ValueError(f"Component '{comp.name}' was not found in fitted JSON: {fit_json_path}")

            fit_isotherm = JsonMultiSiteIsotherm.from_component_json(componentJson)
            if fit_isotherm.isCarrierGas:
                comp.isCarrierGas = True
            fit_isotherms.append(fit_isotherm)

        return fit_isotherms

    def _read_fit_data(self, fileName: Union[str, Path]) -> pd.DataFrame:
        required_cols = max(self.columnPressure, self.columnLoading)
        if self.columnError is not None:
            required_cols = max(required_cols, self.columnError)

        df = self._read_component_data(fileName=fileName, min_columns=required_cols + 1)

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

    def _is_plottable_component(self, comp_index: int) -> bool:
        if self.carrierGasComponent is not None and comp_index == self.carrierGasComponent:
            return False

        fit_isotherm = self.fit_isotherms[comp_index]
        if getattr(fit_isotherm, "isCarrierGas", False):
            return False
        if hasattr(fit_isotherm, "sites") and len(fit_isotherm.sites) == 0:
            return False

        return True

    def fit_figure(self, comp_index: int) -> go.Figure:
        if comp_index < 0 or comp_index >= self.numberOfComponents:
            raise IndexError(f"comp_index {comp_index} out of range [0, {self.numberOfComponents - 1}]")
        if not self._is_plottable_component(comp_index):
            raise ValueError(f"Component {comp_index} is a carrier gas or has no fitted isotherm sites")

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

#fig.add_trace(
#go.Scatter(
#x = pressure_grid,
#y = fit_loading,
#mode = "lines",
#name = "fit f(x)",
#line = dict(width = 3.0),
#hovertemplate = ("fit f(x)<br>" + "P=%{x:.4g} Pa<br>" + "q=%{y:.6g} mol/kg<extra></extra>"),
#)
#)

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

    def all_fit_figure(self, show_fit=True) -> go.Figure:
        plottableComponents = [
            comp
            for comp in sorted(self.components, key=self._component_sort_key)
            if self._is_plottable_component(comp.index)
        ]

        if not plottableComponents:
            raise ValueError("No plottable fitted components were found")

        colorway = pc.qualitative.Plotly
        fig = go.Figure()
        ymax = 0.0

        for colorIndex, comp in enumerate(plottableComponents):
            color = colorway[colorIndex % len(colorway)]

            df = self._read_fit_data(self.filenames[comp.index])

            pressure_data = df["pressure"].to_numpy(dtype=float)
            loading_data = df["loading"].to_numpy(dtype=float)
            pressure_grid = self._pressure_grid(pressure_data)

            start_loading = self._evaluate_isotherm(self.start_isotherms[comp.index], pressure_grid)
            fit_loading = self._evaluate_isotherm(self.fit_isotherms[comp.index], pressure_grid)

            for values in (loading_data, start_loading, fit_loading):
                if np.isfinite(values).any():
                    ymax = max(ymax, float(np.nanmax(values)))

            if show_fit:
                fig.add_trace(
                    go.Scatter(
                        x=pressure_grid,
                        y=start_loading,
                        mode="lines",
                        name=f"{comp.name} start f(x)",
                        legendgroup=comp.name,
                        line=dict(width=2.4, dash="dash", color=color),
                        hovertemplate=(
                            f"{comp.name} start f(x)<br>"
                            + "P=%{x:.4g} Pa<br>"
                            + "q=%{y:.6g} mol/kg<extra></extra>"
                        ),
                    )
                )

                fig.add_trace(
                    go.Scatter(
                        x=pressure_grid,
                        y=fit_loading,
                        mode="lines",
                        name=f"{comp.name} fit f(x)",
                        legendgroup=comp.name,
                        line=dict(width=3.0, color=color),
                        hovertemplate=(
                            f"{comp.name} fit f(x)<br>"
                            + "P=%{x:.4g} Pa<br>"
                            + "q=%{y:.6g} mol/kg<extra></extra>"
                        ),
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
                        color=color,
                    )

            fig.add_trace(
                go.Scatter(
                    x=pressure_data,
                    y=loading_data,
                    mode="markers",
                    name=f"{comp.name} raw data",
                    legendgroup=comp.name,
                    marker=dict(size=7, color=color),
                    hovertemplate=(
                        f"{comp.name} raw data<br>"
                        + "P=%{x:.4g} Pa<br>"
                        + "q=%{y:.6g} mol/kg<extra></extra>"
                    ),
                    **scatter_kwargs,
                )
            )

        y_upper = 1.08 * ymax if ymax > 0.0 else 1.0

        fig = self._publication_layout(
            fig,
            xaxis_title="Pressure, p [Pa]",
            yaxis_title="Absolute loading, q [mol/kg]",
            title_text=self._plot_title("Fitted isotherms"),
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
                title=dict(text="Components"),
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