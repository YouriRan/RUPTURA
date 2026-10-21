"""Notebook Plotly figures plus a bounded-memory browser column explorer."""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import replace
import io
import json
from pathlib import Path
import threading

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from ruptura.plot_breakthrough import BreakthroughPlotly, COLUMN_METRICS, COMPONENT_METRICS, MetricSpec
from ruptura.plot_mixture import MixturePredictionPlotly

_CACHE = OrderedDict()
_LOCK = threading.RLock()


def _blocks(handle):
    pending, offset = b"", 0
    while True:
        chunk = handle.read(1 << 20)
        if not chunk:
            break
        pending += chunk
        start = 0
        while True:
            end = pending.find(b"\n\n", start)
            if end < 0:
                break
            yield offset + start, pending[start:end]
            start = end + 2
        pending = pending[start:]
        offset += start
    if pending.strip():
        yield offset, pending


class BlockIndex:
    """Keep offsets for at most 2,000 complete frames, including the tail."""
    def __init__(self, path):
        self.path = path
        self.frames = []
        self.shape = None
        self.total = 0
        stride = 1
        last = None
        with path.open("rb") as handle:
            for offset, block in _blocks(handle):
                lines = [line for line in block.splitlines() if line.strip() and not line.lstrip().startswith(b"#")]
                if not lines:
                    continue
                columns = len(lines[0].split())
                shape = (len(lines), columns)
                if self.shape is None:
                    self.shape = shape
                if shape != self.shape or any(len(line.split()) != columns for line in lines):
                    continue
                entry = (self.total, offset, len(block))
                last = entry
                if self.total % stride == 0:
                    self.frames.append(entry)
                    if len(self.frames) >= 2000:
                        self.frames = self.frames[::2]
                        stride *= 2
                self.total += 1
        if last and (not self.frames or self.frames[-1] != last):
            self.frames.append(last)
        if not self.frames:
            raise ValueError(f"No complete result frames in {path.name} yet.")

    def read(self, frame):
        if frame < 0:
            frame += len(self.frames)
        if frame < 0 or frame >= len(self.frames):
            raise ValueError("The selected frame is outside the available results.")
        _, offset, length = self.frames[frame]
        with self.path.open("rb") as handle:
            handle.seek(offset)
            return np.loadtxt(io.BytesIO(handle.read(length)), comments="#", ndmin=2)


def index(path):
    path = Path(path)
    stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns)
    with _LOCK:
        if key in _CACHE:
            _CACHE.move_to_end(key)
            return _CACHE[key]
        value = BlockIndex(path)
        _CACHE[key] = value
        while len(_CACHE) > 12:
            _CACHE.popitem(last=False)
        return value


def metrics(config):
    result = {**COLUMN_METRICS, **COMPONENT_METRICS}
    result.pop("Dcdt", None)
    if config.get("FluidPhase", "Gas") == "Liquid":
        for key in ("V", "Tg", "dTgdt", "rho"):
            spec = result[key]
            result[key] = replace(spec, ylabel=spec.ylabel.replace("Gas", "Liquid").replace("gas", "liquid"), header_label=None)
        result["pH"] = MetricSpec("pH [-]", "column", 12)
        for key in ("Y", "P"):
            result.pop(key, None)
        result["Pnorm"] = replace(result["Pnorm"], ylabel="Normalized concentration, cᵢ/cᵢ,feed [-]", header_label=None)
    result["Qtotal"] = MetricSpec("Total sorbed loading [mol/kg]", "component", 6)
    return result


def metadata(directory):
    directory = Path(directory)
    config = json.loads((directory / "simulation.json").read_text())
    components = [{"index": i, "name": component["Name"], "carrier": bool(component.get("CarrierGas")),
                   "adsorbing": bool(component.get("PhysisorptionSites") or component.get("ChemisorptionSites"))}
                  for i, component in enumerate(config["Components"])]
    if config.get("SimulationType") == "MixturePrediction":
        plotter = MixturePredictionPlotly.from_simulation_json(directory / "simulation.json", directory)
        figures = [{"id": "pure", "label": "Pure-component isotherms"}, {"id": "loading", "label": "Mixture loading"}, {"id": "adsorbed", "label": "Adsorbed mole fractions"}]
        if plotter.has_chemisorption():
            figures += [{"id": "chemisorption", "label": "Chemisorption loading"}, {"id": "total", "label": "Total sorbed loading"}, {"id": "puretotal", "label": "Pure total loading"}]
        return {"kind": "MixturePrediction", "components": components, "plots": figures, "frames": 0}
    data = index(directory / "column.data")
    first = data.read(0)
    return {"kind": config.get("SimulationType", "Breakthrough"), "components": components,
            "energyBalance": bool(config.get("energyBalance", False)),
            "plots": [{"id": "breakthrough", "label": "Breakthrough curves"}, {"id": "temperature", "label": "Temperature history"},
                      {"id": "profile", "label": "Column profiles"}],
            "metrics": [{"id": key, "label": spec.ylabel, "source": spec.source} for key, spec in metrics(config).items()
                        if spec.source != "column" or spec.col_0based < first.shape[1]],
            "frames": len(data.frames), "totalFrames": data.total, "sampled": len(data.frames) < data.total,
            "times": [float(data.read(i)[0, 1]) for i in range(len(data.frames))], "positions": first[:, 2].tolist(),
            "liquid": config.get("FluidPhase", "Gas") == "Liquid"}


def _phase_lines(figure, config, units):
    if config.get("SimulationType") != "SwingAdsorption" or units == "dimensionless":
        return
    time = 0.0
    factor = 1 if units == "s" else 3600 if units == "hr" else 60
    for phase in config.get("SwingAdsorptionPhases", [])[:-1]:
        time += phase["NumberOfTimeSteps"] * config.get("TimeStep", 0.0005)
        figure.add_vline(x=time / factor, line_dash="dot", line_color="#9da8a4")


def figure(directory, options=None):
    directory = Path(directory)
    options = options or {}
    config = json.loads((directory / "simulation.json").read_text())
    mode = options.get("plot", "loading" if config.get("SimulationType") == "MixturePrediction" else "breakthrough")
    if mode == "profile":
        axis = options.get("xAxis", "position")
        if axis not in ("time", "position"):
            raise ValueError("Choose time or position for the column x-axis.")
        if axis == "time":
            mode = "history"
    selected = options.get("components")
    # An explicit component selection is the sole visibility control in the Lab.
    show_carrier = selected is not None or bool(options.get("carrier", False))
    markers = bool(options.get("markers", False))
    if config.get("SimulationType") == "MixturePrediction":
        plotter = MixturePredictionPlotly.from_simulation_json(directory / "simulation.json", directory)
        functions = {"pure": plotter.pure_components, "loading": plotter.mixture_loading, "adsorbed": plotter.mixture_adsorbed_molfractions,
                     "chemisorption": plotter.mixture_chemisorption_loading, "total": plotter.mixture_total_loading, "puretotal": plotter.pure_components_total_loading}
        if mode not in functions:
            raise ValueError("Unknown mixture plot.")
        if selected is not None:
            plotter.components = [item for item in plotter.components if item.index in selected]
        result = functions[mode](include_carrier_gas=show_carrier, show_markers=markers)
    else:
        plotter = BreakthroughPlotly.from_simulation_json(directory / "simulation.json", directory)
        if selected is not None:
            plotter.components = [item for item in plotter.components if item.index in selected]
        units = options.get("xUnits", "min")
        grid = int(options.get("grid", -1))
        if mode == "breakthrough":
            y_units = options.get("yUnits", "normalized_concentration")
            if config.get("FluidPhase") == "Liquid" and y_units in ("partial_pressure", "molefraction"):
                raise ValueError("Select concentration or normalized concentration for liquid simulations.")
            result = plotter.breakthrough(x_units=units, y_units=y_units, include_carrier_gas=show_carrier, show_markers=markers, grid_index=grid, max_points=600)
        elif mode == "temperature":
            result = plotter.temperature_triplet_time(x_units=units, y_units=options.get("temperatureUnits", "kelvin"), show_markers=markers, grid_index=grid, max_points=600)
            if config.get("FluidPhase") == "Liquid":
                result.data[0].name = "Liquid temperature"
        elif mode in ("profile", "history"):
            specs = metrics(config)
            key = options.get("metric", "C")
            if key not in specs:
                raise ValueError("Unknown column metric.")
            spec = specs[key]
            result = go.Figure()
            columns = index(directory / "column.data")
            frame = int(options.get("frame", len(columns.frames) - 1))
            sources = [(None, directory / "column.data")] if spec.source == "column" else [
                (component, directory / f"component_{component.index}_{component.name}.data") for component in plotter.components
                if show_carrier or not component.isCarrierGas]
            for component, path in sources:
                blocks = index(path)
                if blocks.shape[1] <= spec.col_0based:
                    raise ValueError("The selected metric is absent from these result files.")
                if mode == "profile":
                    block = blocks.read(frame)
                    x = block[:, 2]
                    y = block[:, spec.col_0based]
                    if key == "Qtotal":
                        y = y + block[:, 8]
                else:
                    rows = []
                    for i in np.linspace(0, len(blocks.frames) - 1, min(600, len(blocks.frames)), dtype=int):
                        block = blocks.read(int(i))
                        position = grid if grid >= 0 else block.shape[0] + grid
                        if not 0 <= position < block.shape[0]:
                            raise ValueError("The selected position is outside the column grid.")
                        rows.append(block[position])
                    rows = np.asarray(rows)
                    x = rows[:, 1] * (60 if units == "s" else 1 / 60 if units == "hr" else 1)
                    y = rows[:, spec.col_0based]
                    if key == "Qtotal":
                        y = y + rows[:, 8]
                result.add_trace(go.Scatter(x=x, y=y, name=component.name if component else spec.ylabel,
                                           mode="lines+markers" if markers else "lines", line={"width": 2.5}))
            result.update_layout(title=f"{config.get('DisplayName', 'Column')} · {spec.ylabel}",
                                 xaxis_title="Position [m]" if mode == "profile" else f"Time [{units}]", yaxis_title=spec.ylabel)
        else:
            raise ValueError("Unknown plot.")
        if mode in ("breakthrough", "temperature", "history"):
            _phase_lines(result, config, units)
    light = options.get("theme") == "light"
    ink, muted, panel, inset, grid = ("#182e3b", "#4d6575", "#ffffff", "#f4f7f9", "#dce5eb") if light else ("#e6edf3", "#a8b9c7", "#171f29", "#111a23", "#273743")
    result.update_layout(template="plotly_white" if light else "plotly_dark", autosize=True, width=None, height=None,
                         font={"family": "Arial, sans-serif", "size": 14, "color": ink},
                         title={"font": {"size": 15}},
                         margin={"l": 76, "r": 28, "t": 70, "b": 100}, paper_bgcolor=panel, plot_bgcolor=inset,
                         legend={"bgcolor": panel, "bordercolor": grid, "font": {"color": ink}},
                         hoverlabel={"bgcolor": panel, "font": {"color": ink}},
                         colorway=["#83d4c1", "#f3a59c", "#8cbcf3", "#c9b0ef", "#e9cb8a", "#82ccd7", "#dbb6a0"])
    if light:
        result.update_layout(colorway=["#176e60", "#b74b40", "#326aaf", "#7651a4", "#907014", "#267b89", "#875d38"])
    result.update_xaxes(color=muted, linecolor=muted, gridcolor=grid, zerolinecolor=grid)
    result.update_yaxes(color=muted, linecolor=muted, gridcolor=grid, zerolinecolor=grid)
    return result


def json_figure(directory, options):
    return json.loads(pio.to_json(figure(directory, options), engine="json"))


def report(directory):
    directory = Path(directory)
    meta = metadata(directory)
    plots = [item for item in meta["plots"] if item["id"] not in ("profile", "history")]
    snippets = []
    for item in plots:
        chart = figure(directory, {"plot": item["id"]})
        snippets.append(chart.to_html(full_html=False, include_plotlyjs=not snippets, config={"responsive": True}))
    path = directory / "plots.html"
    path.write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
                    '<meta name="color-scheme" content="dark"><title>Ruptura · Results</title><style>body{max-width:1200px;margin:30px auto;background:#10171e;font:16px Arial;color:#e6edf3}h1{padding:0 30px}.plotly-graph-div{min-height:550px}</style>'
                    '<h1>Ruptura · Simulation results</h1>' + "\n".join(snippets) + "</html>", encoding="utf-8")
    return path
