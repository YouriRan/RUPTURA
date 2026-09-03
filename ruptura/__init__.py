"""Python library interface for Ruptura.

The package exposes the C++ simulation objects directly and adds a small
Python-friendly layer for loading JSON inputs and working with computed data as
NumPy arrays.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

try:
    from . import _ruptura as _core
except ImportError as _relative_import_error:
    try:
        import _ruptura as _core  # type: ignore[no-redef]
    except ImportError:
        _core = None
        _extension_import_error: Optional[ImportError] = _relative_import_error
    else:
        _extension_import_error = None
else:
    _extension_import_error = None


_CORE_EXPORTS = [
    "Breakthrough",
    "BreakthroughIntegrationScheme",
    "CVODE",
    "Chemisorption",
    "Column",
    "Component",
    "Fitting",
    "InputReader",
    "Isotherm",
    "MixturePrediction",
    "MultibedBreakthrough",
    "MultibedColumn",
    "MultibedSwingAdsorption",
    "MultiSiteChemisorption",
    "MultiSiteIsotherm",
    "RungeKutta3",
    "SemiImplicitRungeKutta3",
    "SwingAdsorptionSubStage",
    "SwingAdsorption",
    "compute_equilibrium_loadings",
    "compute_first_derivatives",
    "compute_weno",
    "read_input",
    "simulation_from_file",
    "update_velocity_and_pressure",
]

if _core is not None:
    for _name in _CORE_EXPORTS:
        if hasattr(_core, _name):
            globals()[_name] = getattr(_core, _name)


MIXTURE_RESULT_COLUMNS = (
    "pressure",
    "pure_component_loading",
    "mixture_loading",
    "gas_mol_fraction",
    "adsorbed_mol_fraction",
    "reduced_grand_potential",
)

MIXTURE_SORPTION_RESULT_COLUMNS = (
    "pressure",
    "pure_physisorption_loading",
    "pure_chemisorption_loading",
    "pure_total_loading",
    "mixture_physisorption_loading",
    "mixture_chemisorption_loading",
    "mixture_total_loading",
)

# Columns of the ``component_<index>_<name>.data`` files written by
# ``MixturePrediction::run()``. The last four are appended only when a component
# defines a nested chemisorption equilibrium isotherm; in that case columns 2 and
# 3 hold the physisorption contribution alone.
MIXTURE_FILE_COLUMNS = (
    "pressure",
    "pure_physisorption_loading",
    "mixture_physisorption_loading",
    "gas_mol_fraction",
    "adsorbed_mol_fraction",
    "hypothetical_pressure",
    "reduced_grand_potential",
)

MIXTURE_FILE_CHEMISORPTION_COLUMNS = (
    "pure_chemisorption_loading",
    "mixture_chemisorption_loading",
    "pure_total_loading",
    "mixture_total_loading",
)

BREAKTHROUGH_BASE_COLUMNS = (
    "dimensionless_time",
    "time_min",
    "position",
    "interstitial_gas_velocity",
    "total_pressure",
)

BREAKTHROUGH_COMPONENT_COLUMNS = (
    "physisorption",
    "equilibrium_adsorption",
    "partial_pressure",
    "normalized_partial_pressure",
    "physisorption_dot",
)


@dataclass(frozen=True)
class SimulationResult:
    """Computed simulation data plus labels for the final array axis."""

    data: np.ndarray
    columns: tuple[str, ...]
    kind: str

    @property
    def shape(self) -> tuple[int, ...]:
        return self.data.shape

    def __array__(self, dtype: Any = None) -> np.ndarray:
        return np.asarray(self.data, dtype=dtype)

    def column_index(self, name: str) -> int:
        return self.columns.index(name)

    def column(self, name: str) -> np.ndarray:
        return self.data[..., self.column_index(name)]


def _require_core() -> Any:
    if _core is None:
        raise ImportError(
            "The Ruptura extension module '_ruptura' is not available. "
            "Build the Python extension with CMake using BUILD_PYTHON=ON."
        ) from _extension_import_error
    return _core


def load_input(file_name: str | Path) -> Any:
    """Parse a Ruptura JSON input file and return an InputReader."""

    core = _require_core()
    return core.InputReader(str(file_name))


def load_simulation(file_name: str | Path) -> Any:
    """Construct the simulation object described by a Ruptura JSON input file."""

    core = _require_core()
    return core.simulation_from_file(str(file_name))


def breakthrough_columns(number_of_components: int) -> tuple[str, ...]:
    columns = list(BREAKTHROUGH_BASE_COLUMNS)
    for component in range(number_of_components):
        columns.extend(
            f"component_{component}_{field}"
            for field in BREAKTHROUGH_COMPONENT_COLUMNS
        )
    return tuple(columns)


def result_columns(simulation: Any) -> tuple[str, ...]:
    core = _require_core()
    if isinstance(simulation, core.MixturePrediction):
        return MIXTURE_RESULT_COLUMNS
    if isinstance(simulation, (core.Breakthrough, core.MultibedBreakthrough)):
        return breakthrough_columns(int(simulation.number_of_components))
    raise TypeError(
        "Only MixturePrediction and Breakthrough have in-memory NumPy compute results."
    )


def compute(simulation: Any, *, labeled: bool = True) -> np.ndarray | SimulationResult:
    """Compute a simulation and return NumPy data.

    By default this returns a SimulationResult, which NumPy treats as an array
    via ``np.asarray(result)``. Pass ``labeled=False`` to get the raw ndarray.
    """

    data = np.asarray(simulation.compute())
    if not labeled:
        return data

    core = _require_core()
    if isinstance(simulation, core.MixturePrediction):
        return SimulationResult(data=data, columns=MIXTURE_RESULT_COLUMNS, kind="mixture_prediction")
    if isinstance(simulation, (core.Breakthrough, core.MultibedBreakthrough)):
        return SimulationResult(
            data=data,
            columns=breakthrough_columns(int(simulation.number_of_components)),
            kind="breakthrough",
        )
    return SimulationResult(data=data, columns=(), kind=type(simulation).__name__)


def compute_sorption(simulation: Any, *, labeled: bool = True) -> np.ndarray | SimulationResult:
    """Compute physical, chemical, and total equilibrium mixture loadings.

    ``simulation`` must be a :class:`MixturePrediction`. The returned array has
    pressure points on axis 0, components on axis 1, and the columns listed in
    :data:`MIXTURE_SORPTION_RESULT_COLUMNS` on axis 2.
    """

    core = _require_core()
    if not isinstance(simulation, core.MixturePrediction):
        raise TypeError("compute_sorption requires a MixturePrediction")
    data = np.asarray(simulation.compute_sorption())
    if not labeled:
        return data
    return SimulationResult(
        data=data,
        columns=MIXTURE_SORPTION_RESULT_COLUMNS,
        kind="mixture_sorption_prediction",
    )


def read_mixture_prediction(file_name: str | Path) -> SimulationResult:
    """Read one MixturePrediction ``component_<index>_<name>.data`` output file.

    The returned SimulationResult has pressure points on axis 0 and the columns
    listed in :data:`MIXTURE_FILE_COLUMNS`, extended by
    :data:`MIXTURE_FILE_CHEMISORPTION_COLUMNS` when the simulation included a
    chemisorption equilibrium isotherm.
    """

    data = np.loadtxt(file_name, comments="#", ndmin=2)
    physical_width = len(MIXTURE_FILE_COLUMNS)
    total_width = physical_width + len(MIXTURE_FILE_CHEMISORPTION_COLUMNS)
    if data.shape[1] == total_width:
        columns = MIXTURE_FILE_COLUMNS + MIXTURE_FILE_CHEMISORPTION_COLUMNS
    elif data.shape[1] == physical_width:
        columns = MIXTURE_FILE_COLUMNS
    else:
        raise ValueError(
            f"{file_name} has {data.shape[1]} columns; expected "
            f"{physical_width} or {total_width} from MixturePrediction"
        )
    return SimulationResult(data=data, columns=columns, kind="mixture_prediction_file")


def run(file_name: str | Path, *, compute_result: bool = True) -> Any:
    """Start a simulation from Python using a JSON input file.

    If the object has an in-memory ``compute`` method, the default return value
    is a labeled SimulationResult. Set ``compute_result=False`` to call the
    object's file-writing ``run`` method and return the simulation object.
    """

    simulation = load_simulation(file_name)
    if compute_result and hasattr(simulation, "compute"):
        return compute(simulation)
    simulation.run()
    return simulation


def asarray(value: Any, dtype: Any = None) -> np.ndarray:
    """Convert Ruptura results or vector-like fields to a NumPy array."""

    if isinstance(value, SimulationResult):
        return np.asarray(value.data, dtype=dtype)
    return np.asarray(value, dtype=dtype)


def __getattr__(name: str) -> Any:
    if name in _CORE_EXPORTS:
        core = _require_core()
        return getattr(core, name)
    raise AttributeError(f"module 'ruptura' has no attribute {name!r}")


__all__ = [
    *_CORE_EXPORTS,
    "BREAKTHROUGH_BASE_COLUMNS",
    "BREAKTHROUGH_COMPONENT_COLUMNS",
    "MIXTURE_FILE_CHEMISORPTION_COLUMNS",
    "MIXTURE_FILE_COLUMNS",
    "MIXTURE_RESULT_COLUMNS",
    "MIXTURE_SORPTION_RESULT_COLUMNS",
    "SimulationResult",
    "asarray",
    "breakthrough_columns",
    "compute",
    "compute_sorption",
    "load_input",
    "load_simulation",
    "read_mixture_prediction",
    "result_columns",
    "run",
]
