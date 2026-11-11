from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import json
from typing import Dict, List, Optional, Union

import numpy as np


@dataclass
class ComponentInfo:
    index: int
    name: str
    initialGasMoleFraction: Optional[float] = None
    isCarrierGas: bool = False


def load_simulation_metadata(simulation_json: Union[str, Path]) -> Dict[str, object]:
    path = Path(simulation_json)
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    components: List[ComponentInfo] = []
    for i, comp in enumerate(data.get("Components", [])):
        components.append(
            ComponentInfo(
                index=i,
                name=comp.get("Name", f"component {i}"),
                initialGasMoleFraction=comp.get("GasPhaseMolFraction"),
                isCarrierGas=bool(comp.get("CarrierGas", False)),
            )
        )

    return {
        "simulationType": data.get("SimulationType"),
        "displayName": data.get("DisplayName", path.stem),
        "temperature": data.get("Temperature", next((phase["Temperature"] for phase in data.get("SwingAdsorptionPhases", [])[:1] if "Temperature" in phase), 433.0)),
        "inletPressure": data.get("InletPressure", next((phase["InletPressure"] for phase in data.get("SwingAdsorptionPhases", [])[:1] if "InletPressure" in phase), 0.0)),
        "fluidPhase": data.get("FluidPhase", "Gas"),
        "pressureStart": data.get("PressureStart"),
        "pressureEnd": data.get("PressureEnd"),
        "numberOfPressurePoints": data.get("NumberOfPressurePoints"),
        "pressureScale": str(data.get("PressureScale", "log")).lower(),
        "components": components,
    }


def read_blocks(path: Union[str, Path]) -> List[np.ndarray]:
    """
    Read whitespace-separated numeric data split into blocks by blank lines.
    Lines starting with '#' are ignored.
    Malformed rows and rows with a non-dominant column count are skipped, which
    keeps readers usable while a simulation is still appending partial rows.
    Returns list of arrays, each (nrows, ncols).
    """
    path = Path(path)
    blocks: List[np.ndarray] = []
    rows: List[List[float]] = []

    def flush_block():
        nonlocal rows
        if rows:
            column_counts = Counter(len(row) for row in rows)
            ncols = max(column_counts, key=lambda count: (column_counts[count], count))
            rectangular_rows = [row for row in rows if len(row) == ncols]
            if rectangular_rows:
                blocks.append(np.asarray(rectangular_rows, dtype=float))
        rows = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                flush_block()
                continue
            try:
                row = [float(value) for value in s.split()]
                rows.append(row)
            except ValueError:
                pass

    flush_block()

    return blocks
