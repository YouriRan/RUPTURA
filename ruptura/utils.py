from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import glob
import io
import json
import os
import re
from typing import Dict, List, Optional, Union

import numpy as np


@dataclass
class ComponentInfo:
    index: int
    name: str
    Yi0: Optional[float] = None
    isCarrierGas: bool = False


@dataclass
class HeaderComponentInfo:
    name: str
    props: Dict[str, int]  # prop -> 0-based column index


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
                Yi0=comp.get("GasPhaseMolFraction"),
            )
        )

    return {
        "simulationType": data.get("SimulationType"),
        "displayName": data.get("DisplayName", path.stem),
        "temperature": data.get("Temperature"),
        "pressureStart": data.get("PressureStart"),
        "pressureEnd": data.get("PressureEnd"),
        "numberOfPressurePoints": data.get("NumberOfPressurePoints"),
        "pressureScale": str(data.get("PressureScale", "log")).lower(),
        "components": components,
    }


def infer_components_from_files(
    data_dir: Union[str, Path] = ".",
) -> List[ComponentInfo]:
    pattern = str(Path(data_dir) / "component_*.data")
    files = sorted(glob.glob(pattern))

    components: List[ComponentInfo] = []
    regex = re.compile(r"component_(\d+)_(.+)\.data$")

    for fileName in files:
        base = os.path.basename(fileName)
        match = regex.match(base)
        if match is None:
            continue
        components.append(
            ComponentInfo(
                index=int(match.group(1)),
                name=match.group(2),
            )
        )

    components.sort(key=lambda x: x.index)
    return components


def _read_header_lines(path: Path) -> List[str]:
    header_lines: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.strip() == "" or line.lstrip().startswith("#"):
                header_lines.append(line)
            else:
                break
    return header_lines


def parse_header_metadata(path: Union[str, Path]) -> Dict:
    """
    Parse header comments to detect:
      - z_col: 0-based index of z column (default 0)
      - components: {comp_idx: HeaderComponentInfo(name, props)}
      - columns: {col_number_1based: description}

    Recognized patterns:
      # column 6: component 0 P (partial pressure)
      # component 0: Helium (y_i=0.93)   (optional)
    """
    path = Path(path)
    header_lines = _read_header_lines(path)

    columns: Dict[int, str] = {}
    components: Dict[int, HeaderComponentInfo] = {}

    col_pat = re.compile(r"^\s*#\s*column\s+(\d+)\s*:\s*(.+?)\s*$", re.IGNORECASE)
    comp_name_pat = re.compile(r"^\s*#\s*component\s+(\d+)\s*[:=]\s*(.+?)\s*$", re.IGNORECASE)

    comp_names: Dict[int, str] = {}
    for line in header_lines:
        m = comp_name_pat.match(line)
        if m:
            comp_names[int(m.group(1))] = m.group(2).strip()

    z_col_0 = 0
    for line in header_lines:
        m = col_pat.match(line)
        if not m:
            continue

        col_num = int(m.group(1))
        desc = m.group(2).strip()
        columns[col_num] = desc

        if re.match(r"^\s*z\b", desc, flags=re.IGNORECASE):
            z_col_0 = col_num - 1

        m2 = re.match(r"^\s*component\s+(\d+)\s+([A-Za-z0-9_]+)\b", desc, flags=re.IGNORECASE)
        if m2:
            ci = int(m2.group(1))
            prop = m2.group(2)
            if ci not in components:
                name = comp_names.get(ci, f"component {ci}")
                components[ci] = HeaderComponentInfo(name=name, props={})
            components[ci].props[prop] = col_num - 1

    for ci, nm in comp_names.items():
        components.setdefault(ci, HeaderComponentInfo(name=nm, props={}))

    return {"columns": columns, "z_col": z_col_0, "components": components}


def read_blocks(path: Union[str, Path]) -> List[np.ndarray]:
    """
    Read whitespace-separated numeric data split into blocks by blank lines.
    Lines starting with '#' are ignored.
    Returns list of arrays, each (nrows, ncols).
    """
    path = Path(path)
    blocks: List[np.ndarray] = []
    buf: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith("#"):
                if buf:
                    arr = np.loadtxt(io.StringIO("".join(buf)))
                    if arr.ndim == 1:
                        arr = arr.reshape(1, -1)
                    blocks.append(arr)
                    buf = []
                continue
            buf.append(line)

    if buf:
        arr = np.loadtxt(io.StringIO("".join(buf)))
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        blocks.append(arr)

    return blocks
