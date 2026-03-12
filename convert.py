#!/usr/bin/env python3
import argparse
import json
import re
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

COMMENT_RE = re.compile(r"//.*$")
WS_RE = re.compile(r"\s+")

KEY_MAP: Dict[str, str] = {
    "SimulationType": "SimulationType",
    "MixturePredictionMethod": "MixturePredictionMethod",
    "IASTMethod": "IASTMethod",
    "BreakthroughIntegrator": "BreakthroughIntegrator",
    "VelocityProfile": "VelocityProfile",
    "ReadColumnFile": "ReadColumnFile",
    "DisplayName": "DisplayName",
    "Temperature": "Temperature",
    "ColumnVoidFraction": "ColumnVoidFraction",
    "DynamicViscosity": "DynamicViscosity",
    "ParticleDiameter": "ParticleDiameter",
    "ParticleDensity": "ParticleDensity",
    "TotalPressure": "TotalPressure",
    "PressureStart": "PressureStart",
    "PressureEnd": "PressureEnd",
    "NumberOfPressurePoints": "NumberOfPressurePoints",
    "PressureScale": "PressureScale",
    "PressureGradient": "PressureGradient",
    "ColumnEntranceVelocity": "ColumnEntranceVelocity",
    "NumberOfTimeSteps": "NumberOfTimeSteps",
    "NumberOfInitTimeSteps": "NumberOfInitTimeSteps",
    "PulseBreakthrough": "PulseBreakthrough",
    "PulseTime": "PulseTime",
    "TimeStep": "TimeStep",
    "PrintEvery": "PrintEvery",
    "WriteEvery": "WriteEvery",
    "ColumnLength": "ColumnLength",
    "NumberOfGridPoints": "NumberOfGridPoints",
    "ColumnPressure": "ColumnPressure",
    "ColumnLoading": "ColumnLoading",
    "ColumnError": "ColumnError",
    "InfluxTemperature": "InfluxTemperature",
    "internalDiameter": "internalDiameter",
    "outerDiameter": "outerDiameter",
    "wallDensity": "wallDensity",
    "gasThermalConductivity": "gasThermalConductivity",
    "wallThermalConductivity": "wallThermalConductivity",
    "heatTransferGasSolid": "heatTransferGasSolid",
    "heatTransferGasWall": "heatTransferGasWall",
    "heatTransferWallExternal": "heatTransferWallExternal",
    "heatCapacityGas": "heatCapacityGas",
    "heatCapacitySolid": "heatCapacitySolid",
    "heatCapacityWall": "heatCapacityWall",
    "energyBalance": "energyBalance",
}

COMP_KEY_MAP: Dict[str, str] = {
    "FileName": "FileName",
    "CarrierGas": "CarrierGas",
    "GasPhaseMolFraction": "GasPhaseMolFraction",
    "MassTransferCoefficient": "MassTransferCoefficient",
    "AxialDispersionCoefficient": "AxialDispersionCoefficient",
    "MolecularWeight": "MolecularWeight",
    "HeatOfAdsorption": "HeatOfAdsorption",
    "NumberOfIsothermSites": "NumberOfIsothermSites",
    # allow Name to appear inside component blocks
    "MoleculeName": "Name",
    "Name": "Name",
}

ISOTHERM_TYPES: Dict[str, int] = {
    "Langmuir": 2,
    "Anti-Langmuir": 2,
    "BET": 3,
    "Henry": 1,
    "Freundlich": 2,
    "Sips": 3,
    "Langmuir-Freundlich": 3,
    "Redlich-Peterson": 3,
    "Toth": 3,
    "Unilan": 3,
    "O'Brian&Myers": 3,
    "Quadratic": 3,
    "Temkin": 3,
    "Bingel&Walton": 3,
}

KEY_MAP_L = {k.lower(): k for k in KEY_MAP}
COMP_KEY_MAP_L = {k.lower(): k for k in COMP_KEY_MAP}
ISOTHERM_TYPES_L = {k.lower(): k for k in ISOTHERM_TYPES}


def strip_comment(line: str) -> str:
    return COMMENT_RE.sub("", line).rstrip("\n")


def trim(s: str) -> str:
    return s.strip()


def split_tokens(s: str) -> List[str]:
    return [tok for tok in WS_RE.split(trim(s)) if tok]


def parse_bool(token: str) -> bool:
    t = token.strip().lower()
    if t in ("true", "yes", "1", "on"):
        return True
    if t in ("false", "no", "0", "off"):
        return False
    raise ValueError(f"Cannot parse boolean from '{token}'")


def parse_number(token: str) -> Any:
    t = token.strip()
    if t.lower() == "auto":
        return "auto"
    if re.fullmatch(r"[+-]?\d+", t):
        return int(t)
    if re.fullmatch(r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?", t):
        return float(t)
    return t


def parse_component_header(tokens: List[str], line_number: int) -> Tuple[int, Optional[str]]:
    # Accept:
    # 1) Component <id> MoleculeName <name>
    # 2) Component <id> <name>
    # 3) Component <id>            (name provided later via MoleculeName/Name)
    if len(tokens) < 2:
        raise RuntimeError(f"Malformed Component header at line {line_number}")

    try:
        cid = int(tokens[1])
    except ValueError:
        raise RuntimeError(f"Component id not an integer at line {line_number}: '{tokens[1]}'")

    name: Optional[str] = None
    # Look for MoleculeName token anywhere
    for j in range(2, len(tokens) - 1):
        if tokens[j].lower() == "moleculename":
            name = tokens[j + 1]
            break

    # Otherwise, if there is a third token, treat it as name
    if name is None and len(tokens) >= 3:
        name = tokens[2]

    # Otherwise leave as None (set later)
    return cid, name


def parse_isotherm(keyword_canonical: str, value_tokens: List[str], line_number: int) -> Dict[str, Any]:
    if not value_tokens:
        raise RuntimeError(f"Missing isotherm parameters for '{keyword_canonical}' at line {line_number}")
    params = [parse_number(t) for t in value_tokens]
    for p in params:
        if isinstance(p, str):
            raise RuntimeError(f"Non-numeric isotherm parameter '{p}' for '{keyword_canonical}' at line {line_number}")
    return {"Type": keyword_canonical, "Parameters": params}


def ruptura_to_json(text: str) -> Dict[str, Any]:
    # Keep top-level and components separate so we can append Components last.
    top: Dict[str, Any] = {}
    component_by_id: Dict[int, Dict[str, Any]] = {}
    current_component: Optional[Dict[str, Any]] = None

    for line_number, raw in enumerate(text.splitlines(True), start=1):
        raw_no_comment = strip_comment(raw)
        if trim(raw_no_comment) == "" or trim(raw_no_comment).startswith("#"):
            continue

        tokens = split_tokens(raw_no_comment)
        if not tokens:
            continue

        indented = raw_no_comment[:1].isspace()
        first = tokens[0]
        first_l = first.lower()

        # Component header
        if first_l == "component":
            cid, name = parse_component_header(tokens, line_number)
            comp: Dict[str, Any] = {"_ComponentId": cid, "IsothermSites": []}
            if name is not None:
                comp["Name"] = name
            component_by_id[cid] = comp
            current_component = comp
            continue

        # Decide if this line should be parsed as component-level:
        # - indented lines under a current component, OR
        # - non-indented lines that match known component keys/isotherms while a component is active.
        is_component_key = first_l in COMP_KEY_MAP_L or first_l in ISOTHERM_TYPES_L
        treat_as_component_line = current_component is not None and (indented or is_component_key)

        if treat_as_component_line:
            if first_l in COMP_KEY_MAP_L:
                canonical = COMP_KEY_MAP_L[first_l]
                json_key = COMP_KEY_MAP[canonical]
                if len(tokens) < 2:
                    raise RuntimeError(f"Missing value for '{canonical}' at line {line_number}")

                if canonical == "CarrierGas":
                    current_component[json_key] = parse_bool(tokens[1])
                elif canonical == "NumberOfIsothermSites":
                    current_component[json_key] = int(parse_number(tokens[1]))
                elif canonical in ("FileName", "MoleculeName", "Name"):
                    current_component[json_key] = tokens[1]
                else:
                    current_component[json_key] = parse_number(tokens[1])
                continue

            if first_l in ISOTHERM_TYPES_L:
                iso_canonical = ISOTHERM_TYPES_L[first_l]
                current_component.setdefault("IsothermSites", []).append(
                    parse_isotherm(iso_canonical, tokens[1:], line_number)
                )
                continue

            raise RuntimeError(f"Unknown component keyword '{first}' at line {line_number}")

        # Top-level key
        if first_l in KEY_MAP_L:
            canonical = KEY_MAP_L[first_l]
            json_key = KEY_MAP[canonical]
            if len(tokens) < 2:
                raise RuntimeError(f"Missing value for '{canonical}' at line {line_number}")

            if canonical in (
                "SimulationType",
                "MixturePredictionMethod",
                "IASTMethod",
                "BreakthroughIntegrator",
                "VelocityProfile",
                "ReadColumnFile",
                "DisplayName",
            ):
                top[json_key] = tokens[1]
            elif canonical in ("NumberOfTimeSteps",):
                top[json_key] = parse_number(tokens[1])  # allow "auto"
            elif canonical in ("PulseBreakthrough", "energyBalance"):
                top[json_key] = parse_bool(tokens[1])
            else:
                top[json_key] = parse_number(tokens[1])
            continue

        if first.startswith("//") or first.startswith("#"):
            continue

        raise RuntimeError(f"Unknown keyword '{first}' at line {line_number}")

    # Components last
    comps = sorted(component_by_id.values(), key=lambda c: c["_ComponentId"])
    for c in comps:
        c.pop("_ComponentId", None)
        if "IsothermSites" in c and not c["IsothermSites"]:
            c.pop("IsothermSites", None)
        if "Name" not in c:
            raise RuntimeError("Component without Name/MoleculeName encountered")

    top["Components"] = comps
    return top


def convert_one(input_path: Path, output_path: Path, indent: int) -> None:
    text = input_path.read_text(encoding="utf-8", errors="replace")
    data = ruptura_to_json(text)
    output_path.write_text(json.dumps(data, indent=indent, sort_keys=False) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert Ruptura keyword input files to JSON.")
    ap.add_argument("input", nargs="?", type=Path, help="Input Ruptura file (ignored if --batch)")
    ap.add_argument("-o", "--output", type=Path, default=None, help="Output JSON file (single-file mode only)")
    ap.add_argument("--indent", type=int, default=2, help="JSON indentation (default: 2)")
    ap.add_argument("--batch", action="store_true", help="Convert all ./*/*/simulation.input to simulation.json")
    ap.add_argument("--traceback", action="store_true", help="Print full traceback on failures")
    args = ap.parse_args()

    if args.batch:
        ok = 0
        fail = 0
        for input_path in Path(".").glob("*/*/simulation.input"):
            out_path = input_path.parent / "simulation.json"
            try:
                convert_one(input_path, out_path, args.indent)
                print(f"Wrote {out_path}")
                ok += 1
            except Exception as ex:
                fail += 1
                print(f"[FAILED] {input_path.parent}: {ex}")
                if args.traceback:
                    traceback.print_exc()
        print(f"Done. Converted: {ok}, Failed: {fail}")
        return

    if args.input is None:
        raise SystemExit("Error: provide an input file, or use --batch")

    input_path = args.input
    out_path = args.output if args.output is not None else (input_path.parent / "simulation.json")
    try:
        convert_one(input_path, out_path, args.indent)
        print(f"Wrote {out_path}")
    except Exception as ex:
        print(f"[FAILED] {input_path.parent}: {ex}")
        if args.traceback:
            traceback.print_exc()
        # exit non-zero so callers can detect failure if they want
        raise SystemExit(1)


if __name__ == "__main__":
    main()