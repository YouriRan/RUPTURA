"""Convert a legacy RUPTURA text input file to the JSON input format.

The legacy text format is deprecated.  This converter covers every setting
accepted by the original text InputReader.  ``PulseBreakthrough`` and
``PulseTime`` have no equivalent in the JSON reader; conversion fails if they
are present unless ``--allow-lossy`` is supplied.

Examples:
    python convert_legacy_input_to_json.py simulation.input
    python convert_legacy_input_to_json.py simulation.input simulation.json
    python convert_legacy_input_to_json.py --allow-lossy simulation.input
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import warnings
from pathlib import Path
from typing import Any, NoReturn


class LegacyInputDeprecationWarning(DeprecationWarning):
    """Warning emitted when the deprecated legacy input format is converted."""


# A CLI user should see this warning; Python normally hides DeprecationWarning.
warnings.simplefilter("default", LegacyInputDeprecationWarning)


ISOTHERM_PARAMETER_COUNTS = {
    "langmuir": ("Langmuir", 2),
    "anti-langmuir": ("Anti-Langmuir", 2),
    "bet": ("BET", 3),
    "henry": ("Henry", 1),
    "freundlich": ("Freundlich", 2),
    "sips": ("Sips", 3),
    "langmuir-freundlich": ("Langmuir-Freundlich", 3),
    "redlich-peterson": ("Redlich-Peterson", 3),
    "toth": ("Toth", 3),
    "unilan": ("Unilan", 3),
    "o'brian&myers": ("OBrien&Myers", 3),
    "o'brien&myers": ("OBrien&Myers", 3),
    "obrien&myers": ("OBrien&Myers", 3),
    "quadratic": ("Quadratic", 3),
    "temkin": ("Temkin", 3),
    "bingel&walton": ("Bingel&Walton", 3),
}

STRING_KEYS = {
    "simulationtype": "SimulationType",
    "mixturepredictionmethod": "MixturePredictionMethod",
    "displayname": "DisplayName",
    # Used by some older example files before DisplayName became canonical.
    "columnname": "DisplayName",
}

FLOAT_KEYS = {
    "temperature": "Temperature",
    "columnvoidfraction": "ColumnVoidFraction",
    "particledensity": "ParticleDensity",
    "pressurestart": "PressureStart",
    "pressureend": "PressureEnd",
    "pressuregradient": "PressureGradient",
    "columnentrancevelocity": "ColumnEntranceVelocity",
    "timestep": "TimeStep",
    "columnlength": "ColumnLength",
}

INTEGER_KEYS = {
    "numberofpressurepoints": "NumberOfPressurePoints",
    "printevery": "PrintEvery",
    "writeevery": "WriteEvery",
    "numberofgridpoints": "NumberOfGridPoints",
    "columnpressure": "ColumnPressure",
    "columnloading": "ColumnLoading",
    "columnerror": "ColumnError",
}

COMPONENT_STRING_KEYS = {"filename": "FileName"}

COMPONENT_FLOAT_KEYS = {
    "gasphasemolfraction": "GasPhaseMolFraction",
    "masstransfercoefficient": "MassTransferCoefficient",
    "axialdispersioncoefficient": "AxialDispersionCoefficient",
}

UNSUPPORTED_KEYS = {"pulsebreakthrough", "pulsetime"}


class ConversionError(ValueError):
    """A legacy input cannot be converted safely."""


def fail(line_number: int, message: str) -> NoReturn:
    raise ConversionError(f"line {line_number}: {message}")


def first_argument(arguments: str, keyword: str, line_number: int) -> str:
    tokens = arguments.split()
    if not tokens:
        fail(line_number, f"{keyword} requires a value")
    return tokens[0]


def parse_float(arguments: str, keyword: str, line_number: int) -> float:
    token = first_argument(arguments, keyword, line_number)
    try:
        value = float(token)
    except ValueError:
        fail(line_number, f"{keyword} requires a number, got {token!r}")
    if not math.isfinite(value):
        fail(line_number, f"{keyword} requires a finite number")
    return value


def parse_integer(arguments: str, keyword: str, line_number: int) -> int:
    token = first_argument(arguments, keyword, line_number)
    if not re.fullmatch(r"[+]?[0-9]+", token):
        fail(line_number, f"{keyword} requires a non-negative integer, got {token!r}")
    return int(token)


def parse_boolean(arguments: str, keyword: str, line_number: int) -> bool:
    token = first_argument(arguments, keyword, line_number).casefold()
    if token in {"true", "yes"}:
        return True
    if token in {"false", "no"}:
        return False
    fail(line_number, f"{keyword} requires true/false or yes/no, got {token!r}")


def parse_float_list(arguments: str, keyword: str, line_number: int) -> list[float]:
    values = []
    for token in arguments.split():
        try:
            value = float(token)
        except ValueError:
            break
        if not math.isfinite(value):
            fail(line_number, f"{keyword} parameters must be finite numbers")
        values.append(value)
    if not values:
        fail(line_number, f"{keyword} requires numeric parameters")
    return values


def canonical_choice(
    arguments: str,
    keyword: str,
    line_number: int,
    choices: dict[str, str],
) -> str:
    token = first_argument(arguments, keyword, line_number)
    value = choices.get(token.casefold())
    if value is None:
        expected = ", ".join(sorted(set(choices.values())))
        fail(line_number, f"invalid {keyword} {token!r}; expected one of: {expected}")
    return value


def current_component(
    components: list[dict[str, Any]], keyword: str, line_number: int
) -> dict[str, Any]:
    if not components:
        fail(line_number, f"{keyword} appears before the first Component")
    return components[-1]


def parse_component(arguments: str, components: list[dict[str, Any]], line_number: int) -> None:
    tokens = arguments.split()
    if len(tokens) < 2:
        fail(line_number, "Component requires an index and name")
    try:
        component_id = int(tokens[0])
    except ValueError:
        fail(line_number, f"invalid component index {tokens[0]!r}")
    # The text reader assigned IDs by appearance and did not use the stated
    # index.  It also had historical examples using ``Component 0 Helium``.
    if tokens[1].casefold() == "moleculename":
        if len(tokens) < 3:
            fail(line_number, "MoleculeName must be followed by a name")
        component_name = tokens[2]
    else:
        component_name = tokens[1]
    components.append({"Name": component_name, "PhysisorptionSites": []})


def convert_text(text: str, *, allow_lossy: bool = False) -> tuple[dict[str, Any], list[str]]:
    """Convert legacy input text and return ``(JSON object, warnings)``."""
    output: dict[str, Any] = {}
    components: list[dict[str, Any]] = []
    declared_site_counts: list[tuple[int, int, int]] = []
    lossy_messages: list[str] = []

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        # Both comment styles are accepted by legacy examples.  The legacy
        # reader consumes only the first scalar token, so removing an inline
        # comment preserves its behavior.
        line = re.split(r"//|#", raw_line, maxsplit=1)[0].strip()
        if not line:
            continue

        parts = line.split(maxsplit=1)
        keyword = parts[0]
        arguments = parts[1].strip() if len(parts) == 2 else ""
        folded = keyword.casefold()

        if folded in UNSUPPORTED_KEYS:
            message = f"line {line_number}: {keyword} has no JSON-reader equivalent"
            if not allow_lossy:
                raise ConversionError(message + "; use --allow-lossy to omit it")
            lossy_messages.append(message + " and was omitted")
            continue

        if folded == "component":
            parse_component(arguments, components, line_number)
            continue

        if folded in COMPONENT_STRING_KEYS:
            component = current_component(components, keyword, line_number)
            component[COMPONENT_STRING_KEYS[folded]] = first_argument(arguments, keyword, line_number)
            continue

        if folded == "carriergas":
            component = current_component(components, keyword, line_number)
            component["CarrierGas"] = parse_boolean(arguments, keyword, line_number)
            continue

        if folded in COMPONENT_FLOAT_KEYS:
            component = current_component(components, keyword, line_number)
            component[COMPONENT_FLOAT_KEYS[folded]] = parse_float(arguments, keyword, line_number)
            continue

        if folded == "numberofisothermsites":
            current_component(components, keyword, line_number)
            declared_site_counts.append(
                (len(components) - 1, parse_integer(arguments, keyword, line_number), line_number)
            )
            continue

        if folded in ISOTHERM_PARAMETER_COUNTS:
            component = current_component(components, keyword, line_number)
            isotherm_type, parameter_count = ISOTHERM_PARAMETER_COUNTS[folded]
            values = parse_float_list(arguments, keyword, line_number)
            if len(values) < parameter_count:
                fail(
                    line_number,
                    f"{keyword} requires {parameter_count} parameters, got {len(values)}",
                )
            component["PhysisorptionSites"].append(
                {"Type": isotherm_type, "Parameters": values[:parameter_count]}
            )
            continue

        if folded in STRING_KEYS:
            if folded == "simulationtype":
                value = canonical_choice(
                    arguments,
                    keyword,
                    line_number,
                    {
                        "breakthrough": "Breakthrough",
                        "mixtureprediction": "MixturePrediction",
                        "fitting": "Fitting",
                        "test": "Test",
                    },
                )
            elif folded == "mixturepredictionmethod":
                value = canonical_choice(
                    arguments,
                    keyword,
                    line_number,
                    {"iast": "IAST", "siast": "SIAST", "ei": "EI", "sei": "SEI"},
                )
            elif folded in {"displayname", "columnname"}:
                value = first_argument(arguments, keyword, line_number)
            output[STRING_KEYS[folded]] = value
            continue

        if folded == "iastmethod":
            output["IASTMethod"] = canonical_choice(
                arguments,
                keyword,
                line_number,
                {
                    "fastias": "FastIAST",
                    "fastiast": "FastIAST",
                    "bisection": "NestedLoopBisection",
                    "nestedloopbisection": "NestedLoopBisection",
                },
            )
            continue

        if folded == "pressurescale":
            output["PressureScale"] = canonical_choice(
                arguments,
                keyword,
                line_number,
                {"log": "Log", "linear": "Linear", "normal": "Linear"},
            )
            continue

        if folded == "totalpressure":
            output["InletPressure"] = parse_float(arguments, keyword, line_number)
            # The legacy model specified inlet pressure and inlet velocity.
            output.setdefault("BoundaryCondition", "InletPressureInletVelocity")
            continue

        if folded == "numberoftimesteps":
            token = first_argument(arguments, keyword, line_number)
            if token.casefold() == "auto":
                output["NumberOfTimeSteps"] = "auto"
            else:
                output["NumberOfTimeSteps"] = parse_integer(arguments, keyword, line_number)
            continue

        if folded in FLOAT_KEYS:
            output[FLOAT_KEYS[folded]] = parse_float(arguments, keyword, line_number)
            continue

        if folded in INTEGER_KEYS:
            output[INTEGER_KEYS[folded]] = parse_integer(arguments, keyword, line_number)
            continue

        fail(line_number, f"unknown legacy keyword {keyword!r}")

    for component_id, expected, line_number in declared_site_counts:
        actual = len(components[component_id]["PhysisorptionSites"])
        if actual < expected:
            fail(
                line_number,
                f"Component {component_id} declares {expected} isotherm sites but contains {actual}",
            )
        if actual > expected:
            components[component_id]["PhysisorptionSites"] = components[component_id][
                "PhysisorptionSites"
            ][:expected]
            lossy_messages.append(
                f"line {line_number}: Component {component_id} declares {expected} "
                f"isotherm sites; ignored {actual - expected} surplus site(s), matching "
                "the legacy reader's effective behavior"
            )

    for component in components:
        if not component["PhysisorptionSites"]:
            del component["PhysisorptionSites"]

    if components:
        output["Components"] = components

    if output.get("SimulationType", "Breakthrough") == "Breakthrough":
        geometry: dict[str, Any] = {"Type": "PackedBed"}
        if "ColumnVoidFraction" in output:
            geometry["ColumnVoidFraction"] = output["ColumnVoidFraction"]
        output["Geometry"] = geometry

    return output, lossy_messages


def default_output_path(input_path: Path) -> Path:
    return input_path.with_suffix(".json")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a deprecated RUPTURA text input file to JSON.",
    )
    parser.add_argument("input", type=Path, help="legacy text input file")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="output JSON file (default: input name with .json suffix)",
    )
    parser.add_argument(
        "--allow-lossy",
        action="store_true",
        help="omit unsupported PulseBreakthrough/PulseTime settings instead of failing",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing output file",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level (default: 2)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_path = args.output or default_output_path(args.input)

    warnings.warn(
        "The legacy RUPTURA text input format is deprecated; use JSON input files.",
        LegacyInputDeprecationWarning,
        stacklevel=2,
    )

    try:
        if output_path.exists() and not args.force:
            raise ConversionError(
                f"output file {output_path} already exists; use --force to overwrite it"
            )
        text = args.input.read_text(encoding="utf-8")
        converted, lossy_messages = convert_text(text, allow_lossy=args.allow_lossy)
        output_path.write_text(
            json.dumps(converted, indent=args.indent, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    except (OSError, ConversionError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    for message in lossy_messages:
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    print(f"Converted {args.input} -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
