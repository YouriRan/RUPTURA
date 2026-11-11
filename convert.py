"""Convert a legacy RUPTURA input file to the current JSON input format.

Two kinds of input are converted:

* The deprecated text format.  This converter covers every setting accepted by
  the original text InputReader.  ``PulseBreakthrough`` and ``PulseTime`` have
  no equivalent in the JSON reader; conversion fails if they are present unless
  ``--allow-lossy`` is supplied.
* JSON inputs with the step-based time keys of RUPTURA 2 (``NumberOfTimeSteps``,
  ``WriteEvery``, ``PrintEvery``, and swing-phase ``NumberOfTimeSteps``).  They
  become times in s (``SimulationTime``, ``WriteInterval``, phase ``Duration``);
  CVODE and ROS34PW2 inputs also lose ``TimeStep`` and ``NumberOfInitTimeSteps``.
* JSON inputs with the temperature keys that temperature laws replaced
  (``nonIsothermal``, the five-parameter ``GAB``, the General chemisorption
  ``AdsorptionActivationEnergy``/``DesorptionActivationEnergy``, and the reaction
  ``forwardActivationEnergy``/``gibbsFreeEnergy``).  The parameters become
  temperature laws, ``{"Value": ..., "Law": ...}``.

Without an output path a JSON file is converted in place and the original is kept
next to it as ``_legacy_<name>``, or as ``_legacy_2_<name>``, ... if that exists.

The rules (see ``convert_time_keys`` and ``convert_temperature_keys``) are the ones
the C++ InputReader prints when it rejects such an input.

Examples:
    python convert.py simulation.input
    python convert.py simulation.input simulation.json
    python convert.py --allow-lossy simulation.input
    python convert.py simulation.json
"""

from __future__ import annotations

import argparse
import copy
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


# The step-based time keys of RUPTURA 2 counted steps of TimeStep; these were its defaults.
LEGACY_TIME_STEP = 5.0e-4  # s
LEGACY_WRITE_EVERY = 10000
STEP_KEYS = ("numberoftimesteps", "writeevery", "printevery")


def find_key(settings: dict[str, Any], key: str) -> str | None:
    """Return the spelling of ``key`` in ``settings`` that the C++ reader uses, or None.

    Input keys are case-insensitive.  Like the reader, an exact match wins; otherwise the first match in sorted
    key order (the order of nlohmann::json objects) does, so a duplicate in other capitalization is ignored.
    """
    if key in settings:
        return key
    folded = key.casefold()
    return next((name for name in sorted(settings) if name.casefold() == folded), None)


def seconds(count: Any, time_step: float) -> Any:
    """Return ``count`` steps of ``time_step`` in s, rounded to 12 significant digits as the C++ reader prints it."""
    if isinstance(count, bool) or not isinstance(count, (int, float)):
        return count
    return float(f"{count * time_step:.12g}")


def convert_time_keys(settings: dict[str, Any]) -> dict[str, Any] | None:
    """Return ``settings`` with the step-based time keys replaced by times in s, or None if it has none.

    SimulationTime = NumberOfTimeSteps * TimeStep ("auto" stays "auto"), WriteInterval = WriteEvery * TimeStep
    (WriteEvery defaults to 10000), and each swing phase's Duration = NumberOfTimeSteps * TimeStep, with TimeStep
    defaulting to 5e-4 s.  PrintEvery is dropped (progress is printed at every output), and so is the root
    NumberOfTimeSteps of a SwingAdsorption input (the phase durations set its length).  CVODE and ROS34PW2 inputs
    also lose TimeStep and NumberOfInitTimeSteps, which only the fixed-step RungeKutta3 uses.
    """
    integrator_key = find_key(settings, "BreakthroughIntegrator")
    adaptive = integrator_key is not None and str(settings[integrator_key]).casefold() in ("cvode", "ros34pw2")
    time_step_key = find_key(settings, "TimeStep")
    time_step_value = settings.get(time_step_key) if time_step_key is not None else None
    time_step = (
        float(time_step_value)
        if isinstance(time_step_value, (int, float)) and not isinstance(time_step_value, bool)
        else LEGACY_TIME_STEP
    )
    phases_key = find_key(settings, "SwingAdsorptionPhases")
    phases = settings[phases_key] if phases_key is not None and isinstance(settings[phases_key], list) else None

    has_step_keys = any(find_key(settings, key) is not None for key in STEP_KEYS) or any(
        isinstance(phase, dict) and find_key(phase, "NumberOfTimeSteps") is not None for phase in phases or []
    )
    dropped_for_adaptive = {"timestep", "numberofinittimesteps"} if adaptive else set()
    if not has_step_keys and not any(key.casefold() in dropped_for_adaptive for key in settings):
        return None

    steps_key = find_key(settings, "NumberOfTimeSteps")
    write_every_key = find_key(settings, "WriteEvery")
    has_write_interval = write_every_key is not None or find_key(settings, "WriteInterval") is not None
    converted: dict[str, Any] = {}
    for key, value in settings.items():
        folded = key.casefold()
        if folded in dropped_for_adaptive or folded == "printevery":
            continue
        if folded == "numberoftimesteps":
            # Only the spelling the reader used counts; duplicates in other capitalization were ignored.
            if phases is None and key == steps_key:
                is_auto = isinstance(value, str) and value.casefold() == "auto"
                converted["SimulationTime"] = "auto" if is_auto else seconds(value, time_step)
                if not has_write_interval:
                    converted["WriteInterval"] = seconds(LEGACY_WRITE_EVERY, time_step)
                    has_write_interval = True
            continue
        if folded == "writeevery":
            if key == write_every_key:
                converted["WriteInterval"] = seconds(value, time_step)
            continue
        if folded == "swingadsorptionphases" and phases is not None:
            # Each phase keeps its key order, with Duration in place of NumberOfTimeSteps.
            converted[key] = [
                {
                    ("Duration" if name.casefold() == "numberoftimesteps" else name): (
                        seconds(item, time_step) if name.casefold() == "numberoftimesteps" else item
                    )
                    for name, item in phase.items()
                }
                if isinstance(phase, dict)
                else phase
                for phase in value
            ]
            continue
        converted[key] = value
    if has_step_keys and not has_write_interval:
        converted["WriteInterval"] = seconds(LEGACY_WRITE_EVERY, time_step)
    return converted


# Gas constant in J/(mol K), as in the C++ code.
GAS_CONSTANT = 8.31446261815324

# The isotherm types that "nonIsothermal" scaled, and the position of their affinity in Parameters.
AFFINITY_INDEX = {"langmuir": 1, "ph-langmuir": 1, "sips": 1, "langmuir-freundlich": 1, "toth": 1}


def is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def renamed(mapping: dict[str, Any], replacements: dict[str, tuple[str, Any] | None]) -> dict[str, Any]:
    """Return ``mapping`` in its key order, with each key in ``replacements`` (case-insensitive) renamed to a new
    (key, value) pair, or dropped when its replacement is None."""
    folded = {key.casefold(): replacement for key, replacement in replacements.items()}
    result: dict[str, Any] = {}
    for key, value in mapping.items():
        if key.casefold() not in folded:
            result[key] = value
        elif folded[key.casefold()] is not None:
            new_key, new_value = folded[key.casefold()]
            result[new_key] = new_value
    return result


def exponential_law(value: Any, law: str, energy: Any, reference_temperature: float | None) -> Any:
    """Return the law object of ``value * exp(+-energy/(R T))`` (VantHoff: +, Arrhenius: -), or ``value`` itself for
    a zero energy.

    These laws were uncentred. A law inherits the ReferenceTemperature of its component, so when the component has
    one, the value is written at that temperature, V = value exp(+-energy/(R T_ref)), which is the same function.
    """
    if not is_number(value) or not is_number(energy) or energy == 0:
        return value
    coefficient_key = "Energy" if law == "VantHoff" else "ActivationEnergy"
    if reference_temperature is not None:
        sign = 1.0 if law == "VantHoff" else -1.0
        value = value * math.exp(sign * energy / (GAS_CONSTANT * reference_temperature))
    return {"Value": value, "Law": law, coefficient_key: energy}


def convert_component_temperature_keys(component: dict[str, Any], where: str, global_component: dict[str, Any] | None,
                                       notes: list[str]) -> dict[str, Any] | None:
    """Return a component object (or a ComponentParameters override) with temperature laws, or None if it needs none.

    ``global_component`` is the component an override applies to; an override without sites of its own may only
    repeat its nonIsothermal setting.
    """
    def temperature_settings(item: dict[str, Any] | None) -> tuple[bool | None, float | None]:
        if item is None:
            return None, None
        flag = find_key(item, "nonIsothermal")
        reference = find_key(item, "ReferenceTemperature")
        return (item[flag] is True if flag is not None else None,
                float(item[reference]) if reference is not None and is_number(item[reference]) else None)

    flag_key = find_key(component, "nonIsothermal")
    reference_key = find_key(component, "ReferenceTemperature")
    heat_key = find_key(component, "HeatOfAdsorption")
    own_flag, own_reference = temperature_settings(component)
    global_flag, global_reference = temperature_settings(global_component)
    # An override starts from the component, so it scales when either sets nonIsothermal, with either's T_ref.
    non_isothermal = own_flag if own_flag is not None else bool(global_flag)
    effective_reference = own_reference if own_reference is not None else global_reference
    site_kinds = [kind for kind in ("PhysisorptionSites", "ChemisorptionSites") if find_key(component, kind)]

    if global_component is not None and own_flag is not None and own_flag != bool(global_flag) and not site_kinds:
        raise ConversionError(f"{where}: nonIsothermal differs from the component's but the override has no sites "
                              "of its own; write the temperature laws on the sites")

    # The reference temperature only centred the nonIsothermal scaling. It stays (as ReferenceTemperature) where it
    # did, and the laws written here inherit it; an override inherits the component's when it gives none.
    write_reference = non_isothermal and effective_reference is not None
    global_kept_reference = global_reference if global_flag else None
    reference_temperature = effective_reference if write_reference else global_kept_reference

    converted = copy.deepcopy(component)
    changed = flag_key is not None
    component_heat = component[heat_key] if heat_key is not None else None
    for kind in site_kinds:
        sites = converted[find_key(converted, kind)]
        if not isinstance(sites, list):
            continue
        affinity_heats = []
        for index, site in enumerate(sites):
            if not isinstance(site, dict):
                continue
            type_key = find_key(site, "Type")
            parameters_key = find_key(site, "Parameters")
            isotherm_type = str(site[type_key]).casefold() if type_key is not None else ""
            parameters = site[parameters_key] if parameters_key is not None else None
            if isinstance(parameters, list) and isotherm_type == "gab" and len(parameters) == 5:
                # C = C_0 exp(H_C/(R T)) and K = K_0 exp(H_K/(R T)) were scaled by exp(1/(R T)) to the power H_C, H_K.
                q, c0, k0, hc, hk = parameters
                site[parameters_key] = [q, exponential_law(c0, "VantHoff", hc, reference_temperature),
                                        exponential_law(k0, "VantHoff", hk, reference_temperature)]
                changed = True
            elif non_isothermal and isinstance(parameters, list) and isotherm_type in AFFINITY_INDEX:
                position = AFFINITY_INDEX[isotherm_type]
                enabled = parameters and is_number(parameters[0]) and parameters[0] != 0
                if enabled and len(parameters) > position and is_number(parameters[position]) and parameters[position]:
                    parameters[position] = {"Value": parameters[position], "Law": "VantHoff"}
                    site_heat_key = find_key(site, "HeatOfAdsorption")
                    affinity_heats.append(site[site_heat_key] if site_heat_key is not None else component_heat)
            rate_key = find_key(site, "RateEquationParameters")
            rate_parameters = site.get(rate_key) if rate_key is not None else None
            if kind == "ChemisorptionSites" and isinstance(rate_parameters, dict):
                replacements: dict[str, tuple[str, Any] | None] = {}
                for direction in ("Adsorption", "Desorption"):
                    energy_key = find_key(rate_parameters, f"{direction}ActivationEnergy")
                    if energy_key is None:
                        continue
                    replacements[energy_key] = None
                    coefficient_key = find_key(rate_parameters, f"{direction}RateCoefficient")
                    if coefficient_key is not None:
                        replacements[coefficient_key] = (
                            coefficient_key,
                            exponential_law(rate_parameters[coefficient_key], "Arrhenius",
                                            rate_parameters[energy_key], reference_temperature))
                if replacements:
                    site[rate_key] = renamed(rate_parameters, replacements)
                    changed = True
        if kind == "ChemisorptionSites" and len(set(map(repr, affinity_heats))) > 1:
            # The chemisorption mixture scaled every site with the heat of the first one.
            raise ConversionError(f"{where}: the chemisorption sites have different heats ({affinity_heats}), but "
                                  "nonIsothermal scaled all of them with the first; write their VantHoff laws by hand")

    if not changed:
        return None
    replacements = {}
    if flag_key is not None:
        replacements[flag_key] = None
    if reference_key is not None:
        replacements[reference_key] = ("ReferenceTemperature", component[reference_key]) if write_reference else None
    if heat_key is not None and component_heat == 1.0 and is_gab_heat_placeholder(converted):
        replacements[heat_key] = None
        notes.append(f"{where}: dropped HeatOfAdsorption 1.0, which only made GAB see exp(1/(RT)) as its scale")
    result = renamed(converted, replacements)
    if write_reference and reference_key is None and global_kept_reference is None:
        # The component's T_ref centred this override's scaling, but the component does not keep it.
        result["ReferenceTemperature"] = effective_reference
    return result


def is_gab_heat_placeholder(component: dict[str, Any]) -> bool:
    """Return whether a HeatOfAdsorption of 1.0 served only GAB: the component has a GAB site, every enabled
    physisorption site is GAB, and every chemisorption site gives its own HeatOfAdsorption."""
    sites = {kind: component.get(find_key(component, kind) or "", []) for kind in ("PhysisorptionSites",
                                                                                   "ChemisorptionSites")}
    all_sites = [site for kind_sites in sites.values() if isinstance(kind_sites, list) for site in kind_sites
                 if isinstance(site, dict)]

    def site_type(site: dict[str, Any]) -> str:
        type_key = find_key(site, "Type")
        return str(site[type_key]).casefold() if type_key is not None else ""

    def enabled(site: dict[str, Any]) -> bool:
        parameters = site.get(find_key(site, "Parameters") or "", [])
        first = parameters[0] if isinstance(parameters, list) and parameters else 0
        return not is_number(first) or first != 0

    return (any(site_type(site) == "gab" for site in all_sites)
            and all(site_type(site) == "gab" or not enabled(site) for site in sites["PhysisorptionSites"] or [])
            and all(find_key(site, "HeatOfAdsorption") is not None for site in sites["ChemisorptionSites"] or []))


def convert_reaction_temperature_keys(kinetics: dict[str, Any]) -> dict[str, Any] | None:
    """Return reaction Kinetics with temperature laws, or None if it has no forwardActivationEnergy or
    gibbsFreeEnergy.

    k+ = k0 exp(-Ea/(R T)) becomes an Arrhenius ForwardRateCoefficient. K = K0 exp(-dG/(R T)) becomes a VantHoff
    EquilibriumConstant whose Energy is the HeatOfReaction -dG, which was also the heat released.
    """
    energy_key = find_key(kinetics, "forwardActivationEnergy")
    gibbs_key = find_key(kinetics, "gibbsFreeEnergy")
    if energy_key is None and gibbs_key is None:
        return None
    energy = kinetics[energy_key] if energy_key is not None else 0
    gibbs = kinetics[gibbs_key] if gibbs_key is not None else 0
    replacements: dict[str, tuple[str, Any] | None] = {}
    if energy_key is not None:
        replacements[energy_key] = None
    if (rate_key := find_key(kinetics, "ForwardRateCoefficient")) is not None:
        replacements[rate_key] = ("ForwardRateCoefficient",
                                  exponential_law(kinetics[rate_key], "Arrhenius", energy, None))
    if (constant_key := find_key(kinetics, "EquilibriumConstant")) is not None:
        constant = kinetics[constant_key]
        replacements[constant_key] = ("EquilibriumConstant",
                                      {"Value": constant, "Law": "VantHoff"} if is_number(gibbs) and gibbs else constant)
    if gibbs_key is not None:
        replacements[gibbs_key] = ("HeatOfReaction", -gibbs) if is_number(gibbs) and gibbs else None
    return renamed(kinetics, replacements)


def convert_temperature_keys(settings: dict[str, Any], notes: list[str] | None = None) -> dict[str, Any] | None:
    """Return ``settings`` with the temperature keys replaced by temperature laws, or None if it has none.

    ``"nonIsothermal": true`` becomes a VantHoff law on the affinity of each enabled Langmuir, pH-Langmuir, Sips,
    Langmuir-Freundlich, and Toth site, which inherits the HeatOfAdsorption (and the ReferenceTemperature, now
    spelled so). ``"nonIsothermal": false`` and its referenceTemperature are dropped. The five-parameter GAB
    [q_m, C_0, K_0, H_C, H_K] becomes [q_m, VantHoff C_0 with Energy H_C, VantHoff K_0 with Energy H_K], and a
    HeatOfAdsorption of 1.0 that only served GAB is dropped. The General chemisorption activation energies become
    Arrhenius laws of their rate coefficients. Reactions are converted by convert_reaction_temperature_keys. These are
    the rules the C++ InputReader prints when it rejects such an input. Notes on dropped keys are appended to
    ``notes``.
    """
    notes = notes if notes is not None else []
    converted = copy.deepcopy(settings)
    changed = False

    components_key = find_key(converted, "Components")
    components = converted.get(components_key) if components_key is not None else None
    by_name: dict[str, dict[str, Any]] = {}
    if isinstance(components, (list, dict)):
        entries = list(components.items()) if isinstance(components, dict) else list(enumerate(components))
        for index, component in entries:
            if not isinstance(component, dict):
                continue
            new = convert_component_temperature_keys(component, f"Components[{index}]", None, notes)
            name_key = find_key(component, "Name")
            if name_key is not None:
                by_name[str(component[name_key])] = component
            if new is not None:
                components[index] = new
                changed = True

    adsorbents_key = find_key(converted, "Adsorbents")
    adsorbents = converted.get(adsorbents_key) if adsorbents_key is not None else None
    for adsorbent_index, adsorbent in enumerate(adsorbents if isinstance(adsorbents, list) else []):
        if not isinstance(adsorbent, dict):
            continue
        adsorbent_components = adsorbent.get(find_key(adsorbent, "Components") or "")
        for index, component in enumerate(adsorbent_components if isinstance(adsorbent_components, list) else []):
            if isinstance(component, dict):
                new = convert_component_temperature_keys(
                    component, f"Adsorbents[{adsorbent_index}].Components[{index}]", None, notes)
                if new is not None:
                    adsorbent_components[index] = new
                    changed = True
        overrides = adsorbent.get(find_key(adsorbent, "ComponentParameters") or "")
        for name, override in (overrides.items() if isinstance(overrides, dict) else []):
            if isinstance(override, dict):
                new = convert_component_temperature_keys(
                    override, f"Adsorbents[{adsorbent_index}].ComponentParameters.{name}", by_name.get(name), notes)
                if new is not None:
                    overrides[name] = new
                    changed = True

    reactions_key = find_key(converted, "Reactions")
    reactions = converted.get(reactions_key) if reactions_key is not None else None
    for reaction in reactions if isinstance(reactions, list) else []:
        kinetics_key = find_key(reaction, "Kinetics") if isinstance(reaction, dict) else None
        if kinetics_key is not None and isinstance(reaction[kinetics_key], dict):
            new = convert_reaction_temperature_keys(reaction[kinetics_key])
            if new is not None:
                reaction[kinetics_key] = new
                changed = True

    return converted if changed else None


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

    return convert_time_keys(output) or output, lossy_messages


def default_output_path(input_path: Path) -> Path:
    return input_path.with_suffix(".json")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a deprecated RUPTURA text input file, or a JSON input with step-based time keys or "
        "removed temperature keys, to the current JSON format.",
    )
    parser.add_argument("input", type=Path, help="legacy text input file, or JSON input file (*.json)")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="output JSON file (default: input name with .json suffix; a JSON input is converted in place and "
        "the original kept as _legacy_<name>, or _legacy_2_<name>, ... if that exists)",
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
        default=None,
        help="JSON indentation level (default: 2, or that of the JSON input)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.input.suffix.casefold() == ".json":
        output_path = args.output or args.input
        # The original is kept under the first free name, so that an earlier conversion's archive survives.
        archive_path = args.input.with_name("_legacy_" + args.input.name)
        archive_number = 2
        while archive_path.exists():
            archive_path = args.input.with_name(f"_legacy_{archive_number}_{args.input.name}")
            archive_number += 1
        notes: list[str] = []
        try:
            text = args.input.read_text(encoding="utf-8")
            settings = json.loads(text)
            if not isinstance(settings, dict):
                raise ConversionError(f"{args.input} is not a JSON object")
            time_converted = convert_time_keys(settings)
            converted = convert_temperature_keys(time_converted or settings, notes) or time_converted
            if converted is None:
                print(f"{args.input} has no step-based time keys or temperature keys; nothing to convert")
                return 0
            if output_path != args.input and output_path.exists() and not args.force:
                raise ConversionError(f"output file {output_path} already exists; use --force to overwrite it")
            if output_path == args.input:
                archive_path.write_text(text, encoding="utf-8")
            # Keep the indentation of the input file, so that a converted file differs only in the converted keys.
            indented = re.search(r"\n( +)\S", text)
            indent = args.indent if args.indent is not None else (len(indented.group(1)) if indented else 2)
            output_path.write_text(json.dumps(converted, indent=indent, ensure_ascii=False) + "\n", encoding="utf-8")
        except (OSError, ValueError) as error:
            print(f"error: {error}", file=sys.stderr)
            return 1
        kept = f" (original kept as {archive_path})" if output_path == args.input else ""
        print(f"Converted {args.input} -> {output_path}{kept}")
        for note in notes:
            print(f"  note: {note}")
        return 0

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
            json.dumps(converted, indent=2 if args.indent is None else args.indent, ensure_ascii=False) + "\n",
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
