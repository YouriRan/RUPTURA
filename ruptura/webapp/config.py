"""Lossless simulation documents and migration of saved Qt Lab workspaces."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import uuid

from .field_help import FIELD_HELP

TYPES = ["Breakthrough", "SwingAdsorption", "MixturePrediction"]
METHODS = ["IAST", "SIAST", "EI", "SEI", "SCI", "SPI"]
PRESET_INFO = {
    "breakthrough": {"title": "Gas breakthrough", "description": "Gas feed through a packed adsorption column."},
    "liquid": {"title": "Liquid & pH", "description": "Liquid transport with pH-dependent adsorption."},
    "mixture": {"title": "Mixture prediction", "description": "Pure and mixture equilibrium isotherms."},
    "swing": {"title": "Swing adsorption", "description": "Ordered adsorption and regeneration phases."},
    "chemisorption": {"title": "Chemisorption", "description": "Physical and chemical uptake together."},
    "multibed": {"title": "Multibed column", "description": "Distinct adsorbents in a layered column."},
    "paper-t15": {"title": "T15 · Marx et al.", "description": "CO₂ / N₂ / H₂ 40/30/30 on activated carbon. Energy balance, 298 K, 10 bar, 4,000 s.", "group": "paper"},
    "paper-bea-alkanes": {"title": "BEA · C₇ alkanes", "description": "Seven C₇ isomers at 1% each in helium. Isothermal IAST, 552 K, 1.2 bar, 1,200 s.", "group": "paper"},
}
PRESET_INFO["paper-monolith-400ppm-2mmol"] = {
    "title": "Monolith · 400 ppm · 2 mmol/s",
    "description": "Stampi-Bombelli et al. · CO₂ in N₂, square channels, dual uptake, literature dispersion, 24,000 s.",
    "group": "paper",
}
PAPER_PRESETS = tuple(key for key, info in PRESET_INFO.items() if info.get("group") == "paper")
EXAMPLE_ORDER = ("breakthrough", *PAPER_PRESETS)
for _number, _key in enumerate(EXAMPLE_ORDER, 1):
    PRESET_INFO[_key]["number"] = _number
    if _key in PAPER_PRESETS:
        PRESET_INFO[_key]["title"] = f"Example {_number}: {PRESET_INFO[_key]['title']}"
ISOTHERMS = {
    "Langmuir": [("qmax", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity [1/Pa]", 1e-5)],
    "pH-Langmuir": [("qmax", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity [m³/mol]", 0.5), ("pKa", "pKa", 7.0)],
    "Anti-Langmuir": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity [1/Pa]", 1e-5)],
    "BET": [("monolayerLoading", "Monolayer capacity [mol/kg]", 1.0), ("cParameter", "C", 10.0), ("saturationPressureFactor", "Saturation pressure factor", 1e-5)],
    "Henry": [("henryConstant", "Henry constant [mol/kg/Pa]", 1e-5)],
    "Freundlich": [("freundlichConstant", "Freundlich constant", 1e-5), ("heterogeneityExponent", "Exponent", 1.0)],
    "Sips": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("heterogeneityExponent", "Exponent", 1.0)],
    "Langmuir-Freundlich": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("heterogeneityExponent", "Exponent", 1.0)],
    "Redlich-Peterson": [("linearConstant", "Linear constant", 1e-5), ("affinityConstant", "Affinity", 1e-5), ("exponent", "Exponent", 1.0)],
    "Toth": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("heterogeneityParameter", "Heterogeneity", 1.0)],
    "Unilan": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("energyDistributionWidth", "Energy width", 1.0)],
    "OBrien&Myers": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("variance", "Variance", 1.0)],
    "Quadratic": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("linearAffinity", "Linear affinity", 1e-5), ("quadraticAffinity", "Quadratic affinity", 1e-10)],
    "Temkin": [("saturationLoading", "Capacity [mol/kg]", 1.0), ("affinityConstant", "Affinity", 1e-5), ("interactionParameter", "Interaction", 0.1)],
    "Bingel&Walton": [("qmax", "Capacity [mol/kg]", 1.0), ("a", "a [1/Pa]", 1e-5), ("b", "b [1/Pa]", 1e-5)],
    "GAB": [("qmax", "Monolayer capacity [mol/kg]", 1.0), ("c0", "C₀", 1.0), ("k0", "K₀", 0.5), ("dHC", "ΔH_C [J/mol]", 0.0), ("dHK", "ΔH_K [J/mol]", 0.0)],
}
KINETICS = {
    "FirstOrder": [("MassTransferCoefficient", "Rate coefficient [1/s]", 0.04)],
    "PseudoNth": [("RateCoefficient", "Rate coefficient", 0.04), ("Order", "Order", 1.0)],
    "Avrami": [("RateCoefficient", "Rate coefficient", 0.04), ("Order", "Order", 1.0)],
    "General": [("AdsorptionRateCoefficient", "Adsorption rate", 2.5e-5), ("AdsorptionActivationEnergy", "Adsorption activation energy [J/mol]", 0.0),
                ("DesorptionRateCoefficient", "Desorption rate", 1e-6), ("DesorptionActivationEnergy", "Desorption activation energy [J/mol]", 0.0),
                ("PoreConcentrationOrder", "Pore concentration order", 1), ("CapacityOrder", "Capacity order", 1), ("DesorptionOrder", "Desorption order", 1),
                ("FilmMassTransferCoefficient", "Film transfer [m/s]", 0.0015), ("PoreDiffusivity", "Pore diffusivity [m²/s]", 8e-11), ("UsePoreSurfaceTransport", "Pore surface transport", True)],
    "Elovich": [("Alpha", "Alpha", 0.02), ("Beta", "Beta", 3.0), ("FilmMassTransferCoefficient", "Film transfer [m/s]", 0.001),
                ("PoreDiffusivity", "Pore diffusivity [m²/s]", 5e-11), ("UsePoreSurfaceTransport", "Pore surface transport", True)],
}


def field(key, label, default=None, options=None, kind="number"):
    return {"key": key, "label": label, "default": default, "type": "select" if options else kind, "options": options,
            "description": FIELD_HELP[key]}


SCHEMA = {
    "help": FIELD_HELP,
    "types": TYPES, "methods": METHODS, "isotherms": ISOTHERMS, "kinetics": KINETICS,
    "operating": [field("Temperature", "Temperature [K]", 433), field("FluidPhase", "Fluid phase", "Gas", ["Gas", "Liquid"]),
                  field("MixturePredictionMethod", "Mixture method", "IAST", METHODS), field("IASTMethod", "IAST solver", "FastIAST", ["FastIAST", "NestedLoopBisection"])],
    "integration": [field("BreakthroughIntegrator", "Integrator", "RungeKutta3", ["RungeKutta3", "CVODE", "SIRK3"]),
                    field("ConvectionScheme", "Convection scheme", "WENO3", ["WENO3", "Upwind"]),
                    field("SimulationTime", "Simulation time [s] (or auto)", "auto", kind="text"),
                    field("WriteInterval", "Output interval [s]", 5.0), field("NumberOfGridPoints", "Grid intervals", 100),
                    field("TimeStep", "Time step [s]", 0.0005), field("NumberOfInitTimeSteps", "Initialization steps", 0)],
    "cvode": [field("CVODERelativeTolerance", "Relative tolerance", 1e-4), field("CVODEAbsoluteToleranceConcentration", "Concentration tolerance [mol/m³]", 1e-5),
              field("CVODEAbsoluteToleranceLoading", "Loading tolerance [mol/kg]", 1e-6), field("CVODEAbsoluteToleranceTemperature", "Temperature tolerance [K]", 1e-3),
              field("CVODEMaximumTimeStep", "Maximum internal step [s] (−1 auto, 0 unlimited)", -1),
              field("CVODEMaximumOrder", "Maximum BDF order", 5),
              field("CVODELinearSolver", "Linear solver", "Dense", ["Dense", "SPGMR"]), field("CVODEKrylovDimension", "Krylov dimension", 30)],
    "pressure": [field("PressureStart", "Start pressure [Pa]", 1), field("PressureEnd", "End pressure [Pa]", 1e6),
                 field("NumberOfPressurePoints", "Pressure points", 100), field("PressureScale", "Pressure scale", "Log", ["Log", "Linear"])],
    "column": [field("ColumnLength", "Column length [m]", 0.3), field("ColumnVoidFraction", "Void fraction", 0.4), field("ParticleDensity", "Particle density [kg/m³]", 1000),
               field("ParticleDiameter", "Particle diameter [m]", 0.001), field("DynamicViscosity", "Dynamic viscosity [Pa s]", 1e-5),
               field("BoundaryCondition", "Boundary condition", "InletPressureInletVelocity", ["InletPressureInletVelocity", "InletPressureOutletPressure", "InletVelocityOutletPressure", "FixedVelocity", "FixedPressureInletVelocity"]),
               field("InletPressure", "Inlet pressure [Pa]"), field("OutletPressure", "Outlet pressure [Pa]"),
               field("ColumnEntranceVelocity", "Inlet velocity [m/s]"), field("PressureGradient", "Pressure gradient [Pa/m]", 0)],
    "liquid": [field("LiquidDensity", "Liquid density [kg/m³]", 1000), field("pHMode", "pH mode", "Fixed", ["Fixed", "HPlus", "OHMinus"]),
               field("pHValue", "Fixed pH", 7), field("pKw", "pKw", 14), field("pHComponent", "Acid/base component index", 0)],
    "thermal": [field("energyBalance", "Enable energy balance", False, kind="boolean"), field("InfluxTemperature", "Feed temperature [K] (blank inherits)"),
                field("wallDensity", "Wall density [kg/m³]", 0), field("gasThermalConductivity", "Gas conductivity [W/m/K]", 0),
                field("wallThermalConductivity", "Wall conductivity [W/m/K]", 0), field("heatTransferGasSolid", "Gas–solid transfer [W/m²/K]", 0),
                field("heatTransferGasWall", "Gas–wall transfer [W/m²/K]", 0), field("heatTransferWallExternal", "Wall–ambient transfer [W/m²/K]", 0),
                field("heatCapacityGas", "Gas heat capacity [J/kg/K]", 1), field("heatCapacitySolid", "Solid heat capacity [J/kg/K]", 1), field("heatCapacityWall", "Wall heat capacity [J/kg/K]", 1),
                field("liquidThermalConductivity", "Liquid conductivity [W/m/K]", 0.6), field("heatTransferLiquidSolid", "Liquid–solid transfer [W/m²/K]", 0),
                field("heatTransferLiquidWall", "Liquid–wall transfer [W/m²/K]", 0), field("heatCapacityLiquid", "Liquid heat capacity [J/kg/K]", 4180)],
    "component": [field("GasPhaseMolFraction", "Feed mole fraction", 0), field("LiquidPhaseConcentration", "Inlet concentration [mol/m³]", 0),
                  field("InitialLiquidPhaseConcentration", "Initial concentration [mol/m³]", 0), field("MassTransferCoefficient", "Mass transfer [1/s]"),
                  field("AxialDispersionCoefficient", "Axial dispersion [m²/s]"), field("MolecularWeight", "Molecular weight [kg/mol]"),
                  field("HeatOfAdsorption", "Heat of adsorption [J/mol]"), field("referenceTemperature", "Reference temperature [K]"),
                  field("nonIsothermal", "Non-isothermal isotherm", False, kind="boolean")],
    "geometry": [field("Type", "Geometry", "PackedBed", ["PackedBed", "Monolith"]), field("ColumnVoidFraction", "Geometry void fraction"),
                 field("ParticleDiameter", "Geometry particle diameter [m]"), field("InternalDiameter", "Internal diameter [m]"), field("OuterDiameter", "Outer diameter [m]"),
                 field("ChannelShape", "Channel shape", "circular", ["circular", "square", "triangular", "hexagonal"]),
                 field("InternalChannelDimension", "Channel dimension [m]", 0.001), field("NumberOfChannels", "Number of channels", 10),
                 field("WashcoatThickness", "Washcoat thickness [m]", 0.0001), field("WashcoatVolumePerChannelVolume", "Washcoat / channel volume"), field("ForchheimerCoefficient", "Forchheimer coefficient [1/m]", 0)],
}

_KEYS = "SimulationType DisplayName Components Adsorbents ColumnSections Geometry PhysisorptionSites ChemisorptionSites Name Type Parameters RateEquation RateEquationParameters MaximumLoading HeatOfAdsorption Reactions Phase Style Site Reactants Products Stoichiometry Kinetics ForwardOrders BackwardOrders RateLimitTime CarrierGas FileName ReadColumnFile ColumnDistances DebugForceMultibed MPDSettings SwingAdsorptionPhases Adsorbent AdsorbentLength Length InterfaceLength ComponentParameters MixFraction ColumnPressure ColumnLoading ColumnError forwardRateCoefficient forwardActivationEnergy equilibriumConstant gibbsFreeEnergy".split()
_KEYS += [item["key"] for group in SCHEMA.values() if isinstance(group, list) for item in group if isinstance(item, dict)]
KEYS = {key.lower(): key for key in _KEYS}
ENUMS = {item["key"]: item["options"] for group in SCHEMA.values() if isinstance(group, list) for item in group if isinstance(item, dict) and item["options"]}
ENUMS.update(SimulationType=TYPES + ["Fitting", "Test"], Type=list(ISOTHERMS) + ["PackedBed", "Monolith"],
             RateEquation=list(KINETICS), Phase=["Physisorbed", "Chemisorbed", "PoreConcentration"],
             Style=["GeneralPowerLaw", "LangmuirHinshelwood", "LangmuirHinshelwoodHougenWatson"])




def canonicalize(value, names=False):
    if isinstance(value, list):
        return [canonicalize(item) for item in value]
    if not isinstance(value, dict):
        return value
    result = {}
    for old_key, item in value.items():
        key = old_key if names else KEYS.get(old_key.lower(), old_key)
        if key in result:
            raise ValueError(f"Duplicate input key: {key}")
        item = copy.deepcopy(item) if key == "RateEquationParameters" else canonicalize(item, key in ("Components", "ComponentParameters"))
        if isinstance(item, str):
            item = next((option for option in ENUMS.get(key, []) if option.lower() == item.lower()), item)
        result[key] = item
    return result


def normalize_input(document, directory=None):
    if not isinstance(document, dict):
        raise ValueError("A simulation must be a JSON object.")
    result = canonicalize(document)
    kind = result.get("SimulationType", "Breakthrough")
    if kind not in TYPES or str(result.get("MixturePredictionMethod", "IAST")).upper() == "MPD" or "MPDSettings" in result:
        raise ValueError("Ruptura Lab supports breakthrough, swing adsorption, and analytical mixture prediction. MPD and fitting are excluded.")
    if isinstance(result.get("Components"), dict):
        result["Components"] = list(result["Components"].values())
    if not isinstance(result.get("Components"), list) or not result["Components"]:
        raise ValueError("Add at least one component.")
    geometry = result.get("Geometry")
    if isinstance(geometry, dict) and "Type" not in geometry:
        # Older form resets removed this required discriminator while showing a
        # default geometry. Preserve channel geometries when repairing that input.
        channel_keys = {"ChannelShape", "InternalChannelDimension", "NumberOfChannels",
                        "WashcoatThickness", "WashcoatVolumePerChannelVolume"}
        geometry["Type"] = "Monolith" if channel_keys.intersection(geometry) else "PackedBed"
    if directory:
        result = resolve_files(result, Path(directory))
    return result


def resolve_files(value, directory):
    if isinstance(value, list):
        return [resolve_files(item, directory) for item in value]
    if not isinstance(value, dict):
        return value
    return {key: str((directory / item).resolve()) if key in ("ReadColumnFile", "FileName") and isinstance(item, str) and item else resolve_files(item, directory)
            for key, item in value.items()}


def apply_changes(original, before, after):
    if before == after:
        return copy.deepcopy(original)
    if isinstance(before, dict) and isinstance(after, dict):
        result = copy.deepcopy(original) if isinstance(original, dict) else {}
        for key in before.keys() - after.keys():
            result.pop(key, None)
        for key, value in after.items():
            if key not in before or before[key] != value:
                result[key] = apply_changes(result.get(key), before.get(key), value)
        return result
    if isinstance(original, list) and isinstance(before, list) and isinstance(after, list) and len(original) == len(before) == len(after):
        return [apply_changes(a, b, c) for a, b, c in zip(original, before, after)]
    return copy.deepcopy(after)


def presets():
    gas = {"SimulationType": "Breakthrough", "DisplayName": "CO₂ / N₂ breakthrough", "Temperature": 300.0,
           "MixturePredictionMethod": "IAST", "BoundaryCondition": "InletPressureInletVelocity", "InletPressure": 100000.0,
           "ColumnEntranceVelocity": 0.08, "ColumnLength": 0.1, "ColumnVoidFraction": 0.4, "ParticleDensity": 1000.0,
           "ParticleDiameter": 0.001, "DynamicViscosity": 1.8e-5, "BreakthroughIntegrator": "RungeKutta3",
           "SimulationTime": 2.0, "TimeStep": 0.002, "WriteInterval": 0.02, "NumberOfGridPoints": 30,
           "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": 0.4, "ParticleDiameter": 0.001},
           "Components": [{"Name": "N2", "CarrierGas": True, "GasPhaseMolFraction": 0.85, "MolecularWeight": 0.028},
                          {"Name": "CO2", "GasPhaseMolFraction": 0.15, "MassTransferCoefficient": 0.08, "AxialDispersionCoefficient": 1e-6,
                           "MolecularWeight": 0.044, "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [2.0, 2e-5]}]}]}
    liquid = copy.deepcopy(gas)
    liquid.update(DisplayName="Liquid pH-dependent adsorption", FluidPhase="Liquid", LiquidDensity=997.0, pHMode="Fixed", pHValue=6.5,
                  MixturePredictionMethod="SPI", ColumnEntranceVelocity=0.001, DynamicViscosity=0.00089, TimeStep=0.05, SimulationTime=20.0,
                  WriteInterval=0.5)
    liquid["Components"] = [{"Name": "Solute", "LiquidPhaseConcentration": 2.0, "InitialLiquidPhaseConcentration": 0.0,
                              "MassTransferCoefficient": 0.08, "AxialDispersionCoefficient": 1e-7,
                              "PhysisorptionSites": [{"Type": "pH-Langmuir", "Parameters": [2.0, 0.5, 7.0]}]}]
    mixture = {key: copy.deepcopy(gas[key]) for key in ("Temperature", "MixturePredictionMethod", "Components")}
    mixture.update(SimulationType="MixturePrediction", DisplayName="CO₂ equilibrium isotherm", PressureStart=1.0,
                   PressureEnd=1e6, NumberOfPressurePoints=100, PressureScale="Log")
    swing = copy.deepcopy(gas)
    del swing["SimulationTime"]  # the phase durations set the simulated time
    swing.update(SimulationType="SwingAdsorption", DisplayName="Pressure swing adsorption", SwingAdsorptionPhases=[
        {"Name": "Adsorption", "InletPressure": 100000.0, "Duration": 1.0},
        {"Name": "Regeneration", "InletPressure": 20000.0, "Duration": 1.0}])
    chem = copy.deepcopy(gas)
    chem["DisplayName"] = "CO₂ chemisorption"
    chem["Components"][1]["ChemisorptionSites"] = [{"Type": "Langmuir", "Parameters": [0.8, 2.5e-5], "RateEquation": "FirstOrder",
        "RateEquationParameters": {"MassTransferCoefficient": 0.04}, "MaximumLoading": 0.8, "HeatOfAdsorption": 42000.0}]
    multibed = copy.deepcopy(gas)
    multibed["DisplayName"] = "Layered adsorption column"
    multibed["Adsorbents"] = [{"Name": "Front bed", "ComponentParameters": {"CO2": {"PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [2.0, 2e-5]}]}}},
                              {"Name": "Polishing bed", "ComponentParameters": {"CO2": {"PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [4.0, 4e-5]}]}}}]
    multibed["ColumnSections"] = [{"Adsorbent": "Front bed", "Length": 0.05}, {"InterfaceLength": 0.0}, {"Adsorbent": "Polishing bed", "Length": 0.05}]
    result = {"breakthrough": gas, "liquid": liquid, "mixture": mixture, "swing": swing, "chemisorption": chem, "multibed": multibed}
    for key in PAPER_PRESETS:
        result[key] = json.loads((Path(__file__).with_name("presets") / f"{key}.json").read_text())
    for key in EXAMPLE_ORDER:
        value = result[key]
        value["DisplayName"] = f"Example {PRESET_INFO[key]['number']}: {value['DisplayName']}"
    return result


def migrate_example_catalog(state):
    """Update untouched defaults once while retaining user edits and run associations."""
    if state.get("exampleCatalogVersion", 0) >= 2:
        return False
    catalog = presets()
    simulations = state["simulations"]
    for path in Path(__file__).with_name("presets").joinpath("legacy").glob("*.json"):
        key = path.stem
        original = json.loads(path.read_text())
        for item in list(simulations):
            if item["id"] != f"default-{key}":
                continue
            before = {k: v for k, v in original.items() if k != "DisplayName"}
            after = {k: v for k, v in item["input"].items() if k != "DisplayName"}
            name = item["input"].get("DisplayName", "")
            if before == after and name == original["DisplayName"] and not item.get("lastRunId"):
                simulations.remove(item)
            else:
                item.pop("presetId", None)
                if name == original["DisplayName"]:
                    used = {entry["input"].get("DisplayName") for entry in simulations}
                    number = state.get("nextUntitledIndex", 1)
                    while f"untitled_{number}" in used:
                        number += 1
                    item["input"]["DisplayName"] = f"untitled_{number}"
                    state["nextUntitledIndex"] = number + 1
    for item in simulations:
        name = item["input"].get("DisplayName", "")
        for key, value in catalog.items():
            legacy_name = value["DisplayName"].split(": ", 1)[-1]
            if name == legacy_name:
                item["input"]["DisplayName"] = value["DisplayName"]
                break
    if state["selectedId"] not in {item["id"] for item in simulations}:
        state["selectedId"] = next((item["id"] for item in simulations if item.get("presetId") == "paper-monolith-400ppm-2mmol"), simulations[0]["id"])
    state["exampleCatalogVersion"] = 2
    return True


def _phy(site):
    kind = site.get("Type", "Langmuir")
    parameters = site.get("Parameters", {})
    if isinstance(parameters, dict):
        parameters = [parameters.get(key, default) for key, _, default in ISOTHERMS[kind]]
    return {**site.get("Settings", {}), "Type": kind, "Parameters": parameters}


def _chem(site):
    if isinstance(site.get("Parameters"), list):
        return {key: copy.deepcopy(value) for key, value in site.items() if key != "Id"}
    kind = site.get("Type", "FirstOrder")
    parameters = site.get("Parameters", {})
    result = _phy(site.get("Isotherm", parameters.get("Isotherm", {"Type": "Langmuir"})))
    kinetics = {}
    for key, _, default in KINETICS[kind]:
        old = "rateCoefficient" if kind == "FirstOrder" and key == "MassTransferCoefficient" else key[0].lower() + key[1:]
        kinetics[key] = parameters.get(old, default)
    result.update(RateEquation=kind, RateEquationParameters=kinetics, MaximumLoading=parameters.get("maximumLoading", 0.8),
                  HeatOfAdsorption=parameters.get("heatOfChemisorption", 0.0))
    return result


def migrate_qt(state):
    components = {item["Id"]: item for item in state.get("Components", [])}
    columns = {item["Id"]: item for item in state.get("Columns", [])}
    experiments = []
    for simulation in state.get("Simulations", []):
        kind = simulation.get("Type", "Breakthrough")
        if kind not in TYPES or simulation.get("MixturePredictionMethod", "IAST").upper() == "MPD":
            raise ValueError("This saved state contains an excluded MPD or fitting simulation.")
        column = columns.get(simulation.get("SelectedColumnId"), next(iter(columns.values()), {}))
        feeds = {item["ComponentId"]: item.get("GasPhaseMolFraction", 0.0) for item in simulation.get("ComponentFeeds", [])}
        identifiers = list(feeds)
        reactions = state.get("Reactions", []) if kind != "MixturePrediction" else []
        for reaction in reactions:
            for participant in reaction.get("Reactants", []) + reaction.get("Products", []):
                if isinstance(participant, dict) and participant.get("ComponentId") in components and participant["ComponentId"] not in identifiers:
                    identifiers.append(participant["ComponentId"])
        output_components = []
        for identifier in identifiers:
            component = components[identifier]
            item = {**component.get("Settings", {}), "Name": component["Name"], "GasPhaseMolFraction": feeds.get(identifier, 0.0)}
            for key in ("MassTransferCoefficient", "AxialDispersionCoefficient", "MolecularWeight", "HeatOfAdsorption"):
                if component.get(key) is not None:
                    item[key] = component[key]
            if component.get("ReferenceTemperature") is not None:
                item["referenceTemperature"] = component["ReferenceTemperature"]
            if component.get("NonIsothermal"):
                item["nonIsothermal"] = True
            if component.get("CarrierGas"):
                item["CarrierGas"] = True
            else:
                for key, convert in (("PhysisorptionSites", _phy), ("ChemisorptionSites", _chem)):
                    if component.get(key):
                        item[key] = [convert(site) for site in component[key]]
            output_components.append(item)
        result = {"SimulationType": kind, "DisplayName": simulation.get("DisplayName", "Imported simulation"), "Temperature": column.get("Temperature", 300),
                  "MixturePredictionMethod": simulation.get("MixturePredictionMethod", "SIAST"), "Components": output_components}
        if kind == "MixturePrediction":
            result.update(PressureStart=simulation.get("PressureStart", 0.0001), PressureEnd=simulation.get("PressureEnd", 1e7),
                          NumberOfPressurePoints=simulation.get("PressurePoints", 100), PressureScale=simulation.get("PressureScale", "log"))
        else:
            mapping = {"BoundaryCondition": "BoundaryCondition", "DynamicViscosity": "DynamicViscosity", "VoidFraction": "ColumnVoidFraction",
                       "ParticleDensity": "ParticleDensity", "ParticleDiameter": "ParticleDiameter", "InletPressure": "InletPressure",
                       "PressureGradient": "PressureGradient", "Velocity": "ColumnEntranceVelocity", "Length": "ColumnLength"}
            result.update({target: column[source] for source, target in mapping.items() if source in column})
            integrator = simulation.get("Integrator", "CVODE")
            result.update(BreakthroughIntegrator=integrator, WriteInterval=simulation.get("WriteInterval", 5.0),
                          NumberOfGridPoints=simulation.get("GridPoints", 100))
            if kind != "SwingAdsorption":  # the phase durations set the simulated time of a swing run
                result["SimulationTime"] = "auto" if simulation.get("AutoTime", True) else simulation.get("SimulationTime", 3600.0)
            if integrator != "CVODE":  # CVODE chooses its own internal steps
                result.update(TimeStep=simulation.get("TimeStep", 0.0005), NumberOfInitTimeSteps=simulation.get("InitSteps", 0))
            if column.get("OutletPressure") is not None:
                result["OutletPressure"] = column["OutletPressure"]
            geometry_type = column.get("GeometryType", "PackedBed")
            geometry = {"Type": geometry_type, "OuterDiameter": column.get("OuterDiameter", 0.045)}
            if geometry_type == "Monolith":
                geometry.update({key: column.get(key, default) for key, default in (("ChannelShape", "square"), ("InternalChannelDimension", 0.001), ("NumberOfChannels", 10), ("WashcoatThickness", 0.0001))})
                if column.get("WashcoatVolumePerChannelVolume") is not None:
                    geometry["WashcoatVolumePerChannelVolume"] = column["WashcoatVolumePerChannelVolume"]
            else:
                geometry.update(ColumnVoidFraction=column.get("VoidFraction", 0.4), ParticleDiameter=column.get("ParticleDiameter", 0.001), InternalDiameter=column.get("InternalDiameter", 0.04))
            result["Geometry"] = {**geometry, **column.get("GeometrySettings", {})}
            beds = column.get("Beds", [])
            if column.get("Mode", "single") != "single" and beds:
                adsorbents, sections, length = [], [], 0.0
                for index, bed in enumerate(beds):
                    name = bed.get("Name") or f"Bed {index + 1}"
                    adsorbent = {"Name": name, "ColumnVoidFraction": bed.get("VoidFraction", 0.4), "ParticleDensity": bed.get("ParticleDensity", 1000),
                                 "ParticleDiameter": bed.get("ParticleDiameter", 0.001), **bed.get("Settings", {})}
                    if column["Mode"] == "mixed":
                        adsorbent["MixFraction"] = bed.get("MixFraction", 0.5)
                    else:
                        section = {"Adsorbent": name, "Length": bed.get("Length", 0.15)}
                        length += section["Length"]
                        if bed.get("GridPoints", 0) > 0:
                            section["NumberOfGridPoints"] = bed["GridPoints"]
                        sections.append(section)
                        if index + 1 < len(beds):
                            interface = {"InterfaceLength": bed.get("InterfaceLengthAfter", 0)}
                            length += interface["InterfaceLength"]
                            if interface["InterfaceLength"] > 0 and bed.get("InterfaceGridPointsAfter", 0) > 0:
                                interface["NumberOfGridPoints"] = bed["InterfaceGridPointsAfter"]
                            sections.append(interface)
                    adsorbents.append(adsorbent)
                result["Adsorbents"] = adsorbents
                if column["Mode"] == "multibed":
                    result.update(ColumnSections=sections, ColumnLength=length)
            if column.get("EnergyBalance"):
                result["energyBalance"] = True
                for key in ("InfluxTemperature", "WallDensity", "GasThermalConductivity", "WallThermalConductivity", "HeatTransferGasSolid", "HeatTransferGasWall", "HeatTransferWallExternal", "HeatCapacityGas", "HeatCapacitySolid", "HeatCapacityWall"):
                    if key in column:
                        result[key if key == "InfluxTemperature" else key[0].lower() + key[1:]] = column[key]
            if kind == "SwingAdsorption":
                result["SwingAdsorptionPhases"] = [{key: value for key, value in phase.items() if key != "Id"} for phase in simulation.get("SwingAdsorptionPhases", [])]
            if reactions:
                output_reactions = []
                for reaction in reactions:
                    item = {key: reaction.get(key, default) for key, default in (("Phase", "Physisorbed"), ("Style", "GeneralPowerLaw"), ("Site", 0), ("RateLimitTime", 1e-4))}
                    item["Stoichiometry"] = {}
                    for side, order_key in (("Reactants", "ForwardOrders"), ("Products", "BackwardOrders")):
                        participants = reaction.get(side, [])
                        item[side] = [identifiers.index(part["ComponentId"]) for part in participants]
                        item["Stoichiometry"][side] = [part.get("Stoichiometry", 1.0) for part in participants]
                        if any(part.get("Order") is not None for part in participants):
                            item[order_key] = [part.get("Order") if part.get("Order") is not None else part.get("Stoichiometry", 1) for part in participants]
                    item["Kinetics"] = {key[0].lower() + key[1:]: reaction.get(key, default) for key, default in (("ForwardRateCoefficient", 1), ("ForwardActivationEnergy", 0), ("EquilibriumConstant", 1), ("GibbsFreeEnergy", 0))}
                    output_reactions.append(item)
                result["Reactions"] = output_reactions
        result.update(column.get("Settings", {}))
        result.update(simulation.get("Settings", {}))
        if simulation.get("PHComponentId") in identifiers:
            result["pHComponent"] = identifiers.index(simulation["PHComponentId"])
        if simulation.get("OriginalInput") and simulation.get("OriginalModel"):
            original = simulation["OriginalInput"]
            result = apply_changes(original, simulation["OriginalModel"], result)
            renamed = {item["Name"]: components[f"component-{index + 1}"]["Name"]
                       for index, item in enumerate(original.get("Components", []))
                       if f"component-{index + 1}" in components and item["Name"] != components[f"component-{index + 1}"]["Name"]}
            for bed in result.get("Adsorbents", []):
                if "ComponentParameters" in bed:
                    bed["ComponentParameters"] = {renamed.get(name, name): value for name, value in bed["ComponentParameters"].items()}
                for component in bed.get("Components", []):
                    component["Name"] = renamed.get(component["Name"], component["Name"])
            for reaction in result.get("Reactions", []):
                for side in ("Reactants", "Products"):
                    reaction[side] = [next((new for old, new in renamed.items() if old.lower() == item.lower()), item)
                                      if isinstance(item, str) else item for item in reaction[side]]
        experiments.append({"id": simulation.get("Id", uuid.uuid4().hex), "input": normalize_input(result), "lastRunId": simulation.get("LastRun", {}).get("id")})
    if not experiments:
        raise ValueError("The saved workspace has no simulations.")
    index = min(max(state.get("SelectedSimulationIndex", 0), 0), len(experiments) - 1)
    return {"Format": "RupturaWebLabState", "Version": 1, "simulations": experiments, "selectedId": experiments[index]["id"], "legacyLibrary": state}


def import_document(document, directory=None):
    if document.get("Format") == "RupturaLabState" or ("Simulations" in document and "Columns" in document):
        state = migrate_qt(document)
    elif document.get("Format") == "RupturaWebLabState":
        state = copy.deepcopy(document)
        if not state.get("simulations"):
            raise ValueError("The workspace has no simulations.")
    else:
        item = {"id": uuid.uuid4().hex, "input": normalize_input(document, directory)}
        state = {"Format": "RupturaWebLabState", "Version": 1, "simulations": [item], "selectedId": item["id"]}
    seen = set()
    for item in state["simulations"]:
        if not isinstance(item.get("id"), str) or item["id"] in seen:
            raise ValueError("Simulation identifiers must be unique strings.")
        seen.add(item["id"])
        item["input"] = normalize_input(item["input"], directory)
        if directory and item.get("libraryDefaults"):
            item["libraryDefaults"] = resolve_files(item["libraryDefaults"], Path(directory))
    if directory:
        library = state.get("library", {})
        if not isinstance(library, dict):
            raise ValueError("Workspace libraries must be a JSON object.")
        for entries in library.values():
            for entry in entries:
                entry["data"] = resolve_files(entry["data"], Path(directory))
    if state.get("selectedId") not in seen:
        state["selectedId"] = state["simulations"][0]["id"]
    # Resolve saved modular workspaces once, then detach every simulation. Their
    # complete inputs remain usable without shared definitions or references.
    if "library" in state:
        from .library import ensure_library
        ensure_library(state)
        state.pop("library")
    from .about import example_notes
    for item in state["simulations"]:
        item["input"] = normalize_input(item["input"], directory)
        for key in ("bindings", "libraryDefaults", "templateId"):
            item.pop(key, None)
        item.setdefault("notes", example_notes(item))
        if not isinstance(item["notes"], str):
            raise ValueError("Simulation notes must be text.")
    state["Version"] = 3
    return state
