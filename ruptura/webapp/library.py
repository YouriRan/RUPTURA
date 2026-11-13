"""Compatibility reader for saved modular workspaces.

The resolved input remains a standalone CLI document. libraryDefaults records the
last inherited values so imports can recover the resolved scientific input before
converting it into independent simulations.
"""
from __future__ import annotations

import copy
import uuid

FEED_KEYS = ["GasPhaseMolFraction", "LiquidPhaseConcentration", "InitialLiquidPhaseConcentration", "CarrierGas"]
ADSORPTION_KEYS = ["PhysisorptionSites", "ChemisorptionSites", "MassTransferCoefficient", "AxialDispersionCoefficient",
                   "HeatOfAdsorption", "referenceTemperature", "nonIsothermal"]
COLUMN_KEYS = ["Geometry", "Adsorbents", "ColumnSections", "ColumnLength", "ColumnVoidFraction", "ParticleDensity",
               "ParticleDiameter", "wallDensity", "wallThermalConductivity", "heatTransferGasSolid", "heatTransferGasWall",
               "heatTransferWallExternal", "heatCapacitySolid", "heatCapacityWall", "heatTransferLiquidSolid", "heatTransferLiquidWall"]
GROUPS = ["components", "columns", "reactions"]


def _entry(library, group, data, name):
    for entry in library[group]:
        if entry["data"] == data:
            return entry["id"]
    identifier = uuid.uuid4().hex
    library[group].append({"id": identifier, "name": name, "data": copy.deepcopy(data)})
    return identifier


def register_experiment(state, experiment):
    """Split a standalone input into shared species, column, and reaction entries."""
    source, library = experiment["input"], state["library"]
    identifiers, parameters = [], {}
    for component in source["Components"]:
        species = {key: value for key, value in component.items() if key not in FEED_KEYS + ADSORPTION_KEYS}
        if not isinstance(species.get("Name"), str) or not species["Name"].strip():
            raise ValueError("A component definition needs a Name.")
        identifier = _entry(library, "components", species, species.get("Name", "Component"))
        identifiers.append(identifier)
        parameters[identifier] = {key: component[key] for key in ADSORPTION_KEYS if key in component}
    column = {"settings": {key: source[key] for key in COLUMN_KEYS if key in source}, "componentParameters": parameters}
    name = source.get("DisplayName", "Column").split(": ", 1)[-1]
    column_id = _entry(library, "columns", column, name)
    reactions = []
    for index, reaction in enumerate(source.get("Reactions", [])):
        data = copy.deepcopy(reaction)
        for side in ("Reactants", "Products"):
            data[side] = [source["Components"][part]["Name"] if isinstance(part, int) else part for part in data.get(side, [])]
        reactions.append(_entry(library, "reactions", data, reaction.get("Name", f"Reaction {index + 1}")))
    experiment["bindings"] = {"components": identifiers, "column": column_id, "reactions": reactions}
    experiment["libraryDefaults"] = library_defaults(state, experiment)


def library_defaults(state, experiment):
    library, links = state["library"], experiment["bindings"]

    def get(group, identifier):
        entry = next((entry for entry in library[group] if entry["id"] == identifier), None)
        if entry is None:
            raise ValueError(f"Missing {group} library entry: {identifier}")
        return copy.deepcopy(entry["data"])

    column = get("columns", links["column"])
    result = copy.deepcopy(column["settings"])
    result["Components"] = [{**get("components", identifier), **copy.deepcopy(column["componentParameters"].get(identifier, {}))}
                            for identifier in links["components"]]
    reactions = []
    names = [part["Name"].lower() for part in result["Components"]]
    if len(set(names)) != len(names):
        raise ValueError("Choose only one library definition per component name.")
    # A reusable material may describe more species than a particular experiment.
    # The CLI accepts bed overrides only for species present in its input.
    for bed in result.get("Adsorbents", []):
        if "ComponentParameters" in bed:
            bed["ComponentParameters"] = {name: data for name, data in bed["ComponentParameters"].items() if name.lower() in names}
        if "Components" in bed:
            bed["Components"] = [part for part in bed["Components"] if part["Name"].lower() in names]
    for index, identifier in enumerate(links["reactions"]):
        reaction = get("reactions", identifier)
        for side in ("Reactants", "Products"):
            parts = reaction.get(side, [])
            if any(not isinstance(part, str) or part.lower() not in names for part in parts):
                raise ValueError("A linked reaction uses a component missing from this experiment.")
            previous = experiment["input"].get("Reactions", [])
            string_names = index < len(previous) and any(isinstance(part, str) for part in previous[index].get(side, []))
            reaction[side] = parts if string_names else [names.index(part.lower()) for part in parts]
        reactions.append(reaction)
    if reactions:
        result["Reactions"] = reactions
    return result


def ensure_library(state):
    from .config import apply_changes
    from .about import example_notes

    library = state.setdefault("library", {group: [] for group in GROUPS})
    if not isinstance(library, dict):
        raise ValueError("Workspace libraries must be a JSON object.")
    for group in GROUPS:
        entries = library.setdefault(group, [])
        if not isinstance(entries, list):
            raise ValueError("Workspace libraries must contain lists.")
        seen = set()
        for entry in entries:
            if not isinstance(entry, dict) or not isinstance(entry.get("id"), str) or entry["id"] in seen:
                raise ValueError("Library identifiers must be unique strings within each library.")
            if not isinstance(entry.get("name"), str) or not isinstance(entry.get("data"), dict):
                raise ValueError("Library entries need a name and a JSON definition.")
            seen.add(entry["id"])
            data = entry["data"]
            if group == "components" and not isinstance(data.get("Name"), str):
                raise ValueError("A component definition needs a Name.")
            if group == "columns" and (not isinstance(data.get("settings"), dict) or not isinstance(data.get("componentParameters"), dict)):
                raise ValueError("A column definition needs settings and componentParameters objects.")
    species_ids = {entry["id"] for entry in library["components"]}
    for column in library["columns"]:
        if set(column["data"]["settings"]) - set(COLUMN_KEYS):
            raise ValueError("Column library settings must contain column properties.")
        for identifier, parameters in column["data"]["componentParameters"].items():
            if identifier not in species_ids or not isinstance(parameters, dict):
                raise ValueError("Column adsorption parameters must refer to a component in the library.")
    for experiment in state["simulations"]:
        experiment.setdefault("notes", example_notes(experiment))
        if not isinstance(experiment["notes"], str):
            raise ValueError("Experiment notes must be text.")
        if not experiment.get("bindings"):
            register_experiment(state, experiment)
        else:
            links = experiment["bindings"]
            if not isinstance(links, dict) or not isinstance(links.get("components"), list) or not links["components"] or not isinstance(links.get("reactions"), list):
                raise ValueError("Experiment bindings need component and reaction lists.")
            base = library_defaults(state, experiment)
            old = experiment.get("libraryDefaults", base)
            experiment["input"] = apply_changes(base, old, experiment["input"])
            experiment["libraryDefaults"] = base
    state["Version"] = 2
    return state
