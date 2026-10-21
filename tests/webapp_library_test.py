"""Archived modular workspaces become lossless, independent simulations."""
import copy
import json

import pytest

from ruptura.webapp.config import PAPER_PRESETS, import_document, presets
from ruptura.webapp.library import ensure_library, register_experiment
from ruptura.webapp.runs import StateStore, atomic_json


def modular_workspace(*names):
    state = {"Format": "RupturaWebLabState", "Version": 2, "selectedId": "0",
             "simulations": [{"id": str(i), "input": presets()[name]}
                             for i, name in enumerate(names)]}
    return ensure_library(state)


def assert_independent(state):
    assert state["Version"] == 3
    assert "library" not in state
    for item in state["simulations"]:
        assert not {"bindings", "libraryDefaults", "templateId"}.intersection(item)


@pytest.mark.parametrize("name", list(presets()))
def test_modular_migration_is_lossless_and_idempotent(name):
    archived = modular_workspace(name)
    original = copy.deepcopy(archived)
    state = import_document(archived)
    assert_independent(state)
    assert state["simulations"][0]["input"] == presets()[name]
    assert import_document(state) == state
    assert archived == original


def test_shared_edits_and_local_overrides_are_materialized_then_detached():
    archived = modular_workspace("breakthrough", "breakthrough", "breakthrough")
    first, second, _ = archived["simulations"]
    first["input"]["Temperature"] = 350
    first["input"]["Components"][1]["GasPhaseMolFraction"] = 0.3
    second["input"]["Components"][1]["MolecularWeight"] = 0.05
    species_id = first["bindings"]["components"][1]
    species = next(item for item in archived["library"]["components"] if item["id"] == species_id)
    species["data"]["MolecularWeight"] = 0.045
    column = archived["library"]["columns"][0]
    column["data"]["settings"]["ColumnLength"] = 0.2
    column["data"]["componentParameters"][species_id]["PhysisorptionSites"][0]["Parameters"][0] = 5
    state = import_document(archived)
    assert_independent(state)
    inputs = [item["input"] for item in state["simulations"]]
    assert [item["Components"][1]["MolecularWeight"] for item in inputs] == [0.045, 0.05, 0.045]
    assert inputs[0]["Temperature"] == 350
    assert inputs[0]["Components"][1]["GasPhaseMolFraction"] == 0.3
    assert all(item["ColumnLength"] == 0.2 for item in inputs)
    assert all(item["Components"][1]["PhysisorptionSites"][0]["Parameters"][0] == 5 for item in inputs)
    inputs[0]["Components"][1]["MolecularWeight"] = 0.06
    inputs[0]["Components"][1]["PhysisorptionSites"][0]["Parameters"][0] = 8
    inputs[0]["ColumnLength"] = 0.4
    assert inputs[2]["Components"][1]["MolecularWeight"] == 0.045
    assert inputs[2]["Components"][1]["PhysisorptionSites"][0]["Parameters"][0] == 5
    assert inputs[2]["ColumnLength"] == 0.2
    assert import_document(state) == state


def test_distinct_column_parameters_remain_distinct():
    archived = modular_workspace("breakthrough", "chemisorption")
    first, second = archived["simulations"]
    species = first["bindings"]["components"][1]
    column = next(item for item in archived["library"]["columns"] if item["id"] == second["bindings"]["column"])
    column["data"]["componentParameters"][species]["MassTransferCoefficient"] = 0.2
    state = import_document(archived)
    assert state["simulations"][0]["input"]["Components"][1]["MassTransferCoefficient"] == 0.08
    assert state["simulations"][1]["input"]["Components"][1]["MassTransferCoefficient"] == 0.2


def test_shared_reactions_are_resolved_to_local_participants_and_kinetics():
    archived = modular_workspace("breakthrough", "breakthrough")
    for item in archived["simulations"]:
        item["input"]["Reactions"] = [{"Reactants": [1], "Products": [0], "Kinetics": {"forwardRateCoefficient": 0.5}}]
        register_experiment(archived, item)
    archived["library"]["reactions"][0]["data"]["Kinetics"]["forwardRateCoefficient"] = 2
    state = import_document(archived)
    for item in state["simulations"]:
        assert item["input"]["Reactions"] == [{"Reactants": [1], "Products": [0], "Kinetics": {"forwardRateCoefficient": 2}}]
    state["simulations"][0]["input"]["Reactions"][0]["Kinetics"]["forwardRateCoefficient"] = 7
    assert state["simulations"][1]["input"]["Reactions"][0]["Kinetics"]["forwardRateCoefficient"] == 2


def test_notes_identity_and_history_survive_persisted_upgrade(tmp_path):
    archived = modular_workspace("paper-t15", "paper-monolith-400ppm-2mmol")
    archived.update(selectedId="1", nextUntitledIndex=12, includedPresets=list(PAPER_PRESETS), exampleCatalogVersion=1)
    archived["simulations"][0].update(notes="My observations", lastRunId="saved-run", templateId="paper-t15", presetId="paper-t15")
    atomic_json(tmp_path / ".ruptura-lab" / "state.json", archived)
    state = StateStore(tmp_path).get()
    assert_independent(state)
    assert len(state["simulations"]) == 2
    assert state["selectedId"] == "1"
    assert state["nextUntitledIndex"] == 12
    assert state["simulations"][0]["notes"] == "My observations"
    assert state["simulations"][0]["lastRunId"] == "saved-run"
    assert state["simulations"][0]["presetId"] == "paper-t15"
    assert "10.1021/acs.iecr.4c01392" in state["simulations"][1]["notes"]
    assert "notes" not in state["simulations"][0]["input"]
    assert json.loads((tmp_path / ".ruptura-lab" / "state.json").read_text()) == state
    assert StateStore(tmp_path).get() == state


def test_broken_library_links_are_rejected():
    state = modular_workspace("breakthrough")
    state["library"]["columns"] = []
    with pytest.raises(ValueError, match="Missing columns"):
        import_document(state)


def test_unused_species_do_not_leak_into_resolved_bed_parameters():
    archived = modular_workspace("multibed")
    bed = archived["library"]["columns"][0]["data"]["settings"]["Adsorbents"][0]
    bed["ComponentParameters"]["CH4"] = {"HeatOfAdsorption": 15000}
    bed["Components"] = [{"Name": "CH4", "HeatOfAdsorption": 15000}]
    state = import_document(archived)
    resolved = state["simulations"][0]["input"]["Adsorbents"][0]
    assert "CH4" not in resolved["ComponentParameters"]
    assert resolved["Components"] == []
    assert bed["ComponentParameters"]["CH4"]["HeatOfAdsorption"] == 15000


def test_standalone_imports_remain_independent():
    state = import_document({"Format": "RupturaWebLabState", "simulations": [
        {"id": str(i), "input": presets()["breakthrough"]} for i in range(2)]})
    assert_independent(state)
    assert state["simulations"][0]["input"] == state["simulations"][1]["input"]
    state["simulations"][0]["input"]["Components"][1]["MolecularWeight"] = 0.06
    assert state["simulations"][1]["input"]["Components"][1]["MolecularWeight"] == 0.044
    state["simulations"][0]["notes"] = {}
    with pytest.raises(ValueError, match="notes"):
        import_document(state)
