"""Local web Lab integration tests; run with python -m pytest tests/webapp_test.py."""
from __future__ import annotations

import copy
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import threading
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import numpy as np
import pytest

from ruptura.webapp.config import PAPER_PRESETS, PRESET_INFO, SCHEMA, apply_changes, import_document, migrate_qt, normalize_input, presets
from ruptura.webapp.runs import RunManager, StateStore, atomic_json
from ruptura.webapp.plots import BlockIndex, figure, json_figure, metadata, report

ROOT = Path(__file__).resolve().parents[1]
CLI = next((Path(path) for path in [os.environ.get("RUPTURA_CLI_EXECUTABLE"), shutil.which("ruptura"),
            ROOT / "build" / "ruptura", *sorted((ROOT / "build").glob("*/ruptura"))]
            if path and Path(path).is_file()), ROOT / "build" / "ruptura")


@pytest.fixture
def manager(tmp_path):
    if not CLI.is_file():
        pytest.skip("Build the ruptura CLI before running web integration tests.")
    value = RunManager(tmp_path, CLI)
    try:
        yield value
    finally:
        value.close()


def wait_for(manager, identifier, statuses=("completed", "failed", "canceled")):
    deadline = time.monotonic() + 20
    while time.monotonic() < deadline:
        record = manager.get(identifier)
        if record["status"] in statuses:
            return record
        time.sleep(0.02)
    pytest.fail(f"Run did not reach {statuses}: {manager.get(identifier)}")


def small_input(name):
    config = presets()[name]
    if name == "mixture":
        config["NumberOfPressurePoints"] = 5
    else:
        config.update(NumberOfTimeSteps=3, TimeStep=1e-6, WriteEvery=1, PrintEvery=1)
        if name == "swing":
            for phase in config["SwingAdsorptionPhases"]:
                phase["NumberOfTimeSteps"] = 3
    return config


@pytest.mark.parametrize("name", list(presets()))
def test_all_presets_validate_in_cpp(manager, name, tmp_path):
    config = manager.validate(presets()[name])
    assert config["SimulationType"] in ("Breakthrough", "SwingAdsorption", "MixturePrediction")
    assert not list((tmp_path / "simulations").iterdir())


@pytest.mark.parametrize("kind,method", [("Fitting", "IAST"), ("fitting", "SPI"), ("MixturePrediction", "MPD"), ("Breakthrough", "mpd")])
def test_excluded_features_are_rejected(kind, method):
    config = presets()["mixture"]
    config.update(SimulationType=kind, MixturePredictionMethod=method)
    with pytest.raises(ValueError, match="excluded"):
        normalize_input(config)


@pytest.mark.parametrize("path", ["examples/Liquid/pH-dependent/simulation.json", "examples/Liquid/multibed-pH-front/simulation.json",
    "examples/Chemisorption-CO2/breakthrough/simulation.json", "examples/Multibed/CO2-C3H8-gradient/simulation.json",
    "examples/Multibed/CO2-N2-polishing/simulation.json", "examples/Reactions/pore-langmuir-hinshelwood/simulation.json"])
def test_import_preserves_entire_scientific_input(path):
    config = json.loads((ROOT / path).read_text())
    before = copy.deepcopy(config)
    imported = import_document(config)
    assert imported["simulations"][0]["input"] == config
    assert config == before


def test_relative_resources_and_case_insensitive_keys(tmp_path):
    config = presets()["liquid"]
    config["fluidphase"] = config.pop("FluidPhase").lower()
    config["components"] = {"0": config.pop("Components")[0]}
    config["ReadColumnFile"] = "initial.json"
    result = normalize_input(config, tmp_path)
    assert result["FluidPhase"] == "Liquid"
    assert isinstance(result["Components"], list)
    assert result["ReadColumnFile"] == str(tmp_path / "initial.json")


@pytest.mark.parametrize("name", [name for name, value in presets().items() if "Geometry" in value])
def test_missing_geometry_discriminator_is_repaired_without_changing_settings(manager, name):
    config = presets()[name]
    original = copy.deepcopy(config)
    del config["Geometry"]["Type"]
    repaired = manager.validate(config)
    assert repaired == original
    assert "Type" not in config["Geometry"]


def test_saved_example_with_reset_geometry_is_repaired_on_open(tmp_path):
    original = {"id": "default-paper-t15", "presetId": "paper-t15", "notes": "My observations",
                "lastRunId": "saved-run", "input": presets()["paper-t15"]}
    broken = copy.deepcopy(original)
    del broken["input"]["Geometry"]["Type"]
    path = tmp_path / ".ruptura-lab" / "state.json"
    atomic_json(path, {"Format": "RupturaWebLabState", "Version": 3, "selectedId": original["id"],
                      "simulations": [broken], "includedPresets": list(PAPER_PRESETS), "exampleCatalogVersion": 2})
    state = StateStore(tmp_path).get()
    assert state["simulations"] == [original]
    assert json.loads(path.read_text()) == state
    assert StateStore(tmp_path).get() == state


def test_state_is_atomic_persistent_and_preserves_unknown_settings(tmp_path):
    store = StateStore(tmp_path)
    value = store.get()
    value["simulations"][0]["input"]["CVODEMaximumTimeStep"] = 0.013
    value["simulations"][0]["input"]["IASTMethod"] = "NestedLoopBisection"
    store.save(value)
    restored = StateStore(tmp_path)
    assert restored.get() == value
    invalid = copy.deepcopy(value)
    invalid["simulations"][0]["input"]["SimulationType"] = "Fitting"
    with pytest.raises(ValueError):
        store.save(invalid)
    assert StateStore(tmp_path).get() == value


def test_paper_defaults_upgrade_existing_workspaces_once(tmp_path):
    path = tmp_path / ".ruptura-lab" / "state.json"
    original = {"id": "my-experiment", "input": presets()["liquid"], "lastRunId": "saved-run"}
    original["input"]["DisplayName"] = "My edited experiment"
    atomic_json(path, {"Format": "RupturaWebLabState", "Version": 1,
                      "simulations": [original], "selectedId": original["id"]})
    store = StateStore(tmp_path)
    value = store.get()
    assert value["selectedId"] == original["id"]
    assert {key: value["simulations"][0][key] for key in original} == original
    assert "bindings" not in value["simulations"][0]
    assert "library" not in value
    assert value["Version"] == 3
    assert {item["presetId"] for item in value["simulations"][1:]} == set(PAPER_PRESETS)
    assert StateStore(tmp_path).get() == value
    value["simulations"][1]["input"]["Temperature"] = 305
    removed = value["simulations"].pop()["id"]
    store.save(value)
    assert StateStore(tmp_path).get() == value
    assert removed not in {item["id"] for item in StateStore(tmp_path).get()["simulations"]}


def test_paper_presets_have_the_published_operating_conditions():
    catalog = presets()
    t15 = catalog["paper-t15"]
    assert t15["energyBalance"] and t15["Temperature"] == 298
    assert t15["InletPressure"] == 1e6 and t15["ColumnLength"] == 1.2
    assert [item["GasPhaseMolFraction"] for item in t15["Components"]] == [0, 0.4, 0.3, 0.3]
    assert t15["heatTransferGasWall"] == 54.26
    bea = catalog["paper-bea-alkanes"]
    assert bea["Temperature"] == 552 and bea["InletPressure"] == 120000
    assert len(bea["Components"]) == 8 and bea["Components"][0]["GasPhaseMolFraction"] == pytest.approx(0.93)
    assert all(item["GasPhaseMolFraction"] == pytest.approx(0.01) for item in bea["Components"][1:])
    for key in PAPER_PRESETS:
        value = catalog[key]
        assert value["BreakthroughIntegrator"] == "CVODE"
        if "monolith" in key:
            assert value["Geometry"]["Type"] == "Monolith"
            assert value["Geometry"]["NumberOfChannels"] == 144
            assert value["Components"][1]["GasPhaseMolFraction"] == pytest.approx(400e-6)
            assert value["TimeStep"] * value["NumberOfTimeSteps"] == 24000
    # Preset edits must not mutate the bundled defaults.
    t15["Temperature"] = 400
    assert presets()["paper-t15"]["Temperature"] == 298


def test_example_catalog_names_and_single_monolith():
    catalog = presets()
    assert [key for key in catalog if "monolith" in key] == ["paper-monolith-400ppm-2mmol"]
    assert sorted(info["number"] for info in PRESET_INFO.values() if "number" in info) == [1, 2, 3, 4]
    for key, value in catalog.items():
        if "number" in PRESET_INFO[key]:
            assert value["DisplayName"].startswith(f"Example {PRESET_INFO[key]['number']}: ")
        else:
            assert not value["DisplayName"].startswith("Example ")
        if key in PAPER_PRESETS:
            assert PRESET_INFO[key]["title"].startswith(f"Example {PRESET_INFO[key]['number']}: ")
        else:
            assert not PRESET_INFO[key]["title"].startswith("Example ")
    for group in SCHEMA.values():
        if isinstance(group, list):
            for item in group:
                if isinstance(item, dict):
                    assert item["description"]


def test_retired_monolith_defaults_preserve_custom_work(tmp_path):
    legacy = ROOT / "ruptura/webapp/presets/legacy"
    untouched = json.loads((legacy / "paper-monolith-400ppm-3mmol.json").read_text())
    edited = json.loads((legacy / "paper-monolith-400ppm-1mmol.json").read_text())
    edited["Temperature"] = 310
    state = {"Format": "RupturaWebLabState", "Version": 1,
             "selectedId": "default-paper-monolith-400ppm-3mmol", "simulations": [
                 {"id": "original-gas", "input": {**presets()["breakthrough"], "DisplayName": "CO₂ / N₂ breakthrough"}},
                 {"id": "default-paper-monolith-400ppm-3mmol", "presetId": "paper-monolith-400ppm-3mmol", "input": untouched},
                 {"id": "default-paper-monolith-400ppm-1mmol", "presetId": "paper-monolith-400ppm-1mmol", "input": edited}]}
    atomic_json(tmp_path / ".ruptura-lab/state.json", state)
    value = StateStore(tmp_path).get()
    assert "default-paper-monolith-400ppm-3mmol" not in {item["id"] for item in value["simulations"]}
    custom = next(item for item in value["simulations"] if item["id"] == "default-paper-monolith-400ppm-1mmol")
    assert custom["input"]["Temperature"] == 310 and custom["input"]["DisplayName"] == "untitled_1"
    assert "presetId" not in custom
    assert value["selectedId"] == "default-paper-monolith-400ppm-2mmol"
    assert value["simulations"][0]["input"]["DisplayName"].startswith("Example 1: ")
    assert StateStore(tmp_path).get() == value


def qt_state():
    return {"Format": "RupturaLabState", "Version": 4, "SelectedSimulationIndex": 0,
            "Components": [{"Id": "component-1", "Name": "N2", "CarrierGas": True},
                           {"Id": "component-2", "Name": "CO2", "MassTransferCoefficient": 0.04,
                            "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": {"qmax": 1, "affinityConstant": 1e-5}}],
                            "ChemisorptionSites": [{"Type": "FirstOrder", "Parameters": {"rateCoefficient": 0.03, "maximumLoading": 0.8, "heatOfChemisorption": 42000},
                                                     "Isotherm": {"Type": "Langmuir", "Parameters": {"qmax": 0.8, "affinityConstant": 2.5e-5}}}]}],
            "Columns": [{"Id": "column-1", "Mode": "single", "Temperature": 300, "InletPressure": 100000, "Velocity": 0.08,
                         "Length": 0.1, "VoidFraction": 0.4, "ParticleDensity": 1000, "ParticleDiameter": 0.001, "DynamicViscosity": 1.8e-5,
                         "BoundaryCondition": "InletPressureInletVelocity", "GeometryType": "PackedBed"}],
            "Simulations": [{"Id": "simulation-1", "Type": "Breakthrough", "DisplayName": "Legacy experiment", "SelectedColumnId": "column-1",
                             "Integrator": "RungeKutta3", "AutoSteps": False, "TimeSteps": 5,
                             "ComponentFeeds": [{"ComponentId": "component-1", "GasPhaseMolFraction": 0.85}, {"ComponentId": "component-2", "GasPhaseMolFraction": 0.15}],
                             "Settings": {"CVODERelativeTolerance": 2e-8}}]}


def test_qt_state_migrates_legacy_chemisorption_and_saved_snapshots(manager):
    old = qt_state()
    migrated = migrate_qt(old)
    config = migrated["simulations"][0]["input"]
    manager.validate(config)
    chemical = config["Components"][1]["ChemisorptionSites"][0]
    assert chemical["Type"] == "Langmuir"
    assert chemical["RateEquation"] == "FirstOrder"
    assert chemical["RateEquationParameters"]["MassTransferCoefficient"] == 0.03
    assert config["CVODERelativeTolerance"] == 2e-8
    assert migrated["legacyLibrary"] == old
    original = copy.deepcopy(config)
    original.pop("PressureGradient", None)
    original["CVODELinearSolver"] = "SPGMR"
    old["Simulations"][0].update(OriginalInput=original, OriginalModel=copy.deepcopy(config))
    assert migrate_qt(old)["simulations"][0]["input"] == original


def test_snapshot_patch_changes_only_edited_fields():
    source = {"Components": [{"Name": "A", "Extra": 2}], "MissingDefault": None}
    before = {"Components": [{"Name": "A", "Default": 1}]}
    after = {"Components": [{"Name": "B", "Default": 1}]}
    result = apply_changes(source, before, after)
    assert result == {"Components": [{"Name": "B", "Extra": 2}], "MissingDefault": None}


@pytest.mark.parametrize("name", list(presets()))
def test_run_snapshot_outputs_and_interactive_plots(manager, name):
    config = small_input(name)
    submitted = manager.start(config, "experiment-1")
    config["DisplayName"] = "Edited after starting"
    record = wait_for(manager, submitted["id"])
    assert record["status"] == "completed", manager.log(record["id"])["text"]
    directory = Path(record["directory"])
    snapshot = json.loads((directory / "simulation.json").read_text())
    assert snapshot["DisplayName"] != "Edited after starting"
    assert (directory / "analysis.ipynb").is_file()
    meta = metadata(directory)
    assert meta["components"]
    if name != "mixture":
        assert meta["energyBalance"] == bool(snapshot.get("energyBalance", False))
    assert figure(directory, {"plot": meta["plots"][0]["id"]}).data
    carriers = [item["Name"] for item in snapshot["Components"] if item.get("CarrierGas")]
    assert all(not any(carrier in trace.name for carrier in carriers)
               for trace in figure(directory, {"plot": meta["plots"][0]["id"]}).data)
    assert figure(directory).layout.paper_bgcolor == "#171f29"
    light = figure(directory, {"theme": "light"})
    assert light.layout.paper_bgcolor == "#ffffff"
    assert light.layout.font.color == "#182e3b"
    if name != "mixture":
        assert meta["frames"] >= 2
        for mode in ("profile", "history"):
            chart = figure(directory, {"plot": mode, "metric": "C", "grid": 0, "frame": 0, "components": [len(meta["components"])-1]})
            assert len(chart.data) == 1
            assert len(chart.data[0].x) > 1
        assert json_figure(directory, {"plot": "temperature"})["data"]
    if name == "liquid":
        assert meta["liquid"]
        assert "pH" in [item["id"] for item in meta["metrics"]]
        assert figure(directory, {"plot": "profile", "metric": "pH", "frame": 0}).data
    assert "plotly" in report(directory).read_text().lower()


def test_cpp_validation_creates_no_notebook_or_results(tmp_path):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(presets()["liquid"]))
    process = subprocess.run([str(CLI), "--validate", str(path)], cwd=tmp_path, capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    assert sorted(item.name for item in tmp_path.iterdir()) == ["input.json"]


@pytest.mark.parametrize("kind", ["breakthrough", "mixture"])
def test_explicit_component_selection_controls_carrier_visibility(manager, kind):
    record = wait_for(manager, manager.start(small_input(kind))["id"])
    assert record["status"] == "completed"
    modes = ["breakthrough", "profile", "history"] if kind == "breakthrough" else ["pure", "loading", "adsorbed"]
    for mode in modes:
        options = {"plot": mode, "metric": "C", "components": [0], "carrier": False}
        selected = figure(record["directory"], options)
        assert len(selected.data) == 1 and "N2" in selected.data[0].name
        options["components"] = [1]
        assert all("N2" not in trace.name for trace in figure(record["directory"], options).data)
        options["components"] = []
        assert not figure(record["directory"], options).data


def test_invalid_input_never_starts_job(manager):
    value = presets()["breakthrough"]
    del value["Geometry"]
    with pytest.raises(ValueError, match="Geometry"):
        manager.start(value)
    assert manager.list() == []


def test_cancellation_queue_and_shutdown(tmp_path):
    executable = tmp_path / "fake-ruptura"
    executable.write_text(f"#!{sys.executable}\nimport json,sys,time\n"
                          "if '--capabilities' in sys.argv: print(json.dumps({'input_validation':True,'sundials':False}))\n"
                          "elif '--validate' in sys.argv: print('valid')\n"
                          "else:\n print('running',flush=True)\n time.sleep(30)\n")
    executable.chmod(0o700)
    manager = RunManager(tmp_path, executable)
    try:
        first = manager.start(presets()["liquid"])
        wait_for(manager, first["id"], ("running",))
        queued = manager.start(presets()["liquid"])
        assert queued["status"] == "queued"
        assert manager.cancel(queued["id"])["status"] == "canceled"
        manager.cancel(first["id"])
        assert wait_for(manager, first["id"])["status"] == "canceled"
        assert manager.log(first["id"])["offset"] > 0
    finally:
        manager.close()


def test_discovery_does_not_trust_manifest_paths_or_restart_processes(tmp_path):
    run = tmp_path / "simulations" / "old-run"
    run.mkdir(parents=True)
    atomic_json(run / "manifest.json", {"id": "wrong-id", "directory": "/private", "status": "running"})
    manager = RunManager(tmp_path, CLI)
    try:
        record = manager.get("old-run")
        assert record["directory"] == str(run)
        assert record["status"] == "interrupted"
        with pytest.raises(FileNotFoundError):
            manager.run_file("old-run", "../../state.json")
    finally:
        manager.close()


def test_stream_index_bounds_memory_and_keeps_last_frame(tmp_path):
    path = tmp_path / "column.data"
    with path.open("w") as handle:
        handle.write("# test output\n")
        for frame in range(5001):
            handle.write(f"{frame} {frame/60} 0 0\n{frame} {frame/60} 1 1\n\n")
    data = BlockIndex(path)
    assert data.total == 5001
    assert len(data.frames) <= 2000
    assert data.read(-1)[0, 0] == 5000
    assert data.read(0).shape == (2, 4)


def test_http_api_session_import_run_and_download(tmp_path):
    from ruptura.webapp.server import LabServer
    server = LabServer(("127.0.0.1", 0), tmp_path, CLI)
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    base = f"http://127.0.0.1:{server.server_port}"

    def request(path, value=None, session=True, headers=None):
        headers = dict(headers or {})
        if session:
            headers["X-Ruptura-Session"] = server.session
        if value is not None:
            headers["Content-Type"] = "application/json"
        return urlopen(Request(base + path, data=json.dumps(value).encode() if value is not None else None, headers=headers), timeout=15)
    try:
        with request("/", session=False) as response:
            assert b"Ruptura Lab" in response.read()
        with request("/api/bootstrap", session=False) as response:
            boot = json.load(response)
            assert boot["session"] == server.session
            assert boot["about"]["citation"]["doi"] == "10.1080/08927022.2023.2202757"
            assert boot["about"]["authors"][0]["name"] == "Youri Ran"
            assert len(boot["about"]["authors"]) == 11
            assert boot["about"]["licenseText"] == (ROOT / "LICENSE.md").read_text()
        with pytest.raises(HTTPError) as error:
            request("/api/runs", {}, session=False)
        assert error.value.code == 403
        with pytest.raises(HTTPError) as error:
            request("/api/runs", {}, headers={"Origin": "https://example.org"})
        assert error.value.code == 403
        with pytest.raises(HTTPError) as error:
            request("/api/import", {"document": {"Format": "RupturaWebLabState", "library": None,
                    "simulations": [{"id": "bad-library", "input": presets()["liquid"]}]}})
        assert error.value.code == 400
        with request("/api/import", {"document": presets()["liquid"]}) as response:
            imported = json.load(response)
        with request("/api/state", imported) as response:
            assert json.load(response)["selectedId"] == imported["selectedId"]
        with request("/api/runs", {"input": small_input("liquid"), "simulationId": imported["selectedId"]}) as response:
            record = json.load(response)
        assert wait_for(server.manager, record["id"])["status"] == "completed"
        with request(f"/api/runs/{record['id']}/plots") as response:
            assert json.load(response)["frames"] >= 2
        with request(f"/api/runs/{record['id']}/plot", {"plot": "profile", "metric": "pH", "frame": 0}) as response:
            assert json.load(response)["data"]
        with request(f"/api/runs/{record['id']}/file/analysis.ipynb") as response:
            assert json.load(response)["nbformat"] == 4
        with pytest.raises(HTTPError) as error:
            request(f"/api/runs/{record['id']}/file/%2E%2E%2Fmanifest.json")
        assert error.value.code == 404
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)

@pytest.mark.parametrize("kind", ["breakthrough", "multibed"])
def test_carrier_can_adsorb_physically_and_chemically(manager, kind):
    config = small_input(kind)
    config["Components"][0].update(MassTransferCoefficient=0.1,
        PhysisorptionSites=[{"Type": "Langmuir", "Parameters": [1, 1e-5]}],
        ChemisorptionSites=[{"Type": "Langmuir", "Parameters": [0.5, 1e-5], "RateEquation": "FirstOrder",
                            "RateEquationParameters": {"MassTransferCoefficient": 0.1}, "MaximumLoading": 0.5}])
    record = wait_for(manager, manager.start(config)["id"])
    assert record["status"] == "completed", manager.log(record["id"])["text"]
    data = np.loadtxt(Path(record["directory"]) / "component_0_N2.data")
    assert metadata(record["directory"])["components"][0]["adsorbing"]
    assert np.max(data[:, 6]) > 0  # Actual physical uptake, rather than ignored UI settings.
    assert np.max(data[:, 8]) > 0  # Chemical uptake also includes the initial carrier.
    assert np.max(data[:, 5]) > 0.99  # The initial column still contains the carrier.


def test_column_axis_switch_uses_position_or_time_with_matching_slider(manager):
    record = wait_for(manager, manager.start(small_input("breakthrough"))["id"])
    directory = record["directory"]
    meta = metadata(directory)
    spatial = figure(directory, {"plot": "profile", "xAxis": "position", "frame": 0, "metric": "C", "components": [1]})
    temporal = figure(directory, {"plot": "profile", "xAxis": "time", "grid": 0, "metric": "C", "components": [1]})
    np.testing.assert_allclose(spatial.data[0].x, meta["positions"])
    np.testing.assert_allclose(temporal.data[0].x, meta["times"])
    assert spatial.layout.xaxis.title.text == "Position [m]"
    assert temporal.layout.xaxis.title.text == "Time [min]"
    with pytest.raises(ValueError, match="x-axis"):
        figure(directory, {"plot": "profile", "xAxis": "invalid"})


@pytest.mark.parametrize("axis", ["position", "time"])
def test_mp4_contains_selected_animation_frames_at_requested_fps(manager, axis):
    from ruptura.webapp.video import export_mp4
    if not shutil.which("ffprobe"):
        pytest.skip("FFprobe is required to inspect the encoded video.")
    config = small_input("breakthrough")
    config["NumberOfGridPoints"] = 3
    record = wait_for(manager, manager.start(config)["id"])
    meta = metadata(record["directory"])
    output = export_mp4(record["directory"], {"xAxis": axis, "fps": 12, "metric": "C", "components": [1], "theme": "light"})
    inspected = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_streams", "-of", "json", str(output)],
                               capture_output=True, text=True, check=True)
    stream = json.loads(inspected.stdout)["streams"][0]
    assert stream["codec_name"] == "h264"
    assert stream["pix_fmt"] == "yuv420p"
    assert (stream["width"], stream["height"]) == (1280, 720)
    assert stream["r_frame_rate"] == "12/1"
    assert int(stream["nb_frames"]) == (meta["frames"] if axis == "position" else len(meta["positions"]))
    assert not list(Path(record["directory"]).glob(".video-*"))
    with pytest.raises(ValueError, match="Select at least one component"):
        export_mp4(record["directory"], {"components": []})
    for fps in (0, 61, 1.5, True, "20", None):
        with pytest.raises(ValueError, match="Frames per second"):
            export_mp4(record["directory"], {"fps": fps})
