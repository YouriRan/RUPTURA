"""convert.py: step-based time keys become times in s.

The cases are those of TimeInput.RejectsStepBasedKeysWithTheConvertedTimes in example_inputreader_tests.cpp,
whose error message gives the same replacement, so that the converter and the reader agree.
"""

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("convert", ROOT / "convert.py")
convert = importlib.util.module_from_spec(spec)
spec.loader.exec_module(convert)


def test_runge_kutta_keeps_its_time_step():
    settings = {"BreakthroughIntegrator": "RungeKutta3", "TimeStep": 0.002, "NumberOfTimeSteps": 750000,
                "WriteEvery": 1000, "PrintEvery": 50000}
    assert convert.convert_time_keys(settings) == {"BreakthroughIntegrator": "RungeKutta3", "TimeStep": 0.002,
                                                   "SimulationTime": 1500.0, "WriteInterval": 2.0}


def test_cvode_drops_time_step_and_ramp():
    settings = {"BreakthroughIntegrator": "CVODE", "TimeStep": 0.002, "NumberOfInitTimeSteps": 10,
                "NumberOfTimeSteps": "auto", "WriteEvery": 1000, "PrintEvery": 50000}
    assert convert.convert_time_keys(settings) == {"BreakthroughIntegrator": "CVODE", "SimulationTime": "auto",
                                                   "WriteInterval": 2.0}


def test_cvode_time_step_alone_is_dropped():
    settings = {"BreakthroughIntegrator": "CVODE", "TimeStep": 0.5, "SimulationTime": 10.0, "WriteInterval": 1.0}
    assert convert.convert_time_keys(settings) == {"BreakthroughIntegrator": "CVODE", "SimulationTime": 10.0,
                                                   "WriteInterval": 1.0}


def test_swing_phases_get_durations_with_the_former_default_time_step():
    settings = {"SimulationType": "SwingAdsorption", "NumberOfTimeSteps": 3000000,
                "SwingAdsorptionPhases": [{"Name": "Adsorption", "NumberOfTimeSteps": 1000000},
                                          {"Name": "Regeneration", "InletPressure": 1.0e5,
                                           "NumberOfTimeSteps": 2000000}]}
    assert convert.convert_time_keys(settings) == {
        "SimulationType": "SwingAdsorption",
        "SwingAdsorptionPhases": [{"Name": "Adsorption", "Duration": 500.0},
                                  {"Name": "Regeneration", "InletPressure": 1.0e5, "Duration": 1000.0}],
        "WriteInterval": 5.0,
    }


def test_duplicate_capitalization_uses_the_key_the_reader_used():
    # The reader takes the exact spelling first, so "NumberOfTimesteps" was ignored.
    settings = {"NumberOfTimeSteps": 112501, "TimeStep": 0.0005, "NumberOfTimesteps": "auto"}
    assert convert.convert_time_keys(settings) == {"SimulationTime": 56.2505, "WriteInterval": 5.0,
                                                   "TimeStep": 0.0005}


def test_time_based_input_is_left_alone():
    assert convert.convert_time_keys({"SimulationTime": "auto", "WriteInterval": 5.0, "TimeStep": 1e-3}) is None


def test_json_file_is_converted_in_place_and_archived(tmp_path):
    path = tmp_path / "simulation.json"
    original = '{\n    "BreakthroughIntegrator": "CVODE",\n    "TimeStep": 2.0,\n    "NumberOfTimeSteps": 750\n}'
    path.write_text(original)

    assert convert.main([str(path)]) == 0
    assert (tmp_path / "_legacy_simulation.json").read_text() == original
    assert json.loads(path.read_text()) == {"BreakthroughIntegrator": "CVODE", "SimulationTime": 1500.0,
                                            "WriteInterval": 20000.0}
    assert path.read_text().startswith('{\n    "')  # the indentation of the input is kept

    # A second run finds nothing to convert and keeps the archive.
    assert convert.main([str(path)]) == 0
    assert (tmp_path / "_legacy_simulation.json").read_text() == original


def test_legacy_text_input_gets_time_keys(tmp_path):
    path = tmp_path / "simulation.input"
    path.write_text("SimulationType Breakthrough\nTimeStep 0.001\nNumberOfTimeSteps 100000\nWriteEvery 500\n"
                    "PrintEvery 1000\n")
    output = tmp_path / "simulation.json"
    assert convert.main([str(path), str(output)]) == 0
    converted = json.loads(output.read_text())
    assert converted["SimulationTime"] == 100.0
    assert converted["WriteInterval"] == 0.5
    assert converted["TimeStep"] == 0.001
    assert not {"NumberOfTimeSteps", "WriteEvery", "PrintEvery"} & set(converted)
