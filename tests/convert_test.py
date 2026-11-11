"""convert.py: step-based time keys become times in s.

The cases are those of TimeInput.RejectsStepBasedKeysWithTheConvertedTimes in example_inputreader_tests.cpp,
whose error message gives the same replacement, so that the converter and the reader agree.
"""

import importlib.util
import json
import math
from pathlib import Path

import pytest

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


def test_ros34pw2_time_step_and_ramp_are_dropped():
    settings = {"BreakthroughIntegrator": "ROS34PW2", "TimeStep": 0.5, "NumberOfInitTimeSteps": 10,
                "SimulationTime": 10.0, "WriteInterval": 1.0}
    assert convert.convert_time_keys(settings) == {"BreakthroughIntegrator": "ROS34PW2", "SimulationTime": 10.0,
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


# Temperature laws.  The cases follow the replacements that TemperatureInput.RejectsRemovedKeysWithTheirReplacement in
# example_inputreader_tests.cpp expects from the reader.


def langmuir_component(**keys):
    return {"Name": "CO2", **keys, "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [2.0, 1.0e-5]}]}


def test_non_isothermal_affinity_becomes_a_vant_hoff_law_that_inherits_the_heat():
    settings = {"Components": [langmuir_component(HeatOfAdsorption=25000.0, referenceTemperature=318.0,
                                                  nonIsothermal=True)]}
    assert convert.convert_temperature_keys(settings) == {"Components": [{
        "Name": "CO2", "HeatOfAdsorption": 25000.0, "ReferenceTemperature": 318.0,
        "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [2.0, {"Value": 1.0e-5, "Law": "VantHoff"}]}]}]}


def test_disabled_sites_and_the_carrier_placeholder_get_no_law():
    component = {"Name": "He", "nonIsothermal": True, "HeatOfAdsorption": 1.0,
                 "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [0.0, 1.0]},
                                        {"Type": "Langmuir", "Parameters": [1.0, 0.0]}]}
    converted = convert.convert_temperature_keys({"Components": [component]})
    assert converted["Components"][0]["PhysisorptionSites"] == component["PhysisorptionSites"]
    assert "nonIsothermal" not in converted["Components"][0]


def test_isothermal_flag_and_its_reference_temperature_are_dropped():
    settings = {"Components": [langmuir_component(HeatOfAdsorption=25000.0, nonIsothermal=False,
                                                  referenceTemperature=300.0)]}
    assert convert.convert_temperature_keys(settings) == {"Components": [langmuir_component(HeatOfAdsorption=25000.0)]}


def test_gab_energies_become_vant_hoff_laws_and_the_placeholder_heat_is_dropped():
    notes = []
    settings = {"Components": [{"Name": "H2O", "HeatOfAdsorption": 1.0, "nonIsothermal": True,
                                "PhysisorptionSites": [{"Type": "GAB", "Parameters": [35.0, 0.15, 0.51, 6600.0, 0.0]}]}]}
    assert convert.convert_temperature_keys(settings, notes) == {"Components": [{
        "Name": "H2O",
        "PhysisorptionSites": [{"Type": "GAB", "Parameters": [35.0, {"Value": 0.15, "Law": "VantHoff", "Energy": 6600.0},
                                                              0.51]}]}]}
    assert len(notes) == 1 and "HeatOfAdsorption 1.0" in notes[0]


def test_chemisorbed_gab_keeps_a_component_heat_that_physisorption_uses():
    component = {"Name": "H2O", "HeatOfAdsorption": 1.0, "nonIsothermal": True,
                 "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [1.0, 1.0e-5]}],
                 "ChemisorptionSites": [{"Type": "GAB", "Parameters": [2.58, 0.155, 0.871, 6600.0, 0.0],
                                         "HeatOfAdsorption": 50700.0, "MaximumLoading": 2.58}]}
    converted = convert.convert_temperature_keys({"Components": [component]})["Components"][0]
    assert converted["HeatOfAdsorption"] == 1.0
    assert converted["ChemisorptionSites"][0]["Parameters"][1] == {"Value": 0.155, "Law": "VantHoff", "Energy": 6600.0}


def test_general_activation_energies_become_arrhenius_laws():
    site = {"Type": "Langmuir", "Parameters": [2.0, 1.0e-5], "RateEquation": "General",
            "RateEquationParameters": {"AdsorptionRateCoefficient": 2.5e-5, "AdsorptionActivationEnergy": 30000.0,
                                       "DesorptionRateCoefficient": 1.0e-6, "DesorptionActivationEnergy": 0.0,
                                       "CapacityOrder": 1},
            "MaximumLoading": 2.0}
    converted = convert.convert_temperature_keys({"Components": [{"Name": "CO2", "ChemisorptionSites": [site]}]})
    assert converted["Components"][0]["ChemisorptionSites"][0]["RateEquationParameters"] == {
        "AdsorptionRateCoefficient": {"Value": 2.5e-5, "Law": "Arrhenius", "ActivationEnergy": 30000.0},
        "DesorptionRateCoefficient": 1.0e-6, "CapacityOrder": 1}


def test_an_uncentred_law_is_written_at_the_reference_temperature_the_component_keeps():
    # k0 exp(-E/(R T)) = V exp[-E/R (1/T - 1/T_ref)] with V = k0 exp(-E/(R T_ref)).
    site = {"Type": "Langmuir", "Parameters": [2.0, 1.0e-5], "RateEquation": "General",
            "RateEquationParameters": {"AdsorptionRateCoefficient": 2.5e-5, "AdsorptionActivationEnergy": 30000.0},
            "MaximumLoading": 2.0}
    component = {"Name": "CO2", "HeatOfAdsorption": 40000.0, "referenceTemperature": 300.0, "nonIsothermal": True,
                 "ChemisorptionSites": [site]}
    converted = convert.convert_temperature_keys({"Components": [component]})["Components"][0]
    law = converted["ChemisorptionSites"][0]["RateEquationParameters"]["AdsorptionRateCoefficient"]
    assert law["Value"] == 2.5e-5 * math.exp(-30000.0 / (convert.GAS_CONSTANT * 300.0))
    assert converted["ChemisorptionSites"][0]["Parameters"][1] == {"Value": 1.0e-5, "Law": "VantHoff"}


def test_chemisorption_sites_with_different_heats_are_not_converted():
    sites = [{"Type": "Langmuir", "Parameters": [1.0, 1.0e-5], "HeatOfAdsorption": heat, "MaximumLoading": 1.0}
             for heat in (40000.0, 60000.0)]
    with pytest.raises(convert.ConversionError, match="different heats"):
        convert.convert_temperature_keys({"Components": [{"Name": "CO2", "nonIsothermal": True,
                                                          "ChemisorptionSites": sites}]})


def test_override_inherits_the_flag_and_the_reference_temperature():
    settings = {"Components": [langmuir_component(HeatOfAdsorption=25000.0, nonIsothermal=False,
                                                  referenceTemperature=300.0)],
                "Adsorbents": [{"Name": "A", "ComponentParameters": {"CO2": {
                    "nonIsothermal": True, "PhysisorptionSites": [{"Type": "Sips", "Parameters": [1.0, 2.0e-5, 1.2]}]}}}]}
    converted = convert.convert_temperature_keys(settings)
    assert converted["Components"][0] == langmuir_component(HeatOfAdsorption=25000.0)
    # The component no longer keeps its T_ref, so the override states it.
    assert converted["Adsorbents"][0]["ComponentParameters"]["CO2"] == {
        "PhysisorptionSites": [{"Type": "Sips", "Parameters": [1.0, {"Value": 2.0e-5, "Law": "VantHoff"}, 1.2]}],
        "ReferenceTemperature": 300.0}


def test_override_without_sites_cannot_change_the_flag():
    settings = {"Components": [langmuir_component(nonIsothermal=False)],
                "Adsorbents": [{"Name": "A", "ComponentParameters": {"CO2": {"nonIsothermal": True}}}]}
    with pytest.raises(convert.ConversionError, match="no sites"):
        convert.convert_temperature_keys(settings)


def test_reaction_gibbs_energy_becomes_the_heat_of_reaction():
    kinetics = {"forwardRateCoefficient": 2.0, "forwardActivationEnergy": 55000.0, "equilibriumConstant": 1.0,
                "gibbsFreeEnergy": -250000.0}
    converted = convert.convert_temperature_keys({"Reactions": [{"Kinetics": kinetics}]})
    assert converted["Reactions"][0]["Kinetics"] == {
        "ForwardRateCoefficient": {"Value": 2.0, "Law": "Arrhenius", "ActivationEnergy": 55000.0},
        "EquilibriumConstant": {"Value": 1.0, "Law": "VantHoff"}, "HeatOfReaction": 250000.0}

    constant = {"forwardRateCoefficient": 2.0, "forwardActivationEnergy": 0.0, "equilibriumConstant": 80.0,
                "gibbsFreeEnergy": 0.0}
    assert convert.convert_temperature_keys({"Reactions": [{"Kinetics": constant}]}) == {
        "Reactions": [{"Kinetics": {"ForwardRateCoefficient": 2.0, "EquilibriumConstant": 80.0}}]}


def test_converted_input_has_nothing_left_to_convert():
    settings = {"Components": [langmuir_component(HeatOfAdsorption=25000.0, nonIsothermal=True)],
                "Reactions": [{"Kinetics": {"forwardRateCoefficient": 2.0, "forwardActivationEnergy": 0.0,
                                            "equilibriumConstant": 1.0, "gibbsFreeEnergy": -1000.0}}]}
    assert convert.convert_temperature_keys(convert.convert_temperature_keys(settings)) is None


def test_archive_takes_the_first_free_name(tmp_path):
    path = tmp_path / "simulation.json"
    (tmp_path / "_legacy_simulation.json").write_text("time-key original")
    original = json.dumps({"Components": [langmuir_component(HeatOfAdsorption=25000.0, nonIsothermal=True)]})
    path.write_text(original)

    assert convert.main([str(path)]) == 0
    assert (tmp_path / "_legacy_simulation.json").read_text() == "time-key original"
    assert (tmp_path / "_legacy_2_simulation.json").read_text() == original
    assert "nonIsothermal" not in json.loads(path.read_text())["Components"][0]
