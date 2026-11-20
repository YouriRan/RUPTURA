"""Generate the RUPTURA mixture-prediction and breakthrough inputs of the MPD validation (Sharma et al. 2023, MOR).

Templates: examples/MOR-CO2-C3H8/siast and examples/MOR-CO2-C3H8/breakthrough of the repository, changed to the paper:
  - pure-component dual-site Langmuir fits of Table 3 (site 1 small pockets, site 2 large pockets), mol/kg and 1/Pa
  - breakthrough of Sec. 4.2 / Tables 1-2: 10 % CO2, 10 % propane, 80 % helium at 1e5 Pa and 300 K,
    rho_p = 1711.06 kg/m3, eps_b = 0.128, v_in = 0.1 m/s, L = 0.1 m, k = 0.06 1/s; integrated with ROS34PW2
    (Rosenbrock, default tolerances)
MPD runs read the reference distribution written by ../mpd/build_distribution.py; runs whose distribution does not
exist yet are skipped. Rerun to regenerate every simulation directory.

  mixture_equimolar_{MPD,SIAST,IAST}       300 K, y = 0.5/0.5, 1e2-1e8 Pa (Fig. 6 of the paper)
  mixture_pure_{CO2,C3H8}_MPD              300 K, pure components, against the Table 3 fits
  mixture_equimolar_MPD_T350K_ref300K      350 K from the 300 K distribution (temperature reweighting in RUPTURA)
  mixture_equimolar_MPD_T350K_ref350K      350 K from the 350 K distribution (reference)
  breakthrough_{MPD,SIAST,IAST}            Fig. 10a of the paper
  pure_{CO2,C3H8}_MPD_T<T>K_ref<T_ref>K    one-component distributions of the pure windows
                                           (../mpd/build_pure_distribution.py): against the Table 3 fits and,
                                           at 350 K, the temperature reweighting
  mixture_equimolar_MPD_low_loading        300 K, 1e2-400 Pa, from the well-sampled low-loading states only
                                           (../mpd/build_low_loading_distribution.py)
"""

import copy
import json
from pathlib import Path

HERE = Path(__file__).parent
EXAMPLE = HERE.resolve().parents[3] / "examples" / "MOR-CO2-C3H8"  # this file: validation/3.2.2_mpd/sharma/runs
MPD = HERE.parent / "mpd"

TABLE_3 = {
    "CO2": [{"Type": "Langmuir", "Parameters": [1.60, 2.617e-4]}, {"Type": "Langmuir", "Parameters": [3.62, 4.859e-7]}],
    "C3H8": [{"Type": "Langmuir", "Parameters": [1.73, 6.506e-9]}, {"Type": "Langmuir", "Parameters": [1.09, 2.376e-4]}],
}

if __name__ == "__main__":
    mixture = json.loads((EXAMPLE / "siast" / "simulation.json").read_text())
    mixture.update({"Temperature": 300.0, "PressureStart": 1.0e2, "PressureEnd": 1.0e8, "NumberOfPressurePoints": 61})
    pure = {}
    for name, fractions in {"CO2": {"CO2": 1.0}, "C3H8": {"C3H8": 1.0}}.items():
        pure[name] = copy.deepcopy(mixture)
        for component in pure[name]["Components"]:
            component["GasPhaseMolFraction"] = fractions.get(component["Name"], 0.0)
    breakthrough = json.loads((EXAMPLE / "breakthrough" / "simulation.json").read_text())
    breakthrough.update({"Temperature": 300.0, "ColumnVoidFraction": 0.128, "ParticleDensity": 1711.06,
                         "InletPressure": 1.0e5, "ColumnEntranceVelocity": 0.1, "ColumnLength": 0.1,
                         "BreakthroughIntegrator": "ROS34PW2"})
    breakthrough["Geometry"]["ColumnVoidFraction"] = 0.128
    for component in breakthrough["Components"]:
        component["GasPhaseMolFraction"] = 0.8 if component.get("CarrierGas") else 0.1

    runs = {
        "mixture_equimolar_MPD": (mixture, "MPD", 300.0),
        "mixture_equimolar_SIAST": (mixture, "SIAST", 300.0),
        "mixture_equimolar_IAST": (mixture, "IAST", 300.0),
        "mixture_pure_CO2_MPD": (pure["CO2"], "MPD", 300.0),
        "mixture_pure_C3H8_MPD": (pure["C3H8"], "MPD", 300.0),
        "mixture_equimolar_MPD_T350K_ref300K": ({**mixture, "Temperature": 350.0}, "MPD", 300.0),
        "mixture_equimolar_MPD_T350K_ref350K": ({**mixture, "Temperature": 350.0}, "MPD", 350.0),
        "breakthrough_MPD": (breakthrough, "MPD", 300.0),
        "breakthrough_SIAST": (breakthrough, "SIAST", 300.0),
        "breakthrough_IAST": (breakthrough, "IAST", 300.0),
    }
    for name, (template, method, distribution_temperature) in runs.items():
        settings_file = MPD / f"particle_distribution_T{distribution_temperature:.0f}K.json"
        if method == "MPD" and not settings_file.is_file():
            print(f"skipped {name}: no {distribution_temperature:.0f} K distribution yet")
            continue
        simulation = copy.deepcopy(template)
        simulation["DisplayName"] = name
        simulation["MixturePredictionMethod"] = method
        for component in simulation["Components"]:
            if component["Name"] in TABLE_3:
                # MPD has no pure-component isotherms; IAST and SIAST use the Table 3 fits
                component.pop("PhysisorptionSites", None)
                if method != "MPD":
                    component["PhysisorptionSites"] = copy.deepcopy(TABLE_3[component["Name"]])
        if method == "MPD":
            settings = json.loads(settings_file.read_text())
            settings["FileName"] = f"../../mpd/{settings['FileName']}"
            simulation["MPDSettings"] = settings
        directory = HERE / name
        directory.mkdir(exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=2))
        print(f"wrote {name}")

    for name in TABLE_3:
        for temperature, distribution_temperature in [(300.0, 300.0), (350.0, 300.0), (350.0, 350.0)]:
            run = f"pure_{name}_MPD_T{temperature:.0f}K_ref{distribution_temperature:.0f}K"
            settings_file = MPD / f"particle_distribution_T{distribution_temperature:.0f}K_{name}.json"
            if not settings_file.is_file():
                print(f"skipped {run}: no {distribution_temperature:.0f} K distribution yet")
                continue
            simulation = copy.deepcopy(pure[name])
            simulation.update({"DisplayName": run, "MixturePredictionMethod": "MPD", "Temperature": temperature})
            simulation["Components"] = [component for component in simulation["Components"]
                                        if component.get("CarrierGas") or component["Name"] == name]
            for component in simulation["Components"]:
                component.pop("PhysisorptionSites", None)
            settings = json.loads(settings_file.read_text())
            settings["FileName"] = f"../../mpd/{settings['FileName']}"
            simulation["MPDSettings"] = settings
            directory = HERE / run
            directory.mkdir(exist_ok=True)
            (directory / "simulation.json").write_text(json.dumps(simulation, indent=2))
            print(f"wrote {run}")

    settings_file = MPD / "particle_distribution_T300K_low_loading.json"
    if settings_file.is_file():
        simulation = copy.deepcopy(mixture)
        simulation.update({"DisplayName": "mixture_equimolar_MPD_low_loading", "MixturePredictionMethod": "MPD",
                           "PressureEnd": 400.0, "NumberOfPressurePoints": 16})
        for component in simulation["Components"]:
            component.pop("PhysisorptionSites", None)
        settings = json.loads(settings_file.read_text())
        settings["FileName"] = f"../../mpd/{settings['FileName']}"
        simulation["MPDSettings"] = settings
        directory = HERE / "mixture_equimolar_MPD_low_loading"
        directory.mkdir(exist_ok=True)
        (directory / "simulation.json").write_text(json.dumps(simulation, indent=2))
        print("wrote mixture_equimolar_MPD_low_loading")
