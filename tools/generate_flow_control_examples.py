#!/usr/bin/env python3
"""Generate short main-JSON examples for the agreed flow and thermal controls."""
import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "examples/new_examples/integrated/flow_controls"

def base(liquid=True):
    kind = "liquid_breakthrough" if liquid else "gas_breakthrough"
    d = json.loads((ROOT / f"examples/new_examples/unified/fvm/{kind}/simulation.json").read_text())
    d["ConvectionScheme"] = "WENO3"
    d.update(NumberOfGridPoints=12, SimulationTime=2.0, WriteInterval=0.1)
    d["Solver"].update(MaximumTimeStep=0.1)
    d["InitialConditions"] = {"Type": "Empty" if liquid else "CarrierGas", "Loading": "zero"}
    d["Geometry"]["ParticleSphericity"] = 1.0
    return d

def save(name, d):
    folder = TARGET / name
    folder.mkdir(parents=True, exist_ok=True)
    d["DisplayName"] = name.replace("_", " ")
    d["Output"]["Plots"] = True
    (folder / "simulation.json").write_text(json.dumps(d, indent=2) + "\n")

def main():
    d = base(); d["BoundaryCondition"] = "FixedVelocity"
    save("liquid_constant_velocity", d)
    d = base(); d.update(BoundaryCondition="FixedPressureGradient", PressureGradient=-10000.0,
                         InletPressure=110000.0)
    d["Geometry"]["ParticleSphericity"] = 0.8
    d["ParticleTransport"] = {"Model": "HSDM", "Geometry": "CylinderRadial", "Radius": 0.0005,
                              "RadialCells": 8, "SurfaceDiffusivity": {"A": 1e-10, "B": 0.0}}
    save("liquid_gradient_cylindrical_hsdm", d)
    d = base(); d.update(BoundaryCondition="FixedVelocity", ThermalModel="NonAdiabatic", energyBalance=True,
                         FluidSolidThermalEquilibrium=True, InfluxTemperature=305.0,
                         heatCapacitySolid=900.0, wallDensity=7800.0, heatCapacityWall=450.0,
                         heatTransferLiquidWall=50.0, liquidThermalConductivity=0.6)
    d["Coolant"] = {"Model": "Transient", "Volume": 1e-5, "Density": 1000.0, "HeatCapacity": 4180.0,
                    "InitialTemperature": 298.15, "InletTemperature": 293.15, "MassFlowRate": 0.001,
                    "HeatTransferCoefficient": 50.0, "FlowDirection": "CounterCurrent"}
    save("liquid_common_temperature_coolant", d)
    d = base(False); d["BoundaryCondition"] = "FixedVelocity"
    d["Notes"] = ["Prescribed gas velocity throughout the column; Ergun does not determine velocity. Pressure follows ideal-gas concentration storage."]
    save("gas_constant_velocity", d)
    for preset in ("CarrierGas", "FeedPressurized", "GasPressurized", "Evacuated"):
        d = base(False); d["InitialConditions"] = {"Type": preset, "Loading": "zero"}
        if preset == "GasPressurized": d["InitialConditions"]["GasPhaseMolFractions"] = {"A": 0.0, "B": 1.0}
        if preset == "Evacuated": d["InitialConditions"]["Pressure"] = 1000.0
        save("initial_" + preset.lower(), d)
    d = base(); d.update(SimulationType="SwingAdsorption", NumberOfCycles=2)
    d.pop("SimulationTime")
    d["Columns"] = {"A": {}}
    d["Sources"] = {"wash": {"Pressure": 100000.0, "Temperature": 298.15,
                               "LiquidPhaseConcentrations": {"A": 0.0, "B": 0.0}}}
    d["SwingAdsorptionPhases"] = [
        {"Name": "Adsorption", "Duration": 0.5, "Columns": {"A": {
            "FlowControl": {"Model": "FixedVelocity"},
            "LeftPort": {"Mode": "Velocity", "Source": "feed", "Velocity": 0.01},
            "RightPort": {"Mode": "Velocity", "Source": "outlet", "Velocity": 0.0}}}},
        {"Name": "Reverse wash", "Duration": 0.5, "Columns": {"A": {
            "FlowControl": {"Model": "FixedPressureGradient", "PressureGradient": 1000.0},
            "LeftPort": {"Mode": "Pressure", "Source": "outlet"},
            "RightPort": {"Mode": "Pressure", "Source": "wash"}}}},
        {"Name": "Hold", "Duration": 0.2, "Columns": {"A": {
            "LeftPort": {"Mode": "Closed"}, "RightPort": {"Mode": "Closed"}}}},
    ]
    save("liquid_cycle_velocity_gradient_wash", d)
    print(f"Wrote {len(list(TARGET.glob('*/simulation.json')))} examples")

if __name__ == "__main__":
    main()
