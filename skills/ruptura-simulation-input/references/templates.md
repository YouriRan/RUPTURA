# Templates

STRUCTURE ONLY: the numbers are placeholders. Never copy them as data; every physical value comes
from the user.

Contents: gas breakthrough (packed bed) | liquid breakthrough with a fixed pH | mixture isotherm |
fitting | fragments (swing phases, layered column, chemisorption sites, reaction, energy balance,
MPD)

## Gas breakthrough, packed bed

```json
{
  "SimulationType": "Breakthrough",
  "DisplayName": "CO2/N2 breakthrough",
  "FluidPhase": "Gas",
  "Temperature": 298.15,
  "MixturePredictionMethod": "IAST",
  "BreakthroughIntegrator": "CVODE",
  "BoundaryCondition": "InletPressureInletVelocity",
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.05,
  "DynamicViscosity": 1.8e-05,
  "ParticleDensity": 1100.0,
  "ColumnLength": 0.2,
  "NumberOfGridPoints": 100,
  "SimulationTime": "auto",
  "WriteInterval": 1.0,
  "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": 0.4, "ParticleDiameter": 0.002,
               "InternalDiameter": 0.02, "OuterDiameter": 0.024},
  "Components": [
    {"Name": "He", "CarrierGas": true, "GasPhaseMolFraction": 0.8, "MolecularWeight": 0.0040026},
    {"Name": "CO2", "GasPhaseMolFraction": 0.15, "MolecularWeight": 0.04401,
     "MassTransferCoefficient": 0.1, "AxialDispersionCoefficient": 1.0e-05,
     "HeatOfAdsorption": 25000.0, "ReferenceTemperature": 298.15,
     "PhysisorptionSites": [{"Type": "Langmuir",
                             "Parameters": [5.0, {"Value": 1.0e-04, "Law": "VantHoff"}]}]},
    {"Name": "N2", "GasPhaseMolFraction": 0.05, "MolecularWeight": 0.028014,
     "MassTransferCoefficient": 0.5, "AxialDispersionCoefficient": 1.0e-05,
     "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [3.0, 1.0e-06]}]}
  ]
}
```

## Liquid breakthrough with a fixed pH

```json
{
  "SimulationType": "Breakthrough",
  "DisplayName": "Solute in water",
  "FluidPhase": "Liquid",
  "Temperature": 298.15,
  "LiquidDensity": 997.0,
  "DynamicViscosity": 8.9e-04,
  "pHMode": "Fixed",
  "pHValue": 6.5,
  "MixturePredictionMethod": "SPI",
  "BreakthroughIntegrator": "CVODE",
  "BoundaryCondition": "InletPressureInletVelocity",
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.001,
  "ParticleDensity": 1200.0,
  "ColumnLength": 0.1,
  "NumberOfGridPoints": 100,
  "SimulationTime": "auto",
  "WriteInterval": 10.0,
  "Geometry": {"Type": "PackedBed", "ColumnVoidFraction": 0.4, "ParticleDiameter": 0.001,
               "InternalDiameter": 0.01, "OuterDiameter": 0.012},
  "Components": [
    {"Name": "Solute", "LiquidPhaseConcentration": 2.0, "InitialLiquidPhaseConcentration": 0.0,
     "MassTransferCoefficient": 0.01, "AxialDispersionCoefficient": 1.0e-08,
     "PhysisorptionSites": [{"Type": "pH-Langmuir", "Parameters": [2.0, 0.5, 7.0]}]}
  ]
}
```

## Mixture isotherm

```json
{
  "SimulationType": "MixturePrediction",
  "DisplayName": "CO2/N2 IAST",
  "Temperature": 298.15,
  "PressureStart": 100.0,
  "PressureEnd": 1.0e+06,
  "NumberOfPressurePoints": 100,
  "PressureScale": "Log",
  "MixturePredictionMethod": "IAST",
  "Components": [
    {"Name": "CO2", "GasPhaseMolFraction": 0.15,
     "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [5.0, 1.0e-04]}]},
    {"Name": "N2", "GasPhaseMolFraction": 0.85,
     "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [3.0, 1.0e-06]}]}
  ]
}
```

## Fitting

```json
{
  "SimulationType": "Fitting",
  "DisplayName": "CO2 fit",
  "ColumnPressure": 1,
  "ColumnLoading": 2,
  "PressureScale": "Log",
  "Components": [
    {"Name": "CO2", "FileName": "co2_298K.dat",
     "PhysisorptionSites": [{"Type": "Sips", "Parameters": [0.0, 0.0, 0.0]}]}
  ]
}
```

## Fragments

Merge these into a full file.

Swing phases (root; omit root `SimulationTime`):

```json
"SwingAdsorptionPhases": [
  {"Name": "Adsorption", "Temperature": 298.15, "InletPressure": 100000.0, "Duration": 100.0},
  {"Name": "Blowdown", "InletPressure": 10000.0, "Duration": 250.0}
]
```

Layered column (root; omit `ColumnLength`):

```json
"Adsorbents": [
  {"Name": "front", "ColumnVoidFraction": 0.38, "ParticleDensity": 1150.0,
   "ParticleDiameter": 8.0e-04},
  {"Name": "polisher", "ColumnVoidFraction": 0.42, "ParticleDensity": 900.0,
   "ParticleDiameter": 1.2e-03,
   "ComponentParameters": {"CO2": {"PhysisorptionSites": [
     {"Type": "Langmuir", "Parameters": [2.0, 5.0e-05]}]}}}
],
"ColumnSections": [
  {"Adsorbent": "front", "Length": 0.06, "NumberOfGridPoints": 30},
  {"InterfaceLength": 0.01, "NumberOfGridPoints": 5},
  {"Adsorbent": "polisher", "Length": 0.04, "NumberOfGridPoints": 20}
]
```

Chemisorption sites (inside a component that also has `PhysisorptionSites`):

```json
"ChemisorptionSites": [
  {"Type": "Langmuir", "Parameters": [0.8, 2.5e-05], "MaximumLoading": 0.8,
   "HeatOfAdsorption": 42000.0, "RateEquation": "PseudoNth",
   "RateEquationParameters": {"RateCoefficient": 0.04, "Order": 2.0}},
  {"Type": "Langmuir", "Parameters": [0.55, 1.8e-05], "MaximumLoading": 0.55,
   "HeatOfAdsorption": 52000.0, "RateEquation": "General",
   "RateEquationParameters": {
     "AdsorptionRateCoefficient": {"Value": 2.5e-05, "Law": "Arrhenius", "ActivationEnergy": 30000.0},
     "DesorptionRateCoefficient": 1.0e-06,
     "PoreConcentrationOrder": 1, "CapacityOrder": 1, "DesorptionOrder": 1,
     "FilmMassTransferCoefficient": 0.0015, "PoreDiffusivity": 8.0e-11,
     "UsePoreSurfaceTransport": true}}
]
```

Reaction (root):

```json
"Reactions": [
  {"Phase": "Physisorbed", "Style": "GeneralPowerLaw", "Site": 0,
   "Reactants": ["CO", "O2"], "Products": ["CO2"],
   "Stoichiometry": {"Reactants": [1.0, 0.5], "Products": [1.0]},
   "Kinetics": {"ForwardRateCoefficient": {"Value": 1.0e+05, "Law": "Arrhenius", "ActivationEnergy": 55000.0},
                "EquilibriumConstant": {"Value": 1.0, "Law": "VantHoff"}, "HeatOfReaction": 250000.0}}
]
```

Energy balance (root, gas):

```json
"energyBalance": true, "InfluxTemperature": 298.15,
"heatCapacityGas": 1010.0, "heatCapacitySolid": 900.0, "heatCapacityWall": 500.0,
"wallDensity": 7800.0, "gasThermalConductivity": 0.026, "wallThermalConductivity": 16.0,
"heatTransferGasSolid": 100.0, "heatTransferGasWall": 20.0, "heatTransferWallExternal": 5.0
```

MPD (root, with `"MixturePredictionMethod": "MPD"`):

```json
"MPDSettings": {
  "FileName": "distribution.data", "ReferenceTemperature": 298.0,
  "ReferenceFugacity": 100000.0, "ReferenceFrameworkMass": 1.3284312537e-22,
  "ComponentBounds": [{"Component": "CO2", "NMin": 0, "NMax": 80, "DeltaN": 1},
                      {"Component": "N2", "NMin": 0, "NMax": 80, "DeltaN": 1}]
}
```
