# Optional blocks: energy balance, MPD, restart

Add these only when requested.

## Energy balance

Set `"energyBalance": true` (default false = isothermal). Keys:

- `InfluxTemperature` [K] (feed temperature; default `Temperature`), `heatCapacitySolid`,
  `heatCapacityWall` [J/(kg K)], `wallDensity` [kg/m3], `wallThermalConductivity` [W/(m K)],
  `heatTransferWallExternal` [W/(m2 K)] (0 = adiabatic), and, per fluid phase:
  - gas: `heatCapacityGas`, `gasThermalConductivity`, `heatTransferGasSolid`, `heatTransferGasWall`
  - liquid: `heatCapacityLiquid`, `liquidThermalConductivity`, `heatTransferLiquidSolid`,
    `heatTransferLiquidWall`
- Also give Geometry `InternalDiameter` and `OuterDiameter`, and a `HeatOfAdsorption` for every
  adsorbing component.

Ask for missing thermal data; do not fill in generic values. For the isotherms to follow the
temperature, write their affinities as VantHoff laws ([isotherms.md](isotherms.md)). For a coolant
along the wall and per-phase wall control, see [swing-adsorption.md](swing-adsorption.md).

## MPD

`MixturePredictionMethod` `"MPD"`, gas only. Root `"MPDSettings"`:

```json
{"FileName": "<path relative to simulation.json>", "ReferenceTemperature": 298.0,
 "ReferenceFugacity": 100000.0, "ReferenceFrameworkMass": 1.3284312537e-22,
 "ComponentBounds": [{"Component": "<name>", "NMin": 0, "NMax": 80, "DeltaN": 1}]}
```

- `ReferenceTemperature` [K], `ReferenceFugacity` [Pa], `ReferenceFrameworkMass` [kg, mass of the
  simulated framework].
- Give one bounds entry per non-carrier component, in component order. NMax >= NMin, DeltaN > 0,
  and (NMax - NMin) must be divisible by DeltaN.
- The file has one row per macrostate in C order (last component varies fastest). Lines starting
  with `#` are comments. Each row holds Pi in [0, 1], optionally followed by the conditional mean
  energy <H> [J], with the same column count on every row.
- MPD components need no `PhysisorptionSites`.

## Restart

`"ReadColumnFile": "<column-state JSON>"` initializes the column from a saved Ruptura state
(`column.json`, which every output rewrites, or `failed_state.json`), which must match the grid,
components, sites and adsorbents. The run goes on from the time of that state: `SimulationTime` (or
the phase durations) further, with the outputs appended to the files of the run directory.
`"ResetClock": true` starts a new experiment from the state instead, at t = 0 and with new output
files. Use them only on request.
