# Swing adsorption

Contents: phases | phase keys | cycles and cyclic steady state | balances in phases.data |
MomentumBalance ports | several columns, valves, and buffers | coolant | PSA recipe

A SwingAdsorption run needs everything a Breakthrough needs (SKILL.md), plus
`"SwingAdsorptionPhases"`, a non-empty ordered array of phase objects. The phases run in array order
and the simulated time is the sum of their durations, so omit the root `SimulationTime` (it is an
error for SwingAdsorption).

## Phase keys

- `Name`: optional string.
- `Temperature`: optional [K]; omitted = the value of the previous setting.
- `InletPressure`: optional [Pa]; omitted = the base `InletPressure`.
- `Duration`: required phase length [s] (`"auto"` is not allowed here).
- `GasPhaseMolFractions`: optional object `{component name: mole fraction}`, the feed of this phase
  (gas only); components not listed are absent from it; omitted = the `GasPhaseMolFraction` values
  of `Components`.
- `ColumnEntranceVelocity`: optional [m/s]; omitted = the base value.
- `FlowDirection`: optional `"Forward"` (feed at z = 0) or `"Reverse"` (feed at z = L, e.g. a
  counter-current purge); omitted = the root `FlowDirection`.
- `ProductEndGasPhaseMolFractions`: optional, forward phases only (gas): the gas that flows back in
  at z = L when the uptake exceeds what the feed brings; omitted = the gas of the last cell.
- `ClosedEnd`: optional `"FeedEnd"` (z = 0) or `"ProductEnd"` (z = L): zero velocity at that end;
  the other end is open at `InletPressure` and gas leaves or enters there (`GasPhaseMolFractions`).
  Do not combine with `FlowDirection`, `ColumnEntranceVelocity`, or
  `ProductEndGasPhaseMolFractions`; needs `BoundaryCondition` `InletPressureInletVelocity` and
  `InletCondition` `Danckwerts`.
- `PressureRamp`: optional `{"Shape": "Exponential" or "Linear", "TimeConstant": s}`: the inlet
  pressure (at the open end of a closed-end phase) goes from its value at the end of the previous
  phase to `InletPressure` (Linear: within `TimeConstant`). Needs `BoundaryCondition`
  `InletPressureInletVelocity`.
- `StopCriterion`: optional `{"Metric": ..., "Operator": ">=" or "<=", "Value": number, "End":
  "FeedEnd" or "ProductEnd" (default), "Component": name}`: the phase ends at the time at which the
  metric reaches the value, and `Duration` is then the longest it lasts. `Metric` is `Pressure`
  [Pa, of the cell at End], `AveragePressure` [Pa], `MaximumTemperature` or `MinimumTemperature`
  [K, solid], `MoleFraction` (needs `Component`, in the cell at End), or `AmountIn` or `AmountOut`
  [mol/m2 net through End in this phase; of `Component`, or in total]. Several criteria:
  `{"All": [criterion, ...]}` or `{"Any": [criterion, ...]}`.
- `MinimumDuration`: optional [s], with `StopCriterion`: the criteria apply from then on.
- `TimeoutPolicy`: optional, with `StopCriterion`: `"Advance"` (default: go on to the next phase at
  `Duration`) or `"Error"` (stop the run).
- `FeedEnd`, `ProductEnd`: only with `BoundaryCondition` `MomentumBalance`: the ports at z = 0 and
  z = L ([below](#momentumbalance-ports)).
- `ThermalControl`: optional, `energyBalance` only: `{"Temperature": K beyond the outer wall,
  "HeatTransferCoefficient": W/(m2 K) of that wall}` for this phase; omitted = the phase
  `Temperature` and `heatTransferWallExternal`. With the root `"Coolant"` also `"MassFlowRate"`:
  kg/s of the fluid in this phase.
- `Columns`, `Connections`: only with the root `Columns` ([below](#several-columns-valves-and-buffers)).

`GasPhaseMolFractions` and `ProductEndGasPhaseMolFractions` also take the `Name` of a phase (a
string) in place of the object: the mean gas that left the column in that phase the last time it
ran, e.g. a purge with the product of the adsorption phase.

Without a `PressureRamp` the pressure steps to the phase `InletPressure` at the start of a phase,
and the gas is rescaled to it at its composition.

## Cycles and cyclic steady state

- `"NumberOfCycles"` (root, positive integer) repeats the whole phase list.
- `"CyclicSteadyState": {"Tolerance": 1e-3, "ConsecutiveCycles": 1}` (root, optional) ends the run
  before `NumberOfCycles`, which it needs (at least 2, the largest number of cycles), once that many
  cycles in a row end within the relative `Tolerance` of the cycle before, in the state, the amounts
  in the column, and the streams; `cycles.data` lists these changes per cycle.

## Balances in phases.data

`phases.data` holds the gas and adsorbed amounts of each component at the start and end of every
phase [mol/m2 of column cross-section], from which purity and recovery follow: what entered minus
what left in a phase is the change of these amounts.

With `StopCriterion`, `CyclicSteadyState`, a phase name as gas, or `MomentumBalance`, `phases.data`
also holds what entered and left at z = 0 and at z = L per component [mol/m2, net per end], then the
heat that entered through the outer wall and the work that brings gas that left through a pressure
port below the root `"ExhaustPressure"` [Pa] to it, over `"CompressionEfficiency"` (default 1) [both
J/m2].

## MomentumBalance ports

With `BoundaryCondition` `MomentumBalance` (SKILL.md), every phase sets its ends with ports and does
not take `InletPressure`, `ColumnEntranceVelocity`, `FlowDirection`, `ClosedEnd`, or
`PressureRamp`:

- `"FeedEnd"`: port at z = 0, `"ProductEnd"`: port at z = L; an end without a port is closed.
- A port is one of:
  - `{"Mode": "Closed"}`
  - `{"Mode": "Pressure", "Pressure": Pa, optional "PressureRamp": {"Shape": "Linear" or
    "Exponential", "TimeConstant": s}}`: gas flows either way; use a ramp for pressurization and
    blowdown.
  - `{"Mode": "Velocity", "Velocity": m/s into the column, "Pressure": Pa of that gas}` (or
    `"MolarFlux"`: mol/(m2 s) in place of both).
- Gas entering at z = 0 is the phase `GasPhaseMolFractions`, at z = L the
  `ProductEndGasPhaseMolFractions`.
- With `energyBalance` give every component `GasHeatCapacity` [J/(mol K)].

## Several columns, valves, and buffers

For pressure equalization, purge with the gas of another bed, or a product tank:

- The root `"Columns": ["A", "B"]` names identical columns; needs `MomentumBalance`, `ROS34PW2`, and
  the Geometry `InternalDiameter`.
- The keys of a phase hold for every column, and its `"Columns": {"A": {...}}` gives one column its
  own `FeedEnd`, `ProductEnd`, `GasPhaseMolFractions`, `ProductEndGasPhaseMolFractions`,
  `Temperature`, or `ThermalControl`.
- A phase `"Connections"` array holds valves: `{"Name": optional, "From": endpoint, "To": endpoint,
  "FlowModel": ...}` with `FlowModel`:
  - `"Linear"` (default): `"Conductance"` [mol/(s Pa)]
  - `"SquareRoot"`: `"Conductance"` [mol/(s Pa^0.5)]
  - `"Orifice"`: `"Area"` [m2]
  - `"Compressor"`: `"MolarFlow"` [mol/s], optional `"Efficiency"`

  and optional `"Opening"` 0 to 1, `"CheckValve": true` (flow only From to To),
  `"MaximumMolarFlow"` [mol/s].
- An endpoint is `"<column>.FeedEnd"`, `"<column>.ProductEnd"`, or a buffer name; a connected end
  takes no port.
- The root `"Buffers": {"Tank": {"Volume": m3, "Pressure": Pa, "GasPhaseMolFractions": {...},
  optional "Temperature", optional "MaximumPressure": gas above it leaves as the product}}`.
- A `StopCriterion` takes `"Column": name`.
- Each column writes to a directory of its name; `connections.data` and `buffers.data` hold the
  amounts per phase [mol].

## Coolant

`"Coolant": {"MassFlowRate": kg/s, "Volume": m3 of the fluid along the column, optional "Density"
(1000) and "HeatCapacity" (4180), optional "FlowDirection": "CoCurrent" (default) or
"CounterCurrent"}` (root key, `energyBalance` only, also for a Breakthrough): a heat-transfer fluid
along the outer wall, which enters at the temperature beyond the wall and exchanges heat with the
wall cell by cell.

## PSA recipe

For a PSA cycle use closed-end phases with a ramp for pressurization and blowdown, and open phases at
constant pressure for adsorption and purge; the pressure drop must stay small against the lowest
pressure.
