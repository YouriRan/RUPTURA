---
name: ruptura-simulation-input
description: Writes, converts, and checks simulation.json input files for the RUPTURA adsorption simulator - fixed-bed and monolith breakthrough curves, pressure, vacuum, and temperature swing cycles (PSA, VSA, TSA), multicolumn setups, IAST and other mixture-isotherm predictions, and isotherm fitting. Use it whenever the user wants to set up or debug a RUPTURA (Ruptura) run, turn a paper, table, or experiment description into a RUPTURA input, convert isotherm, flow, or column data to RUPTURA's units and conventions, or asks what a simulation.json key means, even if they only mention a breakthrough or PSA simulation in a project that uses RUPTURA.
---

# Ruptura `simulation.json`

Ruptura simulates fixed-bed and monolith adsorption: breakthrough curves, pressure and temperature
swing cycles, mixture isotherms, and isotherm fits. It reads one `simulation.json` from the working
directory and writes its output next to it. This skill turns the user's description, paper, table,
or notes into one valid `simulation.json`.

Two properties of the parser shape everything below:

- Unknown keys are errors. Use only the keys and enum strings in this skill and its reference
  files, and no keys or conventions from other simulators.
- Omitted keys silently take parser defaults, and several defaults are unphysical (a molecular
  weight of 1 kg/mol, heat capacities of 1 J/(kg K)). A file that runs is therefore not necessarily
  right. Write every physically meaningful value explicitly, from the user's data.

## Reference files

Read the files a request needs; most requests need only this file plus units and isotherms.

| File | Read it when |
|---|---|
| [references/units-and-conversions.md](references/units-and-conversions.md) | Any source value is not already in SI on Ruptura's basis: isotherm parameters in bar, flows in mL/min, bulk densities, superficial velocities, Peclet numbers |
| [references/isotherms.md](references/isotherms.md) | Writing `PhysisorptionSites`: isotherm types and parameter order, temperature laws, mixture methods, humidity and pH |
| [references/geometry.md](references/geometry.md) | Writing `Geometry` (packed bed or monolith), an axial dispersion law, or a layered or mixed bed |
| [references/numerics.md](references/numerics.md) | Choosing the integrator, grid, run length, output interval, or solver tolerances |
| [references/swing-adsorption.md](references/swing-adsorption.md) | `SwingAdsorption`: phases, stop criteria, cycles, `MomentumBalance` ports, several columns, buffers, coolant |
| [references/chemisorption-and-reactions.md](references/chemisorption-and-reactions.md) | The model has chemisorption sites, rate equations, or reactions |
| [references/optional-blocks.md](references/optional-blocks.md) | Energy balance, MPD mixture prediction, or a restart from a saved column state |
| [references/templates.md](references/templates.md) | You want the overall structure of a file or a fragment. Never copy their numbers as data |

In a Ruptura source tree, `examples/` holds complete inputs of every simulation type. Use one of the
same kind as a structural model if it helps, never as a source of numbers.

## Workflow

1. Identify the `SimulationType`: `Breakthrough`, `SwingAdsorption`, `MixturePrediction`, or
   `Fitting` ([required keys](#simulation-types-and-their-required-keys)).
2. List every input that type needs. Mark each one as given, derivable by a unit conversion, or
   missing.
3. Missing PHYSICAL data must come from the user. Physical data are compositions, temperatures,
   pressures, flows, dimensions, densities, void fractions, isotherm, kinetic, transport and
   thermal parameters. Ask for all missing items in ONE numbered list grouped by the
   [intake topics](#intake-topics). For each question, say in a few words why it is needed and what
   unit you expect. One complete list respects the user's time more than several rounds.
   Tabulated constants of pure substances, such as molar masses, are not experimental data: fill
   them in and list them under "Assumptions and conversions".
4. You MAY choose NUMERICAL settings yourself: integrator, grid, time step, run length, and output
   interval ([numerics](references/numerics.md)). State every choice under "Assumptions and
   conversions".
5. Estimate physical data only when the user explicitly asks for estimates. Label each estimate and
   its basis, and never present one as measured or published.
6. Convert every value to Ruptura's units ([units](references/units-and-conversions.md)). Keep full
   precision during the conversion and keep the source's significant figures in the result. A value
   derived from several inputs or by an inversion (1/n, particle from bulk density) keeps enough
   digits to give the source values back.
7. Build the smallest file that represents the requested physics. Leave out optional blocks that
   the user did not ask for.
8. Run the [final checklist](#final-checklist). If any item fails or cannot be checked, ask instead
   of answering with a file.
9. If you can run commands, [validate the file with Ruptura](#validating-with-ruptura).

Never hand over a supposedly runnable file that contains guessed values, `null`, `NaN`,
`Infinity`, `"TODO"`, ellipses, comments, trailing commas, duplicate keys, or numbers with unit
suffixes.

### Intake topics

Group the questions of step 3 under these topics, and skip the ones that are complete:

- Goal and SimulationType
- Fluid phase, feed composition (and its basis: mole fraction, partial pressure, concentration)
- Temperature(s), pressure(s), and which pressure or flow the experiment controlled
- Flow rate or velocity (superficial or interstitial; reference state)
- Column: length, inner and outer diameter, packed bed or monolith (channel shape, size, count)
- Adsorbent(s): particle density (or bulk density), particle diameter, bed void fraction, layering
- Isotherm model and parameters, with their equation, units, and reference temperature
- Mass transfer (LDF) and axial dispersion
- Heat of adsorption and thermal data (if non-isothermal)
- Chemisorption, reactions, pH, swing steps (if any)
- Run length, output interval, preferred integrator (optional)

### Validating with Ruptura

Look for the binary with `command -v ruptura` (`pip install .` in the Ruptura sources installs it),
or use `build/ruptura` in a Ruptura source tree. Then, in the directory of the file:

```sh
ruptura --validate simulation.json
```

It parses the file exactly as a run would and prints `Input is valid.` or the parser's error. Fix
each error and validate again. The validator checks keys, types, required inputs, and their
consistency. It cannot tell whether your units, densities, velocities, or isotherm forms are right,
and it does not check the `RungeKutta3` time step, so a valid file is necessary, not sufficient. If
no binary is available, say that the file was not validated.

Do not start the simulation itself unless the user asks: runs take from seconds to hours.

## Output

If you can write files, write `simulation.json` into a new directory with a descriptive name, unless
the user gave a location. Ruptura writes its output next to its input, so every run needs its own
directory. Do not overwrite an existing `simulation.json` without asking.

Then answer with:

- "Assumptions and conversions": a short bullet list covering every conversion, derived value,
  numerical choice, and normalization. Omit the section if there is nothing to report.
- The file: its path if you wrote it, otherwise exactly one fenced `json` block containing the
  complete file. JSON cannot hold comments, so put all explanations outside the file.
- Whether it passed `ruptura --validate`, or that it was not validated.
- Optionally, one line on running it: `cd <directory> && ruptura`.

If the user asks for JSON only, output only the JSON object.

## JSON rules

- The root is one JSON object. Keys are matched case-insensitively, but always write them exactly
  as shown here. EXCEPTION: keys inside `"RateEquationParameters"` are case-sensitive.
- Unknown keys are errors. Use only the keys and enum strings in this skill.
- Booleans are JSON `true`/`false`, not strings. Names and enum values are strings.
- Integers (grid counts, step counts, particle numbers, site indices) must be written without a
  decimal point.
- Any omitted key silently takes the [parser default](#parser-defaults). Defaults are NOT
  experimental data, so write every physically meaningful value explicitly.
- Component order is significant: it fixes the output column order and the MPD dimensions.

## Components

`"Components"` is an array: its order sets the component indices. Each component object has a
required `"Name"` (unique string) and optionally:

- `CarrierGas`: bool. Inert, non-adsorbing gas carrier; it needs no sites.
- `Water`: bool. The water component (at most one, not the carrier, gas only): its partial pressure
  sets the relative humidity of the humidity-dependent isotherms, and its isotherm (usually GAB)
  their water loading q_w ([isotherms](references/isotherms.md)).
- `GasPhaseMolFraction`: feed mole fraction [-] (gas).
- `LiquidPhaseConcentration`: feed concentration [mol/m3] (liquid).
- `InitialLiquidPhaseConcentration`: initial column concentration [mol/m3] (liquid).
- `MolecularWeight`: [kg/mol]. Give it for every gas component of a column run.
- `GasHeatCapacity`: optional [J/(mol K)], the molar heat capacity of the gas at constant pressure.
  Give it for every component or for none: with it the gas takes the heat capacity of its
  composition in place of `heatCapacityGas` [J/(kg K)]. Required for `MomentumBalance` with
  `energyBalance`.
- `MassTransferCoefficient`: k_LDF [1/s].
- `AxialDispersionCoefficient`: D_ax [m2/s]; with a dispersion law
  ([geometry](references/geometry.md)), its velocity-independent (molecular) part D_0.
- `HeatOfAdsorption`: positive magnitude [J/mol], the physisorption heat released (if
  `energyBalance`). An affinity `VantHoff` law without `"Energy"` takes it
  ([isotherms](references/isotherms.md)).
- `ReferenceTemperature`: [K]. The T_ref of the temperature laws of this component that give none.
- `PhysisorptionSites`: array of isotherm sites ([isotherms](references/isotherms.md)).
- `ChemisorptionSites`: array of chemisorption sites
  ([chemisorption](references/chemisorption-and-reactions.md)).
- `FileName`: data file. Fitting only.

Rules:

- Gas column runs (`Breakthrough`, `SwingAdsorption`) need EXACTLY ONE component with
  `"CarrierGas": true`. The column starts filled with pure carrier gas at the inlet pressure; only
  the inlet node holds the feed. No key changes this initial state, except a restart file
  (`ReadColumnFile`, [optional blocks](references/optional-blocks.md)).
- Liquid runs have no carrier. The solvent is implicit and must not be listed.
- Every non-carrier component needs at least one `PhysisorptionSites` entry for every adsorbent,
  unless `MixturePredictionMethod` is `MPD`. This also holds for a purely chemisorbing component or
  a non-adsorbing tracer: give it the disabled site `{"Type": "Henry", "Parameters": [0.0]}` and
  explain that outside the JSON. A site whose first parameter is 0 is disabled.
- Give `MassTransferCoefficient` and `AxialDispersionCoefficient` for every adsorbing component in a
  column run. A value of 0 switches the effect off, so use 0 only deliberately.

## Simulation types and their required keys

**Breakthrough** (default type). Required: `SimulationType`, `DisplayName`, `FluidPhase`,
`Temperature`, `MixturePredictionMethod`, `BreakthroughIntegrator`, `BoundaryCondition` plus that
condition's inputs ([below](#boundary-conditions-flow-direction-and-fluid-phase)),
`DynamicViscosity`, `ParticleDensity`, `ColumnLength` (single bed) or the multistack layout
([geometry](references/geometry.md)), `NumberOfGridPoints`, `SimulationTime`, `WriteInterval`,
`TimeStep` (`RungeKutta3` only), `Geometry`, `Components`. Liquids also need `LiquidDensity`.
`Temperature` is both the initial column temperature and the isotherm temperature.

**SwingAdsorption**. Everything a Breakthrough needs, plus `"SwingAdsorptionPhases"`, a non-empty
ordered array of phases, each with a required `Duration` [s]. The phases run in array order and the
simulated time is the sum of their durations, so omit the root `SimulationTime` (it is an error for
SwingAdsorption). `"NumberOfCycles"` (positive integer) repeats the whole phase list. Read
[references/swing-adsorption.md](references/swing-adsorption.md) before writing any phase.

**MixturePrediction**. Gas only. Required: `SimulationType`, `DisplayName`, `Temperature`,
`PressureStart` and `PressureEnd` [Pa, total pressure], `NumberOfPressurePoints`, `PressureScale`
(`"Log"` or `"Linear"`; Log needs both pressures > 0), `MixturePredictionMethod`, `Components` with
`GasPhaseMolFraction` summing to 1. A carrier is optional and is treated as a real inert diluent.
Use `PressureEnd` >= `PressureStart`.

**Fitting**. Required: `SimulationType` `"Fitting"`, `ColumnPressure`, `ColumnLoading`,
`PressureScale`, `Components`. Optional: `FittingWeighting` (`"Absolute"`, default: minimize the RMSE
of the loading in mol/kg; `"Relative"`: minimize the RMSE of (q_data - q_fit)/q_data, which fits the
low-loading Henry region as closely as the plateau and skips points with zero loading).

- Each component needs `Name`, `FileName` (whitespace-separated columns; lines starting with `#`
  are ignored; pressure in Pa, loading in mol/kg), and `PhysisorptionSites` whose `Parameters`
  arrays have the model's exact length. Only their number matters: the fit starts from random
  values, so zeros are fine. The data file needs more points than the model has parameters.
- `ColumnPressure` and `ColumnLoading` are ONE-based column numbers (defaults 1 and 2). Always set
  them explicitly. `ColumnError` is accepted but not used.
- Do not add column or breakthrough keys.

## Boundary conditions, flow direction, and fluid phase

`"BoundaryCondition"` (default `InletPressureInletVelocity`) and the keys each value requires:

| Value | Requires |
|---|---|
| `InletPressureInletVelocity` | `InletPressure` [Pa], `ColumnEntranceVelocity` [m/s] |
| `InletPressureOutletPressure` | `InletPressure`, `OutletPressure` [Pa] |
| `InletVelocityOutletPressure` | `ColumnEntranceVelocity`, `OutletPressure` |
| `FixedVelocity` | `ColumnEntranceVelocity`, plus `InletPressure` or `OutletPressure` |
| `FixedPressureInletVelocity` | `InletPressure`, `PressureGradient` [Pa/m, signed; negative for forward flow]. The velocity follows from the pressure-drop law. `InletPressure + PressureGradient x ColumnLength` must stay > 0. |
| `MomentumBalance` | `OutletPressure` (the pressure the column starts at, at rest), and for a Breakthrough `ColumnEntranceVelocity`. See below. |

`MomentumBalance` is gas only, needs `InletCondition` `Danckwerts`, and takes no `FlowDirection`.
The pressure of each cell is that of its gas and the velocity follows from the pressure-drop law
between cells, so that the column can be filled and emptied through either end. Use
`BreakthroughIntegrator` `ROS34PW2` (CVODE is accepted only with `CVODELinearSolver` `"KLU"`). In
SwingAdsorption every phase sets its ends with ports instead
([swing adsorption](references/swing-adsorption.md)). With `energyBalance` give every component
`GasHeatCapacity` [J/(mol K)].

Choose the value that matches what the experiment controlled. A flow-controlled feed with a known
inlet pressure is `InletPressureInletVelocity`; a known outlet (back-pressure) is
`InletVelocityOutletPressure`. `FixedVelocity` ignores the change in velocity caused by uptake, so
use it only when the user wants that simplification.

`"FlowDirection"` (default `"Forward"`): the end at which the feed enters. `"Reverse"` feeds the
column at z = L: the inlet keys above then apply at z = L and the outlet keys at z = 0, and the
velocities in `column.data` are negative. The column state is not mirrored: the gas that leaves a
reverse step is that of the first cell, the second node in the output files (the node at z = 0
keeps the state the first cell had when the step started). Use it per phase in SwingAdsorption for
counter-current steps; a Breakthrough rarely needs it. If the uptake exceeds what the feed brings,
the flow can also turn back at the outlet end of a forward step: gas then flows back in there, that
of the last cell unless `ProductEndGasPhaseMolFractions` sets it (swing phases).

The pressure drop uses `DynamicViscosity` [Pa s] of the fluid (gas or liquid) together with the
`Geometry`. Pressures must be > 0 and velocities >= 0.

`"FluidPhase"`: `"Gas"` (default) or `"Liquid"`. For a gas, isotherms use partial pressure in Pa.
For a liquid, isotherms use concentration in mol/m3, and the run needs `LiquidDensity` [kg/m3], the
liquid-phase component concentrations, and no carrier. pH (liquid) and humidity (gas) inputs are in
[references/isotherms.md](references/isotherms.md).

Other root key: `"DisplayName"` (label; default `"Column"`).

## Parser defaults

What silently applies when a key is omitted:

Temperature 298 K | InfluxTemperature = Temperature | ColumnVoidFraction 0.4 | ParticleDiameter
1e-3 m | ParticleDensity 1000 kg/m3 | DynamicViscosity 1e-5 Pa s | LiquidDensity 1000 kg/m3 |
ColumnLength 0.3 m | NumberOfGridPoints 100 | TimeStep 5e-4 s | SimulationTime "auto" |
WriteInterval 5 s | NumberOfPressurePoints 100 | PressureScale Log | PressureGradient 0 |
MolecularWeight 1 kg/mol | all transport, heat, and thermal coefficients 0 |
heatCapacityGas/Solid/Wall 1 J/(kg K) | heatCapacityLiquid 4180 | liquidThermalConductivity 0.6 |
SimulationType Breakthrough | FluidPhase Gas | MixturePredictionMethod IAST |
BreakthroughIntegrator RungeKutta3 | BoundaryCondition InletPressureInletVelocity.

Several of these defaults are unphysical, such as MolecularWeight 1 kg/mol and heat capacities of
1 J/(kg K). That is why every material value must be written explicitly.

## Final checklist

Verify silently before answering:

- [ ] The JSON parses: no comments, null, NaN, placeholders, or duplicate keys; integers have no
  decimal point; booleans are true/false.
- [ ] Only keys from this skill, spelled as shown; `RateEquationParameters` keys are exact.
- [ ] The SimulationType has every required key and no keys from other types.
- [ ] Every value is in SI units on the correct basis: interstitial velocity, particle density, bed
  void fraction, isotherm parameters rewritten into Ruptura's equations, positive heats.
- [ ] Gas column run: exactly one carrier, `MolecularWeight` on every component, feed fractions
  >= 0 summing to 1. Liquid: no carrier, `LiquidDensity` set, concentrations >= 0.
- [ ] Every non-carrier component has an isotherm for every adsorbent (unless MPD), each with the
  exact parameter count. No Unilan. Temperature laws only at the positions allowed in
  [isotherms](references/isotherms.md), each with its coefficient key; an affinity `VantHoff` law
  without `Energy` needs a `HeatOfAdsorption`.
- [ ] The mixture method is compatible: SPI for GAB, pH-Langmuir, and the humidity-dependent
  isotherms (with a Water component or `RelativeHumidity`); Langmuir for EI and SEI; a single
  model per site for SCI; no chemisorption with EI; MPD and MixturePrediction are gas only.
- [ ] The `BoundaryCondition` inputs are present; the fixed-gradient pressure stays positive.
- [ ] The `Geometry` is complete and feasible (outer > inner; channel area < monolith area).
- [ ] Multistack: names match, the order matches, one interface between each pair of neighbours,
  grid counts all or none, no `ColumnLength`.
- [ ] Numerical settings are stable (RK3 time step; no `TimeStep` with CVODE or ROS34PW2), the
  output interval resolves the breakthrough, and every swing phase has a `Duration` in s.
- [ ] Optional blocks (energy, chemisorption, reactions, MPD, pH, solver tuning, restart) are
  present only when requested and are complete.
- [ ] "Assumptions and conversions" lists every conversion, estimate, and numerical choice.
