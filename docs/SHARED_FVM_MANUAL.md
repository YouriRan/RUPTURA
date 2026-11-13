# Main JSON inputs for breakthrough and adsorption cycles

This integration adds a shared conservative column model to the November 2026
main source. Select it with `"SpatialDiscretization": "FDM"` or `"FVM"`. Breakthrough is
one operation of that model; a cycle schedules repeated operations of the same
column. Both accept gas or incompressible liquid and packed beds or monoliths.

The original main calculation remains the default for existing inputs without
this selection. The shared path also adapts the original pure-isotherm catalogue,
transported pH and `ChemisorptionSites`; general reaction networks are available
for breakthrough and rejected for cycles. Gas-phase MPD is also available for
both shared axial methods, using the original distribution kernel.
See [current status](UNIFIED_COLUMN_STATUS.md) for the verification scope before
moving an existing scientific campaign between paths.

The [flow, thermal and IC update](FLOW_THERMAL_UPDATE.md) defines the four current
boundary modes, pellet sphericity, common-temperature option and initialization
presets. It supersedes the earlier flow/default descriptions where noted.

## Build and first run

```bash
./install.sh --python
build/ruptura examples/new_examples/unified/fvm/gas_breakthrough/simulation.json --validate
build/ruptura examples/new_examples/unified/fvm/gas_breakthrough/simulation.json
```

The installer builds Release code with `-march=native` disabled. This avoids
building for a newer CPU than the one used on HAL. It checks C++23 standard
library support and can install a pinned compiler in a project environment.
Its Zig wrapper explicitly sets `-mcpu=baseline`; disabling the CMake native
flag alone does not override Zig’s automatic host CPU selection.
LAPACK/BLAS, a C compiler and zlib development files must already be available.
An existing SUNDIALS prefix can be supplied with `--sundials-prefix`.

`--validate` checks parsing, geometry, initial states, every phase's initial RHS,
ports, connections and ledger dimensions without running the simulation or
creating results. It is not a proof that a complete cycle reaches CSS.

```bash
build/ruptura /path/to/simulation.json --no-plots
build/ruptura /path/to/simulation.json --cycles 5
build/ruptura /path/to/simulation.json --restart --cycles 10
build/ruptura /path/to/simulation.json --plot-only
```

A fresh run replaces its own result files. Use separate case directories for
different experiments. `Output.Directory` must be a dedicated results folder,
not the input folder or its parent. Checkpoints include all column states,
particle histories, transient coolant states, buffer states, stream ledgers and
schedule position. Restart requires compatible physics, component order,
geometry, solver and schedule. Increasing the cycle count is supported; changing
a completed breakthrough horizon is a different schedule and requires a new run.

## Common input conventions

Names and structure follow main: `Components`, `Geometry`, `Adsorbents`,
`ColumnSections`, `PhysisorptionSites`, `Parameters`, `energyBalance`,
`BreakthroughIntegrator`, and `SwingAdsorptionPhases`. New optional blocks hold
particle transport, coolant, network equipment and cycle controls. Public examples
are independent `simulation.json` files; no user-facing YAML is needed.

| Setting or quantity | Units / meaning |
|---|---|
| `Temperature`, initial and source temperatures | K |
| `InletPressure`, `OutletPressure`, source/port pressure | Absolute Pa |
| `ColumnLength`, diameters, particle radius, coating thickness | m |
| `NumberOfGridPoints` | Number of axial **intervals**, as in main; states = intervals + 1 |
| `ColumnEntranceVelocity`, port `Velocity` | Interstitial m/s; positive left to right |
| `MolarFlowRate` / root `InletMolarFlowRate` | mol/s; use instead of velocity, gas only |
| `VolumetricFlowRate` / root `InletVolumetricFlowRate` | m³/s; use instead of velocity |
| `LiquidPhaseConcentration` | mol/m³ liquid; solvent is implicit |
| Adsorbed loading, site capacity | mol/kg active adsorbent |
| `ParticleDensity` | kg/m³ adsorbent solid |
| Component `GasHeatCapacity` | J/(mol K), must exceed R for gas, including carrier |
| Root `heatCapacityGas`, `heatCapacityLiquid`, `heatCapacitySolid` | J/(kg K) |
| Component `AdsorbedHeatCapacity` | J/(mol K) |
| `HeatOfAdsorption`, site `EquilibriumHeat` | Positive heat release, J/mol |
| `MassTransferCoefficient` | LDF 1/s |
| `AxialDispersionCoefficient`, `SurfaceDiffusivity` | m²/s |
| `Duration`, `SimulationTime`, `WriteInterval`, ramp times | s |

Gas fractions may be fractions or percentages; the adapter normalizes their sum,
matching main's 40/30/30 T15 input. Liquid concentrations are not normalized.
Duplicate JSON keys and unknown top-level/process/particle keys are errors.
Never specify both velocity and a molar or volumetric flow for the same inlet.

Component order is shared by the connected columns. Component maps may be keyed
by the exact component names. Initial composition defaults to the carrier, or to
feed if no carrier is marked. Explicit initial conditions take precedence:

```json
"InitialConditions": {
  "Type": "CarrierGas",
  "Pressure": 100000.0,
  "GasTemperature": 298.15,
  "SolidTemperature": 298.15,
  "WallTemperature": 298.15,
  "Loading": "equilibrium"
}
```

Use `Loading: "zero"` for a clean bed. Scalars or compatible axial profiles are
accepted for temperature and pressure. For liquid use
`LiquidPhaseConcentrations` instead of gas fractions. A gas initial pressure does
not determine a liquid concentration.

## Materials, mixtures and geometry

`PhysisorptionSites` retain main parameter order. For Sips, main's third parameter
is the reciprocal exponent; conversion to the internal form is explicit.
`EquilibriumHeat` and site `referenceTemperature` describe affinity variation.
Without a reference temperature, the affinity is the absolute Arrhenius
prefactor. With it, the affinity is its value at that temperature. Component
`HeatOfAdsorption` independently controls the energy source. This distinction
preserves the MFI two-site equilibrium heats and its constant effective energy
heat. No equilibrium refitting is performed during conversion.

Shared mixture methods are SPI, SCI, IAST, SIAST, SEI, original-main EI and gas-phase MPD. SEI requires
Langmuir/zero sites. SCI requires matching site models; a nonadsorbing carrier is
padded with zero-capacity sites. Available shared pure models include Langmuir,
Langmuir–Freundlich, Sips, Toth, Henry, Freundlich, anti-Langmuir,
modified-Redlich, fixed-pH Langmuir and the imported humidity/water models.
The original main adapter supplies additional pure models and mixed site types,
including BET, Redlich-Peterson, Unilan, O'Brien-Myers, Quadratic, Temkin and
Bingel-Walton. It also supplies transported HPlus/OHMinus pH. Their original
model/fluid restrictions still apply. Liquid tracked species use explicit
isotherms (a zero Henry site represents a nonadsorbing solute), with no gas
carrier flag.

For MPD, set `MixturePredictionMethod: "MPD"` and retain the original
`MPDSettings`: distribution `FileName`, `ReferenceTemperature` (K),
`ReferenceFugacity` (Pa), `ReferenceFrameworkMass` (kg) and `ComponentBounds`.
Bounds follow the non-carrier component order; the last dimension varies
fastest in the file. Each row contains a probability in [0,1], optionally
followed by its conditional mean Hamiltonian in joules. Loadings are converted
from particle counts to mol/kg by dividing by Avogadro's constant and the
represented framework mass. No analytic isotherm is required for MPD.

The shared model passes each cell's local temperature and partial pressures to
the original kernel. It retains the original temperature/fugacity reweighting;
the optional energy column adds its Hamiltonian correction. The distribution
is cached when the material is constructed, and its content digest participates
in restart compatibility. Relative paths use the simulation's directory. Lab
exposes the same `MPDSettings`; cluster bundles include the distribution and
rewrite its `FileName` to a portable path. See [MPD integration notes](MPD_SHARED_MODEL.md)
and the matched `gas_mpd_breakthrough` examples.

`CapacityTemperature` supports exponential capacity variation with a reference
temperature and coefficient. `ExtraParameters.local_c1: true` supplies a smooth
small-driving-force continuation for powered isotherms, useful for low-exponent
Freundlich/HSDM calculations. It changes the local small-concentration law and
must be recorded when comparing a fit with data.

Main's adsorbent override styles are accepted: named `ComponentParameters`
objects or `Components` arrays. Use the original alternating section format:

```json
"ColumnSections": [
  {"Adsorbent": "Front", "Length": 0.04, "NumberOfGridPoints": 20},
  {"InterfaceLength": 0.01, "NumberOfGridPoints": 10},
  {"Adsorbent": "Rear", "Length": 0.05, "NumberOfGridPoints": 20}
]
```

An interface is a finite transition between adsorbents; zero length gives a sharp
interface. Per-section counts construct the axial coordinates. Alternatively use
`AdsorbentLength`, or a `MixFraction` for every material for a uniform mix.
Mixture fractions must sum to one and cannot be combined with sections.

Monolith geometry uses `Geometry.Type: "Monolith"`, `ChannelShape`,
`ChannelDiameter`/`ChannelSize`, `NumberOfChannels` or an open fraction,
`WashcoatThickness`, `WashcoatSolidFraction`, `SubstrateDensity` and
`SubstrateHeatCapacity`. Circular, square, triangular and hexagonal channels are
available. Channel void volume, active coating mass and substrate thermal mass
are accounted separately. HSDM is currently restricted to nonreactive fixed beds.

## Shared gas/liquid chemisorption

Use `Components[].ChemisorptionSites` for FirstOrder, PseudoNth, Avrami, General,
Elovich and DualNth on either axial method and in breakthrough or cycles.
DualNth has two kinetic branches sharing one physical site's equilibrium.
The old ChemicalKinetics block is converted during import; there is no separate
gas-only public chemistry interface. See [the common chemisorption manual](SHARED_CHEMISORPTION.md)
for signed regeneration rates, temperature conventions, loading/fractional rate
bases, Avrami history policies, GUI controls, examples and checkpoint migration.

## Particle transport

```json
"ParticleTransport": {
  "Model": "HSDM",
  "Geometry": "Sphere",
  "Radius": 0.0005,
  "RadialCells": 12,
  "SurfaceDiffusivity": {"A": 5e-11, "B": 0.0},
  "DiffusivityActivationEnergy": {"A": 0.0, "B": 0.0},
  "referenceTemperature": 298.15,
  "FilmResistance": true,
  "FilmCoefficient": {"A": 5e-5, "B": 0.0}
}
```

`Sphere` and `CylinderRadial` use conservative radial shell volumes; cylinder
ends are excluded. The mean loading source is the shell-volume-weighted radial
derivative. The finite-film calculation solves the coupled multicomponent
surface equilibrium. Particle shell states are retained in checkpoints and NPZ
results and can be inspected in the Lab particle plot.

`Model: "LDF"` can use the supplied rate or a diffusivity-derived rate through
`LDFCoefficientSource`. The particle branch supports Arrhenius coefficients.
Powder transport is LDF only. These particle options are physisorption models;
combining them with chemical uptake or reactions is rejected.

## Flow, boundaries and discretization

Each phase has left and right ports with `Mode` set to `Pressure`, `Velocity` or
`Closed`. Gas uses ideal-gas storage and Ergun momentum, or `MomentumModel:
"Darcy"`; monoliths use channel geometry. Pressure sources can reverse flow, so
their temperature and composition are also specified for backflow. Species use
conservative advective boundary fluxes and Danckwerts inlet / zero dispersive
outlet fluxes. A closed boundary has zero normal mass flux.

For a prescribed gas velocity, `InletDensityBasis: "Reservoir"` retains the
cycle donor convention, c = Preservoir/(RTsource). The explicit alternative
`"LocalPressure"` recovers inlet face pressure from the local Darcy/Ergun
relation and uses that pressure for inflowing density. It is useful for the
Haghpanah boundary correspondence check. It cannot be combined with imposed
molar flow. The choice changes the molar feed and must accompany validation
results. Default reservoir and explicit common-molar-feed campaigns remain
distinct, intentional feed specifications.

Liquid solves incompressible hydraulic continuity. Its solute concentrations are
independent of pressure; solvent flow is the same through the column. A closed
end cannot coexist with nonzero imposed through-flow. A liquid connection uses
`FlowModel: "Incompressible"` and a volumetric conductance in m³/(s Pa).
Gas valves use molar conductance in mol/(s Pa). The liquid network also supports
imposed volumetric transfers and finite, variable-volume tanks.

Axial FDM uses nodal reconstruction; FVM reconstructs control-volume averages.
Both use conservative flux differences and boundary half volumes, with the same
column balances, material laws and radial particle model. `ConvectionScheme`
selects `Upwind`, minmod-limited `TVD`, `WENO3` or `WENO5`. WENO uses a three-
or five-point interior stencil and falls back near boundaries/interfaces.
Nonuniform coordinates are supported;
thermal conduction uses local node spacing. Large grid-size jumps can produce
ill-conditioned WENO stencils and are rejected. Scheme availability does not
establish grid convergence for a scientific case.

## Heat and coolant

With `energyBalance: true`, fluid, solid and wall temperatures are solved.
`GasSolidThermalEquilibrium: true` constrains the two phase temperatures for
applicable gas reference cases. Wall end conditions can be adiabatic or fixed
temperature. Phases may apply ambient exchange, jacket/HTF control, heat flux or
electric heat through `ThermalControl`.

```json
"Coolant": {
  "Model": "Transient", "Volume": 0.12, "Density": 1000.0,
  "HeatCapacity": 4180.0, "InitialTemperature": 298.15,
  "InletTemperature": 293.15, "MassFlowRate": 0.5,
  "HeatTransferCoefficient": 500.0, "FlowDirection": "CoCurrent"
}
```

The transient coolant has an axial thermal state and finite heat storage.
`CounterCurrent` reverses its direction. Per-phase thermal control overrides
the global coolant defaults, including in multicolumn inputs. A hold keeps its
stated thermal conditions, so waiting time and utility use remain in the balance.

Gas storage uses molar Cv = Cp − R and advected enthalpy CpT. Adsorbed energy
uses a stated component heat at the adsorption reference temperature and its
adsorbed heat capacity. Liquid uses volumetric sensible storage rho Cp T;
liquid pressure work and viscous heating are excluded. Pump work is accounted
as a separate utility. These choices must be compared with a reference paper's
equations; conserved internal ledgers alone do not prove equation equivalence.

## Cycle scheduling, valves and tanks

`Columns` maps one to four names to main-style overrides. All columns share the
component order and connected-fluid properties. A phase may supply common ports
for one column, or explicit actions for each named column:

```json
"SwingAdsorptionPhases": [{
  "Name": "Adsorption", "Duration": 300.0,
  "Columns": {
    "A": {
      "LeftPort": {"Mode": "Velocity", "Source": "feed", "MolarFlowRate": 4.055373437876738},
      "RightPort": {"Mode": "Pressure", "Source": "light_receiver"},
      "StreamTags": {"LeftIn": "feed_adsorption", "RightOut": "light_product"}
    },
    "B": {"LeftPort": {"Mode": "Closed"}, "RightPort": {"Mode": "Closed"}}
  }
}]
```

`Sources` hold reservoir pressure, temperature and component composition or
concentration. `Connections` within a phase link `A.left`, `A.right`,
`buffer:light` or `source:feed`. They carry a name, `From`, `To`, `FlowModel`,
conductance, optional `Opening`, `CheckValve`, `MaximumMolarFlow` and stream tag.
Gas flow models include linear, square-root and compressible orifice laws.
Compressor connections can specify discharge pressure and efficiency.

`Buffers` describe finite gas tanks with volume, initial pressure, temperature,
composition and thermal mode. Gas `PressureControl` may export above its
setpoint, using `ControllerConductance`, `Deadband` and an `ExportTag`.
Liquid tanks use initial concentration, initial `Volume`, fixed pressure and an
optional `MaximumVolume`; the volume and solute inventory evolve. Empty or
overfilled liquid tanks fail with an explanation. Liquid tanks do not use an
ideal-gas pressure controller.

An unmentioned column holds with closed ports. Independent action stop rules
can wait until all enabled column actions complete. Finished actions close their
ports and retain thermal control; an explicit `WaitAction` can override those
waiting controls. Connected phases require a shared phase stop or duration, so
a pressure-equalization connection cannot continue after an independent column
closure. Product-composition recycle and finite tanks are different choices:
remembered composition does not create a finite product inventory.

## Stops, CSS and accounting

Timed stages use `Duration`. Conditional stages use `StopCriterion`, optional
`MinimumDuration`, `MaximumDuration` and `TimeoutPolicy` (`Error` or `Advance`).
A rule uses `Metric`, `Column`, optional `Component`/`Port`, `Operator`, and
`Value`. `All` and `Any` combine rules. Metrics include outlet composition or
concentration ratio, pressure, temperature extrema, inventories and buffer
states. See the supplied converted phases for pressure and thermal examples.

`SimulationTime: "auto"` in breakthrough requires outlet arrival **and** loading
equilibrium throughout the selected fed components, including HSDM shells, before
the optional observation interval. `AutomaticStopping` controls `Components`,
`Threshold`, `LoadingRelativeTolerance`, `MinimumTime`, `MaximumTime`, and
`ExtraFraction`. This prevents a carrier or fast component from ending a run
before a slower component saturates. A maximum time without a satisfied criterion
does not establish saturation.

`CyclicSteadyState` controls state, inventory and metric relative tolerances,
state absolute tolerance, consecutive cycles, maximum cycles and diagnostic
period/window. Finite tanks and particle/coolant states participate in the state
comparison. Fixed-cycle completion is distinct from CSS convergence. For multiple
columns compare a full schedule including role exchange, not one half-cycle.

`Performance.FreshFeedTags` includes every external fresh-feed stream, including
repressurization and cooling. Product withdrawal for purge is subtracted from net
export. Internal column/tank transfers are not fresh feed or external product.
Productivity uses total active adsorbent mass and full schedule duration,
including waiting. Liquid purity is the fraction of **tracked solutes**, excluding
implicit solvent; its exported concentrations and solution volumes are recorded
separately. Utility balances include compressor, vacuum, pump, heating, cooling,
buffer and waiting duties. These metrics are simulation outputs, not a completed
plant optimization or calibrated TEA.

## Solvers, outputs, Python and Lab

`BreakthroughIntegrator` applies equally to breakthrough and cycles, and
independently of the axial discretization. Select `BDF1` (backward Euler), `BDF`
(native adaptive BDF/NDF), `CVODE`, `SSPRK3` or `SemiImplicitSSPRK3`. The latter
uses the third-order SSP3(4,3,3) additive scheme with explicit transport and
implicit remaining terms; it is distinct from the original compatibility SIRK3.
Runge-Kutta methods have adaptive error control and accept a `TimeStep` cap.
CVODE requires a SUNDIALS-enabled build and offers Dense and SPGMR;
the shared SPGMR backend uses a right band preconditioner and an internal node-major
state order; public files and checkpoints keep the documented state order. Small examples select Dense.
`Solver` contains `RelativeTolerance`, `AbsoluteTolerance`, `MaximumTimeStep`,
`FirstTimeStep`, `MaximumSteps`, `JacobianMode`, `CacheRHS`, `NonnegativeConstraints` and
`EventTimeTolerance`. `IsothermPowerRegularizer` controls the powered-isotherm
small-driving-force continuation (Pa for gas pressure inputs, mol/m³ for liquid).
Zero retains the exact unregularized pure model and enables analytic Sips spreading
pressure/inversion; T15 explicitly uses zero to match its main Sips law. This
is distinct from the positive `KineticConcentrationRegularizer` used in pore-LDF
ratios. The native structured finite-difference pattern is used for
disconnected columns; general connections retain a conservative dense treatment.
Accepted states are checked; invalid states are not silently committed as
successful calculations. CVODE internally solves z = y + 10 atol and constrains z to be nonnegative.
This permits physical roundoff within ten absolute tolerances during stage restarts
without clipping the public state or changing its differential equation. The
accepted-state checks apply independently. `Solver.NonnegativeConstraints: false`
retains those checks while disabling the CVODE inequality constraint. The T15
comparison retains constraints and uses `CVODEMaximumOrder: 2` to handle its
sharp trace-component startup. When its dense-output polynomial overshoots below zero,
CVODE samples use a convex interpolation between accepted states. This does not
clip or replace the accepted trajectory or integrated ledgers. Film closures
use the nonnegative trial-edge continuation imported from the donor, while radial
transport uses the original trial states. This is not a universal positivity proof
for either integrator.

`WriteInterval` controls saved samples, not integration steps. An explicit
`Output.SamplesPerPhase` can cap temporal resolution. Compact output retains a
rolling `HistoryCycles` profile window while CSV cycle/phase/stream totals cover
the full run. Reduce output knowingly; plots cannot recover discarded profiles.

Results include `summary.json`, `metadata.json`, compressed `states.npz`,
`cycle_metrics.csv`, `phase_metrics.csv`, `streams.csv`, component balances,
convergence diagnostics and a checkpoint if enabled. `Output.MainFormatFiles`
also writes main-compatible `column.data` and component files for existing
analysis workflows. For multiple columns each has its own directory. The NPZ
archive includes particle histories and coolant profiles beyond those text files.
Automatic Matplotlib plots are PNG/PDF; a built-in SVG overview is a fallback.
`--plot-only` uses the saved run's input metadata.

```python
from ruptura import ProcessSimulation, load_simulation
simulation = ProcessSimulation("simulation.json")
simulation.validate()
summary = simulation.run(plots=False)
saved = simulation.results()
# load_simulation("simulation.json") dispatches explicit FDM/FVM inputs here.
```

The native solver does the integration; Python supplies configuration and results.
The Lab includes shared-model presets and a Cycle designer with movable column,
reservoir and buffer symbols, phase selection, valve switching/highlighting,
equipment forms, thermal controls, stop conditions and CSS settings. Its forms,
Input JSON tab, native CLI and cluster exporter use the same simulation object.
Drawing positions are separate workspace metadata. The cluster bundle contains
the canonical `simulation.json`, `run.sh` and `submit_hal.sh`; it does not submit
jobs. Additional JSON editors cover particle/coolant/solver settings and recorded
NPZ arrays supply column and particle plots.
The particle view shows HSDM shell loading. Existing main controls remain usable
for original inputs. A full visual browser regression was not part of this build;
native validation, DOM component checks and automated Lab/Python tests are
recorded under `verification/unified/`.

## HAL campaigns

```bash
tools/submit_integrated_hal.sh cycles
tools/submit_integrated_hal.sh breakthrough
# Optional fixed-cycle diagnostic:
tools/submit_integrated_hal.sh cycles --cycles 3 --no-plots
```

The cycle campaign submits **162 cases in one array without a concurrency cap**;
the breakthrough campaign submits 27. Each task uses one CPU and its own case
directory. Per-case `run` and `submit_sim_hal.sh` launchers are also included.
Set `RUPTURA_EXECUTABLE` when the executable is outside `build/ruptura`.
The archived campaign inputs were converted and checked, not all rerun to CSS.
