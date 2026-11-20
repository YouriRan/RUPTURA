Breakthrough, IAST, and isotherm fitting code
=============================================

This software is a simulation package to compute breakthrough curves and 
IAST mixture predictions. It has been developed at the Delft University of 
Technology (Delft, The Netherlands), during 2022 in active collaboration
with the University of Amsterdam (Amsterdam, The Netherlands), Eindhoven 
University of Technology (Eindhoven, Netherlands), Pablo de Olavide 
University (Seville, Spain), and Shell Global Solutions International B.V.
Amsterdam.

Features
========
* Unlimited amount of components
* Fast (sub-second) IAST mixture computation
* Stable breakthrough computation, including
  - step breakthrough
  - Linearized driving model (LDF)
  - Axial disperion
  - Pressure gradient
* Automatic picture/movie generation
* Isotherm models
  - Langmuir
  - Anti-Langmuir
  - BET
  - Henry
  - Freundlich
  - Sips
  - Langmuir-Freundlich
  - Redlich-Peterson
  - Toth
  - Unilan
  - O’Brien & Myers
  - Quadratic
  - Temkin
* Fitting raw data to isotherm models

Terms of use
============
If you use this software for scientific publications, please cite:<br>
"RUPTURA: Simulation Code for Breakthrough, Ideal Adsorption Solution
Theory Computations, and Fitting of Isotherm Models"<br>
S. Sharma, S. Balestra, R. Baur, U. Agarwal, E. Zuidema, M. Rigutto,
S. Calero, T.J.H. Vlugt, and D. Dubbeldam, 
Molecular Simulation Journal, 49(9), 2023
https://www.tandfonline.com/doi/full/10.1080/08927022.2023.2202757

Authors
=======
Shrinjay Sharma,        Delft University of Technology, The Netherlands<br>
Youri Ran,              University of Amsterdam, The Netherlands<br>
Salvador R.G. Balestra, Pablo de Olavide University, Spain<br>
Richard Baur,           Shell Global Solutions International B.V., The Netherlands<br>
Umang Agarwal,          Shell Global Solutions International B.V., The Netherlands<br>
Eric Zuidema,           Shell Global Solutions International B.V., The Netherlands<br>
Marcello Rigutto,       Shell Global Solutions International B.V., The Netherlands<br>
Sofia Calero,           Eindhoven University of Technology, The Netherlands<br>
Thijs J.H. Vlugt,       Delft University of Technology, The Netherlands<br>
David Dubbeldam,        University of Amsterdam, The Netherlands<br>

Compilation
===========
Ruptura needs a C++ compiler with `<print>` and `<mdspan>`: GCC 16 or newer,
Clang with libc++ 17 or newer, or Apple Clang 17 (Xcode 16.3 or newer).

The simplest route is the conda environment in `env.yml`, which also provides
the compiler. These commands build the `ruptura` executable and the `ruptura`
Python package and install both into the environment:

```
conda env create -f env.yml
conda activate ruptura
cmake --preset conda
cmake --build --preset conda --target install
```

On macOS the conda compiler uses the SDK of the Xcode Command Line Tools
(`xcode-select --install`).

Every supported setup is a CMake preset; `cmake --list-presets` shows the ones
that apply to your machine, and

```
cmake --workflow --preset <preset>
```

configures, builds and runs the unit tests in `build/<preset>`.

| Preset | Setup |
| ------ | ----- |
| `conda` | the `env.yml` environment on Linux or macOS; installs into `$CONDA_PREFIX` |
| `macos` | Apple Clang with `brew install cmake ninja suitesparse` |
| `linux-gcc` | any Linux distribution with GCC 16 or newer |
| `linux-clang` | any Linux distribution with Clang and libc++ 17 or newer |
| `linux-ubuntu-24`, `linux-ubuntu-26`, `linux-debian-13`, `linux-fedora-44`, `linux-redhat-9`, `linux-redhat-10`, `linux-opensuse-leap-15.6`, `linux-opensuse-tumbleweed`, `linux-archlinux` | the packages to install are listed in `docker/<distribution>/install.sh` |

The presets other than `conda` download SUNDIALS and GoogleTest, and build the
Python module for the active Python environment, which needs nanobind
(`python3 -m venv .venv && . .venv/bin/activate && pip install nanobind numpy`).

Without presets:

```
cmake . -B build
cmake --build build
```

Relevant options: `BUILD_CLI` (command-line executable, ON), `BUILD_SUNDIALS`
(CVODE integrator, ON), `BUILD_PYTHON` (nanobind extension, OFF), `BUILD_APP`
(Ruptura Lab launcher, OFF), `BUILD_DOCS` (Doxygen target, ON), `BUILD_TESTING`
(C++ unit tests, ON), `RUPTURA_NATIVE_ARCH` (`-march=native`, ON; OFF for `pip install`
and the conda packages). A build with it runs only on CPUs like the one it was built
on ("Illegal instruction" elsewhere, e.g. on older cluster nodes): pass
`-DRUPTURA_NATIVE_ARCH=OFF` for a binary that is used on other machines.

to clean:<br>
```
rm -rf build
```

to build and host the documentation
```
cmake --build build -- documentation
cd build/html
python -m http.server 8000
```

Testing the build
-----------------
The workflows in `.github/workflows` build and test the presets:

* `pull-request-checks.yml` (pull requests): Ubuntu 24.04, macOS, and the conda
  environment.
* `test-matrix.yml` (started by hand): every Linux image in `docker/` on x86_64
  and aarch64, macOS on Apple silicon and Intel, the conda environment on four
  platforms, and `pip install .`.
* `conda-packages.yml` (started by hand, and on every GitHub release): builds
  the conda packages `ruptura` (executable) and `rupturalib` (Python module)
  from the recipes in `conda/`, and uploads them to anaconda.org once the
  repository has an `ANACONDA_API_TOKEN` secret and an `ANACONDA_OWNER`
  variable.

The Linux jobs also run locally in Docker with act; see `docker/README.md`. The
conda packages build locally with
[rattler-build](https://rattler.build):

```
rattler-build build -r conda/ruptura -m conda/variants.yaml -c conda-forge --output-dir build/conda-packages
rattler-build build -r conda/rupturalib -m conda/variants.yaml -c conda-forge --output-dir build/conda-packages
```

Running
=======
`ruptura` runs the `simulation.json` of the current directory and writes its output there:

```
cd examples/BEA-alkanes-C7/breakthrough
ruptura
```

`ruptura --validate` checks the input without running it.

Ruptura Lab (local web application)
===================================
Ruptura Lab is a browser workspace for configuring simulations, running the local
C++ engine, and exploring its results with the same Plotly functions used by the
analysis notebooks. It supports gas and liquid breakthrough, swing adsorption,
analytical mixture prediction, chemisorption, reactions, and multibed columns.
MPD and fitting remain CLI-only.

The workspace opens in dark mode with the Ruptura molecular logo. Built-in paper
examples include Marx T15 (CO₂/N₂/H₂ with the energy balance), seven C₇ alkane
isomers on BEA, and one 400 ppm monolith at 2 mmol/s. They are
added once to new and existing workspaces, preserving your edits and selection;
removed examples stay removed. More copies are available under New simulation →
From the paper. Inputs ship with the app and use CVODE; see the
[preset notes](ruptura/webapp/presets/README.md) for source inputs and output settings.

Examples use numbered names (`Example 1: …`), while new simulations and copies
receive sequential names such as `untitled_1`. The generic Quick start choices
keep their original names. Quick start includes a JSON paste
box that validates an input before creating a simulation. The question-mark
buttons beside the field reset controls explain each setting on hover, keyboard
focus, or tap.

Import accepts files or pasted workspace JSON and adds their simulations to the
current workspace. Each experiment in the workspace has a three-dot menu with
Delete, Duplicate, and Copy simulation JSON actions. Duplicate creates a new
simulation with an independent copy of all settings and notes. Copy simulation
JSON puts the complete input on the clipboard.

The left sidebar lists simulations. Selecting one opens its settings in the main
panel, with tabs for setup, components, column, reactions, notes, and input JSON.
Each simulation contains its own component properties, isotherms, adsorption heats,
column geometry, adsorbent beds, and reaction network. Editing one simulation does
not change another. New components, beds, and reactions are configured directly
inside the selected simulation.

Previously saved modular workspaces are converted automatically: shared values
and local overrides are resolved into complete, independent simulation inputs.
Simulation IDs, notes, selection, and run associations are retained. Workspace
exports save these independent simulations; simulation exports and run snapshots
remain standalone CLI JSON.

About opens on launch and when clicking the logo. It includes the original
publication, author credits, contact information, a full MIT license dialog, and
a downloadable BibTeX citation. Example references live in the experiment's Notes
tab. Settings offers dark and light themes, remembered in the browser and applied
to interactive plots.

Notes use a large, resizable editor. Carrier gas selects the initial column gas;
it may also have physisorption and chemisorption sites. When the saved run enables
the energy balance, Results opens with breakthrough and temperature curves stacked
vertically at the selected column position.

In Column profiles, choose Position or Time for the x-axis. The other coordinate
becomes a slider above the plot. Set Frames per second (1–60) for playback and use
Download MP4 to export the animation with the selected metric, components, units,
and theme. Exports keep fixed axis limits across frames and use the viewer's saved
frames, including its sampling for large runs. MP4 export uses Matplotlib and
FFmpeg; `ruptura[web]` includes Matplotlib and an FFmpeg package, and an existing
`ffmpeg` on PATH is also supported.

Install the plotting dependencies with the Python package:

```sh
pip install ".[web]"
ruptura-lab
```

For development from this checkout, use an environment containing NumPy, pandas,
and Plotly, and point it at a built CLI:

```sh
python -m ruptura.webapp --cli ./build/ruptura
```

The browser opens `http://127.0.0.1:8765/?session=<key>`; the key is new at every start and
keeps other local users and programs out, so open the Lab through the printed address. Use
`--no-browser` to print the address,
`--port 0` to choose an available port, or `--workspace /path/to/project` to select
a workspace. The service accepts local connections only and serves its Plotly
JavaScript locally; viewing results does not require internet access or Jupyter.

A native launcher can also be built:

```sh
cmake -S . -B build -DBUILD_APP=ON
cmake --build build --target ruptura_lab
./build/app/ruptura_lab
```

The launcher uses `python3` on `PATH`, or the interpreter configured with
`-DRUPTURA_LAB_PYTHON_EXECUTABLE=/path/to/environment/bin/python` (the Python module's
interpreter is reused when `BUILD_PYTHON=ON`). Set `RUPTURA_PYTHON_EXECUTABLE` to
override it at launch. Set `RUPTURA_CLI_EXECUTABLE`, pass
`--cli`, or choose the engine in Settings to select a different `ruptura` binary.
The CLI must include the `--capabilities` and `--validate` commands provided by
this version; rebuild an older executable before using it with the web app.

The editor preserves the complete simulation JSON, including settings that are
edited through the Input JSON tab. Import accepts simulation JSON and saved Qt
or web Lab workspaces. The server-side file browser resolves relative resources
against the imported file's directory; browser uploads resolve them against the
workspace. Workspace edits are saved to `.ruptura-lab/state.json` and can also be
downloaded using Save workspace. Saved Qt libraries are retained in the exported
workspace, and each migrated simulation has its own editable input.

Each queued run receives an immutable `simulation.json` snapshot under
`simulations/`, its own log, and an editable `analysis.ipynb`. Results include
breakthrough curves, pure and mixture isotherms, chemisorption and total loadings,
temperature histories, and a column explorer with metric, component, position,
and time controls. Large column files are streamed and frame indices are sampled
to bound memory. Interactive HTML reports, notebooks, inputs, and raw output files
can be downloaded from Results.

Submitting a simulation opens Run log and keeps it visible when the run finishes.
Open Results when ready to explore the output; each component checkbox controls
its curve directly, including a component marked as carrier gas.

Closing the browser tab leaves the local server and runs active. Press Ctrl+C in
the server terminal to stop it and cancel queued or active jobs. On restart,
saved outputs remain available and previously active jobs are marked interrupted.

Web integration tests:

```sh
python -m pytest tests/webapp_test.py tests/plot_breakthrough_test.py tests/plot_mixture_test.py
```

For an optional browser check, install Playwright and its Chromium browser, then
run `python tests/webapp_browser_smoke.py`. It launches a disposable workspace,
checks editing, simulation runs, plots, downloads, autosave, mobile layout, and
cancellation, then stops its server. Use `--launcher /path/to/ruptura_lab` to
check an installed launcher.

Running the `ruptura` command-line binary also creates the matching
`analysis.ipynb` beside `simulation.json`. If the notebook already exists, the
binary leaves it unchanged.

Input
=====
`docs/GENERATIVE_AI_SIMULATION_PROMPT.md` describes the keys of `simulation.json` with their units, and the
directories under `examples/` hold complete inputs. Keys are case-insensitive.

Simulated time and output are set in seconds: `SimulationTime` (a time, or
`"auto"` to stop 10% after the breakthrough or reaction steady state is
reached) and `WriteInterval` (one column profile and progress line every so
many seconds). `"WriteOutletOnly": true` writes only the outlet row of every
output (z = L, or the first cell in a reverse-flow phase) to `column.data` and the
component files, which keeps long runs small; the breakthrough and temperature
plots read these files as usual, while profiles over position need the full
output. `TimeStep` is the fixed step of `RungeKutta3`. It is not checked
for stability: keep it below the time in which the fluid crosses the smallest
grid cell (dz / v), below about 2 / k for the largest `MassTransferCoefficient`
k, and below dz² / (2 D) for the largest axial dispersion coefficient D and,
with the energy balance, for the thermal diffusivity of the gas, λ / (ρ c_p)
(`gasThermalConductivity` over the gas density times `heatCapacityGas`); a
larger step gives wrong curves without an error, or a temperature that
overflows and stops the run with an isotherm error. `CVODE` and
`ROS34PW2` (an adaptive Rosenbrock-W method with the analytical sorption
Jacobian, controlled by `ROS34PW2RelativeTolerance` and the
`ROS34PW2AbsoluteTolerance...` keys) choose their own steps and reject
`TimeStep`. The step-based keys of
RUPTURA 2 (`NumberOfTimeSteps`, `WriteEvery`, `PrintEvery`) are rejected with
the equivalent times; `python convert.py simulation.json` converts a file in
place and keeps the original as `_legacy_simulation.json`.

Every output rewrites `column.json`, the state of the column and its time; a
failed run leaves `failed_state.json`. `"ReadColumnFile": "column.json"` starts a
run from such a state and continues it: the clock goes on from the time of the
state, for `SimulationTime` (or the phase durations) more, and the outputs are
appended to the files of the run directory. With `"ResetClock": true` the state
starts a new experiment instead, at t = 0 and with new output files.

Swing adsorption inputs use `SimulationType: "SwingAdsorption"` and define
ordered phases. Each phase runs for its `Duration` in seconds after the previous one and can
override temperature, inlet pressure, or both:

```json
"SwingAdsorptionPhases": [
  {"Name": "Adsorption", "Temperature": 552.0, "InletPressure": 100000.0, "Duration": 500.0},
  {"Name": "Regeneration", "Temperature": 673.0, "InletPressure": 100000.0, "Duration": 1000.0}
]
```

Swing adsorption inputs must use `SwingAdsorptionPhases`; array-based phase
definitions are not supported in this release. See
`examples/BEA-alkanes-C7/psa` and `examples/BEA-alkanes-C7/tsa`.

A phase can close one end of the column and ramp the pressure, which gives the
pressurization and blowdown steps of a PSA cycle (`examples/PSA-H2-activated-carbon`):

```json
"SwingAdsorptionPhases": [
  {"Name": "Pressurization", "ClosedEnd": "ProductEnd", "InletPressure": 1.0e6, "Duration": 30.0,
   "PressureRamp": {"Shape": "Exponential", "TimeConstant": 5.0}},
  {"Name": "Adsorption", "InletPressure": 1.0e6, "ColumnEntranceVelocity": 0.03, "Duration": 60.0},
  {"Name": "Blowdown", "ClosedEnd": "ProductEnd", "InletPressure": 1.0e5, "Duration": 30.0,
   "PressureRamp": {"Shape": "Exponential", "TimeConstant": 5.0}},
  {"Name": "Purge", "InletPressure": 1.0e5, "ColumnEntranceVelocity": 0.1, "FlowDirection": "Reverse",
   "GasPhaseMolFractions": {"H2": 1.0}, "Duration": 60.0}
]
```

- `ClosedEnd` is `FeedEnd` (z = 0) or `ProductEnd` (z = L). The velocity is zero there, and
  the other end is open at the phase `InletPressure`. Gas leaves through the open end, or
  enters there with the phase `GasPhaseMolFractions`, as the pressure and the uptake
  require. It takes the place of `FlowDirection` and `ColumnEntranceVelocity`, and needs
  `BoundaryCondition` `InletPressureInletVelocity` and the Danckwerts inlet.
- `PressureRamp` takes the inlet pressure (at the open end of a closed-end phase) from the
  pressure at that end at the end of the previous phase, P₀, to the phase `InletPressure`:
  `Exponential` follows P(t) = P + (P₀ − P) exp(−t/`TimeConstant`), and `Linear` reaches P
  after `TimeConstant` and stays there. It needs `InletPressureInletVelocity`. The whole
  column follows the ramp, with the accumulation c dP/dt / P in the continuity closure and,
  with the energy balance, the pressure work (dP/dt)/(ρ c_p) in the gas. The time derivative
  of the pressure drop is left out, so that the pressure drop must stay small against the
  lowest pressure.
- The gas takes the pressure profile of a new phase at its composition,
  c_i ← c_i P/(R T Σc), after a pressure step without a ramp, a temperature step without
  the energy balance, and when the end where the pressure is set changes (a reverse phase
  after a forward one, or the other way), which moves the profile by the pressure drop. This
  holds for `InletPressureInletVelocity` and `InletPressureOutletPressure`; `FixedVelocity`
  and `FixedPressureInletVelocity` let the feed displace the gas, and
  `InletVelocityOutletPressure` does not follow the phase pressure. A large pressure step can
  stop the run with a negative outlet pressure, when what desorbs at once cannot leave
  through the bed at the new pressure; a `PressureRamp` or a `ClosedEnd` then avoids it.
- `phases.data` holds, per phase, the start and end time, the pressures at both ends, the
  largest |R T Σc/P − 1| at the end, and the gas and adsorbed amount of each component at
  the start and end, per unit column cross-section in mol/m². What entered minus what left
  in a phase is the change of these amounts, which gives product purity and recovery
  without integrating the outlet. In Python, `SwingAdsorption.phase_records` holds the same
  after `run()`.
- ROS34PW2 and CVODE both hold Σc at P/(R T) over many cycles; ROS34PW2 is usually several
  times faster. Its stage matrix leaves out the continuity closure, so that with
  `InletPressureInletVelocity` and `InletVelocityOutletPressure` it is projected on the
  closure: the stages change Σc of a cell as the right-hand side does. In the T15
  temperature-swing cycles of the paper (30 cycles of 2 × 1000 s), max |R T Σc/P − 1| stays
  below 1e-4.

A phase can end on a condition, take the gas that left an earlier phase, and the run can end at
the cyclic steady state:

```json
"NumberOfCycles": 40,
"CyclicSteadyState": {"Tolerance": 1.0e-3, "ConsecutiveCycles": 1},
"SwingAdsorptionPhases": [
  {"Name": "Adsorption", "Duration": 300.0, "MinimumDuration": 10.0,
   "StopCriterion": {"Metric": "MoleFraction", "Component": "CO2", "End": "ProductEnd",
                     "Operator": ">=", "Value": 0.05}},
  {"Name": "Purge", "FlowDirection": "Reverse", "GasPhaseMolFractions": "Adsorption", "Duration": 60.0}
]
```

- `StopCriterion` ends a phase at the time at which its `Metric` reaches the `Value`
  (`Operator` `>=` or `<=`); the `Duration` is then the longest the phase lasts. The metrics
  are `Pressure` (gas pressure of the cell at the `End`, Pa), `AveragePressure`,
  `MaximumTemperature` and `MinimumTemperature` (solid temperature, K), `MoleFraction` (of the
  `Component`, in the cell at the `End`), and `AmountIn` and `AmountOut` (net amount that
  entered or left at the `End` in this phase, of the `Component` or in total, mol/m²). `End` is
  `FeedEnd` (z = 0) or `ProductEnd` (z = L, the default). `{"All": [...]}` and
  `{"Any": [...]}` combine criteria. The criteria apply from `MinimumDuration` on.
  `TimeoutPolicy` `Error` stops the run when the `Duration` passes first; `Advance` (default)
  goes on to the next phase. CVODE finds the time with its root finding; RungeKutta3 and
  ROS34PW2 repeat the step in which the criteria turn true and bisect its end. The outputs
  stay at multiples of `WriteInterval`.
- `GasPhaseMolFractions` and `ProductEndGasPhaseMolFractions` take the name of a phase in
  place of mole fractions: the mean gas that left the column in that phase the last time it
  ran, e.g. a purge with the product of the adsorption step. Until that phase has run, the
  feed is that of `Components`, and the product-end gas that of the last cell.
- `CyclicSteadyState` ends the run before `NumberOfCycles`, which is then the largest number
  of cycles, once `ConsecutiveCycles` cycles in a row end within `Tolerance` of the cycle
  before, in three relative measures: every state value, the amount of every component in
  the column, and every stream amount of every phase. Each change is taken relative to the
  size of the value plus the largest of its kind, so that values near zero count against that
  largest one. `cycles.data` lists the three per cycle, and next to them the amount that the
  streams and the column of the cycle do not account for (the mass balance).
- With any of these keys, and with `MomentumBalance` below, the column integrates the amount
  of each component that passes each end, and `phases.data` and
  `SwingAdsorption.phase_records[i].stream_amounts` hold, per phase, what entered and left at
  z = 0 and at z = L (four values per component, mol/m²). They are net per end: a component
  that entered and left at one end in a phase shows the difference.

`BoundaryCondition` `MomentumBalance` drops the continuity closure. The pressure of every cell
is that of its gas, R T Σc, and the velocity between two cells follows from the pressure-drop
law of the bed for their pressure difference. The column can then be filled and emptied through
either end, as in the steps of a PSA cycle (`examples/PSA-H2-activated-carbon-ports`):

```json
"BoundaryCondition": "MomentumBalance",
"BreakthroughIntegrator": "ROS34PW2",
"OutletPressure": 1.0e5,
"SwingAdsorptionPhases": [
  {"Name": "Pressurization", "Duration": 40.0,
   "FeedEnd": {"Mode": "Pressure", "Pressure": 1.0e6,
               "PressureRamp": {"Shape": "Linear", "TimeConstant": 40.0}},
   "ProductEnd": {"Mode": "Closed"}},
  {"Name": "Adsorption", "Duration": 60.0,
   "FeedEnd": {"Mode": "Velocity", "Velocity": 0.03, "Pressure": 1.0e6},
   "ProductEnd": {"Mode": "Pressure", "Pressure": 1.0e6}},
  {"Name": "Blowdown", "Duration": 40.0,
   "FeedEnd": {"Mode": "Pressure", "Pressure": 1.0e5,
               "PressureRamp": {"Shape": "Linear", "TimeConstant": 40.0}},
   "ProductEnd": {"Mode": "Closed"}},
  {"Name": "Purge", "Duration": 60.0, "ProductEndGasPhaseMolFractions": "Adsorption",
   "FeedEnd": {"Mode": "Pressure", "Pressure": 1.0e5},
   "ProductEnd": {"Mode": "Velocity", "Velocity": 0.1, "Pressure": 1.0e5}}
]
```

- The column starts at rest at `OutletPressure`, filled with the carrier gas, and keeps its
  gas from one phase to the next. It needs a gas and the Danckwerts inlet.
- Each phase sets the two ends with the ports `FeedEnd` (z = 0) and `ProductEnd` (z = L);
  an end without a port is closed. `Mode` is `Closed`, `Pressure`, or `Velocity`:
  - `Pressure`: gas at `Pressure` beyond the end. The velocity follows from the pressure-drop
    law over the half cell between that gas and the centre of the end cell, so that gas
    flows in or out as the pressures require. That half cell is the only resistance: a
    `PressureRamp`, from the pressure of the end cell at the start of the phase, stands for
    the valve or the pump that sets how fast the column follows.
  - `Velocity`: the interstitial `Velocity` into the column (negative: out of it) of gas at
    the `Pressure` given with it, or a `MolarFlux` into the column in mol/(m² s) of column
    cross-section.
  Gas that enters at z = 0 is the phase `GasPhaseMolFractions` (default: the `Components`
  feed), gas that enters at z = L the `ProductEndGasPhaseMolFractions` (default: the gas of
  the last cell), both at the phase `Temperature`. The phase keys `InletPressure`,
  `ColumnEntranceVelocity`, `FlowDirection`, `ClosedEnd`, and `PressureRamp` are not used.
- Without swing phases the feed enters at `ColumnEntranceVelocity`, as gas at `OutletPressure`,
  and leaves to `OutletPressure` at z = L.
- With the energy balance the pressure work of the gas follows its own pressure, which leaves
  its heat capacity at constant volume. That needs a molar heat capacity above R: give every
  component its `GasHeatCapacity` in J/(mol K), which then takes the place of
  `heatCapacityGas` (per kg, one value for any composition; too small for helium- or
  hydrogen-rich gas). `GasHeatCapacity` can be used with every boundary condition.
- The relaxation of the pressure is stiff, and its Jacobian is that of the velocity. Use
  ROS34PW2, which holds it, with the gas temperature and the time derivative of a ramped
  port. CVODE with `CVODELinearSolver` `KLU` works without the energy balance; the reader
  rejects CVODE with `Dense` or `SPGMR`, which take many more steps, and RungeKutta3, which
  needs steps of about a microsecond. With CVODE the stream amounts through a pressure port hold its tolerance on a
  flux that follows a small pressure difference, which is less than that of the column
  amounts; `cycles.data` shows it as the mass balance.

Several identical columns run one cycle together when `Columns` names them. The `Connections`
of a phase are valves and compressors between the ends of the columns and `Buffers`, as in a
pressure equalization, a purge with the product of the other bed, or a product tank
(`examples/PSA-H2-activated-carbon-two-beds`):

```json
"Columns": ["A", "B"],
"Buffers": {"Tank": {"Volume": 1.0e-3, "Pressure": 1.0e6, "MaximumPressure": 1.0e6,
                     "GasPhaseMolFractions": {"H2": 1.0}}},
"SwingAdsorptionPhases": [
  {"Name": "Adsorption A, purge B", "Duration": 40.0,
   "Columns": {"A": {"FeedEnd": {"Mode": "Velocity", "Velocity": 0.03, "Pressure": 1.0e6}},
               "B": {"FeedEnd": {"Mode": "Pressure", "Pressure": 1.0e5}}},
   "Connections": [
     {"Name": "Product", "From": "A.ProductEnd", "To": "Tank", "Conductance": 2.0e-8,
      "CheckValve": true},
     {"Name": "Purge", "From": "A.ProductEnd", "To": "B.ProductEnd", "Conductance": 8.0e-10}]},
  {"Name": "Equalization A to B", "Duration": 20.0,
   "Connections": [{"From": "A.ProductEnd", "To": "B.ProductEnd", "Conductance": 5.0e-8}]}
]
```

- `Columns` needs `MomentumBalance`, ROS34PW2, and the `InternalDiameter` of the `Geometry`
  (the connections carry mol/s). All columns have the input of the one column and start
  alike, so that the first cycles are a start-up. Each writes its output files, with its own
  `phases.data`, to a directory of its name.
- The keys of a phase hold for every column. Its `Columns` object gives single columns their
  own `FeedEnd`, `ProductEnd`, `GasPhaseMolFractions`, `ProductEndGasPhaseMolFractions`,
  `Temperature`, or `ThermalControl`. A `StopCriterion` names its `Column` (default: the
  first), and a phase name as gas is the gas that left that same column.
- A connection runs `From` one endpoint `To` another: `<column>.FeedEnd`,
  `<column>.ProductEnd`, or the name of a buffer. An end with a connection takes no port. Gas
  flows from the higher to the lower pressure, `Opening` (0 to 1) times the flow of its
  `FlowModel`, at most `MaximumMolarFlow` mol/s, and with `"CheckValve": true` only from
  `From` to `To`:
  - `Linear` (default): `Conductance` × ΔP, in mol/(s Pa).
  - `SquareRoot`: `Conductance` × √|ΔP|, in mol/(s Pa^0.5).
  - `Orifice`: the isentropic flow of an ideal gas through the effective `Area` in m²
    (discharge coefficient included), choked below the critical pressure ratio.
  - `Compressor`: `MolarFlow` in mol/s from `From` to `To` whatever the pressures. Its work,
    the isentropic compression of that gas over `Efficiency` (default 1), is in
    `connections.data`.
  `Orifice` and `Compressor` need the `GasHeatCapacity` of the components. The gas keeps the
  temperature of its donor (that of a compressor is cooled back to it).
- A buffer is a mixed vessel of `Volume` m³ at a constant `Temperature` (default: the root
  one), which starts at `Pressure` with `GasPhaseMolFractions`. Gas above `MaximumPressure`
  leaves it: the product of a tank held at its pressure.
- `connections.data` holds, per phase and connection, the amount of each component that
  passed in mol and the compressor work in J; `buffers.data`, per phase and buffer, the
  pressure, the amounts, and what left above `MaximumPressure`.
- `python -m ruptura.flowsheet simulation.json > flowsheet.tex` writes a TikZ picture of the
  columns, buffers, valves, and compressors (`\usepackage{tikz}` suffices). With a phase
  name or number after the file name, the connections of that phase carry their direction,
  the others are gray, and its ports are at the column ends.
- Each column keeps its integrator, and the connections exchange gas in coupling steps: the
  flows follow from the pressures at the start of a step, a twentieth of the time in which
  the fastest connection evens out its two pressures, and a column that gives gas takes the
  step before the column that receives it. What leaves one column is then what enters the
  next, for every component. A valve so wide that it evens out the pressures within 1e-4 of
  the phase is cut to that rate. One end cannot take gas from one connection and give it to
  another at once (put a buffer between them), and gas cannot go around a ring of columns.

A phase can hold the wall at a temperature, and the run lists what the cycle costs:

```json
"ExhaustPressure": 1.0e5,
"CompressionEfficiency": 0.75,
"SwingAdsorptionPhases": [
  {"Name": "Heating", "Duration": 600.0,
   "ThermalControl": {"Temperature": 400.0, "HeatTransferCoefficient": 100.0}}
]
```

- `ThermalControl` (energy balance) sets, for its phase, the `Temperature` beyond the outer
  wall and the `HeatTransferCoefficient` of that wall in W/(m² K); without them the phase
  `Temperature` and `heatTransferWallExternal`.
- With the stream amounts, `phases.data` ends with the heat that entered through the outer
  wall in the phase (negative: cooling) and the work that brings the gas that left through
  a pressure port below `ExhaustPressure` to that pressure, both in J/m² of column
  cross-section (`wall_heat` and `exhaust_work` of the phase records in Python). The work
  is the isentropic compression from the port pressure over `CompressionEfficiency`
  (default 1): that of a vacuum pump.
- `"Coolant": {"MassFlowRate": 0.002, "Volume": 2.0e-4}` (root key, energy balance) puts a
  heat-transfer fluid along the outer wall in place of the constant temperature beyond it. It
  enters at that temperature (the `Temperature` of `ThermalControl`) at z = 0, or at z = L
  with `"FlowDirection": "CounterCurrent"`, and each cell of it exchanges heat with the wall
  of its cell through the `HeatTransferCoefficient`, so that the fluid cools or warms on its
  way. `MassFlowRate` is in kg/s (a `ThermalControl` can set it per phase; 0 leaves the fluid
  where it is), `Volume` in m³ is that of the fluid along the column, and `Density` and
  `HeatCapacity` default to water (1000 kg/m³, 4180 J/(kg K)). The wall heat in `phases.data`
  is then the heat that the fluid gave to the wall.

Optional numerics and particle models
-------------------------------------

```json
"ConvectionScheme": "WENO5",
"ParticleTransport": {"Model": "HSDM", "RadialCells": 12}
```

- `ConvectionScheme` is `WENO3` (default), `WENO5`, `TVD`, or `Upwind`. `WENO5` is of fifth
  order on a grid of one cell width and `WENO3` next to the column ends and on a non-uniform
  grid; `TVD` is the second-order MUSCL scheme with the minmod
  limiter. On T15 (40 cells) the CO2 front is 180 s wide (10 to 90 %) with `Upwind`, 166 s
  with `TVD`, 165 s with `WENO3`, and 164 s with `WENO5`; CVODE and ROS34PW2 agree for each.
  The Jacobian of ROS34PW2 and of CVODE with KLU is that of `WENO3` for all three, which
  costs `WENO5` about 60 % more steps with ROS34PW2.
- `ParticleTransport` `HSDM` replaces the linear driving force of physisorption by diffusion
  of the adsorbed phase in a spherical particle (homogeneous surface diffusion) with
  `RadialCells` cells along its radius (2 to 100, default 12). The surface is at equilibrium
  with the gas. The diffusivity is that of the `MassTransferCoefficient` k of each component,
  D_s/R² = k/15 (Glueckauf), so that both models agree at long contact times and the particle
  takes up faster than the linear driving force at short ones; a temperature law of k carries
  over. Not with `Reactions` or CVODE with KLU. ROS34PW2 treats the radial cells explicitly:
  its step stays below about R²/(D_s N²).
- A `PseudoNth` chemisorption site takes `DesorptionRateCoefficient` and `DesorptionOrder`
  in its `RateEquationParameters` for the release (q > q*); both default to those of the
  uptake, and the coefficient can be a temperature law with its own activation energy. Two
  such sites give dual-order kinetics.

Temperature-dependent parameters
--------------------------------

A temperature-dependent parameter is written as a temperature law in place of
its number, wherever the input reader accepts one:

```json
{"Name": "CO2", "HeatOfAdsorption": 24860.0, "ReferenceTemperature": 318.0,
 "PhysisorptionSites": [{"Type": "Sips", "Parameters": [
     {"Value": 12.45, "Law": "Exponential", "Chi": 0.4},
     {"Value": 4.39e-7, "Law": "VantHoff"},
     1.36]}]}
```

`Value` is the parameter at the reference temperature `T_ref`, and the
coefficient key belongs to the law:

| Law | p(T) | Coefficient |
|---|---|---|
| `VantHoff` | V exp[Q/R (1/T − 1/T_ref)] | `Energy` Q [J/mol] |
| `Arrhenius` | V exp[−E/R (1/T − 1/T_ref)] | `ActivationEnergy` E [J/mol] |
| `Exponential` | V exp[χ (1 − T/T_ref)] | `Chi` χ [-] |
| `Linear` | V + a (T − T_ref) | `Slope` a [unit of V per K] |
| `Inverse` | V + α (1 − T_ref/T) | `Alpha` α [unit of V] |

A law without its own `ReferenceTemperature` takes that of its component (or,
for a reaction, of `Kinetics`). Without any, `VantHoff` and `Arrhenius` are
uncentred, V exp[Q/(RT)] and V exp[−E/(RT)]; the other laws need one.

Only some parameters accept a law, and each only the laws of its kind:

- saturation capacities (Langmuir, pH-Langmuir, Sips, Langmuir-Freundlich,
  Toth, GAB): `Exponential`, `Linear`;
- the affinity of Henry, Freundlich, Langmuir, pH-Langmuir, Sips,
  Langmuir-Freundlich, and Toth: `VantHoff`;
- heterogeneity exponents (Freundlich, Sips, Langmuir-Freundlich, Toth,
  Redlich-Peterson): `Linear`, `Inverse`;
- the two Redlich-Peterson affinities and the GAB C and K: `VantHoff` with an
  `Energy`;
- chemisorption rate coefficients (`MassTransferCoefficient`, `RateCoefficient`,
  `AdsorptionRateCoefficient`, `DesorptionRateCoefficient`, `Alpha`) and a
  reaction's `ForwardRateCoefficient`: `Arrhenius`;
- a reaction's `EquilibriumConstant`: `VantHoff`;
- the LH/LHHW `AdsorptionConstants` of a reaction (one per component name):
  `VantHoff` with an `Energy`.

The affinity and the heat of adsorption are linked. A `VantHoff` affinity
without `Energy` takes the component's `HeatOfAdsorption`. Without a
`HeatOfAdsorption`, the heat released in the energy balance is the `Energy`
of the affinity laws, which must then agree. A given `HeatOfAdsorption` is
always the energy-balance heat, and a law with its own `Energy` keeps it. A
chemisorption site uses its own `HeatOfAdsorption`, else the `Energy` of its
affinity law, else the component's. Reactions link `HeatOfReaction` (positive
when exothermic) and a `VantHoff` `EquilibriumConstant` in the same way.

The keys these laws replaced (`nonIsothermal`, the five-parameter GAB, the
`ActivationEnergy` keys of General chemisorption and reactions, and the
reaction `gibbsFreeEnergy`) are rejected with their replacement. `python
convert.py simulation.json` converts a file in place; the original is kept as
`_legacy_simulation.json`, or `_legacy_2_simulation.json` if that exists.

Liquid breakthrough and pH-dependent adsorption
------------------------------------------------

Set `FluidPhase` to `Liquid` for a constant-density liquid breakthrough. The
solvent is implicit; component state values and liquid isotherm driving forces
are solute concentrations in mol/m3. Define the feed and initial bed state for
each component with `LiquidPhaseConcentration` and
`InitialLiquidPhaseConcentration`, and set `LiquidDensity` in kg/m3. Gas-phase
inputs remain unchanged: their mixture-prediction driving forces are mole
fractions and `PartialPressure` continues to be derived from mole fraction and
total pressure.

The pH-dependent Langmuir isotherm is selected with `Type: "pH-Langmuir"` and
parameters `[q_sat, k, pH0]`:

```
q = q_sat * k_pH * c / (1 + k_pH * c)
k_pH = k / (1 + 10^(pH0 - pH))
```

Use `MixturePredictionMethod: "SPI"` with this isotherm. A constant-pH liquid
uses `pHMode: "Fixed"` and `pHValue`. For a transported acid or base front, use
`pHMode: "HPlus"` or `"OHMinus"`, select the concentration component with
`pHComponent`, a zero-based integer component identifier (for example,
`"pHComponent": 0` selects the first entry in `Components`), and optionally
set `pKw` (default 14). The concentration is
converted from mol/m3 when calculating pH.

Complete single-column and multibed inputs are in
`examples/Liquid/pH-dependent` and `examples/Liquid/multibed-pH-front`.

Macrostate particle-distribution mixture prediction
---------------------------------------------------

Set `MixturePredictionMethod` to `MPD` and provide an `MPDSettings` object:

```json
"MPDSettings": {
  "FileName": "particle_distribution.data",
  "ReferenceTemperature": 300.0,
  "ReferenceFugacity": 100000.0,
  "ReferenceFrameworkMass": 1.0e-24,
  "ComponentBounds": [
    {"Component": "A", "NMin": 0, "NMax": 10, "DeltaN": 1},
    {"Component": "B", "NMin": 0, "NMax": 20, "DeltaN": 1}
  ]
}
```

`ComponentBounds` follows the order of non-carrier entries in `Components`.
Bounds are inclusive, and `(NMax - NMin)` must be divisible by `DeltaN`.
The data file contains exactly the product of the component extents in C order,
so the last component varies fastest. Each non-comment row contains either
`Pi_ref` alone or `Pi_ref` followed by the conditional mean Hamiltonian `<H>`
in joules. Full-line comments beginning with `#` and blank lines are allowed.

`ReferenceFugacity` is the common reference component fugacity in Pa and
`ReferenceFrameworkMass` is the represented framework mass in kg. At present,
the pressure grid is used as the target ideal-gas fugacity. Output component
loadings are converted to mol/kg. If the optional energy column is absent, the
first-order Hamiltonian temperature correction is omitted. A missing target
`Temperature` produces a warning when `ReferenceTemperature` is supplied.
See `examples/MPD` for complete CO2/N2 mixture-prediction and helium-carrier
breakthrough inputs based on the same reference distribution. The
`examples/MPD/BEA-alkanes` validation constructs a seven-component reference
distribution from multinomial site families and compares its MPD loadings with
the competitive multisite Langmuir result.

Python Usage
======
Ruptura can be used from Python through the `ruptura` package. The compiled
extension is imported as `ruptura._ruptura`, while the public `ruptura` package
exposes the simulation objects and convenience helpers.

Installation
------------

The quickest route is the bundled conda environment, which carries every
dependency (the Clang compiler, CMake, Ninja, nanobind, scikit-build-core,
BLAS/LAPACK, SuiteSparse/KLU, SUNDIALS, and the notebook stack). Install with
the two `cmake` commands under Compilation, or with pip:

```
micromamba create -f env.yml     # or: conda env create -f env.yml
micromamba activate ruptura
pip install .
```

More generally, from the top level of the source tree, into whatever virtual
environment, conda or micromamba environment is currently active:

```
pip install .
```

That is all that is needed. The build is driven by `scikit-build-core`, which
runs CMake for you, compiles the nanobind extension against the *active*
interpreter, and installs both the `ruptura` Python package and the `ruptura`
command-line executable. `nanobind` and `cmake` are pulled in automatically as
build requirements, so they do not have to be installed by hand.

Requirements on the machine: a compiler as described under Compilation,
LAPACK/BLAS, and -- because the CVODE integrator is enabled by default --
SuiteSparse (KLU). These come from `env.yml` if you used it; otherwise
`brew install suitesparse` on macOS, or
`apt install libsuitesparse-dev liblapack-dev` on Debian/Ubuntu. CMake looks for KLU in the
active conda/micromamba prefix, in Homebrew and in `/usr/local`; point it
somewhere else with `-DKLU_ROOT=<prefix>`.
The first install also clones and builds SUNDIALS v7.5.0, which takes a few
minutes; the build tree is cached under `build/pip-<tag>/`, so later installs
are fast.

To build without CVODE (no SuiteSparse, no SUNDIALS download, roughly a one
minute install):

```
pip install . --config-settings=cmake.define.BUILD_SUNDIALS=OFF
```

Any other CMake variable can be set the same way, for example an already
installed SUNDIALS:

```
pip install . \
  --config-settings=cmake.define.RUPTURA_USE_SYSTEM_SUNDIALS=ON \
  --config-settings=cmake.define.SUNDIALS_ROOT=/opt/sundials
```

For an editable development install (the extension is rebuilt automatically
when the C++ sources change):

```
pip install nanobind scikit-build-core
pip install --no-build-isolation -ve .
```

Alternatively, build the extension in place with plain CMake -- the compiled
module is copied next to `ruptura/__init__.py`, so the package is importable
from the source directory:

```
cmake . -B build -DBUILD_PYTHON=ON
cmake --build build
```

Start a simulation from Python with a JSON input file:

```python
import numpy as np
import ruptura

simulation = ruptura.load_simulation("simulation.json")
result = ruptura.compute(simulation)

data = np.asarray(result)
print(data.shape)
print(result.columns)
```

For a one-line in-memory run:

```python
result = ruptura.run("simulation.json")
```

All bound C++ objects, such as `InputReader`, `Component`, `Isotherm`,
`MixturePrediction`, `Column`, `Breakthrough`, `Fitting`, and
`SwingAdsorption`, are available from the top-level `ruptura` namespace.

`ruptura.selfdriving` (needs scipy: `pip install "ruptura[selfdriving]"`) uses
the interface for a closed-loop optimization in the manner of a self-driving
lab: `PressureSwingCycle` runs a pressure-swing cycle with product backfill to
its cyclic steady state and reads purity, recovery, and productivity from
`SwingAdsorption.phase_records`, and `Campaign` proposes batches by Bayesian
optimization, runs them in parallel processes, and keeps a lab notebook:

```python
from ruptura.selfdriving import Campaign, PressureSwingCycle

cycle = PressureSwingCycle("simulation.json", product="H2", purity_target=0.999)
Campaign(cycle, "campaign", rounds=15, batch_size=8, workers=8).run()
```

`TwoBedPressureSwingCycle` runs the cycle on two beds (`Columns`) with valves
between their product ends, a product tank, and a pressure equalization in place
of the backfill, from a `MomentumBalance` input such as
`examples/PSA-H2-activated-carbon-two-beds`.
`examples/Self-driving-lab-H2-PSA/analysis.ipynb` runs and analyzes both
campaigns for H2 from the T15 feed on activated carbon.

A released version can also be installed from an index, without a source
checkout:

```
pip install ruptura
```
or, from the anaconda.org channel that `conda-packages.yml` uploads to, where
`ruptura` is the executable and `rupturalib` the Python module (`import ruptura`):
```
conda install -c <channel> -c conda-forge ruptura rupturalib
```
