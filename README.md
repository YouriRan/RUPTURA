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
(Qt application, OFF), `BUILD_DOCS` (Doxygen target, ON), `BUILD_TESTING`
(C++ unit tests, ON), `BUILD_BENCHMARKS` (OFF).

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
```
cd examples/Silicalite-CO2-N2/breakthrough/Langmuir<br>
./run
```

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

The browser opens `http://127.0.0.1:8765`. Use `--no-browser` to print the address,
`--port 0` to choose an available port, or `--workspace /path/to/project` to select
a workspace. The service accepts local connections only and serves its Plotly
JavaScript locally; viewing results does not require internet access or Jupyter.

A native launcher can also be built without Qt:

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
The optional legacy Qt application is available with `-DBUILD_QT_APP=ON` as the
`ruptura_lab_qt` target.

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
binary and the Qt application leave it unchanged.

Input
=====
See the cited article.

Swing adsorption inputs use `SimulationType: "SwingAdsorption"` and define
ordered phases. Each phase runs after the previous one and can override
temperature, inlet pressure, or both:

```json
"SwingAdsorptionPhases": [
  {"Name": "Adsorption", "Temperature": 552.0, "InletPressure": 100000.0, "NumberOfTimeSteps": 1000000},
  {"Name": "Regeneration", "Temperature": 673.0, "InletPressure": 100000.0, "NumberOfTimeSteps": 2000000}
]
```

Swing adsorption inputs must use `SwingAdsorptionPhases`; array-based phase
definitions are not supported in this release. See
`examples/BEA-alkanes-C7/psa` and `examples/BEA-alkanes-C7/tsa`.

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
