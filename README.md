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
```
cmake . -B build
cmake --build build
```

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

Running
=======
```
cd examples/Silicalite-CO2-N2/breakthrough/Langmuir<br>
./run
```

Qt Application
==============
The command-line executable is always built from `src/main.cpp`. The optional
native Qt application lives in `app/` and provides widgets for creating components,
columns, simulations, running the local `ruptura` executable, and opening an
analysis notebook.

Build it with Qt 6 Widgets available:

```
cmake . -B build -DBUILD_APP=ON
cmake --build build --target ruptura_lab
```

Run:

```
./build/app/ruptura_lab
```

Each run immediately receives an `analysis.ipynb` tailored to mixture prediction,
breakthrough, or fitting. The Analysis button opens that notebook in Jupyter Lab
without replacing any edits. Simulation run directories are created under
`simulations/` in the current working directory. The Load JSON button can import
an existing single-simulation `ruptura` JSON file into the editable app state.
Save State writes a separate Ruptura Lab state JSON format that preserves all
components, columns, and multiple simulations.

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
`pHMode: "HPlus"` or `"OHMinus"`, name the concentration component with
`pHComponent`, and optionally set `pKw` (default 14). The concentration is
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
breakthrough inputs based on the same reference distribution.

Python Usage
======
Ruptura can be used from Python through the `ruptura` package. The compiled
extension is imported as `ruptura._ruptura`, while the public `ruptura` package
exposes the simulation objects and convenience helpers.

Build the Python extension with:

```
cmake . -B build -DBUILD_PYTHON=ON -DRUPTURA_PYTHON_MODULE_NAME=_ruptura
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

The python package can be installed with:

```
pip install ruptura
```
or
```
conda install ruptura
```
