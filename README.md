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
The command-line executable is built from `main.cpp`. The optional native Qt
application lives in `app/` and provides widgets for creating components,
columns, simulations, running the local `ruptura` executable, and opening an
analysis notebook.

Build it with Qt 6 Widgets available:

```
cmake . -B build -DBUILD_QT_APP=ON
cmake --build build --target ruptura_lab
```

Run:

```
./build/app/ruptura_lab
```

The Analysis button writes `analysis.ipynb` into the experiment run directory
and opens Jupyter Lab with the relevant `ruptura` widgets preloaded. Simulation
run directories are created under `simulations/` in the current working
directory. The Load JSON button can import an existing single-simulation
`ruptura` JSON file into the editable app state. Save State writes a separate
Ruptura Lab state JSON format that preserves all components, columns, and
multiple simulations.

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
