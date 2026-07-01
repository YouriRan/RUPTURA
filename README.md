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

Input
=====
See the cited article.

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
