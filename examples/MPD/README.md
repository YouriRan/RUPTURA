# CO2/N2 macrostate particle distribution

This example uses a two-component reference distribution with particle-number
bounds `n_CO2,n_N2 = 0,...,80`. The single-column data file is in C order, so
`n_N2` varies fastest. It is normalized at a reference temperature of 293 K
and a common reference component fugacity of 100,000 Pa.

The represented framework mass is
`80 / N_A = 1.3284312537390774e-22 kg`. Consequently, a mean population of 80
particles corresponds to a loading of 1 mol/kg.

- `simulation.json` predicts the equimolar CO2/N2 mixture isotherm.
- `breakthrough/simulation.json` predicts breakthrough of 0.1 CO2 and 0.1 N2
  in a 0.8 helium carrier gas.
- `generate_particle_distribution.py` reproduces
  `particle_distribution.data` from the analytical distribution used to make
  the example.

Run either simulation from its directory so the generated output files remain
with that example:

```sh
cd examples/MPD
../../build/ruptura

cd breakthrough
../../../build/ruptura
```
