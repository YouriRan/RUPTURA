# Convection scheme: first-order upwind against WENO3

Compares the two settings of `ConvectionScheme` on the three breakthrough cases of `../breakthrough`. Paper: Results,
"Integrators & benchmark" (Figure `fig:convection`); method in SI "Method of lines" (Eqs.
`eq:discrete_species_balance`–`eq:weno3_gas_flux`); setup and numbers in SI "Integrator benchmark", paragraph
"Convection scheme" (Table `tab:convection_benchmark`).

```sh
python generate.py --clean                                   # simulations/ and manifest.csv (65 runs)
python run.py --jobs 4 --match cost/                         # timed runs first, with little contention
python run.py --jobs 8 --match stability/ grid/ reference/   # about 40 min; the monolith reference takes ~30 min
python run.py --binary /path/to/frozen/ruptura ...           # use a frozen copy if build/ruptura may be rebuilt
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

The notebook prints the numbers quoted in the paper and writes the panels to
`../../../ruptura-2.0-paper/figures/convection/`.

## Cases and studies

| Study | Settings |
|---|---|
| cost | input grid (BEA and T15 100 points, monolith 50); RK3 at the benchmark step (BEA 1e-2 s, T15 2.5e-2 s, monolith 1e-3 s) and CVODE with SPGMR at the default tolerances; two repeats |
| stability | RK3 at BEA 1.25e-2, 1.5e-2, 1.75e-2 s (Courant 1.23, 1.48, 1.72); T15 3e-2, 5e-2 s; monolith 1.5e-3, 2e-3 s |
| grid | CVODE, rtol 1e-6, 25, 50, 100, 200 grid points |
| reference | CVODE, WENO3, rtol 1e-6, 400 grid points |

Cases as in `../breakthrough`: seven C7 isomers at 1% in helium on BEA (no axial dispersion, convection-limited), Marx
et al. (2014) T15 with energy balance (D_ax = 4.3e-6 m²/s, 38 times below the upwind numerical dispersion vΔz/2 at
100 points), and the Stampi-Bombelli et al. (2024) monolith at 400 ppm with the fitted D_ax = 2.6e-3 m²/s (numerical
dispersion 2% of it).

Accuracy: e50 (relative error of the breakthrough times), E1 (time-averaged absolute difference of the normalized
outlet curves), and the 10–90% front width of the most strongly adsorbed component, each against the reference.

Timing: Apple M1 Max, 4 concurrent single-threaded runs (cost study), Apple clang `-O2` (RelWithDebInfo); wall time of
the whole process and CPU time from `/usr/bin/time`.

## Background

- First-order upwinding adds a numerical dispersion of about vΔz/2. WENO3 (Jiang and Shu, 1996) reconstructs the face
  values from three cells with smoothness-dependent weights: third order where the profile is smooth, one-sided and
  free of oscillations across a front.
- In a gas, WENO3 reconstructs the mole fractions and keeps the total molar flux upwind, so that the species fluxes
  sum to the flux that the velocity closures (`InletPressureInletVelocity`, `InletVelocityOutletPressure`) conserve.
  Limiting each species flux separately breaks this: T15 then stops after ~100 s with a negative outlet pressure.
- The inlet value lies on the inlet face, half a cell from the centre of the first cell; the dispersion term of the
  first cell uses the same position. With the vertex-centred inlet difference instead, the monolith t50 was offset
  by about Δz/2 divided by the front speed (+1.1% at 50 points).
- A van Leer-limited MUSCL scheme was tested in a scratch build (2026-10-26) with the same accuracy as WENO3 but a
  smaller RK3 stability range on BEA (unstable at Courant 1.23): on a smooth profile it behaves like Fromm's scheme,
  whose linear RK3 limit is a Courant number of 1.18, against 1.26 for upwind and 1.63 for the third-order WENO3 weights.
