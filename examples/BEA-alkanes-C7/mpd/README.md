# BEA nC7 / C6m2: a macrostate particle distribution from a Langmuir isotherm

A minimal, self-contained MPD benchmark built on the same system as
`examples/BEA-alkanes-C7/breakthrough`, reduced to a binary mixture of
n-heptane (`nC7`) and 2-methylhexane (`C6m2`) in BEA at 552 K.

The point of the example is that a multisite Langmuir isotherm with a shared
saturation capacity per site family **is** the mean of a product of multinomial
distributions. That lets the macrostate particle distribution be written down in
closed form, handed to Ruptura, and checked against the analytical result — no
Monte Carlo needed to exercise the machinery.

Open `mpd_storytelling.ipynb` for the narrated version. It fits the isotherms,
builds the distribution, plots it, runs all four simulations, and reports the
timings.

## Layout

```
fit_isotherms.py                 GCMC fit  ->  fitted_parameters.json
generate_particle_distribution.py             ->  particle_distribution.data
generate_simulations.py                       ->  the four simulation.json inputs
mpd_storytelling.ipynb           the notebook
mixture/mpd,  mixture/iast       MixturePrediction, 64 log-spaced points, 1e3-1e7 Pa
breakthrough/mpd, breakthrough/iast   Breakthrough, 1 % nC7 / 1 % C6m2 in helium
```

Everything downstream of `fitted_parameters.json` is generated, so the MPD and
the Langmuir/IAST inputs cannot drift apart.

## The model

`fit_isotherms.py` fits the GCMC isotherms in `../fitting`
(`Results.dat-BEA-Repeat-552K-nC7` and `-2mC6`) to

```text
q_i(f) = sum_s q_sat,s b_i,s f / (1 + b_i,s f)
```

with **one saturation capacity per site family, shared by both components** —
the constraint a competitive multinomial site model imposes — and one affinity
per component and family. Nelder-Mead on the absolute residuals; NumPy only.

`Delta H` is not fitted: a single-temperature isotherm carries no information
about it. It is read from the GCMC heat of desorption (column 20), averaged over
the Henry-regime points, as `Q_eq = -R <heat of desorption>`, which is the
positive quantity Ruptura's `HeatOfAdsorption` expects.

`generate_particle_distribution.py` then evaluates, per site family,

```text
Pi_s(n) = N_s! / [(N_s - sum_i n_i)! prod_i n_i!]
          * prod_i (b_i,s f_ref)^n_i / (1 + sum_j b_j,s f_ref)^N_s
```

and convolves the two families over the total counts. The particle window is
fixed at `[0, 200]` and split between the families in proportion to the fitted
capacities, so the framework mass follows as
`m_fw = 200 / (N_A sum_s q_sat,s)`. The file is tabulated with `DeltaN = 2`
(101 nodes per component, 5151 populated macrostates of 10201). Section 6 of the
notebook quantifies the discretization error: machine precision over the
plotted range, growing only in the deep Henry and full-saturation limits.

## Running it

Build Ruptura first, then either open the notebook or run each input from its
own directory so the outputs stay with it:

```sh
python3 fit_isotherms.py
python3 generate_particle_distribution.py
python3 generate_simulations.py

cd mixture/mpd    && ../../../../../build/ruptura
cd ../iast        && ../../../../../build/ruptura
cd ../../breakthrough/mpd  && ../../../../../build/ruptura
cd ../iast                 && ../../../../../build/ruptura
```

Set `RUPTURA_BIN` if the executable is elsewhere; the notebook honours it.

## A note on cost

The breakthrough inputs keep every physical setting of
`examples/BEA-alkanes-C7/breakthrough` but reduce `NumberOfGridPoints` from 100
to 20 and the output `TimeStep` from 0.001 s to 0.105 s — the values used by
`examples/MPD/BEA-alkanes/breakthrough`. An MPD equilibrium evaluation walks
every populated macrostate at every grid point and every internal integration
step, so its cost scales with (macrostates) x (grid points) x (steps). With
5151 macrostates the literal defaults of the parent example would take weeks;
as configured, the MPD run takes on the order of an hour and the Langmuir/IAST
run a few seconds. Both models use identical column settings, so the comparison
is unaffected.

The notebook caches each run in that directory's `benchmark.json` and skips
completed runs; pass `force=True` to `run_ruptura` to recompute.
