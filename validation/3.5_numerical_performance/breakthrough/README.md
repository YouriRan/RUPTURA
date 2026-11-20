# Full-breakthrough integrator benchmark

RK3 against CVODE and ROS34PW2 over complete breakthroughs, so that cost and accuracy refer to the same, physically
complete result. Paper: Results, "Integrators & benchmark" (Figure `fig:benchmark`); setup in SI
"Integrator benchmark" (`apx:benchmark_study`, Tables `tab:benchmark_runs`, `tab:monolith_benchmark_settings`,
`tab:benchmark_scaling`). The 0.1 s matrices in `../figures` time single steps only.

```sh
python generate.py --clean    # simulations/ and manifest.csv (147 runs)
python run.py --jobs 4        # wall times to results.csv, run.log per run
python run.py --jobs 4 --repeats 3   # time each run three times: median in elapsed_seconds, all in elapsed_all
python run.py --jobs 4 --binary /path/to/frozen/ruptura   # use a frozen copy if build/ruptura may be rebuilt meanwhile
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

The current results combine two batches. On 2026-11-03 the RK3, CVODE SPGMR, reference, and ROS34PW2 runs of the
precision and stability studies were run with one frozen binary (log: `run_2026-11-03.log`, 51 runs, 40 min):

```sh
python run.py --jobs 4 --binary B --match /rk3/ /cvode/ /reference/ /ros34pw2/
```

The CVODE Dense runs and the components and methods studies are those of 2026-10-28 (`results-2026-10-28.csv` and
`manifest-2026-10-28.csv` keep that batch), made in three batches one after the other, because concurrent `run.py`
processes overwrite each other's `results.csv` (log: `run_2026-10-28.log`):

```sh
python run.py --jobs 4 --binary B --match /rk3/ /cvode/ /reference/ stability components methods   # 50 runs, 38 min
python run.py --jobs 4 --binary B --match marx-T15/cvode-dense monolith-400ppm-fitted-dl/cvode-dense   # 10 runs, 10 min
python run.py --jobs 5 --binary B --match bea-alkanes-C7/cvode-dense   # 5 runs, 3.1 h
```

The notebook writes Figure `fig:benchmark` to `../../../ruptura-2.0-paper/figures/benchmark/`: the top row (a)–(c) as one
full-width figure, `benchmark_precision.pdf` (work–precision for the three cases, shared axes, shared three-row legend:
circles RK3 in blue, squares CVODE SPGMR in orange, triangles ROS34PW2 in green, dark = coarse Δt or rtol, crosses
= unstable RK3 runs and CVODE runs that stop), and the bottom row as `benchmark_evaluations.pdf` (d, wall time against
mixture-prediction calls) and `benchmark_speedup.pdf` (e, speed-up of CVODE SPGMR and ROS34PW2 at their default
tolerances, rtol 1e-5 and 1e-6, against the RK3 Courant number). The SI figures `benchmark_dense_*.pdf` compare the
two CVODE linear solvers.

## Cases and studies

| Study | Case | What limits RK3 | RK3 Δt / s | CVODE |
|---|---|---|---|---|
| precision | `bea-alkanes-C7`: example input, isothermal, 7 × 1% C7 in He, 1200 s | convection (Courant 0.98 at 1e-2 s; 1.25e-2 stable, 1.5e-2 unstable) | shared grid 1e-4 … 1e-2 | SPGMR and Dense, rtol 1e-2 … 1e-6; reference SPGMR 1e-8; ROS34PW2 rtol 1e-2 … 1e-6 |
| precision | `marx-T15`: `marx2/ternary/T15` (energy balance, Marx h_fw, C_p,g per kg), 1500 s | 3e-2 s stable, 5e-2 s unstable (Δz/v = 0.437 s) | shared grid 1e-4 … 1e-2 | SPGMR and Dense, rtol 1e-2 … 1e-6; reference SPGMR 1e-8; ROS34PW2 rtol 1e-2 … 1e-6 |
| precision | `monolith-400ppm-fitted-dl`: `stampi-bombelli/dl_from_stampibombelli/monolith_400ppm_1mmol`, isothermal, 3000 s | dispersion: 0.63 Δz²/D = 1.66e-3 s; 1.5e-3 stable, 2e-3 unstable | shared grid 1e-4 … 1e-2 (3.125e-3 and 1e-2 unstable) | SPGMR and Dense, rtol 1e-2 … 1e-6; reference SPGMR 1e-8; ROS34PW2 rtol 1e-2 … 1e-6 |
| components | BEA with the first 2, 4, 6 isomers (removed feed goes to He) | | 1e-2 | SPGMR, default tolerances |
| methods | T15 with Sips → single-site Langmuir (same q_sat, k); IAST, SIAST, SEI, SPI | | 2.5e-2 | SPGMR, default tolerances |

Added 2026-11-19 for a fair comparison with ROS34PW2, not yet run (`../../../TODO.md` V-ITG-3): three CVODE series at
rtol 1e-2 … 1e-6 without the 4 Δz/v cap (`CVODEMaximumTimeStep` 0), `cvode-nocap` (SPGMR), `cvode-precond` (SPGMR with
the analytical sorption Jacobian as left preconditioner), and `cvode-klu` (KLU with that Jacobian; BEA only, since KLU
rejects the energy balance of T15 and the chemisorption site of the monolith), and the case
`monolith-400ppm-fitted-dl-danckwerts`, the monolith with the Danckwerts inlet of its validation runs, with every series.

Every RK3 step divides the 2 s output interval, CVODE returns its solution every 2 s, and ROS34PW2 ends a step on every
output time, so all runs sample the outlet at the same times. All runs use the default WENO3 convection and the
Dirichlet inlet (`InletCondition` `Dirichlet`, pinned in `generate.py` since the default became Danckwerts). ROS34PW2
uses the tolerances of CVODE (absolute tolerances scaled likewise) and no step cap; its default, rtol 1e-6 with
1e-7 mol/m³, 1e-8 mol/kg, 1e-5 K, is the rtol 1e-6 run. SPGMR is matrix-free and unpreconditioned
(Krylov dimension 30); Dense is the default `CVODELinearSolver`, with a finite-difference Jacobian and a dense LU.
The CVODE absolute tolerances scale with rtol from their defaults (1e-6 mol/m³, 1e-7 mol/kg, 1e-4 K at rtol 1e-5).
CVODE's internal step is capped at 4 Δz/v by default: 4.06e-2 s (BEA), 1.749 s (T15), 0.247 s (monolith), with Δz/v
= 1.015e-2 s, 0.437 s, and 6.18e-2 s.

Error: e50 = max_j |t50,j − t50,j^ref| / t50,j^ref over the adsorbing components. Runs made with the
step-based binary (before 2026-10-26) labelled each output with the pre-step time; their output times are
shifted by one `TimeStep`, which is kept in each run's `_legacy_simulation.json`. Newer runs are labelled
correctly and are not shifted (`analysis.ipynb`, `output_time_offset`).

Stability: an RK3 run is unstable when it stops, when its outlet exceeds twice the largest roll-up of the reference,
or when no adsorbing component reaches half its feed value at the outlet. A completed CVODE run is plotted at its
error, whatever it is; a CVODE run that stops is drawn at e50 = 1 with its wall time until the stop.

Mixture-prediction calls (column-wide equilibrium evaluations): RK3 3 per step; CVODE 2 per right-hand-side
evaluation (two fixed-point passes, `src/integrators/cvode.cpp`), counting the linear-solver evaluations from
`run.log` (Krylov iterations with SPGMR, finite-difference Jacobians with Dense), plus 2 per output; ROS34PW2 4 per
step and 4 per rejected step, as `run.log` reports them ("mixture predictions").

Timing: Apple M1 Max, Apple clang 17, RelWithDebInfo, whose `-O2` the later `-O3 -march=native -ffast-math
-fno-finite-math-only` of `RUPTURA_COMPILE_OPTIONS` overrides; SUNDIALS 7.5.0; wall time of the whole process, each
run timed once. 2026-11-03 batch: a private RelWithDebInfo build of the working tree (13:5x, after the ROS34PW2
rework), frozen as `ruptura_bench_20261103`, 4 runs at a time; the load besides was a VS Code renderer and Spotlight
indexing (about 4 cores). 2026-10-28 batch: Binary: `build/ruptura` of 2026-10-28 14:02 (working tree on top of `596a3ad`), frozen for the batch.
RK3, CVODE SPGMR, and the T15 and monolith Dense runs ran 4 at a time; the five BEA Dense runs ran 5 at a time. They
slowed each other down. A probe of the rtol 1e-5 run alone took about 12 min for the first 100 s; in the batch,
while all five ran, it advanced about 2.4× slower. The likely cause is the memory bandwidth of five 2727 × 2727
dense LUs. The BEA Dense wall times are therefore upper bounds.
A VS Code renderer used about 2 of the 10 cores throughout.

## Results (2026-11-03; ROS34PW2 added; RK3, CVODE SPGMR, and references rerun)

| Case | Integrator | Setting | Wall time / s | e50 | Mixture-prediction calls |
|---|---|---|---|---|---|
| BEA | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 22.1 / 64.7 / 193 / 614 / 1889 | 1.7e-5 / 2.0e-7 / 2.0e-7 / 2.0e-7 / 2.0e-7 | 3.6e5 … 3.6e7 |
| BEA | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 20.7 / 27.7 / 25.1 / 15.5 / 14.8 | 3.1e-3 / 2.1e-4 / 7.1e-5 / 2.3e-7 / 4.7e-8 | 3.6e5 … 5.6e5 |
| BEA | ROS34PW2 | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 1.78 / 1.92 / 1.61 / 1.31 / 1.44 | 5.1e-3 / 4.9e-4 / 1.8e-6 / 1.9e-6 / 9.7e-8 | 3.1e3 … 6.6e3 |
| T15 | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 29.1 / 88.1 / 264 / 821 / 2529 | 5.1e-8 (all) | 4.5e5 … 4.5e7 |
| T15 | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 3.84 / 2.28 / 2.48 / 2.83 / 3.07 | 7.4e-4 / 6.5e-5 / 1.0e-5 / 1.0e-6 / 2.0e-7 | 2.5e4 … 4.5e4 |
| T15 | ROS34PW2 | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 0.94 / 1.04 / 1.24 / 1.90 / 3.27 | 2.6e-4 / 5.7e-5 / 2.0e-5 / 3.0e-6 / 4.9e-7 | 3.3e3 … 2.4e4 |
| Monolith | RK3 | Δt 1e-2, 3.125e-3 (unstable) / 1e-3 / 3.125e-4 / 1e-4 s | – / 74.9 / 232 / 728 | – / 4.7e-2 / 3.2e-5 / 2.9e-8 | 9.0e6 … 9.0e7 |
| Monolith | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 5.95 / 6.49 / 9.38 / 7.92 / 10.7 | 2.1e-3 / 4.6e-4 / 5.5e-4 / 4.0e-4 / 8.9e-8 | 7.8e5 … 1.4e6 |
| Monolith | ROS34PW2 | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 0.56 / 0.52 / 0.56 / 0.64 / 0.53 | 4.5e-7 / 3.9e-7 / 3.8e-7 / 3.8e-7 / 3.8e-7 | 6.0e3 … 6.3e3 |

The references (SPGMR, rtol 1e-8) take 15.1 s (BEA), 5.3 s (T15), and 30.5 s (monolith). They agree with those of
2026-10-28 to e50 1.0e-7 (BEA), 1.7e-9 (T15), and 1.6e-9 (monolith).

- **ROS34PW2 at its default tolerance (rtol 1e-6):** BEA 1.44 s, e50 9.7e-8 (CVODE at its default 15.5 s,
  2.3e-7: ×10.8); T15 3.27 s, 4.9e-7 (CVODE 2.83 s, 1.0e-6); monolith 0.53 s, 3.8e-7 (CVODE 7.92 s, 4.0e-4; CVODE
  at rtol 1e-6 10.7 s, 8.9e-8: ×20). Steps 1019, 6058, 1564 with 1, 2, 1 rejected; mean step 116, 0.57, 31 Δz/v.
  Speed-up × RK3 Courant number: 13–15 (BEA), 0.18–0.20 (T15), 2.2–2.3 (monolith); CVODE: 1.2–1.4, 0.20–0.24,
  0.15.
- **ROS34PW2 at loose tolerances:** on BEA, rtol 1e-4 … 1e-2 rejects 155–288 steps (an oscillating step size); the
  error then rises to 5.1e-3. On the monolith, one step per 2 s output meets every tolerance (1507–1564 steps),
  so the outputs set the step and the error stays at 3.8e-7 … 4.5e-7. That floor does not shrink with shorter steps
  (3.1e-7 at rtol 1e-9 with steps ≤ 0.02 s); on 12 cells converged ROS34PW2 and RK3 agree to 1e-7 of the feed.
- **Cost per call:** ROS34PW2 2.9–4.2e-4 s on BEA, 1.4–2.8e-4 s on T15, 8.5e-5–1.05e-4 s on the monolith, 2–13×
  RK3 and CVODE (4.1–6.1e-5, 5.6–9.0e-5, 6.0–8.3e-6 s): with few calls, the 600–1500 outputs (the same for every
  integrator) and the linear algebra take most of a ROS34PW2 run.
- **RK3 and CVODE against 2026-10-28:** every CVODE run and reference takes exactly as many steps as before; all
  wall times are 16–38% shorter (most 18–25%), uniformly over the integrators, so the ratios between them are
  unchanged.

## Results (2026-10-28; WENO3, CVODE defaults rtol 1e-5 and 4 Δz/v, SPGMR and Dense)

Changes since 2026-10-21: WENO3 convection (default since 2026-10-27), the CVODE defaults rtol 1e-5 and a 4 Δz/v
step cap (`../../../TODO.md` P-ITG-1), energy-conserving gas–wall heat exchange (C-PHY-2) and `ParticleDensity` 456.6
kg/m³ (B-EX-8) for T15, the speed-ups of `707e64c` (P-GEN-1), and the new Dense series. The 2026-10-21 `results.csv`
and `manifest.csv` are kept as `results-2026-10-21.csv` and `manifest-2026-10-21.csv`.

| Case | Integrator | Setting | Wall time / s | e50 | Mixture-prediction calls |
|---|---|---|---|---|---|
| BEA | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 22.4 / 67.6 / 200 / 623 / 1933 | 1.7e-5 / 2.0e-7 / 2.0e-7 / 2.0e-7 / 2.0e-7 | 3.6e5 … 3.6e7 |
| BEA | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 21.2 / 28.2 / 25.6 / 16.1 / 15.2 | 3.1e-3 / 2.1e-4 / 7.1e-5 / 2.3e-7 / 4.7e-8 | 3.6e5 … 5.6e5 |
| BEA | CVODE Dense | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 11208 / 7186 / 6276 / 5956 / 6024 | 1.4e-3 / 2.8e-4 / 1.0e-4 / 8.9e-5 / 1.4e-5 | 4.2e6 … 2.1e7 |
| T15 | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 26.4 / 79.5 / 241 / 733 / 2261 | 5.5e-5 / 5.4e-5 / 5.4e-5 / 5.4e-5 / 5.4e-5 | 4.5e5 … 4.5e7 |
| T15 | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | stops (2.0) / 2.4 / 2.7 / 2.9 / 2.8 | – / 1.9e-4 / 1.1e-5 / 2.8e-6 / 5.8e-7 | 2.9e4 … 3.9e4 |
| T15 | CVODE Dense | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 135 / 104 / 102 / 136 / 127 | 3.1e-3 / 2.0e-4 / 1.9e-5 / 4.7e-6 / 6.7e-7 | 3.9e5 … 7.4e5 |
| Monolith | RK3 | Δt 1e-2, 3.125e-3 (unstable) / 1e-3 / 3.125e-4 / 1e-4 s | – / 86.8 / 264 / 817 | – / 4.7e-2 / 3.2e-5 / 1.4e-8 | 9.0e6 … 9.0e7 |
| Monolith | CVODE SPGMR | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 6.8 / 7.9 / 11.4 / 9.7 / 11.7 | 8.8e-4 / 5.1e-4 / 4.5e-5 / 1.1e-5 / 3.0e-7 | 8.1e5 … 1.4e6 |
| Monolith | CVODE Dense | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 11.0 / 6.2 / 10.5 / 21.6 / 446 | 0.48 / 1.8e-3 / 6.6e-4 / 9.3e-5 / 6.7e-5 | 4.3e5 … 6.0e7 |

The references (SPGMR, rtol 1e-8) take 16.0 s (BEA), 5.0 s (T15), and 35.8 s (monolith).

- **CVODE with SPGMR is now faster than RK3 at every stable RK3 step, in all three cases.** At the default
  tolerance, speed-up × RK3 Courant number is 1.19–1.37 for BEA, 0.18–0.21 for T15, and 0.136–0.144 for the
  monolith. On BEA, break-even is at a Courant number of 1.2–1.4, at the RK3 stability limit (1.23 stable, 1.48
  unstable); with the old Δz/v cap it was near 0.5. At the default tolerance, SPGMR takes 29,864 steps on BEA (the
  4 Δz/v cap allows 29,550), 2,458 on T15 (cap: 858), and 20,211 on the monolith (cap: 12,131). Each step costs
  6.1, 6.8, and 28.0 right-hand-side evaluations.
- **Wall time ∝ mixture-prediction calls for RK3 and SPGMR.** The cost per call is 4.2–6.2e-5 s on BEA,
  5.0e-5–9.2e-5 s on T15 (where CVODE calls cost more), and 8.2–9.6e-6 s on the monolith.
- **Dense is 245–530× slower than SPGMR on BEA and 38–47× slower on T15.** Its cost is not proportional to the
  mixture-prediction calls (panel d): 5.3e-4–1.4e-3 s per call on BEA (10–33× SPGMR) and 1.8–3.3e-4 s on T15. The BEA
  state has 2727 variables, so each Jacobian costs 2727 right-hand sides plus a dense LU; the runs need 759–3862
  Jacobian evaluations. On the monolith (smaller state), Dense matches SPGMR down to rtol 1e-4 and takes 21.6 s at
  1e-5, but 446 s at 1e-6: 235,386 steps (SPGMR: 26,218), 61,940 Newton convergence failures, and 64,576 Jacobian
  evaluations. Speed-up × Courant number of Dense: 0.0032–0.0037 (BEA), 0.0038–0.0045 (T15), 0.061–0.065 (monolith).
- **Dense is also less accurate at the same tolerance on BEA and the monolith.** At rtol 1e-6, e50 is 1.4e-5
  against 4.7e-8 (BEA) and 6.7e-5 against 3.0e-7 (monolith); on T15 the two agree within a factor of two. The BEA
  and monolith Dense runs end at BDF order 1, and the BEA runs have 294–3552 Newton convergence failures. The cause
  was not investigated.
- **Failures.** T15 with SPGMR at rtol 1e-2 stops at t = 187.6 s ("the corrector convergence test failed repeatedly
  or with |h| = hmin", h = 2.2e-6 s); with upwind and the Δz/v cap (2026-10-21) it completed. Dense at the same
  tolerance completes. The monolith with Dense at rtol 1e-2 completes but overshoots: the outlet reaches 4.3 times
  the feed, and t50 = 659 s against 1272 s (e50 = 0.48).
- **RK3 stability runs now complete.** On 2026-10-21, T15 and the monolith stopped with a negative outlet pressure.
  Now BEA at 1.5e-2 s and the monolith at 2e-3 s run into unbounded outlets (58× and 2500× the feed), and T15 at
  5e-2 s collapses: every outlet mole fraction is 0 after a few seconds. The monolith runs at 3.125e-3 and 1e-2 s
  also complete with unbounded outlets. The largest stable steps (1.25e-2, 3e-2, 1.5e-3 s) were not rerun; WENO3
  does not change them (`../convection`).
- **RK3 convergence floor.** RK3 converges to e50 = 2.0e-7 on BEA (2026-10-21, upwind: 3.0e-9 and decreasing) and
  5.4e-5 on T15 (4.9e-5). On the monolith, where the velocity is prescribed, it converges to 1.4e-8. The likely cause
  is the two-pass coupling of the CVODE right-hand side (B-ITG-7), now also visible on BEA; this was not checked.
- **Components (BEA):** RK3 10.6 / 16.0 / 20.9 / 22.4 s and CVODE (SPGMR) 7.0 / 11.1 / 14.5 / 16.1 s for 2 / 4 / 6 / 7
  isomers: ×2.1 and ×2.3.
- **Methods (T15, Langmuir):** RK3 IAST 6.37, SIAST 6.38, SEI 5.14, SPI 7.71 s; CVODE (SPGMR) 1.66, 1.63, 1.47,
  1.82 s.

## Results (2026-10-21; upwind, CVODE defaults rtol 1e-4 and Δz/v; superseded)

| Case | Integrator | Setting | Wall time / s | e50 | Mixture-prediction calls |
|---|---|---|---|---|---|
| BEA | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 29.3 / 84.3 / 248 / 768 / 2367 | 1.8e-6 / 1.1e-7 / 3.4e-8 / 1.0e-8 / 3.0e-9 | 3.6e5 … 3.6e7 |
| BEA | CVODE | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 37.7 / 42.9 / 53.3 / 52.3 / 47.8 | 7.4e-4 / 5.4e-5 / 2.4e-4 / 2.2e-6 / 1.9e-7 | 4.5e5 … 6.9e5 |
| T15 | RK3 | Δt 1e-2 / 3.125e-3 / 1e-3 / 3.125e-4 / 1e-4 s | 39.3 / 110 / 317 / 949 / 2867 | 4.9e-5 (all) | 4.5e5 … 4.5e7 |
| T15 | CVODE | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 3.5 / 3.6 / 3.3 / 3.7 / 3.8 | 3.1e-4 / 7.5e-6 / 4.6e-6 / 1.6e-7 / 3.0e-8 | 2.5e4 … 3.1e4 |
| Monolith | RK3 | Δt 1e-2, 3.125e-3 (unstable) / 1e-3 / 3.125e-4 / 1e-4 s | – / 132 / 409 / 1277 | – / 4.7e-2 / 3.4e-5 / 8.6e-8 | 9.0e6 … 9.0e7 |
| Monolith | CVODE | rtol 1e-2 / 1e-3 / 1e-4 / 1e-5 / 1e-6 | 12.3 / 11.3 / 14.4 / 15.9 / 37.1 | 2.4e-2 / 6.8e-3 / 7.9e-4 / 2.7e-4 / 1.0e-5 | 7.9e5 … 3.0e6 |

- **Main result: wall time ∝ mixture-prediction calls.** Per call: BEA 6.7–8.3e-5 s, T15 7.0e-5–1.4e-4 s
  (CVODE calls cost more: more IAST iterations from a worse starting guess), monolith 1.2–1.5e-5 s. Within a
  case the cost per call varies by less than a factor of two whichever the integrator; the integrators
  differ in the number of calls. RK3: 3 per step, step set by stability. CVODE: 2 per RHS evaluation, 2.9
  (BEA), 3.4 (T15), 10.3 (monolith) RHS evaluations per step, and the step sits at the Δz/v cap in all three
  cases at the default tolerance. Calls(RK3)/calls(CVODE) ≈ (3 / 2r) (Δz/v) / Δt_RK3.
- **Speed-up vs RK3 Courant number (paper panel d):** speed-up × Courant number is constant per case:
  0.44–0.54 (BEA), 0.20–0.27 (T15), 0.14–0.15 (monolith) = 3/(2r) × cost-per-call ratio (SI Eq.
  `eq:benchmark_speedup`). CVODE breaks even near Courant 0.5 on BEA. RK3 stability: BEA stable at 1.25e-2 s,
  unstable at 1.5e-2 s (Courant 1.23 / 1.48); T15 3e-2 / 5e-2 s (0.069 / 0.114); monolith 1.5e-3 / 2e-3 s
  (0.024 / 0.032).
- **BEA (convection-limited):** Δt_RK3 = Δz/v, so CVODE needs 1.3–1.9× more calls and is slower (38–53 s
  against 29 s). Its tolerance then sets only the Newton–Krylov accuracy (default rtol: e50 = 2.4e-4).
- **T15:** RK3 at 1e-2 s (Courant 0.023): 39 s, against 3.3 s for CVODE at default tolerance (e50 = 4.6e-6);
  extrapolated to the largest stable step (3e-2 s) RK3 would be ~4× slower. RK3 converges to 4.9e-5 from the reference:
  the two-pass coupling of the CVODE right-hand side (kept by decision, `../../../TODO.md`). A test build with eight
  passes agrees with RK3 (Δt = 1e-3 s) to 2.4e-6 at 2.3× the wall time.
- **Monolith (dispersion-limited):** vΔz/D = 0.04. RK3 is unstable at 2e-3 s and inaccurate close to the
  limit (e50 = 4.7e-2 at 1e-3 s; a scan of 1e-3 … 1.25e-4 s converges smoothly). RK3 needs 409 s for 3.4e-5; CVODE
  reaches 1.0e-5 in 37 s. At the default tolerance, CVODE gives 7.9e-4, close to atol / c_CO₂,in = 6e-4.
- **Components (BEA):** RK3 13.6 / 21.1 / 27.4 / 29.3 s and CVODE 22.9 / 34.4 / 44.6 / 53.3 s for 2 / 4 / 6 / 7
  isomers: ×2.2 and ×2.3, through the cost per call (the number of calls hardly changes).
- **Methods (T15, Langmuir):** RK3 IAST 10.3, SIAST 10.4, SEI 7.07, SPI 7.77 s; CVODE 2.54, 2.61, 2.23, 2.33 s.

### History

- A first version (before 2026-10-20 evening) used the molar gas heat capacity of the Marx inputs (29.51
  J mol⁻¹ K⁻¹ in a J kg⁻¹ K⁻¹ field, `../../3.3.1_energy_balance_marx_ternary/marx2/conversions.md` §5). The gas energy balance was then 38×
  too light, and RK3 on T15 was unstable at 1e-3 s: 599 s at 5e-4 s against 6 s with CVODE.
- It also returned the CVODE solution every 1 s, which recorded CVODE at odd and RK3 at even seconds; the
  linear interpolation of t50 then differed by 3.4e-4 on the sharp H₂ front of T15 and hid the RK3
  convergence.
- 2026-10-21: the RK3 grid became shared (logspace(-4, -2, 5) with 3.125 for 3.16), and a stability study was added.
  While these runs were going, another session rebuilt `build/ruptura` from in-progress source changes (adsorbing
  carrier gas, `src/mixture_prediction.cpp`). Runs started on the intermediate build (10:41–11:03) were wrong for
  the monolith (t50 −5% at 3.125e-4 s); all of them were rerun with a frozen copy of the 11:03 build, which
  reproduces the CVODE outputs of 2026-10-20 bit for bit. Use `run.py --binary` with a frozen copy for long runs.
- BDF order: on BEA the CVODE error is not monotonic in rtol (order 5 at rtol ≤ 1e-4); a separate build with
  `CVodeSetMaxOrd(3)` gives 5.4e-5 / 2.9e-6 / 6.0e-7 for rtol 1e-3 / 1e-4 / 1e-5 (`../../../TODO.md`).
- 2026-11-03: ROS34PW2 series added (15 directories `precision/<case>/ros34pw2/`), `InletCondition` `Dirichlet`
  pinned in `generate.py` (all earlier runs used that default), and the RK3, CVODE SPGMR, reference, and stability
  runs rerun with the same binary as ROS34PW2; the Dense, components, and methods runs are those of 2026-10-28.
- 2026-10-28: rerun with WENO3 and the new CVODE defaults, and a Dense series in the precision study (15 new
  directories `precision/<case>/cvode-dense/`). `generate.py` now scales the tolerances from the new defaults
  (rtol 1e-5), which gives the same inputs as before, up to one ulp. The stability criterion of `analysis.ipynb`
  now also catches a collapsed outlet (T15 at 5e-2 s), applies to RK3 only, and the notebook draws CVODE runs
  that stop as crosses.
