# ROS34PW2 and the Danckwerts inlet on the validation cases

Reruns every breakthrough and swing-adsorption case of the paper's validation studies with the ROS34PW2 integrator
and the Danckwerts inlet (the default since 2026-11-03), with a binary built at `-O3`. When these runs were made,
the stored results of these studies had been computed with CVODE (RK3 for Elfving) and the Dirichlet inlet. Later
on 2026-11-03 the studies themselves were switched to ROS34PW2 with the Danckwerts inlet and rerun in place
(`../../TODO.md` I-BC-4, I-BC-5; Elfving excepted, B-EX-9), so the copies here duplicate those runs; this directory
remains the record of the comparison against the CVODE runs.

| study | source | cases | source integrator |
|---|---|---|---|
| `marx2_ternary` | `marx2/ternary` | 24 | CVODE |
| `stampi-bombelli` | `stampi-bombelli/dl_from_literature`, `dl_from_stampibombelli` | 20 | CVODE (Dense) |
| `stadler` | `stadler` (breakthrough and swing inputs) | 61 | CVODE |
| `liu` | `liu/simulations` (breakthrough inputs) | 50 | CVODE |
| `3.4.1_boundary_conditions_marx_10bar` | `validation/3.4.1_boundary_conditions_marx_10bar` | 20 | CVODE |
| `3.4.1_boundary_conditions_bea_dilute` | `validation/3.4.1_boundary_conditions_bea_dilute` | 5 | CVODE |
| `3.4.2_swing_adsorption` | `validation/3.4.2_swing_adsorption` | 5 | CVODE |
| `3.4.3_sensitivity` | `validation/3.4.3_sensitivity` | 38 | CVODE |

The integrator benchmark (`../3.5_numerical_performance/breakthrough`, Dirichlet by design) and the convection study compare
integrators and schemes themselves and are not part of this rerun.

## Settings

- `generate.py` copies each source input with `"BreakthroughIntegrator": "ROS34PW2"` and
  `"InletCondition": "Danckwerts"`, and removes `TimeStep`, `NumberOfInitTimeSteps`, and every CVODE key. ROS34PW2
  runs at its defaults: relative tolerance 1e-6; absolute tolerances 1e-7 mol/m³, 1e-8 mol/kg, and 1e-5 K; no step
  cap. Everything else, including `WriteInterval`, is the source input's.
- Study names and case paths are those of `../inlet_condition_dirichlet_vs_danckwerts`, whose
  `simulations/danckwerts/` copies are the same cases with CVODE and the Danckwerts inlet.
- Binary: `build/ruptura-o3-20261103`, a frozen copy of `build/macos/ruptura`. It was built on 2026-11-03 at 14:17
  from the working tree on `596a3ad` with the `macos` preset (`Release`, so `-O3 -DNDEBUG` for RUPTURA and
  SUNDIALS, plus `-march=native -ffast-math -fno-finite-math-only`), configured with
  `-DFETCHCONTENT_SOURCE_DIR_SUNDIALS=build/_deps/sundials-src -DBUILD_PYTHON=OFF -DBUILD_APP=OFF`. Unit tests:
  174/175 (B-EX-3).

## Files

- `generate.py` writes `simulations/<study>/<case>/simulation.json` and `simulations/index.json`.
- `run.py` (from `../inlet_condition_dirichlet_vs_danckwerts`) runs the cases with 8 workers. Each case writes
  its output and `run.log` in its own directory. Pass `--binary` as an absolute path, because every case runs in
  its own directory.

```
python generate.py --clean
python run.py --binary /Users/youriran/pusan/ruptura/build/ruptura-o3-20261103 --timeout 300 simulations/<study> ...
```

## Runs (2026-11-03)

Pass 1: every study except `elfving`, with 8 workers and `--timeout 300` (log `run_pass1.log`; under 8 min; 17 GB
of output). All 223
runs finish, including the six 3.4.1 runs and the three Stampi-Bombelli runs that stop under CVODE with the
Danckwerts inlet, and the 36 Liu runs that stopped with `mxstep` in the stored study.

RUPTURA `Total time` per study, summed over the cases. These runs used 8 workers, the CVODE + Danckwerts copies 4,
so the comparison is indicative only:

| study | ROS34PW2 | CVODE + Danckwerts | stored |
|---|---|---|---|
| `marx2_ternary` | 145 s | 137 s | 149 s |
| `stampi-bombelli` | 7.9 s | 2957 s (17 of 20 finish) | 1356 s (19 of 20) |
| `stadler` | 407 s | 755 s | 912 s |
| `liu` | 3.1 s | — | 61 s (14 of 50 finish) |
| `3.4.1_boundary_conditions_marx_10bar` | 130 s | 168 s (14 of 20) | 202 s (14 of 20) |
| `3.4.1_boundary_conditions_bea_dilute` | 18 s | — | 234 s |
| `3.4.2_swing_adsorption` | 326 s | 282 s | 301 s |
| `3.4.3_sensitivity` | 177 s | 124 s | 129 s |

Largest |Δt50|/t50 over the fed components, against the CVODE + Danckwerts copies where they exist (scratch check
with the functions of `../inlet_condition_dirichlet_vs_danckwerts/compare.py`): `marx2_ternary` 8e-6,
`stampi-bombelli` 8.4e-4, `stadler` 2.8e-4, 3.4.1 Marx 1e-5, 3.4.2 1e-7, 3.4.3 5e-6. Against the stored runs
(Dirichlet): BEA dilute 1e-6 (no physical dispersion, so the inlet does not matter), Liu 6.4e-3 (the 14 runs that
finish there; includes the inlet change).

**Elfving is left out.** Since 2026-11-20 `generate.py` does not copy its 13 inputs (`../../TODO.md` B-EX-16): the
copies stood in `simulations/` unrun, and a rerun of every case waited 16 h for each. The stored inputs keep the old
VantHoff scaling of the empty-column water Langmuir site (`../../TODO.md` B-EX-9, B-PHY-10): b is 7.6e8 times the paper value at 298 K, which makes the site nearly
rectangular. Near p = 0 its uptake is a rate of about 1e11 s⁻¹ in the gas balance, in a chemisorption block that
J~ of ROS34PW2 does not contain. ROS34PW2 then takes steps of 3e-9 s (244 000 steps for the first millisecond), and
RK3 runs only because the clamp to non-negative states removes the overshoot. Without that site, or with the
constant affinity that `elfving/generate.py` now writes, ROS34PW2 takes ordinary steps: run 3 reaches t = 3120 s
of 18 000 s in 300 s. These runs wait for the decision on the water model (B-EX-9).
