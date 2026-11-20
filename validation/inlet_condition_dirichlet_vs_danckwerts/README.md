# Inlet condition: Dirichlet against Danckwerts on the validation cases

`InletCondition` (2026-11-03) selects the species and gas-temperature inlet. The default was Dirichlet when these runs were made; it is
Danckwerts since later that day.

- `Dirichlet`: the inlet face holds the feed state. Axial dispersion (and heat conduction) across the
  inlet face then adds to the feed flux, by about a fraction D_ax/(vL) of the feed up to the breakthrough (TODO
  B-BC-4).
- `Danckwerts`: the stencils of the first cell use the inlet-face value c_b = (u c_in + g c_1)/(u + g), g = D/δ, so
  the convective and dispersive inflow equals the feed flux u c_in (SI, species and fluid-temperature boundary
  conditions).

This study reruns every validation case of the paper that has physical axial dispersion with both inlets, using the
same binary and settings:

| study | source | cases | vL/D_ax |
|---|---|---|---|
| `marx2_ternary` | `marx2/ternary` | 24 | ≈ 3300 |
| `stampi-bombelli` | `stampi-bombelli/dl_from_literature`, `dl_from_stampibombelli` | 20 | 7–450 (literature), 0.04–2 (fitted) |
| `stadler` | `stadler` (breakthrough and swing inputs) | 61 | 120–490 |
| `3.4.1_boundary_conditions_marx_10bar` | `validation/3.4.1_boundary_conditions_marx_10bar` | 20 | ≈ 3300 |
| `3.4.2_swing_adsorption` | `validation/3.4.2_swing_adsorption` | 5 | ≈ 3300 |
| `3.4.3_sensitivity` | `validation/3.4.3_sensitivity` | 38 | 33–3300 |

The BEA cases have no physical dispersion, so both inlets are identical there; they are not rerun.

The CVODE inputs of `stampi-bombelli` do not set `CVODELinearSolver` and were run with the former default, Dense.
`generate.py` writes Dense into these copies, because with SPGMR, the default since 2026-11-02, the steep 5.6%
fronts stall: `monolith_5.6pct_0p02mmol` stops with `mxstep` at t = 883 s (Dirichlet) and 699 s (Danckwerts), but
finishes in 21 s with Dense for both inlets (TODO V-GEN-4).

## Files

- `generate.py`: writes `simulations/<inlet>/<study>/<case>/simulation.json` and `simulations/index.json`.
- `run.py`: runs them with 8 workers. Run the 3.4.1 cases with `--timeout 900`, as in that study, since the imposed-pressure runs of T3 and T9 stall.
- `compare.py`: computes t10, t50, and t90 of every fed component, and the largest outlet mole-fraction difference,
  for Danckwerts against Dirichlet and for Dirichlet against the stored results. It writes `comparison.csv`.

```
python generate.py --clean
python run.py -j 4 --binary <frozen ruptura> simulations/dirichlet/stampi-bombelli simulations/danckwerts/stampi-bombelli ...
python run.py -j 4 --binary <frozen ruptura> --timeout 900 simulations/*/3.4.1_boundary_conditions_marx_10bar
python compare.py
```

The runs of 2026-11-03 used the binary built from the shared tree at 11:50 that day (frozen copy
`ruptura_danckwerts`).

## Results (2026-11-03)

**Reproduction.** Every Dirichlet copy reproduces the stored output of its source byte for byte (component files
compared with `cmp`; `compare.py`: zero difference in all 168 cases). The new switch leaves the Dirichlet path
unchanged, and the stored results of these studies are current. Runs that stop do so at the same time in the stored
study, in the Dirichlet copy, and in the Danckwerts copy:
- in 3.4.1, the four FixedVelocity runs and FixedPressureInletVelocity T3/T9, both inlets;
- `stampi-bombelli/dl_from_stampibombelli/monolith_5.6pct_0p02mmol`, at t = 72 s with both inlets.

**Danckwerts against Dirichlet** (largest changes; all numbers in `comparison.csv`):

| study | t50 | 10–90% width | note |
|---|---|---|---|
| `marx2_ternary` | ≤ 0.05% | ≤ 0.4% | vL/D ≈ 3300 |
| `3.4.1` | ≤ 0.06% | ≤ 0.8% | |
| `3.4.2` | ≤ 0.01% | ≤ 0.04% | |
| `stadler` | ≤ 0.3% for the WGS conversion runs | | 1.1% in the first seconds of the fresh-bed run |
| `3.4.3` | 2.7% at D_ax = 1e-3 m²/s, ≤ 0.24% otherwise | ≤ 1.5% | see below |
| `stampi-bombelli`, literature D_L | 5.6% monolith: +19%, +6%, +2% at 0.02, 0.06, 0.18 mmol/s; others ≤ 2% | 5.6% monolith: +8%, +32%, −50%; others ≤ 0.3% | vL/D 7–450 |
| `stampi-bombelli`, fitted D_L | 1.7–3.4 times later in the monolith | 1.5–5.6 times wider | vL/D 0.04–2 |

- **T15 dispersion sweep (3.4.3).** At D_ax = 1e-3 m²/s, the CO₂ t50 is 2.9% earlier than without dispersion with
  the Dirichlet inlet (the paper's number), but only 0.24% earlier with Danckwerts. The broadening of the front
  (+46%) remains. `figures/t15_dispersion_sweep_t50.pdf`.
- **Stampi-Bombelli, literature D_L (paper runs).** The CO₂ t50 against the experiment is:

  | contactor | Dirichlet | Danckwerts |
  |---|---|---|
  | 5.6% monolith at 0.02 / 0.06 / 0.18 mmol/s | −11% / −8% / −13% | +6% / −2% / −11% |
  | 400 ppm monolith | +14% to +25% | +14% to +25% |
  | 5.6% packed bed | +6% to +17% | +6% to +18% |

  The fronts stay steeper than measured. `figures/monolith_5p6pct_dl_literature.pdf`.
- **Stampi-Bombelli, fitted D_L.** Dirichlet puts the breakthrough far too early: −29% to −72% in the monolith, and
  within seconds in the packed bed (vL/D ≈ 0.04), where the dispersive inflow exceeds the feed by orders of
  magnitude. With Danckwerts, the authors' fitted coefficients give t50 within +13% to +24% (400 ppm monolith),
  −5% and −15% (5.6% monolith, 0.06 and 0.18 mmol/s), and +18% (packed bed, 2 mmol/s). The 10–90% widths are 60–133%
  of the measured ones, against 11–72% in the monolith with the literature D_L and Danckwerts. `figures/monolith_*_dl_fitted.pdf`,
  `figures/packed_bed_5p6pct_dl_fitted.pdf`.
- **New stops with Danckwerts:** fitted-D_L packed bed at 0.5 and 1 mmol/s (CVODE `mxstep` at 72 s and 96 s;
  Dirichlet finishes).
  - With `CVODEMaximumOrder` 1 (TODO B-ITG-11), 1 mmol/s finishes with t50 = 2976 s (experiment 2612 s) and a
    10–90% width of 1493 s (experiment 1834 s).
  - 0.5 mmol/s still crawls at order 1 (t = 768 s after about 25 min; stopped).
  - These inputs are not in the paper.
