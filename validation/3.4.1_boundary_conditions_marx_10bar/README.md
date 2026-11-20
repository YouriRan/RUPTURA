# 3.4.1 Boundary conditions: Marx et al. (2014) ternary experiments at 1e6 Pa

Status: runs done with ROS34PW2 and the Danckwerts inlet (2026-11-03), figure in the paper (Results, "Effect of
boundary conditions"; setup in SI "Boundary-condition study"). Dilute control case:
`../3.4.1_boundary_conditions_bea_dilute`.

Checks the five gas-phase pressure–velocity boundary modes on a strongly adsorbing, non-isothermal case
with experimental reference data. Each base experiment is run once per mode. The modes are
parameterized so that they coincide once the column is filled with the feed (`conversions.md`), so
differences during breakthrough come only from how each mode couples pressure and velocity.

## Base experiments

All at $1\times10^{6}$ Pa and $1\times10^{-5}$ m³/s. Inputs from `../../3.3.1_energy_balance_marx_ternary/marx2/ternary/<run>/simulation.json`,
digitized experiments in `../../marx2/reference_digitized/<run>/`.

| Run | Feed CO₂/N₂/H₂ | $T_\mathrm{feed}$ / K | Note |
|---|---|---|---|
| T3  | 15/80/5  | 298 | |
| T9  | 15/80/5  | 318 | ambient 298 K (SI note on the 318 K experiments) |
| T15 | 40/30/30 | 298 | strongest sorption effect |
| T21 | 40/30/30 | 318 | marx2 had to use `InletVelocityOutletPressure` and stopped at 1174 s |

Only the boundary keys, `SimulationTime` (4000 s), `WriteInterval` (4 s; `NumberOfTimeSteps` 8e6 and
`WriteEvery` 8000 in the step-based input the runs were made with) and `DisplayName` differ from the base
inputs.

## Boundary modes

| JSON name | SI label | Prescribed | Velocity | Pressure |
|---|---|---|---|---|
| `InletVelocityOutletPressure` | BC1 | $v(0)$, $P(L)$ | continuity with sorption sink | Ergun, marched back from the outlet |
| `InletPressureOutletPressure` | BC2 | $P(0)$, $P(L)$ | $v(0)$ by shooting, then continuity | Ergun, marched from the inlet |
| `InletPressureInletVelocity` | BC3 | $P(0)$, $v(0)$ | continuity with sorption sink | Ergun, marched from the inlet |
| `FixedPressureInletVelocity` | BC4 | linear $P(z)$ with gradient $-G$ | local Ergun velocity from local density, no continuity | imposed |
| `FixedVelocity` | BC5 | $v$ everywhere | uniform, no sorption effect | Ergun, marched from the inlet (diagnostic) |

The experiment uses a mass-flow controller and a back-pressure regulator, which corresponds to BC1. All
marx2 runs except T8 and T21 use BC3.

## Smoke test (T15, first 2 s)

All five modes start and run without error. At $t = 1$ s:

- BC1 and BC3 keep $v(0) = 0.02745$ m/s. Uptake at the front halves the downstream velocity (0.0134 m/s).
- BC2 starts at 0.0350 m/s in the helium-filled column (`conversions.md`, section 4). Its inlet velocity
  then rises to 0.060 m/s to hold the prescribed pressure drop against the uptake. The feed rate is
  therefore not controlled, and BC2 cannot represent this experiment.
- BC4 has 0.0350 m/s in the helium region against 0.0275 m/s at the inlet. It shows no sorption effect.
- BC5 stays at 0.02745 m/s everywhere.

## Results (2026-11-20, the inputs of the model of Marx et al., ROS34PW2, `analysis.ipynb`)

Rerun with `marx2/ternary` (TODO V-PHY-12): the fluid fraction is now ε_t, so the Ergun drop is about 3 Pa instead of
31 Pa (outlet pressures 999996.6–999996.9 Pa). CO2 t50 shift against `InletVelocityOutletPressure`:
`InletPressureInletVelocity` −3e-6; `InletPressureOutletPressure` −12.6 % (T3), −13.5 % (T9), −26.6 % (T15),
−28.4 % (T21); `FixedPressureInletVelocity` −0.12 to −0.85 %; `FixedVelocity` +0.01 to −0.33 %. T15 with both
pressures: v(0) = 3.2, 2.75, 2.1 times nominal at 4, 8, 20 s. Without molar continuity the front falls to
2e-8–2e-7 of the feed concentration (20–110 s) and refills after 130–220 s. BEA unchanged. CVODE copies:
`../3.4.1_boundary_conditions_marx_10bar_cvode`.

## Results (2026-11-03, ROS34PW2, `analysis.ipynb`)

`generate.py` now writes `"BreakthroughIntegrator": "ROS34PW2"` (default tolerances) and `"InletCondition":
"Danckwerts"` and drops the CVODE keys of the marx2 inputs (`../../TODO.md` I-BC-5). Run with the frozen `-O3` binary
`build/ruptura-o3-20261103` and `run.py` (8 workers, no timeout); a binary built from the tree after the bulk-reaction
heat change of 14:33 (F-PHY-5) gives byte-identical output. The CVODE runs with the Dirichlet inlet of 2026-11-02
remain, byte for byte, in `../inlet_condition_dirichlet_vs_danckwerts/simulations/dirichlet/`.

Every run completes 4000 s, BC4 and BC5 included, in 2.8–7.8 s (`Total time`; 20 runs, 8 at a time). CO₂ t50 in s
(relative shift against BC1):

| Run | BC1 | BC2 | BC3 | BC4 | BC5 |
|---|---|---|---|---|---|
| T3  | 901.46 | 809.11 (−10.2%) | 901.43 (−3.2e-5) | 900.71 (−0.08%) | 900.96 (−0.05%) |
| T9  | 650.21 | 580.47 (−10.7%) | 650.19 (−3.1e-5) | 649.40 (−0.12%) | 649.60 (−0.09%) |
| T15 | 634.79 | 492.75 (−22.4%) | 634.77 (−2.9e-5) | 632.77 (−0.32%) | 633.51 (−0.20%) |
| T21 | 477.06 | 366.65 (−23.1%) | 477.05 (−3.2e-5) | 474.91 (−0.45%) | 475.48 (−0.33%) |

- BC2 inlet velocity on T15: 2.73 v_in at 4 s and 2.13 v_in at 20 s.
- BC4 and BC5: the total gas concentration at the front falls to 6.3e-8–2.7e-7 of the feed concentration
  between 52 and 100 s (T15 at 20 s: 7.3e-7 with BC4, 1.2e-4 with BC5), and the whole column is above half the
  feed concentration again after 148–204 s. ROS34PW2 rejects about 15% of its steps in these runs (1037–1332 of
  7656–9204), against 1–2 with BC1–BC3.
- Up to the stop of the CVODE runs (Danckwerts copies in `../inlet_condition_dirichlet_vs_danckwerts`), the total
  gas concentrations of ROS34PW2 and CVODE agree to 7e-5 of the feed concentration.
- Dilute BEA control: BC3–BC5 within 1.7e-3 of BC1, BC2 6.4% (`../3.4.1_boundary_conditions_bea_dilute`).

## Results (2026-11-02, `analysis.ipynb`)

Rerun with the marx2 inputs of 2026-11-02 (`ParticleDensity` 456.6, WENO3, CVODE defaults rtol 1e-5 and 4 Δz/v),
the weekend binary `build/ruptura-rerun`, and `run.py --timeout 900`.

| Run | BC1 | BC2 | BC3 | BC4 | BC5 |
|---|---|---|---|---|---|
| T3  | completed | completed | completed | stopped at 12 s | failed at 20 s |
| T9  | completed | completed | completed | stopped at 16 s | failed at 24 s |
| T15 | completed | completed | completed | completed (70 s wall) | failed at 36 s |
| T21 | completed | completed | completed | completed (46 s wall) | failed at 44 s |

BC5 now fails with CVODE `mxstep` (the step size collapses) instead of a negative outlet pressure. BC1–BC3 take 4–13 s
of wall time.

CO₂ t50 in s (relative shift against BC1):

| Run | BC1 | BC2 | BC3 | BC4 |
|---|---|---|---|---|
| T3  | 901.34 | 809.00 (−10.2%) | 901.31 (−3.1e-5) | |
| T9  | 650.12 | 580.41 (−10.7%) | 650.10 (−2.5e-5) | |
| T15 | 634.71 | 492.65 (−22.4%) | 634.70 (−2.7e-5) | 632.59 (−0.33%) |
| T21 | 477.00 | 366.59 (−23.1%) | 476.98 (−3.1e-5) | 474.80 (−0.46%) |

BC2 inlet velocity on T15: 2.7 v_in at 4 s and 2.1 v_in at 20 s. Dilute BEA control unchanged: BC3–BC5 within
1.7e-3 of BC1, BC2 up to 6.4%.

## Results (2026-10-20, `analysis.ipynb`)

Rerun with the marx2 inputs of 2026-10-20: h_fw from the Marx correlation (T3, T15: 54.26; T9, T21: 53.15
W m⁻² K⁻¹) instead of 25.6, and the gas heat capacity per unit mass (`../../marx2/conversions.md` §5).

| Run | BC1 | BC2 | BC3 | BC4 | BC5 |
|---|---|---|---|---|---|
| T3  | completed | completed | completed | stopped after 20 s | failed after 20 s |
| T9  | completed | completed | completed | stopped after 20 s | failed after 20 s |
| T15 | completed | completed | completed | stopped after 20 s | failed after 30 s |
| T21 | completed | completed | completed | stopped after 30 s | failed after 40 s |

CO₂ t50 in s (relative shift against BC1):

| Run | BC1 | BC2 | BC3 |
|---|---|---|---|
| T3  | 904.68 | 816.26 (−9.8%) | 904.61 (−7e-5) |
| T9  | 657.39 | 589.22 (−10.4%) | 657.35 (−7e-5) |
| T15 | 604.59 | 467.15 (−22.7%) | 604.50 (−1.4e-4) |
| T21 | 459.82 | 352.05 (−23.4%) | 459.77 (−1.1e-4) |

Before the rerun (h_fw = 25.6, molar C_p,g): BC2 −9.5%, −10.1%, −22.1%, −22.3%; BC3 within 1.2e-4.

- BC1 and BC3 are the same model here: the pressure drop is 3e-5 of the operating pressure. T21 now
  completes 4000 s with both (marx2 stopped at 1174 s).
- BC2: the inlet velocity follows the uptake (T15: 0.035 m/s in the helium-filled column, 0.077 m/s
  (2.8 v_in) at 4 s, 0.058 m/s at 20 s; before the rerun 0.076 and 0.056 m/s, back to 0.0271 m/s at 3000 s). Every front arrives early: on T15, N₂ at 120 s
  against 200 s and H₂ at 61 s against 120 s. The feed rate is not controlled, so BC2 does not represent
  the experiment (mass-flow controller).
- BC4 and BC5 do not enforce continuity of the molar flux. The gas taken up at the front is not
  replaced: in T15 at 20 s, c_tot drops from 400 to 0.14 mol/m³ (BC5) and to 1.0e-3 mol/m³ (BC4).
  BC5 then stops with "pressure gradient is too large (negative outlet pressure)" (failed state:
  gas temperature −41 K). In BC4 the helium ahead of the front leaves at its Ergun velocity
  (0.035 m/s), faster than the feed arrives, and the gap is evacuated; the run advances at 1–2 s of
  simulated time per minute and is stopped after 15 min (`run.py --timeout 900`). The 15/80/5 runs behave the same, because
  N₂ is also taken up strongly at 1e6 Pa.
- Dilute control (BEA, 7 × 1% in He): BC3, BC4, and BC5 agree with BC1 to 2e-3 in t50 for all seven
  isomers; BC2 shifts t50 by up to 6.4% (`../3.4.1_boundary_conditions_bea_dilute/README.md`).
- Wall time: BC1 and BC3 17–20 min, BC2 31–34 min (8–10 runs in parallel on 10 cores).

## Checks

- [x] Every mode reaches 4000 s. T21 is the critical case (earlier BC3 failure in marx2). BC1–BC3: yes;
      BC4 and BC5: no (above).
- [ ] After complete breakthrough, every mode shows $P(0) = P_\mathrm{in}$, $P(L) = P_\mathrm{out}$,
      $v(0) = v_\mathrm{in}$. Compare against the table in `conversions.md`. (BC2 on T15 ends at
      v(0) = 0.0271 m/s, slightly below v_in, at 3000 s; not yet checked against the 4000 s state.)
- [x] Outlet CO₂ mole fraction of every completed mode against the experiment; t50 per mode (above).
      N₂/H₂ t50 and temperatures per mode: only for T15 so far.
- [x] Velocity profiles across the CO₂ front (paper figure, panel b). BC1–BC3 include the sorption
      effect; BC4 and BC5 do not.
- [x] BC2 and BC4: inlet and helium-region velocity until the column is filled with the feed (above).
- [ ] Global CO₂ balance per mode (fed − eluted − held up). BC4 and BC5 do not conserve the total molar flux.
- [x] Analysis notebook (data-handling functions at the top, one cell per plot, 8×6 PDF figures).
- [x] Paper: map BC1–BC5 to these JSON names in the SI tables (done earlier, see `../../TODO.md`).

## Running

```sh
python3 generate.py --clean    # rewrites simulations/ and deletes old output
python3 run.py --binary /Users/youriran/pusan/ruptura/build/ruptura-o3-20261103   # 20 runs, 8 workers
```

A 4000 s run takes 3–8 s with ROS34PW2. BC2 is the slowest because it solves for the inlet velocity in every
right-hand-side evaluation. Output is about 90 MB per run. With CVODE (before 2026-11-03), BC4 on T3 and T9 and
BC5 on all four runs did not complete; `run.py --timeout 900` stopped the stalling runs.

```sh
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

`analysis.ipynb` reads this directory and `../3.4.1_boundary_conditions_bea_dilute` and writes the four
panels to `../../ruptura-2.0-paper/figures/boundary_conditions/`.
