# 3.4.1 Boundary conditions, dilute control: BEA with seven C₇ isomers in helium

Control case for `../3.4.1_boundary_conditions_marx_10bar`. The Marx runs take up a large share of the
feed at the adsorption front; here seven C₇ isomers at 1% each in helium are taken up, so sorption
changes the molar flow by a few percent at most. The input is the BEA base case of the integrator
benchmark (`../../benchmark/generate.py`), `../../examples/BEA-alkanes-C7/breakthrough/simulation.json`.
It is isothermal at 552 K and 1.2e5 Pa (L = 0.1 m), with dual- or single-site Langmuir isotherms and
IAST (paper SI, table "Dilute case BEA-alkanes-C7").

## Parameterization

The example runs `FixedPressureInletVelocity` with `PressureGradient` = −1000 Pa/m. That operating point
is kept. The other modes get the Ergun velocity of the feed at this gradient (stable root, as in
`../3.4.1_boundary_conditions_marx_10bar/conversions.md`, section 4) and P_out = P_in − G L:

| Quantity | Value |
|---|---|
| Feed molar mass | 0.010736 kg/mol (0.93 He + 7 × 0.01 C₇) |
| Feed density at 552 K, 1.2e5 Pa | 0.2807 kg/m³ |
| G | 1000 Pa/m |
| v_in | 0.098066 m/s (RUPTURA prints 0.0980655 m/s for the example) |
| P_out | 119900 Pa |

Changes against the example: the boundary keys, `ColumnEntranceVelocity`, `BreakthroughIntegrator` ROS34PW2
(default tolerances) and `InletCondition` Danckwerts (since 2026-11-03, as in the Marx case; before, CVODE with
`CVODELinearSolver` SPGMR, since the dense default takes about 10 h here), `SimulationTime` 2000 s, `WriteInterval`
2 s. Without physical axial dispersion, the Danckwerts and the Dirichlet inlet coincide. The example's
`ColumnEntranceVelocity` (0.019 m/s) was never used, because `FixedPressureInletVelocity` takes the velocity from
the gradient.

## Results

2026-11-03, ROS34PW2 (frozen `-O3` binary `build/ruptura-o3-20261103`): all five modes complete 2000 s in 0.9–2.2 s.
Relative to `InletVelocityOutletPressure`, the largest t50 shift over the seven isomers is (CVODE, 2026-11-02, in
brackets):

| Mode | max \|Δt50\|/t50 |
|---|---|
| `InletPressureInletVelocity` | 9.3e-4 (9.4e-4) |
| `FixedPressureInletVelocity` | 1.7e-3 (1.7e-3) |
| `FixedVelocity` | 1.7e-3 (1.7e-3) |
| `InletPressureOutletPressure` | 6.4e-2 (6.4e-2) |

The first three agree because sorption hardly changes the molar flow. `InletVelocityOutletPressure`
evaluates the feed concentration at P_out, 0.08% below P_in, which explains most of the 1e-3 level.
`InletPressureOutletPressure` raises the inlet velocity to hold the pressure drop against the uptake,
so its fronts arrive earlier even in this dilute case.

## Running

```sh
python generate.py --clean    # rewrites simulations/ and deletes old output
python run.py --binary /Users/youriran/pusan/ruptura/build/ruptura-o3-20261103   # 5 runs
```

The analysis is in `../3.4.1_boundary_conditions_marx_10bar/analysis.ipynb`, which reads both run sets.
