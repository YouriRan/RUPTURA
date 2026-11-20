# 3.4.3 Parameter sensitivity

Status: runs done, figure in the paper (Results, "Sensitivity studies"; setup in SI "Sensitivity study").

One-at-a-time sensitivity of the Marx et al. (2014) ternary experiment T15 (CO₂/N₂/H₂ = 40/30/30,
298 K, 1e6 Pa, 1e-5 m³/s), with the validated input `../../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json`. T15 is
also the Marx case of the integrator benchmark (`../../benchmark/`, case `marx-403030_298_1e6_2e-1`)
and of the boundary-condition study (`../3.4.1_boundary_conditions_marx_10bar`).
The unit conversions of the T15 inputs (flow rate to `ColumnEntranceVelocity`, pressures,
temperatures) are in `../3.4.1_boundary_conditions_marx_10bar/conversions.md`.

## Runs (38)

| Kind | Inputs | Levels |
|---|---|---|
| reference | unchanged T15 | – |
| local | q_sat, k, ν, Q_ads of CO₂; ρ_s; k_LDF,CO₂; D_ax (all components); h_fs, h_fw, h_w∞; C_p,s, C_p,w | ×0.9, ×1.1 |
| sweep | k_LDF,CO₂ | ×0.25, ×0.5, ×2, ×4 |
| sweep | D_ax, every component | 0, 1e-5 (benchmark levels), 1e-4, 1e-3 m²/s; the benchmark level 1e-6 is left out |
| sweep | N_z | 10, 20, 50, 200 (benchmark levels) |
| check | N_z = 200 with D_ax = 8.7e-5 m²/s, every component | numerical dispersion |

The dispersion and grid levels are read from `../../benchmark/generate.py`, so both studies stay on
the same levels. Changes common to every run: 2000 s simulated, output every 2 s. ROS34PW2 (from the T15
input) adapts its internal step and ends a step at every output.

Outcome measures at the outlet: t50 of CO₂, front spreading t90 − t10, and the largest outlet gas
temperature rise ΔT_max. Normalized coefficients S = log(m(1.1p)/m(0.9p)) / log(1.1/0.9).

## Results (2026-11-20, the inputs of the model of Marx et al.)

Rerun with `marx2/ternary` (TODO V-PHY-12; the CO2 Sips parameters are law objects whose `Value` is scaled). Reference
t50 541.3 s, t90 − t10 193.9 s, ΔT_max 23.7 K. S(t50): q_sat 0.86, k 0.51, ν 0.87, ρ_s 0.88, Q_ads −0.18 (heat only),
k_LDF 0.005; S(spread): q_sat 0.87, Q_ads 0.73, ρ_s 0.72, h_fw −0.72; |S| of h_w,inf < 0.03. k_LDF ×0.25…×4: t50
−2.0 … +0.4 %; D_ax 1e-3: t50 +0.5 %, spread +27 %; N = 10: spread +25 %, N = 200: t50 −0.15 %. Earlier results below
are for the reference settings.

## Results

Rerun on 2026-11-03 with ROS34PW2 and the Danckwerts inlet, both written in the marx2 T15 input that every run
copies (ρ_s = 456.6 kg/m³, h_fw from the Marx et al. correlation, mass-based C_p,g, WENO3, 100 grid points).
Against the CVODE runs with the Dirichlet inlet of 2026-10-31, only the dispersion sweep changes visibly (below);
everything else moves in the third decimal.

Reference: t50 = 634.7 s (experiment 571 s), t90 − t10 = 146.0 s, ΔT_max = 17.63 K.

| Input | S(t50) | S(t90 − t10) | S(ΔT_max) |
|---|---|---|---|
| q_sat,CO₂ | 0.926 | 0.740 | 0.748 |
| b_CO₂ | 0.504 | 0.172 | 0.511 |
| ν_CO₂ | 0.777 | −1.173 | 1.535 |
| Q_ads,CO₂ | 0.087 | 1.511 | 0.630 |
| k_LDF,CO₂ | 0.008 | −0.206 | 0.191 |
| ρ_s | 0.897 | 0.681 | 0.272 |
| h_fs | 0.001 | −0.005 | 0.085 |
| h_fw | 0.064 | −0.725 | −0.348 |
| h_w∞ | 0.002 | 0.025 | −0.004 |
| C_p,s | −0.003 | −0.000 | −0.297 |
| C_p,w | 0.040 | 0.116 | −0.064 |
| D_ax | 0.000 | 0.002 | −0.002 |

- t50 is set by the equilibrium and the bed density (q_sat, ρ_s, ν, b). A larger heat of adsorption raises the
  affinity at 298 K but also heats the bed; the two nearly cancel in t50 (S = 0.09), while the spreading follows
  it more than any other input (S = 1.51).
- The front spreading follows the heat of adsorption, the gas–wall heat transfer, and the isotherm shape, not
  k_LDF or D_ax.
- k_LDF,CO₂ over a factor 16 (×0.25 to ×4): t50 −2.4% to +0.6%, spreading +63% to −15%.
- D_ax from 0 to 1e-4 m²/s: t50 within 0.04%, spreading within 4.0%. D_ax = 1e-3 m²/s: t50 −0.25%, spreading
  +47%. With the Dirichlet inlet of the 2026-10-31 runs, 1e-3 m²/s advanced t50 by 2.9%: dispersion across the
  inlet face added to the feed (`../inlet_condition_dirichlet_vs_danckwerts`).
- Numerical dispersion: first-order upwinding would add about vΔz/2 = 1.6e-4 m²/s at N_z = 100; WENO3 removes
  most of it. Check run: N_z = 200 with D_ax = 8.7e-5 m²/s broadens the front by 2.8% and differs from the
  100-point reference by up to 0.039 in y/y_in, so it does not reproduce the 100-point curve.
- Grid: N_z = 100 is within 0.11% (t50), 1.1% (spreading), and 2.0% (ΔT_max) of N_z = 200. N_z = 10 widens the
  front by 47%.
- Panel (a) leaves out h_w∞ (|S| < 0.03) so that the labels stay legible at print size.
- ρ_s = 456.6 kg/m³ comes from the marx2 inputs; with S(t50) ≈ 0.9, this choice matters nearly one to one for the
  breakthrough time.

## Running

```sh
python generate.py            # rewrites simulations/ (use --clean to delete old output)
python run.py                 # 38 runs, 8 workers, about 20 s each
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

`analysis.ipynb` writes the four panels to `../../ruptura-2.0-paper/figures/sensitivity/`.

## Open points

- [ ] The Liu et al. (2024) Figs. 2–3 sweeps (`../../liu/simulations/`, all finished) are a second,
      published sensitivity case. They belong to 3.3.4 and are not in this subsection.
- [ ] Drive the sweeps through the Python API, so that this subsection also demonstrates the API.
