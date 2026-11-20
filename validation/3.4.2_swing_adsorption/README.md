# 3.4.2 Pressure- and temperature-swing cycles

Status: TSA runs done, figure in the paper (Results, "Cycle performance"; setup in SI
"Temperature-swing cycles"). PSA deferred.

## Case

Two-step TSA with the existing phase driver (`SwingAdsorptionPhases`) on the validated Marx T15 column
(`../../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json`, energy balance on, `InletPressureInletVelocity`):

| Step | T_m (feed and surroundings) | P_in | Feed |
|---|---|---|---|
| Adsorption | 298 K | 1e6 Pa | T15, CO₂/N₂/H₂ 40/30/30, 0.0275 m/s |
| Regeneration | 373 K | 1e6 Pa | same (the driver changes only T and P_in) |

Equal step durations t½ = 1000, 2000, 4000, 8000, 16000 s; every run covers at least 60000 s and five
cycles (30, 15, 8, 5, 5 cycles). Output every 20 s; ROS34PW2 ends a step at every output. The first adsorption step starts
from the helium-filled column, so the first cycle begins with the T15 breakthrough. Two MixturePrediction
runs give the IAST CO₂ loading of the feed: 3.6955 mol/kg at 298 K and 1.1024 mol/kg at 373 K, an
equilibrium working capacity of 2.593 mol/kg.
The unit conversions of the T15 inputs (flow rate to `ColumnEntranceVelocity`, pressures,
temperatures) are in `../3.4.1_boundary_conditions_marx_10bar/conversions.md`.

## Why this case

- The BEA examples (`../../examples/BEA-alkanes-C7/psa`, `tsa`) cannot be used: their isotherms have
  no `HeatOfAdsorption`, so the TSA "regeneration" desorbs nothing. Its outlet peak at the switch
  (normalized partial pressure 1.25) is the gas inventory reacting to the instantaneous 552 → 673 K
  step, not desorption.
- An isothermal TSA aborts. With the energy balance off, the phase temperature is imposed on the whole
  column at once. The 373 → 298 K step raises the equilibrium loading everywhere, the uptake exceeds the
  feed supply, and the outlet velocity becomes negative (−0.63 m/s); RUPTURA stops with "pressure
  gradient is too large (negative outlet pressure)". A PSA repressurization step would do the same.
  Both need flow reversal (`../../TODO.md`).
- With the energy balance on, the phase temperature sets only the inlet-gas and external temperatures
  (`Column::setTemperature`). The temperature front travels through the column, and the uptake stays
  limited by the feed.

## Results (2026-11-20, the inputs of the model of Marx et al.)

Rerun with `marx2/ternary` (TODO V-PHY-12; the equilibrium runs use SCI): equilibrium working capacity 3.364 − 1.049 =
2.315 mol/kg; working capacity 0.214, 0.489, 0.981, 1.655, 2.217 mol/kg (9–96 %) for t_1/2 = 1000–16000 s;
productivity 1.07, 1.22, 1.23, 1.03, 0.69e-4 mol/(kg s); adsorbent swing 7.4–70.1 K (28.8 K at 4000 s); cycles to a
residual below 1 %: 14, 7, 4, 3, 2; CO2 balance within 3.1e-4. `gas_inventory.py`: c_tot R T_g / P − 1 within
−2.7e-4 … +5e-6, gas share 7–15 %, inventory change ≤ 1.7e-6. Earlier results below are for the reference settings.

## Results (`analysis.ipynb`)

Rerun on 2026-11-19 (V-ITG-2) with ROS34PW2 after its stage projection on the closure (B-ITG-16), from a frozen
build of the working tree (`run.py --rerun --binary`, log `run_2026-11-19_projection.log`). ROS34PW2 and the
Danckwerts inlet are written in the marx2 T15 input that the runs copy; the inputs are those of 2026-11-03.

| t½ / s | Complete cycles | Δq_cyc / mol kg⁻¹ | Fraction of IAST | Productivity / mol kg⁻¹ s⁻¹ | r_k, last | Adsorbent T swing / K | CO₂ balance | Wall time per cycle / s |
|---|---|---|---|---|---|---|---|---|
| 1000 | 30 | 0.207 | 0.080 | 1.03e-4 | 6.3e-6 | 7.4 | 3.2e-4 | 0.59 |
| 2000 | 15 | 0.460 | 0.177 | 1.15e-4 | 7.9e-6 | 14.9 | 3.7e-4 | 1.08 |
| 4000 | 8 | 0.938 | 0.362 | 1.17e-4 | 5.6e-6 | 28.9 | 3.6e-4 | 2.01 |
| 8000 | 5 | 1.702 | 0.656 | 1.06e-4 | 6.1e-7 | 50.5 | 3.2e-4 | 3.52 |
| 16000 | 5 | 2.424 | 0.935 | 7.6e-5 | 1.9e-7 | 69.6 | 1.8e-4 | 5.58 |

The 2026-11-04 runs (before the projection) gave the same Δq_cyc, temperature swings, and cycles to 1e-2 to three
digits, r_k of 3e-5 to 5e-5 in the last cycle, CO₂ balances of 7.5e-5 to 1.7e-4, and 2.6 to 21.7 s per cycle (the
1 to 12 s of the paper are from 2026-11-03).

The outputs now in `simulations/` are a second rerun, made at 19:58 with the B-ITG-17 fix (ROS34PW2 no longer
starts a step from the clamped state; log `run_2026-11-19_clampfix.log`). It gives every column of the table to the
digits shown, except the wall times: 1.05, 1.89, 3.51, 6.15, and 8.60 s per cycle. Those runs shared the machine
with five CVODE Dense benchmark runs (load about 12 on 10 cores). On the first four 1000 s cycles, run alone, both
binaries take 4.9 s of CPU time (8600 and 8674 steps) and give the same CO₂ balance (1.207e-3; 6.174e-4 at rtol
1e-8), so the clamp is not the floor of the balance. The table keeps the wall times of the first rerun, made at a
lower load (about 6); time them again on a quiet machine.

- r_k (CO₂ accumulated per cycle relative to Δq_cyc) falls below 1e-2 after 14, 8, 4, 3, and 2 cycles. In the
  last cycle it is 2e-7 to 8e-6; with CVODE (rtol 1e-5) it levelled off near 8e-4.
- The column wall dominates the heat capacity (wall about 3000 J m⁻¹ K⁻¹ against about 150 J m⁻¹ K⁻¹ for
  the bed) and is coupled to the surroundings through h_w∞ = 5 W m⁻² K⁻¹. Hence the small temperature swing for
  short steps and the productivity maximum at 4000 s.
- State transfer: the CO₂ balance (fed − eluted − accumulated, trapezoidal over the 20 s output) closes to 1.8e-4
  to 3.7e-4 of the fed CO₂. With the Danckwerts inlet the fed amount is the feed flux v c_0; `co2_mass_balance`
  adds the dispersive term only for a Dirichlet inlet. Check on the first four 1000 s cycles (2026-11-19): the
  closure is 1.2e-3, 0.87e-3, and 1.1e-3 with outputs every 20, 5, and 2 s, so it is not the trapezoidal rule. It
  falls with the ROS34PW2 tolerance, 7.7e-4 at rtol 1e-7 and 6.2e-4 at 1e-8, towards the 6.8e-4 of the 2026-11-04
  binary: the stage projection gives the integration error of Σc to the components by their mole fractions.
  Regeneration phases lose CO₂ and adsorption phases gain it, 1e-4 to 6e-4 of the fed amount per phase. The floor
  near 6e-4 was there before the projection, and it is not the clamp of B-ITG-17 either (unchanged by that fix); it
  is not traced.
- Gas inventory (SI remark; `python gas_inventory.py`): the `InletPressureInletVelocity` closure carries the thermal
  expansion of the gas (B-BC-2), and ROS34PW2 projects its stages on it (B-ITG-16). During the swing, c_tot R T_g
  deviates from the Ergun pressure by −9.7e-5 to +4.0e-5 (2026-11-04, before the projection: −2.9 % to +4.0 %). The
  equilibrium uses y_j P_Ergun and is not affected. The gas holds 6–14 % of the CO₂ in the column; putting the gas
  at the Ergun pressure changes the CO₂ inventory at the step ends by at most 2.8e-6 (before: 4e-4).

## Running

```sh
python generate.py            # rewrites simulations/ (use --clean to delete old output)
python run.py                 # 7 runs, 8 workers, 16–28 s each (2026-11-19)
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb
```

`analysis.ipynb` writes the four panels to `../../ruptura-2.0-paper/figures/cycles/`.

## Points for later

- [x] Check the state transfer: the column inventory must be continuous across step boundaries, and the
      cumulative mass balance must close over the whole sequence.
- [x] Repeat TSA cycles in the phase list and assess cyclic steady state after the run; report working
      capacity and productivity.
- [ ] Purity and recovery: with the feed as regeneration gas they are not meaningful. They need per-step
      feed composition (`../../TODO.md`, cycle driver).
- [ ] PSA. In the PSA example, "Blowdown" is a co-current feed at low pressure, not a counter-current
      depressurization. Realistic cycles need flow reversal, per-step feed changes, and a treatment of the
      gas inventory at a pressure step (`../../TODO.md`).
