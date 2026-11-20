# A self-driving lab on a computer: optimizing an H2 PSA cycle

An autonomous lab such as RoboChem (Noël group, Science 2024) closes the loop between a robot that runs experiments, inline
analytics that measure them, and Bayesian optimization that chooses the next experiments, here yield and throughput of
a photochemical reaction. This example does the same on a computer, with the RUPTURA Python interface as the robot. The
technique is in the package `ruptura.selfdriving`; this folder holds the notebook that sets up and analyzes two
campaigns, one on a single column and one on two beds, their lab notebooks (`campaign/experiments.csv`,
`campaign_two_beds/experiments.csv`), and their figures.

| Lab | `ruptura.selfdriving` |
| --- | --- |
| Robot | `PressureSwingCycle`: one PSA cycle with product backfill on the column of a base input, run through `ruptura.load_simulation(...).run()` until the cyclic steady state; `TwoBedPressureSwingCycle`: the cycle on two beds with valves, a product tank, and pressure equalization |
| Inline analytics | `key_performance`: purity, recovery, and productivity of the last cycle, from `SwingAdsorption.phase_records` |
| Brain | `propose_batch`: Gaussian processes, ParEGO scalarization, expected improvement times the probability of feasibility, and batches of 8 by Kriging believer |
| Loop and lab notebook | `Campaign` and `campaign/experiments.csv` |

The separation is H2 purification from the 40/30/30 CO2/N2/H2 feed of experiment T15 of Marx et al. on activated carbon
(`examples/PSA-H2-activated-carbon`), with its validated isotherms (Van 't Hoff in temperature), mass transfer, and
energy balance. CO2 stays in the column; N2 adsorbs weakly and breaks through soon after H2, so that it sets the purity.
The campaign maximizes **productivity** (mol H2 per kg sorbent and second, net of the H2 used for the purge and the
backfill) and **recovery** (net H2 product over H2 fed) under **H2 purity ≥ 99.9 %**, over eight variables.

The cycle has five steps: backfill (H2 product enters at z = L, the feed end closed, and raises the pressure from P_L by
the share b of P_H − P_L), feed pressurization (the feed raises it to P_H from z = 0, the product end closed),
adsorption (feed at z = 0, product out at z = L), blowdown (to P_L through z = 0, the product end closed), and a
counter-current H2 purge at P_L. With b = 0 it is the Skarstrom cycle of `examples/PSA-H2-activated-carbon`.

| Variable | Range |
| --- | --- |
| Temperature of the feed and the surroundings | 288 to 328 K |
| Adsorption pressure | 4 to 10 bar |
| Purge pressure | 1 to 2 bar |
| Feed velocity (interstitial) | 0.01 to 0.06 m/s |
| Adsorption time (= purge time) | 20 to 200 s |
| Purge-to-feed ratio (mol H2 purge per mol feed) | 0.02 to 0.6 |
| Column length, at the sorbent mass of T15 | 0.6 to 3.6 m |
| Backfill share b of the pressure rise | 0 to 1 |

The module description of `ruptura.selfdriving.psa` has the math from these variables to the input, and from the phase
records to the KPIs.

## Two beds

The second campaign runs the cycle on two beds A and B, half a cycle apart (`examples/PSA-H2-activated-carbon-two-beds`,
`BoundaryCondition` `MomentumBalance`). Valves between the product ends and a product tank carry what the one-column
cycle takes from outside: while A adsorbs, its product end feeds the tank, held at P_H, and purges B, through valves
sized from the molar feed flow (the product valve passes it at 2 % above P_H, the purge valve G times it between P_H and
P_L). A pressure equalization from the product end of A to that of B, until B has risen by the share e of P_H − P_L,
takes the place of the backfill; the other seven variables and their ranges are those of the one-column campaign, and
e runs from 0 to 0.5. The product is what reaches the tank, and the productivity is per kg of the sorbent of both beds.
`figures/two_bed_flowsheet.tex` draws the system (`ruptura.flowsheet`).

| Phase (A adsorbing) | Bed A | Bed B | Valves |
| --- | --- | --- | --- |
| Pressurization A, blowdown B | feed end to P_H, linear in t_r | feed end to P_L, linear in t_r | |
| Adsorption A, purge B | feed at v_F and P_H | feed end at P_L | A to tank (check valve), A to B (purge) |
| Equalization A to B | closed | closed | A to B, until B is at P_L + e (P_H − P_L), at most t_r |

The three phases then run with A and B exchanged.

## Running it

Needs `ruptura` with its compiled extension and scipy (`pip install "ruptura[selfdriving]"`), plus pandas and
matplotlib for the figures. Open `analysis.ipynb` in this folder:

```python
from ruptura.selfdriving import Campaign, PressureSwingCycle, TwoBedPressureSwingCycle

cycle = PressureSwingCycle("../../../examples/PSA-H2-activated-carbon/simulation.json", product="H2", purity_target=0.999)
campaign = Campaign(cycle, "campaign", rounds=15, batch_size=8, workers=8, seed=2026)
campaign.run()            # the rounds that campaign/experiments.csv does not have yet
two_beds = TwoBedPressureSwingCycle("../../../examples/PSA-H2-activated-carbon-two-beds/simulation.json", product="H2",
                                    purity_target=0.999)
Campaign(two_beds, "campaign_two_beds", rounds=15, batch_size=8, workers=8, seed=2026).run()
campaign.write_inputs()   # writes the input of every experiment of the lab notebook again
campaign.rerun()          # runs them all again, 8 at a time
```

Each experiment runs in its own Python process (`python -m ruptura.selfdriving <directory>`), 8 at a time: a
simulation keeps the interpreter lock while it runs, and a failed integration stops only its own experiment. It
writes its directory, `campaign/round_<r>/exp_<n>_<setup>/`, with `setup.json`, `robot.json`, `simulation.json`,
`phases.data`, `kpis.json`, and `run.log`. Rerunning a campaign from scratch gives the same setups and results: on
2026-11-19 the one-column campaign gave the same 136 setups and identical KPIs. A campaign that was interrupted and
resumed seeds the generator of its remaining batches anew, so it differs from one run in one go (TODO.md B-PY-8);
both campaigns here ran in one go.

## Model settings and their check

- 80 cells and ROS34PW2 at a relative tolerance of 1e-5. Against 160 cells with CVODE at 1e-6 (checked 2026-11-19,
  with the cyclic-steady-state stop), recovery and productivity of the baseline cycle agree within 0.01 %, of the
  checked Pareto cycles within 0.03 %, and of exp. 130 within 0.1 %. The impurity is that of CVODE on the same 80 cells: 3.2e-4 against 2.5e-4 for the
  baseline, and 5.3e-6 against 9.9e-7 for exp. 130 (b = 0.33); the grid errs on the safe side, by much less than the
  0.1 % that the constraint allows. ROS34PW2 at 1e-4 changes these impurities by less than 0.4 % and is about twice
  as fast, so 1e-5 is conservative. The baseline takes 8 s (16 cycles) in the campaign.
- Pressurization and blowdown change the pressure linearly over t_r = max(10 s, 40 s (L/1.2 m)² (P_H − P_L)/(9 P_L)),
  which keeps the pressure drop, and with it max |R T Σc/P − 1|, below 1 % (`gas_deviation` in the lab notebook).
- A run ends at its cyclic steady state (`CyclicSteadyState` with a tolerance of 2e-3, at most 100 cycles): the first
  cycle whose state, amounts in the column, and stream amounts are within 0.2 % of those of the cycle before. The KPIs
  are those of the stream amounts of that cycle (`phase_records[i].stream_amounts`). `campaign_python_stop/` holds the
  lab notebook of the campaign of 2026-11-11, run with the stop that this replaced, in Python: blocks of 10 cycles
  until the amounts in the column changed by less than 0.2 % of the feed of a cycle and the KPIs by less than 0.2 %.
- The two beds have the same grid, integrator, tolerance, ramps, and stop. A two-bed experiment takes about ten times
  as long as a one-column one: 67 s for the baseline (16 cycles), against 6 s.

## How the brain decides

- One Gaussian process per objective and one for the constraint, on the unit cube of the setup (log scale for the
  velocity, the times, the purge-to-feed ratio, and the length). A log-normal prior keeps the length scales near half
  the cube: in an earlier campaign (97 % purity, seven variables), with a few tens of runs, the likelihood alone
  explained the purity with temperature alone (length scale 0.07, the other variables at the upper bound), and every
  unexplored temperature then looked 30 % feasible. With the prior the purge-to-feed ratio (0.18) and the length (0.47)
  came out as the variables that matter.
- The constraint is the number of nines of the purity, up to 4 (99.99 %), against a target of 3 (99.9 %).
- A run below the target keeps (nines/3)⁴ of its productivity and recovery in the objective models. In the 97 %
  campaign without this, short columns with hardly any purge (recovery 0.9 at purity 0.35) promised so much that an
  expected improvement of 0.3 times a feasibility of 0.1 beat everything else, and whole batches went there. Counting
  such runs as worthless instead kept every batch next to the known feasible cycles. The soft factor keeps the near
  misses, so that the batches go to the purity boundary, where the trade-off between productivity and recovery lies.
- A cycle whose purge and backfill take more H2 than it produces counts as producing none, as does a failed run (see
  the results for what keeping the negative values did).
- ParEGO scalarizes the two objectives with a random weight per proposal; a batch of 8 is filled with the Kriging
  believer, the surrogates taking their own prediction at each chosen setup as measured.

## Results

### One column

The campaign ran 136 experiments in 16 rounds in 6 min on 8 workers (21 min of simulation; 8.7 s per experiment in the
median, 36 s at most). `analysis.ipynb` makes the figures in `figures/` from `campaign/experiments.csv`.

| Cycle | Productivity, mmol/(kg s) | Recovery | H2 purity | T, P_H/P_L, v_F, t_ads, G, L, b |
| --- | --- | --- | --- | --- |
| Hand-tuned baseline (round 0) | 0.584 | 0.150 | 0.99968 | 298 K, 10/1.00 bar, 0.030 m/s, 30 s, 0.27, 1.20 m, 0 |
| Most productive (exp. 48, round 5) | 1.032 | 0.301 | 0.99922 | 288 K, 8.03/1.03 bar, 0.033 m/s, 60 s, 0.22, 1.16 m, 0.01 |
| Pareto front (exp. 133, round 15) | 1.000 | 0.380 | 0.99965 | 288 K, 9.94/1.35 bar, 0.026 m/s, 105 s, 0.16, 1.39 m, 0.09 |
| Pareto front (exp. 90, round 10) | 0.963 | 0.416 | 0.99932 | 288 K, 10/1.31 bar, 0.029 m/s, 126 s, 0.16, 1.60 m, 0.09 |
| Highest recovery (exp. 100, round 11) | 0.489 | 0.446 | 0.99920 | 288 K, 9.01/1.34 bar, 0.060 m/s, 118 s, 0.16, 3.60 m, 0 |

- At 99.9 % H2 the campaign raised the recovery 2.0 to 3.0 times over the baseline, and the productivity up to 1.8
  times; the hypervolume of the feasible experiments over (0, 0) grew from 0.088 to 0.440, five times that of the
  baseline (`figures/hypervolume.pdf`).
- The optimum is a different cycle, not a tuned baseline: the lowest temperature of the box (288 K, where 41 of the 88
  feasible experiments lie), a purge pressure of 1.0 to 1.35 bar, two to four times longer steps (60 to 126 s), less
  purge (G = 0.16 to 0.22 instead of 0.27), and a small backfill (b up to 0.09, 4 to 5 % of the feed in H2). The
  front runs from short columns with short steps (most productive) to the longest column at the fastest feed
  (highest recovery).
- The backfill keeps N2 out of the product in the middle of the front: without it (b = 0, all else equal), exp. 90
  gives 2.9 % N2 instead of 0.07 %, exp. 133 1.5 % instead of 0.035 %, and exp. 112 (b = 0.035) 0.28 % instead of
  0.07 %. The two ends of the front need almost none (exp. 48: b = 0.009, 0.4 % of the feed; exp. 100: b = 0).
- The fine model confirms the front: on 160 cells with CVODE at 1e-6 the impurity of the most productive cycle is
  6.2e-4 (7.8e-4 on 80 cells), of exp. 90 3.1e-4 (6.8e-4), and of the one with the highest recovery 5.9e-4 (8.0e-4);
  recovery and productivity change by less than 0.03 %. The grid errs on the safe side.
- Purity alone is easy at the cost of recovery: 9 of the 16 start setups, and 26 of the 78 proposals that met the
  target, put more H2 into the column with purge and backfill than the cycle produces, which flushes the column with
  product (open symbols in `figures/impurity.pdf`); 17 of the 28 such proposals had b > 0.4. 52 of the 120 proposals
  met the target with a net product, with 2.6 of the 8 variables at a bound on average.
- One run did not reach the cyclic steady state within 100 cycles (purity 0.85); the others took 13 cycles in the
  median.
- The campaign of 2026-11-11 (`campaign_python_stop/`), with the stop in Python, ended at nearly the same hypervolume
  (0.442) on a different front: 295 K, 10/1.00 bar, 124 to 159 s steps, G = 0.12 to 0.13, b = 0.09 to 0.25 (most
  productive 1.017 mmol/(kg s) at a recovery of 0.411). Its stop ran blocks of 10 cycles, so that no run ended before
  10 cycles and 12 runs did not reach its criterion within 100. In a second campaign of that day, the cycles without a
  net product kept their negative objectives (down to minus the best value so far); the batches then went further to
  the bounds of the box, with 4.3 of the 8 variables at a bound on average instead of 2.1, and the hypervolume ended at
  0.33. The expected improvement grows with the uncertainty, which is largest at the corners of the cube; a trust
  region around the incumbents, as in TuRBO, would be the next step.

### Two beds

The campaign ran 136 experiments in 16 rounds, in one go, in 57 min on 8 workers (242 min of simulation; 90 s per
experiment in the median, 360 s at most). `analysis.ipynb` makes `figures/pareto_front_two_beds.pdf`,
`figures/hypervolume_two_beds.pdf`, and `figures/pareto_setups_two_beds.pdf` from `campaign_two_beds/experiments.csv`.

| Cycle | Productivity, mmol/(kg s) | Recovery | H2 purity | T, P_H/P_L, v_F, t_ads, G, L, e |
| --- | --- | --- | --- | --- |
| Hand-tuned baseline (round 0) | 0.562 | 0.144 | 0.99956 | 298 K, 10/1.00 bar, 0.030 m/s, 30 s, 0.27, 1.20 m, 0 |
| Most productive (exp. 133, round 15) | 1.278 | 0.513 | 0.99915 | 288 K, 10/1.00 bar, 0.037 m/s, 41 s, 0.16, 1.17 m, 0.43 |
| Pareto front (exp. 119, round 13) | 1.070 | 0.520 | 0.99915 | 288 K, 10/1.00 bar, 0.036 m/s, 41 s, 0.16, 1.16 m, 0.50 |
| Pareto front (exp. 121, round 14) | 1.037 | 0.556 | 0.99912 | 288 K, 10/1.00 bar, 0.038 m/s, 41 s, 0.13, 1.25 m, 0.50 |
| Highest recovery (exp. 131, round 15) | 0.206 | 0.592 | 0.99999 | 295 K, 6.60/1.00 bar, 0.060 m/s, 77 s, 0.10, 3.38 m, 0.50 |

- Every cycle of the one-column front is dominated by a two-bed cycle: at 1.04 to 1.28 mmol/(kg s) the two beds
  recover 0.51 to 0.56 of the H2, where the most productive single column recovers 0.30. The hypervolume grew from
  0.081 to 0.708, 8.7 times that of the baseline, and was still growing in the last round (0.668 to 0.708).
- The three most productive cycles share 288 K, 10/1 bar, 41 s steps, G = 0.13 to 0.17, a column of 1.2 m, and a
  nearly complete equalization (e = 0.43 to 0.5, the upper bound); the highest recovery comes from a 3.4 m column at
  the fastest feed.
- The equalization does the work of the backfill and more: without it (e = 0, all else equal), exp. 133 recovers
  0.35 instead of 0.51 and gives 0.76 % N2 instead of 0.085 %; exp. 119 0.35 and 0.89 %, exp. 121 0.38 and 1.6 %.
- None of the 67 proposals that met the target lacked a net product (on one column 26 of 78 did): the equalization,
  unlike the backfill, takes no product. 6 proposals made no net product, all below the target; 2.8 of the 8
  variables were at a bound on average.
- The fine model confirms the front: on 160 cells with ROS34PW2 at 1e-6 (CVODE is slow with `MomentumBalance` and the
  energy balance, TODO.md B-ITG-21), the impurities of exp. 133, 121, and 131 are 7.7e-4, 8.0e-4, and 3.1e-6 (80
  cells: 8.5e-4, 8.8e-4, and 9.8e-6), with recovery and productivity within 0.07 %.
- 10 runs did not reach the cyclic steady state within 100 cycles; two of them were near the target (purity 0.9992
  and 0.9987). The others took 14 cycles in the median.
- A first run of this campaign, resumed after round 0 (TODO.md B-PY-8), searched differently and ended at a
  hypervolume of 0.564, with a front up to 1.016 mmol/(kg s) at a recovery of 0.385 and up to a recovery of 0.62 at
  0.14 mmol/(kg s): two searches of 136 experiments in eight dimensions do not converge to one front.
