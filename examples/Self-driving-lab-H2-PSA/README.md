# A self-driving lab on a computer: optimizing an H2 PSA cycle

An autonomous lab such as RoboChem (Noël group, Science 2024) closes the loop between a robot that runs experiments, inline
analytics that measure them, and Bayesian optimization that chooses the next experiments, here yield and throughput of
a photochemical reaction. This example does the same on a computer, with the RUPTURA Python interface as the robot. The
technique is in the package `ruptura.selfdriving`; this folder holds the notebook that sets up and analyzes one
campaign, its lab notebook (`campaign/experiments.csv`), and its figures.

| Lab | `ruptura.selfdriving` |
| --- | --- |
| Robot | `PressureSwingCycle`: one PSA cycle with product backfill on the column of a base input, run through `ruptura.load_simulation(...).run()` until the cyclic steady state |
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

## Running it

Needs `ruptura` with its compiled extension and scipy (`pip install "ruptura[selfdriving]"`), plus pandas and
matplotlib for the figures. Open `analysis.ipynb` in this folder:

```python
from ruptura.selfdriving import Campaign, PressureSwingCycle

cycle = PressureSwingCycle("../PSA-H2-activated-carbon/simulation.json", product="H2", purity_target=0.999)
campaign = Campaign(cycle, "campaign", rounds=15, batch_size=8, workers=8, seed=2026)
campaign.run()            # the rounds that campaign/experiments.csv does not have yet
campaign.write_inputs()   # writes the input of every experiment of the lab notebook again
campaign.rerun()          # runs them all again, 8 at a time
```

Each experiment runs in its own Python process (`python -m ruptura.selfdriving <directory>`), 8 at a time: a
simulation keeps the interpreter lock while it runs, and a failed integration stops only its own experiment. It
writes its directory, `campaign/round_<r>/exp_<n>_<setup>/`, with `setup.json`, `robot.json`, `simulation.json`,
`phases.data`, `kpis.json`, and `run.log`. Rerunning the campaign from scratch gives the same setups and results.

## Model settings and their check

- 80 cells and ROS34PW2 at a relative tolerance of 1e-5. Against 160 cells with CVODE at 1e-6, recovery and
  productivity agree within 0.2 %. The impurity is that of CVODE on the same 80 cells: 3.1e-4 against 2.5e-4 for the
  baseline cycle, and 2.9e-6 against 4.8e-7 for a cycle with 30 % backfill; the grid errs on the safe side, by much less
  than the 0.1 % that the constraint allows. ROS34PW2 at 1e-5 takes 13 s for the baseline (20 cycles), CVODE at 1e-4 on
  the same grid 34 s. At 1e-4, ROS34PW2 is twice as fast again, but its impurity of the backfill cycle is 2.3e-5.
- Pressurization and blowdown change the pressure linearly over t_r = max(10 s, 40 s (L/1.2 m)² (P_H − P_L)/(9 P_L)),
  which keeps the pressure drop, and with it max |R T Σc/P − 1|, below 1 % (`gas_deviation` in the lab notebook).
- A run goes on in blocks of 10 cycles, at most 100, until the column holds the same amount of every component from
  one cycle to the next within 0.2 % of the feed of a cycle, and the KPIs within 0.2 % (the purity within 0.2 % of its
  0.1 % margin).

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

The campaign ran 136 experiments in 16 rounds in 11 min on 8 workers (40 min of simulation; 12.5 s per experiment in the
median, 85 s at most). `analysis.ipynb` makes the figures in `figures/` from `campaign/experiments.csv`.

| Cycle | Productivity, mmol/(kg s) | Recovery | H2 purity | T, P_H/P_L, v_F, t_ads, G, L, b |
| --- | --- | --- | --- | --- |
| Hand-tuned baseline (round 0) | 0.584 | 0.150 | 0.99970 | 298 K, 10/1.00 bar, 0.030 m/s, 30 s, 0.27, 1.20 m, 0 |
| Most productive (exp. 128, round 15) | 1.017 | 0.411 | 0.99950 | 295 K, 10/1.00 bar, 0.018 m/s, 124 s, 0.13, 1.01 m, 0.14 |
| Pareto front (exp. 127, round 14) | 0.899 | 0.428 | 0.99963 | 295 K, 10/1.00 bar, 0.018 m/s, 159 s, 0.12, 1.05 m, 0.25 |
| Highest recovery (exp. 125, round 14) | 0.837 | 0.438 | 0.99959 | 295 K, 10/1.00 bar, 0.018 m/s, 154 s, 0.12, 1.30 m, 0.09 |

- At 99.9 % H2 the campaign raised the recovery 2.7 to 2.9 times and the productivity 1.4 to 1.7 times over the
  baseline; the hypervolume of the feasible experiments over (0, 0) grew from 0.088 to 0.442, five times that of the
  baseline (`figures/hypervolume.pdf`).
- The baseline already makes 99.97 % H2, so the campaign raises recovery and productivity at a purity that is reachable
  from the start. Purity alone is easy at the cost of recovery: 9 of the 15 space-filling setups, and 33 of the 88
  proposals that met the target, put more H2 into the column with purge and backfill than the cycle produces, which
  flushes the column with product (open symbols in `figures/impurity.pdf`). 55 of the 120 proposals met the target
  with a net product.
- The optimum is a different cycle, not a tuned baseline: the full pressure ratio at 295 K, a slower feed (0.018 m/s
  instead of 0.030) for four to five times longer steps (124 to 159 s instead of 30 s), half the purge (G = 0.12 to
  0.13 instead of 0.27), and a small backfill (b = 0.09 to 0.25, 5 to 12 % of the feed in H2). The long steps spread
  the gas lost in blowdown over more feed; the backfill pushes the N2 front back to the feed end before the next feed
  pressurization, which saves purge. Without its backfill the most productive cycle gives 6.9 % N2 in the product
  instead of 0.05 %.
- The fine model confirms the front: on 160 cells with CVODE at 1e-6 the most productive cycle gives an impurity of
  2.9e-4 (5.0e-4 on 80 cells), and the one with the highest recovery 2.4e-4 (4.1e-4); their recovery and productivity
  change by less than 0.02 %.
- The search took its time: the hypervolume stayed at 0.188 from round 2 to round 6, while 30 % of the proposals after
  round 0 were cycles whose purge and backfill take more H2 than they produce (most with b > 0.4). Such cycles give a
  pure product, and with their objectives counted as zero they look like cycles that break even. Keeping their negative
  values instead (down to minus the best value so far) made it worse in a second campaign: the batches went further to
  the bounds of the box, with 4.3 of the 8 variables at a bound on average instead of 2.1, and the hypervolume ended at
  0.33. The expected improvement grows with the uncertainty, which is largest at the corners of the cube; a trust region
  around the incumbents, as in TuRBO, would be the next step.
- 12 runs did not reach the cyclic steady state within 100 cycles; all of them are below 99.4 % H2.
- The front sits at the bounds of the pressures (10 and 1 bar) only. Compared with the 97 % target of an earlier
  campaign (2.1 mmol/(kg s) at a recovery of 0.35, with 25 s steps at 0.058 m/s), the 30-fold lower impurity costs
  half the productivity, and the slow cycles that it needs also raise the recovery.
