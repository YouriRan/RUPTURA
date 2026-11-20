# 3.4.1 Boundary conditions with CVODE: the stalls at the near-vacuum fronts

Status: added 2026-11-20 (TODO V-PHY-12). Replaces the CVODE copies in
`../inlet_condition_dirichlet_vs_danckwerts/simulations/*/3.4.1_boundary_conditions_marx_10bar`, which hold the
reference settings of before 2026-11-19.

The `FixedVelocity` and `FixedPressureInletVelocity` runs of `../3.4.1_boundary_conditions_marx_10bar` (T3, T9, T15,
T21) with CVODE at its default settings instead of ROS34PW2. Each run gets 900 s of wall time; the paper quotes how
far the runs get and how their total gas concentration compares with the ROS34PW2 runs up to there (Results,
integrator performance; SI, boundary-condition study).

```sh
python generate.py                                   # after ../3.4.1_boundary_conditions_marx_10bar/generate.py
python run.py --binary /path/to/frozen/ruptura       # eight runs at a time, up to 15 min
python analysis.py
```

The wall-time limit makes the result depend on the load of the machine: run it alone.
