# Verification on T15: free gas pressure and ends of phases

Status: added 2026-11-20 (TODO V-PHY-12), when T15 moved to the inputs of the model of Marx et al. The numbers of the
verification table (`tab:verification`) and of the SI sections "Free gas pressure" and "Ends of phases" come from
these runs; before, they were ad-hoc runs recorded only in `TODO.md` (B-ITG-21, F-BC-3, the stop criteria).

Base: `../3.3.1_energy_balance_marx_ternary/marx2/ternary/T15/simulation.json` on 30 cells.

| Directory | What |
| --- | --- |
| `free_pressure/momentum_balance/{isothermal,energy_balance}` | one 1500 s phase, `MomentumBalance`, velocity port at z = 0, pressure port (1e6 Pa) at z = L, ROS34PW2 |
| `free_pressure/continuity_closure/{isothermal,energy_balance}` | the same with `InletVelocityOutletPressure` |
| `free_pressure/momentum_balance/isothermal_cvode_klu` | the isothermal free-pressure run with CVODE and KLU, the one CVODE setting the reader accepts with `MomentumBalance` |
| `phase_end/{CVODE,ROS34PW2,RungeKutta3}` | the T15 breakthrough (energy balance) as one phase that ends at an outlet CO2 mole fraction of 0.2 |

```sh
python generate.py
python run.py --binary /path/to/frozen/ruptura
python analysis.py
```

Results (2026-11-20, frozen binary of the paper reruns):

- free gas pressure against the closure: t50 of CO2, N2, H2 650.20, 238.09, 131.14 s isothermal and 545.07/545.06,
  207.56, 119.73 s with the energy balance; outlet mole fractions within 8e-6 and 3e-5; peak gas and adsorbent
  temperatures 329.60 and 332.11 K for both; pressure drop between the first and last interior nodes 3.03 and 3.06 Pa
  for both. The z = 0 row of `column.data` holds the port pressure with `MomentumBalance` and the boundary value with
  the closure, so the drop is taken between interior nodes.
- ROS34PW2 steps: 7871 (`MomentumBalance`, isothermal) and 16762 (energy balance), against 4215 and 4711 for the closure.
- CVODE KLU crawls at the CO2 front: 100 000 steps advance 0.18 s at 714 s, and the run stops.
- Phase end: CVODE 545.0620 s, ROS34PW2 545.0602 s, RK3 545.0603 s; spread 1.85 ms.
