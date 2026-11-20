# Convection schemes on T15

Status: added 2026-11-20 (TODO V-PHY-12); before, these numbers were ad-hoc runs recorded in `TODO.md` (F-ITG-5).
T15 with the inputs of the model of Marx et al. on 40 cells, 1500 s, with `Upwind`, `TVD`, `WENO3`, and `WENO5`, each
with ROS34PW2 and CVODE at their default settings. The SI (Method of lines) and the main text (boundary conditions and
convection schemes) quote the CO2 front width t90 - t10 and the ROS34PW2 steps.

```sh
python generate.py
python run.py --binary /path/to/frozen/ruptura
python analysis.py
```

Results (2026-11-20): t90 - t10 = 212.42 (Upwind), 200.09 (TVD), 199.23 (WENO3), 198.60 s (WENO5), the same with
ROS34PW2 and CVODE to 0.01 s; ROS34PW2 takes 4546, 5076, 4993, and 7764 steps.
