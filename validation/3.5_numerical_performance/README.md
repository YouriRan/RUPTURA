# 3.5 Numerical performance

Source: `breakthrough/`, `convection/`, and `figures/` in this directory (moved here from `../../benchmark/` on 2026-11-19, which keeps `cvode-max-step/` and `multibed/`)
- `breakthrough/`: complete breakthroughs with RK3 and CVODE on three cases that differ in what limits the
  explicit step (BEA C7: convection; Marx T15: non-isothermal; Stampi-Bombelli monolith with fitted D_ax:
  dispersion), plus component count and mixture methods. In the paper: Results, "Integrators & benchmark";
  setup in SI "Integrator benchmark".
- `figures/`: CVODE vs RK3 matrix on BEA and Marx (time step, grid, component count) and a mixture-method
  matrix on Marx, over 0.1 s only. Its Marx inputs still use the molar gas heat capacity (`../../TODO.md`).
- `multibed/`: single-bed vs multibed column, correctness and cost.

Status: full-breakthrough benchmark in the paper; SIRK3 later.

## Points

- [x] Work–precision plot (error in t50 vs wall time) per integrator (`benchmark/breakthrough`).
- [x] Main result: wall time is proportional to the number of mixture-prediction calls; CVODE wins when the
      stable RK3 step lies well below Δz/v (its default step cap), and loses when convection limits RK3.
- [ ] Convergence order in $\Delta t$ and $\Delta z$. RK3 converges in Δt on BEA (e50 1.8e-6 → 1e-8) and on
      the monolith. The Δz dependence is in 3.4.3 (grid sweep and numerical dispersion), not as an order.
- [x] Largest stable RK3 time step: BEA Courant 1 (1e-2 s); T15 between 2.5e-2 and 5e-2 s (Δz/v = 0.44 s);
      monolith with fitted D_ax between 1.5e-3 and 2e-3 s, matching the explicit dispersion limit
      0.63 Δz²/D = 1.66e-3 s. This supports the $\Delta t/(\Delta z)^2$ claim in the introduction for
      dispersion-dominated cases.
- [x] Scaling with the number of components and the mixture method (SI Table `tab:benchmark_scaling`).
- [ ] Scaling with grid points and multibed overhead (data in `benchmark/figures`, `benchmark/multibed`; not
      in the paper).
- [x] CVODE tolerances: the default rtol 1e-4 gives e50 = 2.4e-4 on BEA (cap-limited), 4.6e-6 on T15, and
      7.9e-4 on the 400 ppm monolith (≈ atol / feed concentration). Still open: the dilute DAC artefacts of the
      default absolute tolerance in `../3.3.4_humid_dac_liu/liu/README.md`.
- [ ] SIRK3: to be fixed later. It cannot run the benchmark cases and supports one adsorbent only.
      Include it once fixed, or drop it from the paper.
- [ ] `benchmark/figures/run.py` is deleted in the working tree (staged). `benchmark/breakthrough/run.py`
      was rebuilt from it; restore the original if the 0.1 s matrices are to be regenerated.
- [ ] Fair CVODE comparison (`../../TODO.md` V-ITG-3): CVODE without its step cap, with the analytical sorption
      Jacobian as preconditioner and in KLU, the Danckwerts monolith, and repeated timings, in
      `breakthrough/generate.py` and `run.py`; not yet run. The T15 runs also change with the isotherms at T_s (V-PHY-9).
