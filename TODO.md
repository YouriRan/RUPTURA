# TODO

Open work for RUPTURA 2.0 and its paper (`ruptura-2.0-paper/main.tex`). Validation cases are organized
per results subsection in `validation/`; each subsection README lists its own points.

## Labels

Every item has a label `<type>-<area>-<n>`, e.g. `B-GUI-1` (a bug in the Qt app) or `P-ITG-1` (a performance issue
in the integrators). Labels are permanent: completed items stay as checked boxes with their label, and a new item
takes the next free number of its type and area. Each section lists its open items first and then the completed
ones, both ordered by area as in the second table.

| Type | Section |
| --- | --- |
| `C` | Critical issue: silently wrong results for a whole class of runs |
| `B` | Bug |
| `F` | Feature |
| `I` | Paper inconsistency: the code and the paper or SI disagree |
| `P` | Performance |
| `R` | Refactoring and cleanup |
| `V` | Validation and reruns |
| `W` | Paper writing (no code involved) |

| Area | Scope |
| --- | --- |
| `PHY` | Physics models: isotherms, mixture prediction (IAST, SPI, SCI, MPD), energy balance, chemisorption, reactions, pore transport |
| `ITG` | Integrators and spatial discretization: RK3, SIRK3, CVODE, convection scheme |
| `BC` | Boundary conditions and velocity closure |
| `CYC` | Swing-adsorption phases and cycles |
| `MBD` | Multibed (multistacked) columns |
| `FIT` | Isotherm fitting |
| `IO` | Input reader, output files, checkpoints |
| `GUI` | Qt app (`app/`) |
| `PY` | Python package |
| `BLD` | Build, CI, and unit tests |
| `EX` | Example, validation, and benchmark inputs and their figures |
| `GEN` | Cross-cutting: the whole code base or the whole paper |

The 2026-10-26 code review (read-through of the C++ core, the fitting code, the input reader, the Qt app, and the
Python package; no build, no tests run) contributed C-BLD-1, B-PHY-5–7, B-ITG-4–6, B-ITG-10, B-MBD-1–2, B-FIT-1–4,
B-IO-1–4, B-GUI-1, B-PY-1–2, P-PHY-3–5, P-ITG-5–7, P-MBD-1, P-PY-1, P-EX-1, and every `R` item except R-BLD-1–2.
Its bugs were re-checked against the source; the open items it listed only as likely bugs were not: B-PHY-7,
B-ITG-6, B-MBD-1–2, B-FIT-4, B-IO-2–4, B-GUI-1, B-PY-2. Suggested order: C-BLD-1, then the physics items (B-PHY-5,
B-PHY-6, B-IO-1), then SIRK3 by deleting its duplicate (B-ITG-4), then the fitting items (B-FIT-1–4), then the
CVODE second IAST pass (P-ITG-5), then the `R` items.

## Critical issues (C)

- [x] **C-PHY-1 The film/pore/chemisorption coupling lost mass** (fixed 2026-10-21).
  - The film and pore fluxes conserve mass only if both reservoirs have the solid volume (1−ε). The pore
    sink used γ = ρ_p(1−ε)/ε, so the pores lost (1−ε)/ε times what the solid gained.
  - Fix: γ = ρ_p in `transport.cpp`. In the multibed column, γ is the solid-volume-weighted ρ_p
    (`compute_multibed.cpp`). Both inventory limiters were changed the same way.
  - Check on a copy of `examples/Chemisorption-CO2/breakthrough`: CO₂ closure is now 1.0001 (it was
    0.750).
  - New test: `ChemisorptionTransport.FilmPoreAndChemisorptionExchangeConservesMass`.
  - Remaining: ε_p (`ParticlePorosity`) is not modelled; the code equals the SI model at ε_p = 1. See the
    SI item I-PHY-4.
- [x] **C-PHY-2 Gas–wall heat exchange does not conserve energy.** The gas balance (per unit fluid volume) loses
      h_fw (4/d_i)(T_g − T_w), i.e. ε_b·π d_i h_fw ΔT per unit length, while the wall gains π d_i h_fw ΔT
      (`src/integrators/compute.cpp`, `computeEnergyDerivatives`; the code comment already asks whether a 1/ε is
      missing). With h_fw = h_L, the bed loses ε_b times the heat of the Marx et al. model.
      Fixed 2026-10-26 in the code, which now follows the SI (Eqs. `eq:general_fluid_energy_balance` and
      `eq:wall_energy_balance` share H_fw = 4h_fw/d_in per bed volume): `fluidWallPerFluidVolume` = 4/(ε_b d_in) for
      a packed bed (`src/geometry.cpp`). Monoliths already conserved energy. Test
      `EnergyBalance.GasWallExchangeConservesEnergyInPackedBed`; `marx2/conversions.md` §4 updated. Reruns
      pending (V-GEN-2).
- [x] **C-BLD-1 `-ffast-math` compiled every NaN/Inf guard out** (fixed 2026-10-27).
  - `CMakeLists.txt:262` applied `-ffast-math` to all sources but two; with it, `std::isfinite(x)` became
    `return true` and `std::isnan(x)` `return false`. This disabled `finiteNonnegative`/`safeMultiply` in
    `mixture_prediction.cpp`, `finiteNonNegative`/`limitedAdsorbedRate` in `compute.cpp`, the checks in
    `isotherm.{h,cpp}` and `reaction.cpp`, and all `isfinite` input checks in `inputreader.cpp`.
  - Fix: `-fno-finite-math-only` after `-ffast-math` in `RUPTURA_COMPILE_OPTIONS`, so it covers the library,
    CLI, Python module, Qt app, tests, and benchmarks. The other fast-math optimizations stay. `std::isfinite`
    now compiles to a sign-bit mask and one integer compare (Clang) or `fabs` and one compare (GCC).
  - Hand-written bit tests (`std::bit_cast` and the exponent mask) are no alternative: Apple Clang 17 and
    Clang 20 fold them to `true` under `-ffast-math` too. Only GCC 15 keeps them.
  - Check (Apple Clang 17, Release, `-march=native`): the guards are back in the object files of
    `mixture_prediction`, `isotherm`, `reaction`, `inputreader`, and `compute`. Six shortened examples (BEA C7
    CVODE, RK3, and PSA; MOR SIAST; Chemisorption CO₂; Marx T15 with energy balance) give byte-identical outputs,
    and their runtimes are within run-to-run noise (median ratio 0.91–1.05 over three alternating runs). Unit
    test results are unchanged (on HEAD, 152 of 181 pass in both builds; the failures are B-BLD-1 and B-PHY-3).
- [x] **C-EX-1 Gas heat capacity units.** The gas energy balance divides by `gasDensity` (Σ c_i M_i, kg m⁻³) times
      `heatCapacityGas`, which the inputs and the SI give per mole (29.51 J mol⁻¹ K⁻¹). The gas volumetric heat
      capacity is then ρ_g C_p instead of c_tot C_p, about 1/M̄ ≈ 38 times too small for the Marx feeds. Decision
      (2026-10-20): keep the code and the SI unit (J kg⁻¹ K⁻¹, as Stampi-Bombelli, Elfving, Dantas, and the
      examples use) and fix the Marx inputs: `marx2/heat_capacity.py` sets 29.51/M̄ (1014.1 and 1109.2 J kg⁻¹ K⁻¹;
      `marx2/conversions.md` §5). marx2 T1–T24, 3.4.1, 3.4.2, 3.4.3, and the benchmark rerun. The small gas heat
      capacity had made the gas energy balance stiff: RK3 on T15 failed at Δt = 1e-3 s, and now runs at 3e-2 s.
      The inputs that are still molar are B-EX-6.

## Bugs (B)

- [ ] **B-PHY-1** MPD temperature reweighting adds N·ln(T_ref/T) (`src/macrostate_particle_distribution.cpp:200`),
      but the Langmuir affinity scaling (`src/utils.h:35`) has no such term. Away from T_ref, MPD and
      SCI/IAST therefore differ by a factor T_ref/T in the affinity. Decide which is intended.
- [ ] **B-PHY-2** The reaction key `gibbsFreeEnergy` sets both K(T) = K₀·exp(−ΔG/RT) and the reaction heat (−ΔG·r),
      so it acts as a van 't Hoff reaction enthalpy (`src/reaction.cpp:89`, `src/integrators/compute.cpp:345`).
      Rename it or document this. See also I-PHY-10.
- [ ] **B-PHY-3** Avrami rate: commit `fd613c7` (2026-10-05) dropped the factor n, `k^n t^(n−1) (q* − q)`, and kept
      the old line commented out (`src/chemisorption.cpp:38`). The integrated model, the SI
      (Eq. `eq:avrami_chemisorption`), and `ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel` have
      `n k^n t^(n−1)`. No input uses Avrami. Decide which form is intended.
- [ ] **B-PHY-4** BET isotherm (`src/isotherm.h:151`): q = q_m b p / ((1 − c p)(1 − c + b p)). Textbook BET with
      b = C/p₀ and c = 1/p₀ has (1 − c p + b p) as the second factor; `1 − c` mixes a pure number with 1/Pa.
      The spreading pressure (`src/isotherm.h:270`) matches the code's form, so IAST is self-consistent, but
      the model is not BET. Unchanged since the initial commit; check against RUPTURA 1 and fix or document.
      The SI form is I-PHY-5.
- [ ] **B-PHY-5 Anti-Langmuir IAST inverse has the wrong sign in the exponent.** `isotherm.h:433-437` uses
      `1 − exp(+bψ/a)`; inverting `psiForPressure` (`:268`) gives `1 − exp(−bψ/a)`. Returns negative
      hypothetical pressures.
- [ ] **B-PHY-6 O'Brien–Myers `value()` and `psiForPressure()` are inconsistent with each other**, not only with the
      SI (see I-PHY-6): `isotherm.h:193` has σ²·x(1−x)/(1+x)³ and `:348` ½σ²·x/(1+x)²; the exact
      integral of the implemented `value()` has no ½. Fix both in the same direction.
- [ ] **B-PHY-7** `PseudoNth` chemisorption: `pow(driving, order)` with `driving < 0` is NaN, masked only by the
      argument order of `std::min(0.0, NaN)`; it can never desorb.
- [ ] **B-ITG-1** SIRK3: fix later. It cannot run the benchmark cases and supports a single adsorbent only. See also
      the banded-matrix layout and the missing temperature/chemisorption update (B-ITG-4) and P-ITG-6.
- [ ] **B-ITG-2** RK3 close to its stability limit is inaccurate on the dispersion-dominated monolith (fitted D_ax,
      400 ppm): e50 = 0.098 at Δt = 1.25e-3 s (limit ≈ 1.66e-3 s), 2.8e-4 at 4e-4 s. Possibly the clipping of
      negative concentrations at the front (`updateStateRK` max(0, ·); `stampi-bombelli/generate_runs.py`
      mentions it); not checked. Matters for users who pick the largest stable step.
- [ ] **B-ITG-4 SIRK3 with the energy balance builds the banded LAPACK matrix in the wrong layout.** `rk3_si.cpp:421`
      uses `std::layout_left` (column-major `ab`); its twin `computeConcentrationUpdateMatrixEnergyBalanceFinal`
      (`:522`) uses CTAD, which gives `layout_right` (`git show b8a2562` fixed one copy only). SIRK3 also never
      advances the temperature or chemisorption blocks (`propagate` updates `physisorption` and `concentration`
      only), but the reader rejects only reactions and pore transport (B-ITG-9). Reject `EnergyBalance` and
      chemisorption with SIRK3, and merge the EB/non-EB matrix builders (they differ only in the unused
      `facePressures`).
- [ ] **B-BC-1** Boundary modes. `FixedPressureInletVelocity` imposes a linear pressure profile and takes the
      velocity from the local Ergun balance without continuity. `FixedVelocity` ignores the sorption
      effect. 3.4.1 results (2026-10-19):
  - Both fail on all four Marx runs at 1e6 Pa (T3, T9, T15, T21). The gas taken up at the front is not
    replaced, and c_tot drops from 400 to 1e-3 mol/m³. `FixedVelocity` stops after 20–40 s with
    "pressure gradient is too large (negative outlet pressure)" (gas temperature −41 K in the failed
    T15 state). `FixedPressureInletVelocity` evacuates the gap between the feed front and the faster
    helium (0.035 against 0.0275 m/s) and stalls: 20–32 s simulated in 10–17 min.
  - In the dilute BEA case (7 × 1% in He) both agree with `InletVelocityOutletPressure` to 2e-3 in t50.
  - Either guard them (abort with a clear message when c_tot falls far below the feed value) or state
    in the SI that they are valid only when sorption hardly changes the molar flow.
  - `InletPressureOutletPressure` does not control the feed rate: on T15 its inlet velocity reaches
    2.8 × v_in and CO₂ breaks through 22% early; even in the dilute BEA case t50 shifts by up to 6%.
- [ ] **B-BC-2** Velocity closure without thermal expansion: `InletPressureInletVelocity` and
      `InletVelocityOutletPressure` conserve the molar flux with c_tot from the transported
      concentrations, with no ∂c_tot/∂t term. With the energy balance, c_tot R T_g departs from the Ergun
      pressure: −14% to +23% in the 3.4.2 TSA runs. The equilibrium (y_j P_Ergun) is not affected, but the
      gas inventory is.
- [ ] **B-CYC-1** Phase driver: a step that raises the equilibrium loading in the whole column at once aborts. An
      isothermal TSA (energy balance off, so T_m is imposed everywhere) stops at the 373 → 298 K step: the
      uptake exceeds the feed, v(L) = −0.63 m/s, "negative outlet pressure" (3.4.2). A PSA
      repressurization step would do the same. Needs flow reversal (F-BC-1) or a one-way outlet.
- [ ] **B-MBD-1** Multibed gas–wall area comes from bed 0 only (`compute_multibed.cpp:845`); the film and gas–solid
      areas blend 1/ε_a per bed while the loading prefactor uses ε_total (`:682-687`, `:890` vs `:656-664`).
- [ ] **B-MBD-2** Single-bed vs multibed drift: multibed `computeChemisorptionTransportDerivatives`
      (`compute_multibed.cpp:629`) writes `poreDot = −ρ·dq/dt` for every site, the single bed
      (`transport.cpp:41`) skips non-transport sites; `cvode.cpp:592` passes bed 0's components to
      `computeReactionDerivatives` for all beds.
- [ ] **B-FIT-1 GA uniform crossover is a no-op.** `fitting.cpp:522`: `children[i].genotype.substr(j, 1) = …` assigns
      to a temporary.
- [ ] **B-FIT-2 Nelder–Mead shrink leaves n − 2 function values stale.** `fitting.cpp:931-953` moves every vertex but
      re-evaluates only `f[vg]` and `f[vh]`; correct only for two-parameter isotherms.
- [ ] **B-FIT-3 Fitting column-index defaults underflow.** `inputreader.h:163` defaults `columnPressure{0}`;
      `fitting.cpp:37` subtracts 1 → `SIZE_MAX` → "no pressure points found" unless `ColumnPressure` is set.
- [ ] **B-FIT-4** `Isotherm::randomize()` (`isotherm.cpp:266`) never draws b < 1 (`pow(U, 20(U′−1)) ≥ 1`), so every
      initial GA citizen is saturated at Pa-scale pressures. The NM stop test (`fitting.cpp:965`) checks the mean
      RMSE < 1e-4 rather than the simplex size, so it always runs `MAX_IT`.
- [ ] **B-IO-1 `wallDensity` defaults to 0 and is not validated** (`inputreader.h:125`, `:1245`), so
      `1 / (heatCapacityWall · wallDensity)` (`compute.cpp:777`, `compute_multibed.cpp:877`) is `inf` → NaN
      temperatures for any energy-balance run that omits it. (The second half, `PrintEvery`/`WriteEvery`
      not checked > 0, went with the time-based output on 2026-10-27: `WriteInterval` must be > 0.)
- [ ] **B-IO-2** `size_t` inputs read with `get<size_t>()` wrap silently (`"NumberOfGridPoints": -1` → `SIZE_MAX`);
      `readOptionalNonNegativeInteger` exists but is used for a few keys only (four copies of the parsing).
- [ ] **B-IO-3** Mole-fraction normalisation divides by zero when no component sets `GasPhaseMolFraction`
      (`inputreader.cpp:1569-1581`); the swing-phase branch guards it, the global one does not.
- [ ] **B-IO-4** `catch (const std::runtime_error&)` in `breakthrough.cpp:193` and `swing_adsorption.cpp:298`: other
      exceptions skip the `failed_state.json` dump.
- [ ] **B-GUI-1** App: `app/runner.cpp:478` `stopJupyter` calls `deleteLater()` through a pointer the synchronous
      `finished` handler has just nulled.
- [ ] **B-PY-1** Python: `ColumnDataExplorer` joins the relative `data_dir` twice (`plot_breakthrough.py:991` and
      `:1037`); `tox.ini` lacks `pandas`/`plotly`, so the test suite cannot import.
- [ ] **B-PY-2** Python: `fit_figure()` never draws the fit (trace commented out, `plot_fitting.py:467`);
      locale-dependent `read_text()` on non-ASCII presets and LICENSE (`config.py:221`, `:236`, `about.py:68`);
      `RunManager.start` reads `self.cli` twice outside the lock.
- [ ] **B-PY-3** `tests/test_breakthrough.py`, `test_components.py`, `test_fitting.py`, and `test_mixtureprediction.py`
      import `Components` and other names that `ruptura/__init__.py` no longer exports (also on HEAD `633f2bd`), so
      pytest cannot collect them. Found 2026-10-27.
- [ ] **B-BLD-1** Test inputs that exist only locally, so CI checkouts lack them: `examples/Multibed/*`,
      `examples/MPD/BEA-alkanes/*`, `examples/Liquid/*`, `examples/Chemisorption-CO2*/*`,
      `examples/BEA-alkanes-C7/{psa,tsa}`, `examples/BEA-alkanes-C7/fitting/Results.dat-BEA-Repeat-552K-nC7`.
- [ ] **B-BLD-2** `tests/test_*.py` are v1 tests (`from ruptura import Components`); `.gitignore` (`test_*`) hides
      them from git, but a local `pytest tests` still collects and fails them.
- [ ] **B-EX-1** `examples/BEA-alkanes-C7/tsa` and `psa`: no `HeatOfAdsorption`, so the isotherms do not depend on
      temperature and the TSA regeneration desorbs nothing. The outlet peak at the switch is the gas
      inventory reacting to the instantaneous 552 → 673 K step.
- [ ] **B-EX-2** `benchmark/figures/run.py` was deleted in commit `fd613c7` (2026-10-05). Restore it
      (`git show fd613c7^:benchmark/figures/run.py`) or replace it.
- [ ] **B-EX-3** `examples/MPD/BEA-alkanes/breakthrough/change.py` rewrote all 24 benchmark inputs to N = 50,
      dt = 5e-4 (and added a misspelled `NumberOfTimesteps` key), so
      `ExampleInputReader.LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix` fails on duplicate combinations.
      Regenerate the matrix (`generate_simulations.py`) or adapt the test if the N = 50 change was intended
      (2026-10-26: not changed, waiting for a decision).
      2026-10-27: `change.py` and `generate_simulations.py` write the time keys, and the 24 inputs were converted
      (`_legacy_simulation.json`); the misspelled `NumberOfTimesteps` (ignored by the reader) is gone. The matrix
      still has 9 distinct combinations.
- [ ] **B-EX-4 `examples/Reactions/physisorbed-power-law` is unstable.** RK3 at its `TimeStep` (1e-3 s) reaches a gas
      temperature of 8.5e48 K at t = 1 s (HEAD 633f2bd, unmodified input). Loadings above 100 mol/kg follow, and
      the gas temperature stays above 1e35 K up to 3000 s. Check the reaction-heat coupling under RK3, a smaller
      step, or CVODE.
- [ ] **B-EX-5 Marx gas–wall heat-transfer coefficient (2026-10-19).** All 24 marx2 runs use h_fw = 25.6 W m⁻² K⁻¹,
      the reference value of T21 only. Hassan & Chung (`marx/hassan_data/Cases_*/*/my_BreakLab_params.json`,
      `ColBed[11]`) use run-dependent values: T3 28.95, T15 34.5, T21 25.6, T13 3.45, T5 68.9, T17 64.22, T11 65.5,
      T23 61.0. With these, CO₂ t50 against the experiment: T15 −8.5% → −2.5%, T17 −22% → −5%, T13 +13% → −3%.
      This explains the error pattern over T1–T24 (40/30/30 early, worse at high flow and pressure; 1 bar runs
      late) and the low-flow discrepancy that the paper attributes to the flow controller.
  - [x] 2026-10-20: all 24 runs now use Marx et al. Eq. 3 (Nu = 41.13 Re₀^0.32 exp(−6d_p/d_i), Re₀ with the
    initial helium density and the superficial feed velocity, K_L = 0.0400 W m⁻¹ K⁻¹ from their Fig. 6), via
    `marx2/heat_transfer.py`; derivation in `marx2/conversions.md`. h_fw = 25.4–77.1 W m⁻² K⁻¹. marx2 and
    3.4.3 rerun; SI tables and Energy-modeling text updated.
  - [x] 3.4.1 (boundary conditions) and 3.4.2 (cycles) regenerated from marx2 and rerun (2026-10-20), with the
        correlation h_fw and the corrected gas heat capacity.
- [ ] **B-EX-6** Gas heat capacity still molar (29.51 or 35.8; see C-EX-1): `marx2/binary/B1–B36` and
      `marx2/base.json`, `benchmark/marx-403030_298_1e6_2e-1` (0.1 s matrices), the old `marx/` runs, and
      `pusan/`. None is in the paper; fix them before any reuse.
- [x] **B-PHY-8** SPI and SCI evaluate the Toth isotherm without the temperature factor `scale`
      (`src/mixture_prediction.cpp`: `safePureSiteLoading` and `computeSegregatedCompetitiveIsotherm`), while
      `Isotherm::value` applies it to b. With two or more adsorbing components and a non-isothermal Toth site,
      SPI/SCI lose the temperature dependence. A single adsorbate takes the `Isotherm::value` shortcut and is not
      affected (e.g. the Stampi-Bombelli monolith). Found 2026-10-23. Fixed 2026-10-26; test
      `MixturePrediction.SPIAndSCIPassTemperatureScaleToToth`. The Liu runs (SPI, Toth) have no
      `HeatOfAdsorption` on CO₂ and are not affected.
- [x] **B-PHY-9** Unilan: `temp1` and `temp2` are identical, so q* = 0 always (`src/isotherm.h`, `value`). Found
      while building SI Tables S10–S12. Fixed 2026-10-26 (e^{−σ} in the denominator, and the same in the
      spreading pressure, which was zero too); test `Isotherm.UnilanLoadingAndSpreadingPressureFollowTheUnilanForm`.
- [x] **B-ITG-7** CVODE right-hand side: uptake and velocity are coupled by two fixed-point passes that start from
      zero loading rates (`src/integrators/cvode.cpp:372`, `:500`). In Marx T15 this leaves t50 4.9e-5 (relative)
      from the converged RK3 solution at every tolerance. With eight passes (test build, not committed) CVODE
      agrees with RK3 to 2.4e-6, at 2.3× the wall time (with the molar gas heat capacity: 6e-8 at ~20×).
      Decision (2026-10-20): keep the two passes; the difference is stated in the paper (§3.7 and SI
      "Integrator benchmark").
- [x] **B-ITG-8** CVODE BDF order at the Δz/v step cap: on BEA (convection-limited, Courant ≈ 1) the t50 error is
      not monotonic in rtol (1e-3: 5.4e-5, 1e-4: 2.4e-4, 1e-5: 2.2e-6). All runs take one Newton iteration per step
      at the cap; the tolerance changes the Krylov accuracy and the BDF order (last order 2, 3, 5, 5, 4 for rtol
      1e-2 … 1e-6). With `CVodeSetMaxOrd(3)` (test build, not committed) the error is monotonic: 5.4e-5, 2.9e-6,
      6.0e-7 for rtol 1e-3, 1e-4, 1e-5. BDF4/5 are not A-stable, and the upwind convection modes at Courant ≈ 1 lie
      close to the imaginary axis, where they are poorly damped. Consider a `CVODEMaximumOrder` input (default 5,
      or 3 when the step is capped at Δz/v).
  - 2026-10-26: `CVODEMaximumOrder` (1–5, default 5 = CVODE's own default, so results are unchanged) is read,
    passed to `CVodeSetMaxOrd`, and listed in the GUI, the AI prompt doc, and SI Table S12. Test:
    `ExampleInputReader.ParsesAndValidatesCvodeMaximumOrder`. Still open: whether the default should be 3
    when the step is capped.
- [x] **B-ITG-9** SIRK3 (`rk3_si.cpp`) evaluates no reactions and no pore transport at all and does not say so. Only
      `BulkConcentration` reactions are rejected; reject every reaction and every pore-transport site, or
      implement them. 2026-10-26: the input reader rejects every reaction and every surface-pore transport site
      with SIRK3; test `ExampleInputReader.Sirk3RejectsReactionsAndSurfacePoreTransport`.
- [x] **B-ITG-10 CVODE non-negativity: clamp vs `CVodeSetConstraints`.** The review first listed the
      `clampNonnegative(column.state)` after each `CVode()` (`cvode.cpp:77`) as a mutation of CVODE's state.
      It is not: CVODE keeps its history in its own Nordsieck array and only *writes* the interpolated solution
      into the user vector, which it never reads back; the right-hand side binds to the trial vector CVODE
      passes in. The clamp is output hygiene only (one `std::max` per state entry, unmeasurable).
  - Checked 2026-10-26 with `CVodeSetConstraints` (1.0 for concentrations and loadings, 0.0 for temperatures;
    test build, not committed) on marx2 T15 against the clamp:
    - `TimeStep` 1 s, `WriteEvery` 1, 4000 s: 10.35 → 11.97 s wall (+16%), 9244 → 9763 internal steps,
      18,875 → 25,831 Krylov iterations, one constraint-triggered step retry; the first 400 s (steep fronts)
      1.38 → 2.87 s (+108%).
    - The committed T15 input (`TimeStep` 5e-4, 8,000,000 steps): 919.6 → 909.3 s, because 99% of that run is
      the per-output-step right-hand-side evaluation, not the integration (same 9248 internal steps as the 1 s
      run, which finishes in 10 s; see P-EX-1).
    - All output fields agree within ≤ 8e-4 relative (rtol 1e-4); the constrained run still writes small
      negative interpolated values (≤ 3.6e-5 mol m⁻³), since constraints hold at step ends only.
    - Dropping both clamp and constraints is also safe (mole fractions, partial pressures, reaction activities
      and pH already clip negatives in the RHS, and the kinetics use q* − q), but gains nothing.
    - Decision: keep the clamp, do not use constraints.
- [x] **B-ITG-3 Reaction auto-stop depends on `TimeStep`.** `reactionAutoStopReached` tests |ẋ|·`TimeStep` against
      1e-8 + 1e-5·max(1, |x|). On `examples/Reactions/pore-langmuir-hinshelwood` (CVODE, `"auto"`, measured
      2026-10-26) it first fires at 4 s with `TimeStep` 5e-4, 7 s with 2e-3 (the example's value), 226 s with 1e-2,
      939 s with 0.1, and 2475 s with 1 s. The bulk phase settles at 500–900 s. A rate-based criterion
      (|ẋ|/max(1, |x|) per second) removes the dependence. In that run the reaction-source terms decay as fast as
      the state derivatives. A catalytic flow reactor at steady state keeps a nonzero source, though, so including
      the sources can stop the auto-stop from ever firing.
      Fixed 2026-10-27 (time-based I/O): the criterion is |ẋ| ≤ 1e-7 + 1e-4·max(1, |x|) per second on the phase
      states (and the bulk for bulk reactions), without the reaction sources, checked at every output in
      `Breakthrough::advanceToNextOutput`. Tests `ReactionsAutoStop.*`.
- [x] **B-ITG-5 Infinite loop when the breakthrough is converged at step 0.** `cvode.cpp:64` and `rk3.cpp:423` set
      `numberOfSteps = size_t(1.1 · step)`; at step 0 this is 0 and `step >= numberOfSteps − 1` compares against
      `SIZE_MAX`. Reachable from a `ReadColumnFile` of a broken-through column with `AutoNumberOfTimeSteps`.
      The reaction branch already guards with `max(…, 2)`.
      Fixed 2026-10-27 (time-based I/O): the automatic stop runs after an output in
      `Breakthrough::advanceToNextOutput`; the end time is 1.1·t rounded up to whole outputs, so at least one more
      output follows.
- [x] **B-ITG-6** RK3 passes `step · timeStep` as the elapsed time to the chemisorption kinetics (`rk3.cpp:427`),
      which is wrong when `NumberOfInitTimeSteps > 0` ramps the time (`breakthrough.cpp:163-172`). Only Avrami
      (B-PHY-3) reads it. CVODE passes `t` correctly.
      Fixed 2026-10-27 (time-based I/O): RK3 and SIRK3 pass the start time of each step (`RungeKutta3::advance`),
      also during the ramp.
- [x] **B-FIT-5** `Isotherm::isUnphysical` reads `parameters[2]` for `Freundlich`, which has two parameters
      (`src/isotherm.cpp:196`): an out-of-bounds read whenever the fitting cost
      (`src/multi_site_isotherm.cpp:180`) evaluates a Freundlich site. Check only `parameters[0]` and `[1]`.
      Fixed 2026-10-26; test `Isotherm.FreundlichPhysicalityChecksOnlyItsTwoParameters`.
- [x] **B-IO-5 Output times were one `TimeStep` early.** Output time stamps lag the state by one `TimeStep`:
      `writeOutput` runs after `propagate` with the pre-step time (`src/breakthrough.cpp:125`, `:183`; same in
      `src/swing_adsorption.cpp`). The final state is never written (a SwingAdsorption run with `WriteEvery` = 200
      ends its output 20 s early).
  - `Breakthrough::run` and `SwingAdsorption::run` take `stepTime = realTime` before `propagate` and write
    the propagated state with that label (`src/breakthrough.cpp`, `src/swing_adsorption.cpp`). The output
    labelled t is the state at t + Δt, and the initial state is never written.
  - Negligible for small output steps, but 1 s with the 1 s output steps of the Stadler cycles; the
    Stadler notebook adds one `TimeStep` when reading.
  - Fix: write with `realTime` after the step, and remove the notebook's offset at the same time.
  - Fixed 2026-10-26: each snapshot is written before the step with the time of the state it holds, starting
    with the initial state at t = 0, and the final state is written after the last step (both run loops). The
    derivative columns of the t = 0 snapshot are those of `initialize()`. Checked on Marx T15: outputs at 0,
    0.005, 0.01, and the final 0.0125 s. Outputs of runs before this date are still one `TimeStep` late; the
    Stadler notebook offset goes with the rerun (V-GEN-2).
  - When comparing runs with different `TimeStep`, shift the times by `TimeStep` *and* sample at the same
    times: CVODE with `TimeStep` 1 s, `WriteEvery` 2 records odd seconds, RK3 even seconds, and linear t50
    interpolation then differs by up to 3e-4 on the sharp H₂ front of Marx T15 (`benchmark/breakthrough`).
- [x] **B-BLD-3** 18 unit tests fail on a clean checkout (2026-10-21), mostly stale chemisorption input formats
      (`"Type": "General"` as an isotherm type, a `"Parameters"` object). Migrate them to the current
      format: `example_inputreader_tests.cpp`, `fluid_phase_tests.cpp`.
  - 2026-10-26: 163 of 165 pass. The chemisorption fixtures use the current `Type`/`Parameters` +
    `RateEquation`/`RateEquationParameters` format. The two CVODE regression tests reproduce their stored
    values bit for bit with the tolerances they were recorded with (rtol 1e-3, atol 1e-6, the defaults on
    2026-08-20), so the tests now set those tolerances; only the default tolerances had moved them. The BEA
    snapshot tests read `tests/data/bea_breakthrough_column.json` (100 s of RK3 on the BEA example, derivatives
    recomputed at the saved state); their hard-coded values were re-baselined from it.
  - Still failing: `ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel` (B-PHY-3) and
    `ExampleInputReader.LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix` (B-EX-3).
- [x] **B-BLD-4** Unit tests write fixed file names into `std::filesystem::temp_directory_path()` (e.g.
      `ruptura_sirk3_integrator.json`), so two `ruptura_tests` processes running at the same time (two build trees
      or two sessions) break each other with empty-input parse errors (seen 2026-10-23 in
      `ReactionsInput.RejectsLHOutsidePoreConcentration` and
      `MultibedInterfaces.RejectsIncompleteOrInvalidSectionGridCounts`). Use unique names, e.g. with the process
      id. Run alone, the current working tree has 18 of 159 tests failing. Fixed 2026-10-26:
      `tests/temporary_directory.h` gives each process its own directory (removed at exit); two concurrent
      `ruptura_tests` runs both pass.
- [x] **B-BLD-5** 19 of 138 C++ unit tests fail on every platform (macOS Apple Clang, conda Clang, GCC 16,
      Clang/libc++ 18–22):
  - 11 tests throw `unknown isotherm type 'General'` (also `Avrami`, `FirstOrder`, `Elovich`) while reading
    chemisorption inputs: the input reader does not accept the kinetic model names the tests use;
  - `ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel`: rate 0.0303 instead of 0.0442;
  - `GasBreakthroughRegression.ExistingColumnCvodePathRetainsItsOutcome` and
    `LiquidPhase.TransportedPHAndSuperficialFluxWorkInMultibedColumn`: results moved away from the stored values
    (1e-9 and 1e-10 tolerances);
  - `ExampleMassBalance.RecomputesBeaColumnStateDerivatives` and
    `ExampleMixturePrediction.RecomputesBeaBreakthroughEquilibriumLoadings` need
    `examples/BEA-alkanes-C7/breakthrough/column.json`, which is not in the repository;
  - `MixturePrediction.BeaMultinomialMPDMatchesCompetitiveMultisiteLangmuir` needs
    `examples/MPD/BEA-alkanes/mpd/particle_distribution.data`, which `.gitignore` (`**/*.data`) keeps out;
  - `ExampleInputReader.LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix` finds duplicate runs in the local
    `examples/MPD/BEA-alkanes/breakthrough` directories.
  - 2026-10-26: all fixed except Avrami (B-PHY-3) and the MPD matrix (B-EX-3). The BEA `column.json` is now
    `tests/data/bea_breakthrough_column.json` (in the repository). `particle_distribution.data` is still ignored
    (part of B-BLD-1).
- [x] **B-BLD-6** `SingleAdsorbentComputeEquivalence.{MassBalanceRoutinesAreIdentical,TwoIdenticalAdsorbentsMatchSingleBedCalculations}/{2,3,6,7}`
      fail with Clang 21 (Ubuntu 26.04 image, conda-forge clang 21 on macOS; Clang 18, 19, 22, 23 and GCC 16 pass):
      values of ±9e-13 against an absolute tolerance of 1e-12, which is too tight under `-ffast-math`.
      2026-10-26: the relative tolerance now scales with the largest entry of each field (near-zero derivatives
      come from cancelling terms of that size). Passes with Apple Clang; not yet rerun with Clang 21.
- [x] **B-BLD-7** `find_package(LAPACK)` is optional (`CMakeLists.txt:23`), but `src/integrators/rk3_si.cpp` calls
      `dgtsv_` and `dgbsv_`: without LAPACK, configure succeeds and linking fails. Make it `REQUIRED`. Done
      2026-10-26.
- [x] **B-BLD-8** The SUNDIALS settings before `FetchContent_MakeAvailable` (`SUNDIALS_BUILD_CVODE`,
      `SUNDIALS_BUILD_NVEC_SERIAL`, `EXAMPLES_ENABLE`, `ENABLE_TESTS`) are not SUNDIALS options, so a plain
      `cmake -B build` builds all SUNDIALS packages and C examples, and `cmake --install` installs SUNDIALS too. The
      presets set the real options (`BUILD_ARKODE`, …, `EXAMPLES_ENABLE_C`, `EXAMPLES_INSTALL`). 2026-10-26:
      `CMakeLists.txt` sets the real options (only CVODE, no examples, no example install). SUNDIALS libraries are
      still installed, which a shared build needs at run time.
- [x] **B-BLD-9** With `BUILD_CLI` and `BUILD_PYTHON` both on, the executable `build/ruptura` collided with the build
      directory of the Python module (`build/ruptura/`), so linking failed; this affected `pip install .`. The
      module now builds in `build/python/`.
- [x] **B-EX-7** Marx T6 and T18 used k_LDF,CO₂ = 0.07 and 0.15; Marx et al. (Table 6) give 0.04 and 0.11 for these
      298 K runs. Corrected in `marx2/ternary/T6` and `T18/simulation.json` and in the SI table. Rerun: V-PHY-6.
- [x] **B-EX-8** `ParticleDensity` = 502 kg m⁻³ is the bulk density in Marx et al. (see
      `validation/3.3.1_energy_balance_marx_ternary`). Marx et al. (Table 3) write the balance with ε_t, the
      superficial velocity, and ρ_b. With v = Q/(A ε_t), as in marx2, matching it needs ρ_p = ρ_b ε_b/((1−ε_b) ε_t)
      = 457 kg m⁻³; 502 gives 10% more adsorbent per unit gas volume (later breakthrough).
      2026-10-26: T1–T24 use 456.6 (`marx2/particle_density.py`, after `heat_capacity.py`; derivation in
      `marx2/conversions.md` §6). Rerun and SI tables pending (V-GEN-2).

## Features (F)

F-PHY-1, F-PHY-2, F-BC-1, and F-CYC-1 are described in the paper but not implemented. They are to be implemented in
due course. Until then, the paper text must not claim them (I-PHY-2, I-CYC-1).

- [ ] **F-PHY-1 Humidity-dependent CO₂ isotherms and relative humidity.**
  - Paper (SI table "Condition-dependent isotherms"): humidity-enhanced Langmuir, moisture-swing
    (Wang–Lackner) Langmuir, GAB-coadsorption Langmuir/Sips/Toth, GAB-modified Toth–Liu, humidity
    Toth–Climeworks, and RH evaluated from IAPWS–Wagner, Antoine, or a constant saturation pressure.
  - Code: only `GAB` (with a fixed Antoine equation for water, F-PHY-3) and `Langmuir_pH` exist
    (`src/isotherm.h:54`).
  - Needed for 3.3.4: Liu Fig. 4 now freezes the humid Toth law at the feed humidity.
- [ ] **F-PHY-2 Temperature-dependent saturation capacity and heterogeneity exponent.**
  - Paper (SI temperature-dependent isotherms): q_sat constant/linear/exponential in T; ν
    constant/inverse/linear in T.
  - Code: scales only the affinity, through one `HeatOfAdsorption` per component (`src/component.h:75`).
- [ ] **F-PHY-3** GAB saturation pressure is a fixed Antoine equation for water, valid from 274 to 373 K
      (`src/isotherm.h:222`).
- [ ] **F-PHY-4 Pore transport for any gas.** Not needed for 3.3.5 any more.
  - The `BulkConcentration` reaction phase was added on 2026-10-21 (`docs/STADLER_FEATURES_PLAN.md` §1).
    Stadler's Fig. 4 fit is a catalyst-weighted plug-flow reactor without pores.
  - The general feature remains open: every gas should be able to exchange with the pore phase, using its
    own transport parameters (film coefficient and pore diffusivity), independent of chemisorption.
  - Today the bulk ↔ surface ↔ pore exchange exists only for components with a `General`/`Elovich`
    chemisorption site that sets `usePoreSurfaceTransport` (`src/integrators/transport.cpp:41`).
  - For any other species, a `PoreConcentration` reaction changes pore concentrations that never reach the
    bulk. Verified: `examples/Reactions/pore-langmuir-hinshelwood` gives zero gas-phase CO₂ and H₂O after
    5 s of feed.
- [ ] **F-PHY-5** `BulkConcentration` reactions supply no reaction heat; the input reader rejects them together with
      `energyBalance`. To fix: convert the heat per fluid volume to the solid basis with the local solid mass
      per fluid volume (`compute.cpp`, `computeReactionDerivatives`).
- [ ] **F-PHY-6 Variable axial dispersion**: D_ax per component as a function of velocity, e.g.
      D_ax = D_m √ε_b + u r_p (Stadler et al. eq. 12) or D_L = p₁u + p₂ (Stampi-Bombelli), so it follows each
      swing phase.
- [ ] **F-PHY-7** Optional: outlet dead-volume model for extra-column mixing (Stadler et al. and Marx et al. both
      correct their data for it).
- [ ] **F-ITG-2** Second-order total molar flux and gas-temperature advection. With WENO3 (F-ITG-1) the hydrogen
      t50 of Marx T15 still converges at first order (e50 1.2e-2 at N_z = 100 against a 400-point reference), because
      v c_tot and the gas-temperature advection stay upwind. A WENO3 temperature flux alone changed the T15 t50 by
      < 0.1 s (2026-10-26 test build); a higher-order total flux needs a velocity closure that conserves it.
- [ ] **F-BC-1 Flow reversal.**
  - Paper (SI species and temperature boundary conditions): direction-aware Danckwerts inflow/outflow
    and one-way pressure ports.
  - Code: the WENO reconstruction always takes z = 0 as the inlet (`src/integrators/compute.cpp:872`).
  - Needed for counter-current PSA steps.
- [ ] **F-CYC-1 Cycle driver.**
  - Paper (introduction): cyclic-steady-state stopping, process-dependent step termination, and
    per-step feed composition, velocity, outlet pressure, and flow direction.
  - Code: phases accept `Name`, `Temperature`, `InletPressure`, `NumberOfTimeSteps`, and, since 2026-10-21,
    `GasPhaseMolFractions` (feed composition) and `ColumnEntranceVelocity`. A top-level `NumberOfCycles`
    repeats the phase list. Still missing: cyclic-steady-state stopping, step termination on a condition,
    outlet pressure, and flow direction.
  - Section 2.3 of the paper describes the old parser; update it (I-GEN-1).
- [ ] **F-CYC-2** Optional swing-phase extras:
  - per-phase axial dispersion (not needed for Stadler: numerical dispersion v Δz/2 exceeds Eq. 12);
  - rescale the gas inventory at a pressure switch (FixedVelocity already imposes the new pressure on the
    equilibria);
  - end a phase on an outlet mole-fraction threshold (Stadler et al.'s pilot plant).
- [x] **F-PHY-8 Catalyst per bed section** (2026-10-21): a reaction key `"Adsorbents": [names]` restricts it to
      those multibed sections, weighted by their local fraction. The catalyst amount goes into
      `forwardRateCoefficient`.
- [x] **F-PHY-9 Stadler sorption kinetics** (2026-10-21).
  - `FirstOrder` chemisorption has optional keys `DesorptionMassTransferCoefficient`,
    `DesorptionCoverageEnergy`, and `CoverageReferenceLoading`. For q > q*, the coefficient is
    k1 exp[β (q/q_ref)/RT].
  - Tests: `ChemisorptionRate.FirstOrderUsesCoverageDependentDesorptionAboveEquilibrium` and
    `ChemisorptionInput.ParsesFirstOrderDesorptionBranchAndRequiresReferenceLoading`.
- [x] **F-ITG-1 Numerical dispersion of the convective term** (WENO3, the default since 2026-10-27). First-order
      upwind convection adds about vΔz/2 of dispersion: 1.6e-4 m² s⁻¹ for the Marx column at N_z = 100, 38 times
      the physical D_ax, so fitted or reported D_ax values below vΔz/2 had no effect.
  - Input `ConvectionScheme`: `WENO3` (default) or `Upwind`; SIRK3 (implicit upwind matrices) defaults to and only
    accepts `Upwind`. `computeWENO` (uniform) and a non-uniform overload for multibed grids; gas: mole fractions at
    the faces, normalised, times the upwind total flux v c_tot (the flux that the velocity closures conserve);
    liquid: concentrations. The inlet value sits on the inlet face, half a cell from the first cell centre (WENO
    weights 1/2, 1/2 there; dispersion term (c_2 − 3c_1 + 2c_0)/Δz²). Temperature advection stays upwind (F-ITG-2).
  - Tests: `tests/convection_tests.cpp` (exactness for linear profiles, third order, no overshoot at a step,
    non-uniform = uniform, gas fluxes sum to the upwind total, input key); regression tests that hold upwind
    reference values set `Upwind` explicitly.
  - Study `benchmark/convection` (paper Figure `fig:convection`, SI Table `tab:convection_benchmark`): at a fixed
    grid, RK3 −1 to +6%, CVODE +4 to +20% (T15: 8% more steps, 37% more Krylov iterations per step); RK3
    stability limits unchanged. Outlet error E1 with 100 points: BEA 8.1e-3 → 1.3e-4, T15 9.3e-3 → 3.1e-3; monolith
    t50 error with 50 points 4.9e-3 → 8.9e-5. BEA: WENO3 with 25 points beats upwind with 200 at 1/50 of the time.
  - Tested and rejected (scratch builds): van Leer MUSCL (as accurate, RK3 unstable on BEA at Courant 1.23);
    per-species flux reconstruction (T15 stops at ~100 s, "negative outlet pressure"); WENO3 without the inlet-face
    placement in the dispersion term (monolith t50 offset ≈ Δz/2 / front speed, +1.1% at 50 points).
  - Web app fields (`ConvectionScheme`, `CVODEMaximumOrder` in `ruptura/webapp/config.py`): with the time-key
    refactor of the webapp.
- [x] **F-CYC-3 Swing-adsorption extensions for 3.3.5** (2026-10-21).
  - New phase keys `GasPhaseMolFractions` (by component name, normalised) and `ColumnEntranceVelocity`.
  - Top-level `NumberOfCycles` repeats the phase list; phase names get "(cycle k)".
  - `CVODE::reinitialize(column)` re-derives the Δz/v step cap at each switch.
  - Test: `SwingAdsorptionInput.RepeatsCyclesWithPhaseFeedCompositionAndVelocity`.
- [x] **F-MBD-1** Interface refactor implemented from `Refactor.md`: L = Σ L_ads + Σ L_int, separate linear interfaces
      with their own grid counts, strict input validation, Python bindings, GUI persistence and total length,
      migrated examples, and matching docs and paper. All ten new geometry tests and 56 parameterized
      single-adsorbent equivalence tests pass; the full suite has the same 18 pre-existing failures.
      All four multibed examples run, with later breakthrough for finite interfaces and identical sharp/mixed
      curves. Details: `docs/INTERFACE_REFACTOR_VALIDATION.md` (2026-10-06).

## Paper inconsistencies (I)

I-PHY-5 to I-PHY-10 (and the fixed B-PHY-9) were found while building SI Tables S10–S12.

- [ ] **I-PHY-1** `HeatOfAdsorption` sets both the isotherm affinity scaling and the heat released in the energy
      balance. The SI calls these distinct parameters. Either split the input or change the text.
- [ ] **I-PHY-2** Keep the SI humidity table, q_sat(T)/ν(T), and the flow-reversal boundary conditions consistent
      with what gets implemented (F-PHY-1, F-PHY-2, F-BC-1).
- [ ] **I-PHY-3** SI table of conventional isotherms: it writes the Sips exponent as ν, but the code uses 1/ν
      (`src/isotherm.h:166`). The Marx component table now states the code's form.
- [ ] **I-PHY-4** SI pore transport (Eqs. `surface_concentration_balance`, `pore_concentration_balance`,
      `reaction_pore_reservoir`) does not conserve mass as written.
  - β_s, β_p, and Γ are per moving-fluid volume, but F = a_loc k_f (c − c_s) and P = k_p (c_s − c_p) are
    per particle (or pore) volume.
  - The bulk loses (1−ε_b) a_loc k_f Δc per total volume, while the surface reservoir gains ε_b a_loc k_f Δc.
  - Corrected set, consistent with the code after the γ fix (C-PHY-1), which is the case ε_p = 1, β_s = β_p:
    F = ((1−ε_b)/ε_b) a_loc k_f (c − c_s) and P = β_p k_p (c_s − c_p).
    Details in `docs/STADLER_FEATURES_PLAN.md`, "Status after implementation".
- [ ] **I-PHY-5** BET: the code denominator is (1 − k₂ξ)(1 − k₂ + k₁ξ); SI Table S3 has (1 − k₂ξ)[1 + (k₁ − k₂)ξ].
      See B-PHY-4.
- [ ] **I-PHY-6** O'Brien–Myers: the code lacks the factor 1/2 of the SI form; `value()` and `psiForPressure()` also
      disagree with each other (B-PHY-6).
- [ ] **I-PHY-7** Freundlich (like Sips, I-PHY-3): the code uses ξ^(1/ν); SI Table S3 writes ξ^ν. Table S11 lists the
      input as 1/ν.
- [ ] **I-PHY-8** "Modified-Redlich" in Table S3 is not an input type; Redlich-Peterson is.
- [ ] **I-PHY-9** Elovich: the code uses α c exp(−βq) without activation energy or pore-concentration order; the SI
      has α(T_s) with E_α and c^m.
- [ ] **I-PHY-10** Reactions: the code sets k⁻ = k⁺/K with K = K⁰ exp(−ΔG/RT); the SI gives k⁻ its own Arrhenius
      law and a reversibility flag δ_r. LH/LHHW coverages use that same K for every component instead of
      K_{a,j}(T_s). See also B-PHY-2.
- [ ] **I-BC-1** Species inlet condition: the SI describes a Danckwerts (total-flux) inlet, but the code holds the
      inlet node at the feed state, i.e. a Dirichlet inlet (`src/integrators/compute.cpp:652`). Decide whether to
      implement Danckwerts or to describe Dirichlet.
- [ ] **I-BC-2** Monoliths paragraph: the Hagen–Poiseuille vs Ergun explanation does not apply, because every run uses
      `FixedVelocity` (`validation/3.3.2_monoliths_stampi_bombelli`).
- [ ] **I-CYC-1** Introduction: the claims on cyclic operation go beyond Section 2.3 and the code (F-CYC-1).
- [ ] **I-GEN-1** Paper/SI updates for the 2026-10-21 features:
  - "Homogeneous reactions … are not included" no longer holds: describe the catalyst-weighted
    `BulkConcentration` coupling and the reaction `Adsorbents` restriction (F-PHY-8);
  - the `FirstOrder` desorption branch in the chemisorption SI (F-PHY-9);
  - per-phase feed, velocity, and `NumberOfCycles` in Section 2.3 (F-CYC-3).
- [x] **I-BC-3** SI boundary-mode tables now give the RUPTURA setting for BC1–BC5 (gas and liquid), and their
      descriptions follow the code. BC3 is valid for liquids and the non-existent BC6 is removed. The JSON
      section maps the settings to BC1–BC5, and the example no longer sets an unused `PressureGradient`.
- [x] **I-MBD-1** L_int has N_ads − 1 entries (code: `src/column_multibed.cpp:89`); Geometries text corrected.
- [x] **I-MBD-2** Algorithm S1 with user-defined interface fractions removed. Supplementary Section "Geometries" now
      gives the implemented separate linear interfaces, L = Σ L_ads + Σ L_int, section grids, and the added
      inventory L_ads,a + (L_int,a−1 + L_int,a)/2 (`src/column_multibed.cpp`).
- [x] **I-MBD-3** Interface-overlap TODO is obsolete: interfaces are separate sections, so overlap cannot occur by
      construction. Long interfaces are covered by the refactor tests.
- [x] **I-MBD-4** SI "Spatial discretization" now describes per-section and explicit grids with local spacings for
      multistacked columns (`src/integrators/compute_multibed.cpp`).
- [x] **I-IO-1** SI "Json input" rewritten as three multi-page tables (S10 enumerations, S11 isotherm `Parameters`
      order, S12 all keys with units, defaults, and equation references), taken from `src/inputreader.cpp`. The
      old ChemisorptionSites examples used a format the parser rejects.
- [x] **I-EX-1** SI Marx tables now match the marx2 inputs; S1 (selected low-flow runs) removed.

## Performance (P)

- [ ] **P-PHY-1** `Isotherm::psiForPressure`/`value` are not inlined and save all callee-saved registers on every
      call, because of the heavy cases (Toth series, Redlich–Peterson, Bingel–Walton, GAB). An out-of-line slow
      path for those saves 8–10% on Langmuir IAST (BEA), nothing elsewhere; it needs a function called from one
      place. Measured 2026-10-23 (outputs byte-identical), not applied.
- [ ] **P-PHY-2** The four `std::fill`s at the top of `computeFastIAST` are redundant: ~4.5% on BEA RK3. Measured
      2026-10-23 (outputs byte-identical), not applied.
- [ ] **P-PHY-3** Fast IAST convergence (`mixture_prediction.cpp:739`) uses an absolute 1e-13 on the standard
      deviation of ψ; for ψ ≳ 500 this is below double resolution and every call runs all 50 Newton iterations.
      Make it relative.
- [ ] **P-PHY-4** Hot-loop allocations: `computeReactionDerivatives` (`compute.cpp:224`) `activity`/`maximumLoadings`
      per RHS; `reaction.cpp:39-69` `theta` per grid node; Bingel–Walton `psi` (`isotherm.h:369`) two 100-entry
      vectors per call inside the IAST bisection; `computeSegregatedCompetitiveIsotherm` three vectors per grid
      point; `meanParticleNumbers` three per call; `transport.cpp:120` builds an `mdspan` in the innermost loop.
- [ ] **P-PHY-5** `std::pow` with integer exponents: `chemisorption.cpp:49-51`, `isotherm.h:388` (`pow(4, j)`),
      `fitting.cpp:621`, `:626`, `:993`.
- [ ] **P-ITG-1** CVODE step cap `CVODEMaximumTimeStep` = Δz/v (default): in convection-limited cases CVODE takes the
      same steps as RK3 at Courant 1 (BEA C7: ~118,000 steps at every tolerance) plus a Newton–Krylov solve each,
      so it is 1.3–1.8× slower than RK3. The tolerance then sets only the nonlinear-solve accuracy: the default
      rtol 1e-4 gives a t50 error of 2.4e-4, rtol 1e-3 gives 5.4e-5. CVODE also sits at the cap in T15 (3,541
      steps) and the monolith (50,217), so the cap, not the tolerance, sets its cost in all three benchmark
      cases. Check whether a cap of a few Δz/v still avoids the breakthrough "hook" the cap was introduced for.
- [ ] **P-ITG-2** CVODE's default dense linear solver builds the Jacobian by finite differences, one right-hand side
      per state variable (~1700 for 8 components on 100 cells). The `examples/BEA-alkanes-C7/breakthrough`
      standard settings then take ~10 h (IAST); with `"CVODELinearSolver": "SPGMR"` they take ~12 min.
      Consider SPGMR or a banded Jacobian as the default.
  - Profile (2026-10-23): `examples/Multibed/C6-alkanes-layered` (912 states, dense default) spends 63% of its 174 s
    in the Jacobian setup: 25% difference quotients, 38% dense LU (`SUNDlsMat_denseGETRF`). The state is ordered
    block by block (all c, then all q, …), so a banded solver would need a bandwidth of about a whole block; a
    sparse (KLU) Jacobian or SPGMR avoids that.
- [ ] **P-ITG-3** `Timing::measure` calls `steady_clock::now()` 20 times per RK3 step: 3% of the monolith RK3 run.
      Measured 2026-10-23 (outputs byte-identical), not applied.
- [ ] **P-ITG-5 CVODE right-hand side runs the full IAST pass twice** (`cvode.cpp:397`, `for iteration < 2`; see
      B-ITG-7), to make the RHS a function of the state only (sink → velocity → pressure → IAST). For
      `FixedVelocity`, `FixedPressureInletVelocity`, and the liquid boundary modes the velocity and pressure do not
      depend on the sink, so the second pass repeats the first. Gate it on the boundary mode: byte-identical
      results, about half the RHS cost for those cases.
- [ ] **P-ITG-6** SIRK3 copies the whole `Column` twice per step (`rk3_si.cpp:49`, `:274`), as RK3 did before
      `707e64c`, and allocates `solved`, `implicitInvKLs`, and the LAPACK work arrays per call.
- [ ] **P-MBD-1** Static multibed blends recomputed in every RHS (`compute_multibed.cpp:476-481`, `720-729`,
      `849-872`, `909-942`): `fractionOfAdsorbent` is fixed after `initialize()`; precompute per-node tables (the
      unused `prefactorMassTransfer` looks like it was meant for this).
- [ ] **P-IO-1** The `column.json` restart checkpoint (`ReadColumnFile`) is serialized with nlohmann `dump(4)` at every
      output step (`src/breakthrough.cpp`, `src/swing_adsorption.cpp`). For marx-T15 with CVODE (`WriteEvery` 1,
      750 outputs) this is 0.47 s of the 3.1 s wall time (15%). Decide whether every N outputs and the end
      suffice. The dump also includes every cache (about 3× the restart data).
- [ ] **P-PY-1** Python: `plots.metadata()` runs `np.loadtxt` on up to 2000 frames to read one number; `export_mp4`
      renders every frame twice; `all_column_figures` re-parses all files per metric (~22×); `engine="python"` in
      `base_plotly.py:53` is unnecessary for `\s+`.
- [x] **P-PHY-6** MPD reweighting scanned every dense macrostate and re-decoded n per call (~100 µs per evaluation
      for the BEA C7 MPD). It now keeps only the populated macrostates (295 of 78,125) with n, N, and log Π
      precomputed: ~3 µs, bit-identical output (`src/macrostate_particle_distribution.cpp`).
- [x] **P-ITG-4 CVODE run time scales with 1/`TimeStep`.** `propagateCVODE` calls `CVode(..., CV_NORMAL)` once per
      `TimeStep` and then does one full right-hand-side evaluation to refresh the derived column fields. CVODE's own
      steps do not depend on `TimeStep` (it steps past the target time and interpolates). Measured 2026-10-26 on
      marx-T15 (rtol 1e-4, same 2 s output times): `TimeStep` 2 / 0.2 / 0.02 / 0.002 s takes 2.4 / 3.3 / 10.3 / 77 s
      wall time, with 3561 / 3536 / 3530 / 3535 internal steps; t50 agrees to 5e-6 (relative). The default
      `TimeStep` (5e-4 s) applies to every CVODE run that does not set it. Planned fix: time-based output
      (`WriteInterval` in s, one CVODE call per output) and no `TimeStep` for CVODE. See also P-EX-1.
      Fixed 2026-10-27 (time-based I/O): `SimulationTime` and `WriteInterval` in s; CVODE makes one `CVode` call per
      output (or phase end) and rejects `TimeStep`. marx-T15 runs in 1.5 s, with files byte-identical to the
      `TimeStep` 2 s run.
- [x] **P-ITG-7** RK3 with reactions and `AutoNumberOfTimeSteps` evaluates one extra RHS per step
      (`rk3.cpp:449-451`, +33%).
      Fixed 2026-10-27 (time-based I/O): the reaction stop criterion evaluates the RHS once per output, not per
      step.
- [x] **P-EX-1** `marx2/ternary/T15/simulation.json` (and the other marx2 inputs with `TimeStep` 5e-4, `WriteEvery`
      2000): with CVODE each of the 8,000,000 output steps costs a full right-hand-side/IAST evaluation that is
      written only every 2000th step. `TimeStep` 1, `WriteEvery` 1 gives the same files (same CVODE steps) in 10 s
      instead of 15 min. Change the inputs, or skip the post-step derivative evaluation on steps that are not
      written.
      Fixed 2026-10-27 (time-based I/O): every input was converted (`convert.py`, originals in
      `_legacy_simulation.json`); the marx2 CVODE inputs call CVODE once per output (`WriteInterval` 1 s).
- [x] **P-GEN-1 Performance** (2026-10-26, commit `707e64c`; profiled on the `benchmark/breakthrough` cases). Every
      output file is byte-identical and the unit-test results are unchanged; the runs are 1.16–1.76× faster:
  - RK3 copy-constructed the whole `Column` (both `MixturePrediction` objects, all components) and copy-assigned it
    back every step; the stages now advance in place from a saved copy of the state (4–17% of the RK3 run time).
  - `writeOutput` made one `std::print` call per output line; libc++ runs `isatty` and builds a string for each. Each
    snapshot is now one buffer per file (13% of the marx-T15 CVODE run time, `WriteEvery` 1).
  - The fast IAST Newton loop evaluated ψ about 3n times and the last component's loading n − 1 times per iteration;
    `predictMixture` evaluated every isotherm twice when splitting loadings by site and recomputed the exp in
    `scale(T)` (6–30% of the run time; most for single adsorbates such as the monolith).
  - Per-call `std::vector`s in mixture prediction became scratch members (~3%).

## Refactoring and cleanup (R)

- [ ] **R-PHY-1** `computeReactionDerivatives` has four copies of the stoichiometry loop differing only in the target
      span.
- [ ] **R-ITG-1** One right-hand-side function. The precompute-plus-derivative block is written out five times
      (`rk3.cpp` ×2, `cvode.cpp` ×2, `rk3_si.cpp`, `swing_adsorption.cpp`).
- [ ] **R-ITG-2** Template the duplicated integrator helpers (`updateStateRK`, `reactionAutoStopReached`,
      `breakthroughConverged` are identical for both column types), and share the step loop (init ramp,
      write/print cadence, integrator switch, failure dump) between `Breakthrough::run` and
      `SwingAdsorption::run`.
- [ ] **R-GUI-1** `mainwindow.cpp` has ~950 lines of `addNumberField` calls that a `bind(&Config::member)` helper
      would collapse.
- [ ] **R-PY-1** Python dead code: `_resolve_grid_index`, `_rows_at_grid`, `dt_guess`, `_trace_kwargs`,
      `_safe_minmax`, `infer_components_from_files`, `parse_header_metadata`, the unreachable derived-mole-fraction
      path in `_y_values`, duplicate `Cdot`/`Dcdt` specs, `MANIFEST.in`.
- [ ] **R-BLD-1** Rocky 9/10 images build libc++ from source (`docker/rocky-*/install.sh`); switch to
      gcc-toolset-16 once RHEL ships it.
- [ ] **R-BLD-2** `env.yml` uses Clang because conda-forge's `cxx-compiler` still gives GCC 15 (no `<mdspan>`) on
      Linux; switch back to `cxx-compiler` once it moves to GCC 16.
- [ ] **R-GEN-1 Drop the single-bed code path.** `column.cpp` is 82% identical to `column_multibed.cpp`, and `Column`
      is `MultibedColumn` with one adsorbent (`DebugForceMultibed` and `single_adsorbent_equivalence_tests.cpp`
      show it). Removes `column.{h,cpp}`, `compute.cpp`, `sorption.cpp`, `transport.cpp`, half of `rk3.cpp`,
      `cvode.cpp`, and `swing_adsorption.cpp`, and the class of "fixed in one twin only" bugs (SIRK3 layout
      B-ITG-4, pore-dot drift B-MBD-2, `chemisorptionDot` zeroing): 1,500–2,500 lines.
- [ ] **R-GEN-2** Replace the ~330 hand-written lines of copy constructor/assignment (`column.cpp:97-255`,
      `column_multibed.cpp:231-413`), which exist only to call `bindStateViews()`, by an aggregate data base
      class: `Column(const Column& o) : ColumnData(o) { bindStateViews(); }`.
- [ ] **R-GEN-3** Dead code: members `massFlux`, `facePressures`, `numberOfCalls`, `prefactorMassTransfer` (two
      formulas, never read), `adsorbentScaledVoidFraction`, `adsorbentGridPoints`; ~~`computeWENO`~~ (used since F-ITG-1; was: keep if F-ITG-1
      revives it); the two-argument `makeGeometry` overloads in `transport.cpp`; `Reaction::stoichiometricCoefficient`;
      the undefined `*BitString*` declarations in `special_functions.h`; `hypergeometric()`, `maxVectorDifference`,
      `maxTerms`; the never-set `*Flag` members and dead branches of `fit()` in `fitting.cpp`; the 80-line
      bisection duplicated in `isotherm.h:459-550` and `multi_site_isotherm.cpp:88-172`.

## Validation and reruns (V)

Details in `validation/`.

- [ ] **V-PHY-1** Marx T8/T21: marx2 had to use `InletVelocityOutletPressure` and stop early (1872 s and 1174 s). An
      older `marx2/ternary/T21/failed_state.json` shows a diverged state. With the current code, T21
      completes 4000 s with both `InletVelocityOutletPressure` and `InletPressureInletVelocity` (3.4.1).
      Rerun marx2 T21 (and T8) with `InletPressureInletVelocity` to 4000 s; update the SI run table and figures.
- [ ] **V-PHY-2** Find a validation case for pore transport (film and pore reservoirs). Elfving does not use it, and
      Serna-Guerrero is dropped. Candidates:
  - Rosen, J. Chem. Phys. 20 (1952) 387: analytical breakthrough with film and particle diffusion (3.1);
  - CADET general rate model (Leweke and von Lieres, Comput. Chem. Eng. 113 (2018) 274): code-to-code;
  - Kalyanaraman et al., Chem. Eng. J. 259 (2015) 737: CO₂ on PEI/silica hollow fibers, experiments and
    model;
  - Bollini et al., Ind. Eng. Chem. Res. 51 (2012) 15153: breakthrough on zeolite 13X and amine-SBA-15.
- [ ] **V-PHY-3** 3.3.5 Stadler et al. SE-WGS: set up in `stadler/` following the staged plan in its README.
  - Done (2026-10-21):
    - Figs. 4–8 digitized; `conversions.md`; `generate.py`, `run.py`, and the notebook.
    - Fig. 4: 51 runs, and a grid study showing first-order numerical dispersion.
    - Fig. 5 isotherm.
    - Five sorbent-only cycles (Figs. 6, 8A) and five SE-WGS cycles (Figs. 7, 8B).
  - Flow rates from Stadler et al., Energies 14 (2021) 355, Table 2: 2000 mL/min (STP) for adsorption and
    1000 mL/min for regeneration.
  - Open: the stated SE-WGS conditions do not reproduce Fig. 8B (see the README, "SE-WGS result").
    - RUPTURA follows the Table 3 kinetics at 2000 mL/min: X_CO = 0.753, against 0.896–0.91 measured.
    - The CO₂ front is 19.8 s behind CO, against 10.6 s measured and 10.7 s in the authors' model.
    - A lower flow alone cannot explain both. Ask the authors for the 2022 SE-WGS flow, the catalyst, and
      the sorbent batch.
    - The sorbent-only Figs. 6 and 8A agree: the cycle-3 delay is 17.9 s in RUPTURA and in the authors'
      model, and 16.3 s measured.
- [ ] **V-PHY-4** 3.3.6 Liquid and pH: discuss with Shrinjay.
- [ ] **V-PHY-5** Adiabatic thermal mode: not covered by Marx. Dantas runs 1–8 (`dantas/`) are adiabatic if needed.
- [ ] **V-PHY-6** Rerun Marx T6 and T18 and regenerate `marx_T1_T6.pdf` and `marx_T13_T18.pdf` (inputs already fixed,
      B-EX-7).
- [ ] **V-ITG-1** Rerun `benchmark/breakthrough` (README results table, paper Figure `fig:benchmark` and the SI
      benchmark tables): the wall times predate commit `707e64c` (P-GEN-1), which speeds up RK3 and CVODE by
      different factors.
- [ ] **V-CYC-1** 3.4.2 PSA, and TSA purity/recovery: later. They depend on flow reversal and the cycle driver
      (F-BC-1, F-CYC-1).
- [ ] **V-MBD-1** 3.3.7 Park et al. layered beds (Chem. Eng. Sci. 53 (1998) 3951): obtain the paper and set it up.
      Backup: Yang and Lee, AIChE J. 44 (1998) 1325.
- [ ] **V-GEN-1** 3.1 Verification: all items marked new in `validation/3.1_verification/README.md`.
- [ ] **V-GEN-2 Reruns after the 2026-10-26 fixes.** Results made before this date are out of date:
  - gas–wall heat exchange now conserves energy (C-PHY-2; packed beds with `energyBalance` and a wall): marx2
    T1–T24 (with the pending T6/T18 and T8/T21 reruns, V-PHY-6 and V-PHY-1), 3.4.1 Marx, 3.4.2, 3.4.3, Elfving
    (3.3.3, 13 runs), the Stampi-Bombelli packed-bed runs, Dantas, `benchmark/breakthrough` (also due for
    `707e64c`, V-ITG-1), `benchmark/figures`, `benchmark/multibed`; the monoliths are unchanged;
  - Marx `ParticleDensity` 502 → 456.6 kg m⁻³ (B-EX-8; `marx2/particle_density.py`, `marx2/conversions.md` §6):
    marx2, then regenerate 3.4.1, 3.4.2, 3.4.3, and the benchmark from it; SI tables `tab:simulation_settings` and
    the 3.4.3 parameter table still give 502;
  - output times (B-IO-5): every output is now labelled with the time of its state (one `TimeStep` earlier than
    before); only noticeable for Stadler (1 s steps): rerun and remove `output_time_offset` from
    `stadler/stadler_validation.ipynb` in the same step;
  - convection (F-ITG-1): WENO3 is the default since 2026-10-27, so every breakthrough changes; Figure
    `fig:benchmark` and the 3.4.3 text on numerical dispersion (Results "Sensitivity studies": D_ax ≤ 1e-4 m² s⁻¹
    hardly changes t50, N_z = 200 with D_ax = 8.7e-5 reproduces N_z = 100) describe upwind runs and must be
    redone with WENO3;
  - freeze a copy of the binary for the batch (other sessions rebuild `build/ruptura`).
- [x] **V-BC-1** 3.4.1 Boundary conditions: run `validation/3.4.1_boundary_conditions_marx_10bar` and analyze it,
      with the dilute control `validation/3.4.1_boundary_conditions_bea_dilute` (2026-10-19).
- [x] **V-CYC-2** 3.4.2 TSA: two-step cycles on the Marx T15 column, energy balance on
      (`validation/3.4.2_swing_adsorption`, 2026-10-19).
- [x] **V-GEN-3** 3.4.3 Sensitivity analysis: one-at-a-time study on Marx T15 (`validation/3.4.3_sensitivity`,
      2026-10-19).

## Paper writing (W)

- [ ] **W-EX-1** Marx et al. used the multicomponent Sips equation for the ternary predictions (Sect. 4.2.1); marx2
      uses IAST. For the T15 feed at 298 K, multicomponent Sips gives 3.50 mol kg⁻¹ CO₂ against 3.70 (IAST). Marx
      et al. also fitted the feed velocity per experiment (values not published). Mention in the SI.
- [ ] **W-EX-2** Unused Stampi-Bombelli panels 4a, 4c, and 5c–5f: move to the SI or drop.
- [ ] **W-GEN-1** Results headings: mirror the `validation/` structure; keep "Sensitivity studies" (3.4.3).
- [ ] **W-GEN-2** SI Installation and Python API sections are still placeholders.
- [ ] **W-GEN-3** Abstract and Conclusions are empty.
- [x] **W-BC-1** Gas and liquid boundary-mode tables merged into one (Table S7, `tab:pressure_velocity_boundary_modes`):
      input setting × gas/liquid conditions, equations only; the legacy BC1–BC5 names are gone and the
      explanations are in the text.
- [x] **W-MBD-1** Heading "Construction of the multistacked axial geometry." is orphaned now that Algorithm S1
      floats. (2026-10-26: neither the heading nor the algorithm is in `main.tex` any more.)
- [x] **W-EX-3** Serna-Guerrero SI section removed.
- [x] **W-EX-4** `figures/serna-guerrero-fig7.pdf` is still in the paper folder, unused. It is tracked in the paper's
      git and has no identical copy elsewhere. Delete it if it is not needed. Deleted 2026-10-26 (not yet committed
      in the paper's git).
- [x] **W-EX-5** Unused Marx figures removed (`marx_B*.pdf`, `marx_figure_16.pdf`); identical copies remain in
      `marx2/output/pdf/`.
- [x] **W-GEN-4** Citations fixed:
  - `\cite{ruptura1}` (3×) → `sharmadubbeldam2023a`;
  - `farmahinisarkisov2021a` → `farmahinisarkisov2021`;
  - new verified bib entries for Aspen Adsorption, CADET (Leweke and von Lieres 2018), Stampi-Bombelli
    et al. 2024, and NVT+W (Witman, Mahynski, Smit 2018, where the method is named).
- [x] **W-GEN-5** Algorithm S1 converted to the `algorithmic` standard; `\boxed` removed from the direction-aware
      boundary conditions.
- [x] **W-GEN-6** Conciseness pass: per-property vectors/matrices replaced by Table S5, forward/reverse boundary
      conditions and phase-specific fluid symbols stated inline, gas reservoir state in one display, single-row
      `aligned` removed; MPD written without φ and β; Langmuir validation uses component index j.
- [x] **W-GEN-7** STYLE.md applied: `\si{}` units, every equation labelled and numbered, braced `\mathbf{}`, 1-based
      grid nodes in Algorithm S1, grid/component index clash in the SIRK Jacobians fixed (`G_{i,i'}`).
- [x] **W-GEN-8** BibTeX error on `aspenadsorption` (tfo.bst `misc` pops an empty stack) stopped latexmk after the
      first pass, leaving every reference as ??; the entry is now `@manual`.
