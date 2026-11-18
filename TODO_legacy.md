# TODO legacy

Items of `TODO.md` that were completed by 2026-11-02, moved here on that date. Labels are permanent; open
items and new ones are in `TODO.md`.

## Notes

The 2026-10-26 code review (read-through of the C++ core, the fitting code, the input reader, the Qt app, and the
Python package; no build, no tests run) contributed C-BLD-1, B-PHY-5–7, B-ITG-4–6, B-ITG-10, B-MBD-1–2, B-FIT-1–4,
B-IO-1–4, B-GUI-1, B-PY-1–2, P-PHY-3–5, P-ITG-5–7, P-MBD-1, P-PY-1, P-EX-1, and every `R` item.
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

- [x] **B-PHY-2** The reaction key `gibbsFreeEnergy` sets both K(T) = K₀·exp(−ΔG/RT) and the reaction heat (−ΔG·r),
      so it acts as a van 't Hoff reaction enthalpy (`src/reaction.cpp:89`, `src/integrators/compute.cpp:345`).
      Rename it or document this. See also I-PHY-10.
      Fixed 2026-10-30 (temperature laws, `src/temperature_law.h`): every temperature-dependent parameter is a number
      or a law {"Value", "Law", coefficient, "ReferenceTemperature"} (VantHoff, Arrhenius, Exponential, Linear,
      Inverse) at the positions of `Isotherm::temperatureParameters` and the kinetic/reaction keys. Reactions take
      `ForwardRateCoefficient` (Arrhenius), `EquilibriumConstant` (VantHoff), and `HeatOfReaction` = −ΔH_r: a
      VantHoff K without Energy takes HeatOfReaction, and without HeatOfReaction the heat is K's Energy. The removed
      keys (`nonIsothermal`, `gibbsFreeEnergy`, the activation energies, five-parameter GAB) are rejected with their
      replacement; `convert.py` converted every input in the project (1011 files, originals in `_legacy_*`).
      Before/after (frozen binaries): RK3 runs of T15, Stadler WGS, Chemisorption-CO2, and Liu agree to ≤ 9e-12 of
      each column's scale; BEA is identical; CVODE runs differ at the tolerance level (step-size sensitivity).
- [x] **B-PHY-3** Avrami rate: commit `fd613c7` (2026-10-05) dropped the factor n, `k^n t^(n−1) (q* − q)`, and kept
      the old line commented out (`src/chemisorption.cpp:38`). The integrated model, the SI
      (Eq. `eq:avrami_chemisorption`), and `ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel` have
      `n k^n t^(n−1)`. No input uses Avrami. Decide which form is intended.
      Fixed 2026-10-30: the rate is the SI form `n k^n t^(n−1) (q* − q)` again (commented line removed); the
      Avrami test passes, and the unit tests have 1 of 135 failing (B-EX-3).
- [x] **B-PHY-5 Anti-Langmuir IAST inverse has the wrong sign in the exponent.** `isotherm.h:433-437` uses
      `1 − exp(+bψ/a)`; inverting `psiForPressure` (`:268`) gives `1 − exp(−bψ/a)`. Returns negative
      hypothetical pressures.
      Fixed 2026-10-27: `exp(−bψ/a)`. Test `Isotherm.AntiLangmuirInversePressureInvertsTheSpreadingPressure`.
- [x] **B-PHY-6 O'Brien–Myers `value()` and `psiForPressure()` are inconsistent with each other**, not only with the
      SI (see I-PHY-6): `isotherm.h:193` has σ²·x(1−x)/(1+x)³ and `:348` ½σ²·x/(1+x)²; the exact
      integral of the implemented `value()` has no ½. Fix both in the same direction.
      Fixed 2026-10-27 towards the SI and the literature form: `value()` has the ½, so `psiForPressure` is its exact
      integral; also in the Python copy in `plot_fitting.py`. No input used the model. Test
      `Isotherm.OBrienMyersLoadingAndSpreadingPressureFollowTheOBrienMyersForm`; closes I-PHY-6.
- [x] **B-PHY-7** `PseudoNth` chemisorption: `pow(driving, order)` with `driving < 0` is NaN, masked only by the
      argument order of `std::min(0.0, NaN)`; it can never desorb.
      Fixed 2026-10-27: the sign-preserving SI form (Eq. `eq:pseudo_nth_chemisorption`), written k·sign(Δq)·|Δq|^n so
      that it stays 0 at Δq = 0 also for n < 1. Test
      `ChemisorptionRate.PseudoNthIsSignPreservingAndDesorbsAboveEquilibrium`.
- [x] **B-PHY-8** SPI and SCI evaluate the Toth isotherm without the temperature factor `scale`
      (`src/mixture_prediction.cpp`: `safePureSiteLoading` and `computeSegregatedCompetitiveIsotherm`), while
      `Isotherm::value` applies it to b. With two or more adsorbing components and a non-isothermal Toth site,
      SPI/SCI lose the temperature dependence. A single adsorbate takes the `Isotherm::value` shortcut and is not
      affected (e.g. the Stampi-Bombelli monolith). Found 2026-10-23. Fixed 2026-10-26; test
      `MixturePrediction.SPIAndSCIPassTemperatureScaleToToth` (since 2026-10-30
      `SPIAndSCIEvaluateTothAtTheTemperature`, without `scale`). The Liu runs (SPI, Toth) have no
      `HeatOfAdsorption` on CO₂ and are not affected.
- [x] **B-PHY-9** Unilan: `temp1` and `temp2` are identical, so q* = 0 always (`src/isotherm.h`, `value`). Found
      while building SI Tables S10–S12. Fixed 2026-10-26 (e^{−σ} in the denominator, and the same in the
      spreading pressure, which was zero too); test `Isotherm.UnilanLoadingAndSpreadingPressureFollowTheUnilanForm`.
- [x] **B-PHY-10 GAB on a chemisorption site saw T/50 700.** `makeChemisorptionPrediction` gave the chemical mixture
      the heat of the first chemisorption site (50 700 J/mol in the Elfving inputs), so `scale` = exp(50 700/(RT))
      and GAB, which recovered T from `scale`, used T = 0.006 K: p_sat ≈ 6e53 Pa and a GAB loading of ≈ 0. Every
      chemisorption affinity of such a component was scaled with that heat (the Elfving empty-column Langmuir site
      ×7.6e8 at 298 K). Liu (GAB physisorption, heat 1.0) was correct. Found and fixed 2026-10-30 with the
      temperature laws: GAB takes the column temperature, and each site isotherm carries its own laws. Test
      `TemperatureInput.ChemisorbedGABUsesTheTemperatureOfTheColumn`. The Elfving inputs: B-EX-9.
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
- [x] **B-ITG-4 SIRK3 with the energy balance builds the banded LAPACK matrix in the wrong layout.** `rk3_si.cpp:421`
      uses `std::layout_left` (column-major `ab`); its twin `computeConcentrationUpdateMatrixEnergyBalanceFinal`
      (`:522`) uses CTAD, which gives `layout_right` (`git show b8a2562` fixed one copy only). SIRK3 also never
      advances the temperature or chemisorption blocks (`propagate` updates `physisorption` and `concentration`
      only), but the reader rejects only reactions and pore transport (B-ITG-9). Reject `EnergyBalance` and
      chemisorption with SIRK3, and merge the EB/non-EB matrix builders (they differ only in the unused
      `facePressures`).
      Fixed 2026-10-27 (single/multibed merge): SIRK3 runs on the one column type with one matrix builder
      (column-major band storage for I + Δt²L², `rk3_si.cpp`), and the reader rejects `EnergyBalance`,
      chemisorption, and more than one adsorbent with SIRK3.
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
- [x] **B-ITG-12 CVODE with the Dense linear solver is less accurate than SPGMR at the same tolerance.**
      `benchmark/breakthrough` (2026-10-28): at rtol 1e-6, e50 is 1.4e-5 (Dense) against 4.7e-8 (SPGMR) on BEA, and
      6.7e-5 against 3.0e-7 on the monolith; on T15 the two agree within a factor of two. The BEA and monolith Dense
      runs end at BDF order 1. BEA has 294–3552 Newton convergence failures and 759–3862 Jacobian evaluations. The
      monolith at 1e-6 has 61,940 failures, 64,576 Jacobians, and 235,386 steps, against 26,218 steps for SPGMR. At
      rtol 1e-2 the monolith Dense run overshoots to 4.3 times the feed, with t50 = 659 s against 1272 s. Suspect
      the finite-difference Jacobian (non-smooth right-hand side: WENO3 weights, clamps, iterative IAST); not
      investigated. Dense is the default `CVODELinearSolver` (P-ITG-2).
      2026-11-02 rerun: Dense is 160–300 times slower than SPGMR on BEA and 40–50 times on T15 (72 right-hand-side
      evaluations per step for the finite-difference Jacobians), and at rtol 1e-6 on the monolith it takes 9.1 times
      as many steps with e50 6.6e-5 against 8.9e-8. Described in the SI (`fig:benchmark_dense`, paragraph "Linear
      solver"); figure 8 shows SPGMR only. The default is still Dense.
      Done 2026-11-02 (decision): SPGMR is the default `CVODELinearSolver` (`src/inputreader.h`, `cvode.h`, Qt app,
      Ruptura Lab, field help, `docs/GENERATIVE_AI_SIMULATION_PROMPT.md`; SI tables). Why Dense is less accurate
      was not investigated. Inputs that relied on the old default: V-GEN-4.
- [x] **B-FIT-1 GA uniform crossover is a no-op.** `fitting.cpp:522`: `children[i].genotype.substr(j, 1) = …` assigns
      to a temporary.
      Fixed 2026-10-27: the child genotype is sized and copied bit by bit. The BEA nC7 Toth fit ends at the same
      parameters (to 1e-9) after 179 instead of 126 GA generations.
- [x] **B-FIT-2 Nelder–Mead shrink leaves n − 2 function values stale.** `fitting.cpp:931-953` moves every vertex but
      re-evaluates only `f[vg]` and `f[vh]`; correct only for two-parameter isotherms.
      Fixed 2026-10-27: every vertex but the best is re-evaluated after a shrink.
- [x] **B-FIT-3 Fitting column-index defaults underflow.** `inputreader.h:163` defaults `columnPressure{0}`;
      `fitting.cpp:37` subtracts 1 → `SIZE_MAX` → "no pressure points found" unless `ColumnPressure` is set.
      Fixed 2026-10-27: defaults 1, 2, 3 (one-based); 0 and negative values are rejected. SI Table S12 gave the old
      defaults ("req., 1, 2") and now gives 1, 2, 3. Test
      `ExampleInputReader.FittingColumnsAreOneBasedAndDefaultToOneTwoThree`.
- [x] **B-FIT-5** `Isotherm::isUnphysical` reads `parameters[2]` for `Freundlich`, which has two parameters
      (`src/isotherm.cpp:196`): an out-of-bounds read whenever the fitting cost
      (`src/multi_site_isotherm.cpp:180`) evaluates a Freundlich site. Check only `parameters[0]` and `[1]`.
      Fixed 2026-10-26; test `Isotherm.FreundlichPhysicalityChecksOnlyItsTwoParameters`.
- [x] **B-IO-1 `wallDensity` defaults to 0 and is not validated** (`inputreader.h:125`, `:1245`), so
      `1 / (heatCapacityWall · wallDensity)` (`compute.cpp:777`, `compute_multibed.cpp:877`) is `inf` → NaN
      temperatures for any energy-balance run that omits it. (The second half, `PrintEvery`/`WriteEvery`
      not checked > 0, went with the time-based output on 2026-10-27: `WriteInterval` must be > 0.)
      Fixed 2026-10-27: with `energyBalance` the reader requires a finite `wallDensity` > 0 and `heatCapacityWall` >
      0. All 759 energy-balance inputs in the repository pass.
- [x] **B-IO-2** `size_t` inputs read with `get<size_t>()` wrap silently (`"NumberOfGridPoints": -1` → `SIZE_MAX`);
      `readOptionalNonNegativeInteger` exists but is used for a few keys only (four copies of the parsing).
      Fixed 2026-10-27: the eight top-level `size_t` keys use `readOptionalNonNegativeInteger` (no input in the
      repository has a non-integer or negative value for them); the nested reads already checked the sign. The four
      parsing copies remain.
- [x] **B-IO-3** Mole-fraction normalisation divides by zero when no component sets `GasPhaseMolFraction`
      (`inputreader.cpp:1569-1581`); the swing-phase branch guards it, the global one does not.
      Fixed 2026-10-27: the global branch throws as the swing-phase branch does. Test (with B-IO-1 and B-IO-2)
      `ExampleInputReader.RejectsNegativeCountsAllZeroFeedAndMissingWallProperties`.
- [x] **B-IO-4** `catch (const std::runtime_error&)` in `breakthrough.cpp:193` and `swing_adsorption.cpp:298`: other
      exceptions skip the `failed_state.json` dump.
      Fixed 2026-10-27: both catch `std::exception`; GAB, for example, throws `std::domain_error`.
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
- [x] **B-IO-6 Run scripts still look for "Final timestep".** Since the time-based output (2026-10-27) a finished run
      prints "Final time:". `marx2/run.py`, `liu/run.py`, `stadler/run.py`, and the `run.py` of
      `validation/3.4.1_boundary_conditions_{marx_10bar,bea_dilute}`, `3.4.2_swing_adsorption`, and
      `3.4.3_sensitivity` still test for "Final timestep"; also `benchmark/analysis.ipynb`,
      `liu/liu_paper_figures.ipynb`, and `examples/MPD/BEA-alkanes/presentation/presentation.ipynb`. Fixed in
      `validation/3.4.1_boundary_conditions_marx_10bar/analysis.ipynb` (2026-11-02), which marked every run of the
      t50 panel as not completed.
      Fixed 2026-11-02 in all of them (the two notebooks that parsed the step count now read only the final time);
      the four `validation/3.4.*/run.py` also take `--binary PATH`, like `benchmark/breakthrough/run.py`.
- [x] **B-GUI-1** App: `app/runner.cpp:478` `stopJupyter` calls `deleteLater()` through a pointer the synchronous
      `finished` handler has just nulled.
      Fixed 2026-10-27: `stopJupyter` takes the process out of the member and disconnects its handlers before
      terminating it (the `readyReadStandardOutput` handler dereferenced the member as well).
- [x] **B-GUI-2** App: the reaction model knew only the Physisorbed, Chemisorbed, and PoreConcentration phases, so
      an imported `BulkConcentration` reaction (Stadler WGS) became Physisorbed in the modeled input.
      Fixed 2026-10-30 with F-PHY-10: `BulkConcentration` is a phase option in the app, its input normalizer, and
      the web app.
- [x] **B-PY-1** Python: `ColumnDataExplorer` joins the relative `data_dir` twice (`plot_breakthrough.py:991` and
      `:1037`); `tox.ini` lacks `pandas`/`plotly`, so the test suite cannot import.
      Fixed 2026-10-27: the explorer receives file names relative to `data_dir`; `tox.ini` installs `pandas` and
      `plotly>=5.24`.
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
- [x] **B-EX-2** `benchmark/figures/run.py` was deleted in commit `fd613c7` (2026-10-05). Restore it
      (`git show fd613c7^:benchmark/figures/run.py`) or replace it.
      Restored unchanged 2026-10-27; `--dry-run` selects the 186 runs of the current manifest.
- [x] **B-EX-5 Marx gas–wall heat-transfer coefficient (2026-10-19).** All 24 marx2 runs use h_fw = 25.6 W m⁻² K⁻¹,
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
      Closed 2026-10-27: both parts were done on 2026-10-20.
- [x] **B-EX-7** Marx T6 and T18 used k_LDF,CO₂ = 0.07 and 0.15; Marx et al. (Table 6) give 0.04 and 0.11 for these
      298 K runs. Corrected in `marx2/ternary/T6` and `T18/simulation.json` and in the SI table. Rerun: V-PHY-6.
- [x] **B-EX-8** `ParticleDensity` = 502 kg m⁻³ is the bulk density in Marx et al. (see
      `validation/3.3.1_energy_balance_marx_ternary`). Marx et al. (Table 3) write the balance with ε_t, the
      superficial velocity, and ρ_b. With v = Q/(A ε_t), as in marx2, matching it needs ρ_p = ρ_b ε_b/((1−ε_b) ε_t)
      = 457 kg m⁻³; 502 gives 10% more adsorbent per unit gas volume (later breakthrough).
      2026-10-26: T1–T24 use 456.6 (`marx2/particle_density.py`, after `heat_capacity.py`; derivation in
      `marx2/conversions.md` §6). Rerun and SI tables pending (V-GEN-2).

## Features (F)

- [x] **F-PHY-2 Temperature-dependent saturation capacity and heterogeneity exponent.**
  - Paper (SI temperature-dependent isotherms): q_sat constant/linear/exponential in T; ν
    constant/inverse/linear in T.
  - Code: scales only the affinity, through one `HeatOfAdsorption` per component (`src/component.h:75`).
  - Done 2026-10-30 with the temperature laws (B-PHY-2): q_sat Exponential/Linear, exponents Linear/Inverse,
    affinities VantHoff (also Henry, Freundlich, Redlich-Peterson), GAB C and K VantHoff (GAB now gets T from the
    column, three parameters). The SI forms are centred on T_ref; new SI Tables `tab:temperature_laws` and
    `tab:temperature_scalable_parameters`. `MixturePrediction::setTemperature` evaluates the laws once per call.
  - Note: EI sorts the components by q_sat at the input `Temperature` (once, in `sortComponents`); a q_sat law that
    reorders the capacities with temperature is not re-sorted.
- [x] **F-PHY-8 Catalyst per bed section** (2026-10-21): a reaction key `"Adsorbents": [names]` restricts it to
      those multibed sections, weighted by their local fraction. The catalyst amount goes into
      `forwardRateCoefficient`.
- [x] **F-PHY-9 Stadler sorption kinetics** (2026-10-21).
  - `FirstOrder` chemisorption has optional keys `DesorptionMassTransferCoefficient`,
    `DesorptionCoverageEnergy`, and `CoverageReferenceLoading`. For q > q*, the coefficient is
    k1 exp[β (q/q_ref)/RT].
  - Tests: `ChemisorptionRate.FirstOrderUsesCoverageDependentDesorptionAboveEquilibrium` and
    `ChemisorptionInput.ParsesFirstOrderDesorptionBranchAndRequiresReferenceLoading`.
- [x] **F-PHY-10 LH/LHHW adsorption constants per component** (2026-10-30).
  - `Kinetics.AdsorptionConstants` {component name: K_j}, a number or a VantHoff law with Energy, in 1/activity
    (m³/mol). Required for LH/LHHW; every reactant needs K_j > 0; rejected on GeneralPowerLaw.
  - Both styles use the SI coverage θ_j = K_j c_j / (1 + Σ K_i c_i) (`eq:reaction_surface_fraction`); LH is
    irreversible (`eq:lh_reaction_rate`), LHHW subtracts k⁻ Π θ_p (`eq:lhhw_reaction_rate`). Before, both used the
    reaction's equilibrium constant as the adsorption constant of every component, and LH had independent sites.
  - LH/LHHW also run on `BulkConcentration`. Python: `Reaction.adsorption_constants` and the missing
    `Reaction.Phase.BULK_CONCENTRATION`. The Qt app and the web app edit the constants per component, and keep
    them on component renames. `examples/Reactions/pore-langmuir-hinshelwood` uses K_j = b_j R T at 500 K.
  - Tests: `ReactionsRate.LangmuirHinshelwoodUsesCompetitiveCoveragesOfEachComponent`,
    `ReactionsInput.ParsesAdsorptionConstantsByComponentName`, `ReactionsInput.RejectsInvalidAdsorptionConstants`,
    `test_qt_state_migrates_lh_adsorption_constants_by_name`; `LabInput.RenamingComponents…` covers lower-case keys.
  - One set of constants per reaction: per-adsorbent K_{a,j} need one reaction per adsorbent (`Adsorbents`).
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

- [x] **I-PHY-1** `HeatOfAdsorption` sets both the isotherm affinity scaling and the heat released in the energy
      balance. The SI calls these distinct parameters. Either split the input or change the text.
      Done 2026-10-30: a given `HeatOfAdsorption` is the energy-balance heat; an affinity VantHoff law without
      Energy takes it, a law with its own Energy keeps it; without `HeatOfAdsorption` the heat is the law Energy
      (an error with the energy balance if a component's laws disagree). SI text at `eq:adsorption_heat_source`
      and after `tab:temperature_scalable_parameters`.
- [x] **I-PHY-3** SI table of conventional isotherms: it writes the Sips exponent as ν, but the code uses 1/ν
      (`src/isotherm.h:166`). The Marx component table now states the code's form.
      Done 2026-11-02 (decision: the RUPTURA 1 paper is leading, 1/ν): Table S3 writes (kξ)^(1/ν), Table S11 gives
      the input as ν, and the temperature-law text and `tab:temperature_scalable_parameters` say the law acts on ν.
- [x] **I-PHY-5** BET: the code denominator is (1 − k₂ξ)(1 − k₂ + k₁ξ); SI Table S3 has (1 − k₂ξ)[1 + (k₁ − k₂)ξ].
      See B-PHY-4.
      Done 2026-11-02: Table S3 has the code form, which is Eq. 11 of RUPTURA 1; the units text says that RUPTURA 1
      writes it for the relative pressure p/p₀ with dimensionless k₁ and k₂. B-PHY-4 (1 − k₂ mixes a number with
      1/Pa when ξ = p) stays open.
- [x] **I-PHY-6** O'Brien–Myers: the code lacks the factor 1/2 of the SI form; `value()` and `psiForPressure()` also
      disagree with each other (B-PHY-6).
      Fixed with B-PHY-6 (2026-10-27).
- [x] **I-PHY-7** Freundlich (like Sips, I-PHY-3): the code uses ξ^(1/ν); SI Table S3 writes ξ^ν. Table S11 lists the
      input as 1/ν.
      Done 2026-11-02: Table S3 writes K_F ξ^(1/ν); Table S11 gives ν and [K_F] = mol/kg [ξ]^(−1/ν).
- [x] **I-PHY-8** "Modified-Redlich" in Table S3 is not an input type; Redlich-Peterson is.
      Done 2026-11-02: the row, its units sentence, and the "(including the modified-Redlich parameterization)" in
      the SCI text and `tab:mixture_method_compatibility` are removed.
- [x] **I-PHY-9** Elovich: the code uses α c exp(−βq) without activation energy or pore-concentration order; the SI
      has α(T_s) with E_α and c^m.
      2026-10-30: α(T_s) is an Arrhenius law of `Alpha`; the pore-concentration order c^m remains.
      Done 2026-11-02: `eq:elovich_chemisorption` is α(T_s) c_p exp(−b q), with c_j without pore transport and
      uptake up to min(q*, q_max), as `Chemisorption::rate` and the clamp in `sorption.cpp` do.
- [x] **I-PHY-10** Reactions: the code sets k⁻ = k⁺/K with K = K⁰ exp(−ΔG/RT); the SI gives k⁻ its own Arrhenius
      law and a reversibility flag δ_r. LH/LHHW coverages use that same K for every component instead of
      K_{a,j}(T_s). See also B-PHY-2.
      2026-10-30: the SI now states the code's form, k⁻ = k⁺/K_r with a VantHoff K_r(T) whose energy is the heat of
      reaction (`eq:reaction_arrhenius_coefficients`, `eq:reaction_equilibrium_constant`). The LH/LHHW coverages are
      the other session's change (`AdsorptionConstants`).
      2026-10-30: the code now follows the SI coverages with per-component K_j(T) (F-PHY-10). Still open in the SI:
      the reversibility flag δ_r (`eq:pore_reaction_power_law`, `eq:lhhw_reaction_rate`), which the code lacks;
      "available only for the pore-concentration coupling" (`apx:reaction_lh_lhhw`) and "pore phase only"
      (LH/LHHW rows of the style table), since LH/LHHW also run on `BulkConcentration`; and an
      `AdsorptionConstants` row under "Reactions (Kinetics)" in `tab:temperature_scalable_parameters` and
      `tab:json_parameters`.
      Done 2026-11-02: δ_r removed from the three rate laws (every power-law and LHHW reaction is reversible with
      k⁻ = k⁺/K_r, LH irreversible); LH/LHHW available for the pore and bulk couplings (text and style table);
      coverages θ_{a,j,r} with K_{j,r}(T_s) as a VantHoff law (`eq:reaction_surface_equilibrium_coefficient`);
      `AdsorptionConstants` rows in both tables.
- [x] **I-PHY-11** Anti-Langmuir: Table S3 wrote q_sat kξ/(1 − kξ), one coefficient in numerator and denominator;
      the code and RUPTURA 1 have a ξ/(1 − bξ) with two independent parameters. Found and fixed 2026-11-02: Table S3
      K_H ξ/(1 − kξ), Table S11 parameter 1 K_H.
- [x] **I-ITG-1** The benchmark text describes the old CVODE defaults (rtol 1e-4, step cap Δz/v): `main.tex`
      Results "Integrators & benchmark" ("its step stays at its default limit Δz/v", "At the default relative
      tolerance of 10⁻⁴"), the caption of `fig:benchmark` (panel e, "the default step limit of CVODE"), and SI
      `apx:benchmark_study` ("scaled … from their defaults at a relative tolerance of 10⁻⁴", "limited to Δz/v
      (`CVODEMaximumTimeStep` = −1)", "At the default tolerances"). Revise with the V-GEN-2 rerun; SI Table S12
      already gives the new defaults (2026-10-28).
  - 2026-10-28: rerun done (V-ITG-1), figure regenerated with a Dense series, caption updated (Dense, rtol 1e-5,
    4 Δz/v). The Results paragraphs and SI `apx:benchmark_study`, `tab:benchmark_runs`, and `tab:benchmark_scaling`
    still give the 2026-10-21 numbers (`benchmark/breakthrough/README.md`). What changed:
    - CVODE is now faster than RK3 at every stable step, also on BEA. Speed-up × Courant number is 1.19–1.37,
      0.18–0.21, and 0.136–0.144 (was 0.44–0.54, 0.20–0.27, 0.14–0.15).
    - CVODE sits at the cap only on BEA.
    - The runs use WENO3, not upwind.
    - RK3 now levels off at 2.0e-7 on BEA.
    - The stability runs no longer stop with a negative pressure (B-ITG-13).
    - T15 with SPGMR at rtol 1e-2 stops (B-ITG-14).
    - The Dense results are not yet described (B-ITG-12).
  - Done 2026-11-02 with the numbers of the weekend rerun: Results paragraphs, caption, SI `apx:benchmark_study`
    (defaults rtol 1e-5 and 4 Δz/v, step and evaluation counts, Eq. `eq:benchmark_speedup` with the mean CVODE
    step, stability runs, the T15 and BEA paragraphs), `tab:benchmark_runs`, `tab:benchmark_scaling`, the
    convection paragraph and `tab:convection_benchmark`, and the SI integrator defaults. In this rerun the T15 SPGMR
    run at rtol 1e-2 completes, and RK3 agrees with the reference to 5e-8 on T15.
- [x] **I-BC-2** Monoliths paragraph: the Hagen–Poiseuille vs Ergun explanation does not apply, because every run uses
      `FixedVelocity` (`validation/3.3.2_monoliths_stampi_bombelli`).
      Rewritten 2026-11-02 from the rerun (t50 and 10–90% widths against the digitized experiments), with a caption
      that names the panels shown (4b, 4d, 5a, 5b, all 5.6%).
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
- [x] **I-GEN-1** Paper/SI updates for the 2026-10-21 features:
  - "Homogeneous reactions … are not included" no longer holds: describe the catalyst-weighted
    `BulkConcentration` coupling and the reaction `Adsorbents` restriction (F-PHY-8);
  - the `FirstOrder` desorption branch in the chemisorption SI (F-PHY-9);
  - per-phase feed, velocity, and `NumberOfCycles` in Section 2.3 (F-CYC-3).
  - 2026-11-02: the `FirstOrder` desorption branch is done (SI Eq. `eq:first_order_desorption` and the JSON table).
    Still open: "homogeneous reactions are not included" in the main text (Section 2.2, "Reaction kinetics") and
    five times in the SI (Breakthrough modeling, Mass balance, Pore transport, Reaction kinetics), while Results
    "Reactions" uses a bulk-phase reaction; Section 2.3 still says that phases end after a number of time steps
    (now `Duration`) and that repetition and per-phase feed composition are not included, while Results "Cycle
    performance" repeats the cycle to cyclic steady state.
  - Done 2026-11-02:
    - Bulk reactions: source φ_f R^f_j in the main-text and SI species balances, main-text "Reaction kinetics" (with
      the section weight χ_r), new SI subsection `apx:reaction_bulk` (`eq:bulk_reaction_power_law`, catalyst mass in
      k⁺, `Adsorbents`, no reaction heat yet, F-PHY-5), a coupling-table row, `BulkConcentration` and `Adsorbents`
      in the JSON tables; the six "not included" sentences are gone.
    - Section 2.3: Φ_m = (T_m, P_in, v_in, y_in, τ_m), N_cyc repetitions, durations, state carried between cycles,
      the CVODE step limit recomputed at a switch; the list of missing features is removed (I-CYC-1). JSON table
      rows for `NumberOfCycles`, `ColumnEntranceVelocity`, and `GasPhaseMolFractions`; the TSA SI no longer says
      that the phase driver changes only T and P.

## Performance (P)

- [x] **P-PHY-1** `Isotherm::psiForPressure`/`value` are not inlined and save all callee-saved registers on every
      call, because of the heavy cases (Toth series, Redlich–Peterson, Bingel–Walton, GAB). An out-of-line slow
      path for those saves 8–10% on Langmuir IAST (BEA), nothing elsewhere; it needs a function called from one
      place. Measured 2026-10-23 (outputs byte-identical), not applied.
      Done 2026-10-28: `valueSlowPath` (pH Langmuir, Unilan, GAB) and `psiForPressureSlowPath` (pH Langmuir,
      Redlich–Peterson, Toth, Unilan, Bingel–Walton, GAB) are `[[gnu::noinline]]`. The Bingel–Walton loading stays
      inline: its Romberg ψ evaluates it thousands of times, and out of line that IAST ran 27% slower. Outputs are
      byte-identical (BEA C7 RK3, T15 RK3 and CVODE, MOR SIAST, monolith RK3, Cu(dhbc)2 Bingel–Walton IAST); BEA C7
      RK3 24.4 → 23.1 s (−5%), Bingel–Walton IAST 36.6 → 34.4 s (−6%), the others unchanged.
- [x] **P-PHY-2** The four `std::fill`s at the top of `computeFastIAST` are redundant: ~4.5% on BEA RK3. Measured
      2026-10-23 (outputs byte-identical), not applied.
      Done 2026-10-27, in `computeFastSIAST` too. BEA C7 (CVODE and RK3), MOR SIAST, and Marx T15 outputs are
      byte-identical; a 300 s BEA RK3 run took 3% less time.
- [x] **P-PHY-3** Fast IAST convergence (`mixture_prediction.cpp:739`) uses an absolute 1e-13 on the standard
      deviation of ψ; for ψ ≳ 500 this is below double resolution and every call runs all 50 Newton iterations.
      Make it relative.
      Done 2026-10-27 in fast IAST and fast SIAST: the spread of ψ is divided by max(1, |ψ̄|). Outputs change at
      rounding level (T15 with RK3 ≤ 5e-12, MOR SIAST 4e-13, BEA unchanged); CVODE amplifies this to its tolerance
      (T15 at rtol 1e-4: states ≤ 9e-4 and dT_g/dt 0.11 of their column maximum).
- [x] **P-PHY-5** `std::pow` with integer exponents: `chemisorption.cpp:49-51`, `isotherm.h:388` (`pow(4, j)`),
      `fitting.cpp:621`, `:626`, `:993`.
      Done 2026-10-27: Romberg weights 4^j by multiplication (exact) and R² as r·r. The `chemisorption.cpp` exponents
      are non-integer input orders, not integer powers.
- [x] **P-PHY-6** MPD reweighting scanned every dense macrostate and re-decoded n per call (~100 µs per evaluation
      for the BEA C7 MPD). It now keeps only the populated macrostates (295 of 78,125) with n, N, and log Π
      precomputed: ~3 µs, bit-identical output (`src/macrostate_particle_distribution.cpp`).
- [x] **P-ITG-1** CVODE step cap `CVODEMaximumTimeStep` = Δz/v (default): in convection-limited cases CVODE takes the
      same steps as RK3 at Courant 1 (BEA C7: ~118,000 steps at every tolerance) plus a Newton–Krylov solve each,
      so it is 1.3–1.8× slower than RK3. The tolerance then sets only the nonlinear-solve accuracy: the default
      rtol 1e-4 gives a t50 error of 2.4e-4, rtol 1e-3 gives 5.4e-5. CVODE also sits at the cap in T15 (3,541
      steps) and the monolith (50,217), so the cap, not the tolerance, sets its cost in all three benchmark
      cases. Check whether a cap of a few Δz/v still avoids the breakthrough "hook" the cap was introduced for.
      2026-10-28: `benchmark/cvode-max-step` (caps of k Δz/v for k = 1…32 and uncapped, rtol 1e-3/1e-4/1e-5,
      the three benchmark cases). Uncapped shows the hook (largest outlet deviation 4.7e-2 BEA, 5.9e-2 monolith
      at rtol 1e-4); up to 4 Δz/v the error stays at the default level; the monolith error grows 30× from 8 Δz/v.
      At rtol 1e-4, 2 Δz/v runs BEA 1.24× and T15 1.10× faster, the monolith unchanged. rtol 1e-5 with
      4 Δz/v is 5–240× more accurate than the default and 1.0–1.85× faster in all three cases.
      Done 2026-10-28: the defaults are rtol 1e-5 (absolute tolerances 1e-6 mol/m³, 1e-7 mol/kg, 1e-4 K) and an
      automatic cap of 4 Δz/v (`resolveMaximumTimeStep`; no new input key), also in the Qt app, the webapp, the
      AI prompt, and SI Table S12. Runs that set these keys keep their values. See V-GEN-2 and I-ITG-1.
- [x] **P-ITG-2** CVODE's default dense linear solver builds the Jacobian by finite differences, one right-hand side
      per state variable (~1700 for 8 components on 100 cells). The `examples/BEA-alkanes-C7/breakthrough`
      standard settings then take ~10 h (IAST); with `"CVODELinearSolver": "SPGMR"` they take ~12 min.
      Consider SPGMR or a banded Jacobian as the default.
  - Profile (2026-10-23): `examples/Multibed/C6-alkanes-layered` (912 states, dense default) spends 63% of its 174 s
    in the Jacobian setup: 25% difference quotients, 38% dense LU (`SUNDlsMat_denseGETRF`). The state is ordered
    block by block (all c, then all q, …), so a banded solver would need a bandwidth of about a whole block; a
    sparse (KLU) Jacobian or SPGMR avoids that.
  - `benchmark/breakthrough` (2026-10-28, full breakthroughs, rtol 1e-2 … 1e-6): Dense takes 1.7–3.1 h on BEA C7
    (2727 states), 245–530× the SPGMR wall time, and 102–136 s on T15 (38–47×). On the monolith, Dense matches
    SPGMR down to rtol 1e-4 but takes 446 s against 11.7 s at 1e-6. It is also less accurate (B-ITG-12).
  - Done 2026-11-02: SPGMR is the default (B-ITG-12).
- [x] **P-ITG-3** `Timing::measure` calls `steady_clock::now()` 20 times per RK3 step: 3% of the monolith RK3 run.
      Measured 2026-10-23 (outputs byte-identical), not applied.
      Done 2026-10-28: RK3 and SIRK3 time `Total time` once per `advance` call (per output), and the three-phase
      breakdown is measured and printed only with the CMake option `RUPTURA_TIMING_BREAKDOWN` (off by default; the
      notebooks read only `Total time`). Outputs are byte-identical with and without it; the monolith RK3 run gained
      about 1% (7.34 → 7.27 s, near the noise).
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
- [x] **P-ITG-5 CVODE right-hand side runs the full IAST pass twice** (`cvode.cpp:397`, `for iteration < 2`; see
      B-ITG-7), to make the RHS a function of the state only (sink → velocity → pressure → IAST). For
      `FixedVelocity`, `FixedPressureInletVelocity`, and the liquid boundary modes the velocity and pressure do not
      depend on the sink, so the second pass repeats the first. Gate it on the boundary mode: byte-identical
      results, about half the RHS cost for those cases.
      Fixed 2026-10-27: one pass for `FixedVelocity`, `FixedPressureInletVelocity`, and the liquid phase. Byte-identical
      for the monolith and liquid regression cases; where IAST iterates, the second pass had restarted from the
      first pass's warm start, so BEA C7 (`FixedPressureInletVelocity`, CVODE) moves within the CVODE tolerance
      (≤5e-4 of the column maximum) and runs 25% faster (2.35 → 1.76 s, the same 6033 steps).
- [x] **P-ITG-6** SIRK3 copies the whole `Column` twice per step (`rk3_si.cpp:49`, `:274`), as RK3 did before
      `707e64c`, and allocates `solved`, `implicitInvKLs`, and the LAPACK work arrays per call.
      Fixed 2026-10-27 (single/multibed merge): SIRK3 advances the column in place and allocates its work arrays
      once per step.
- [x] **P-ITG-7** RK3 with reactions and `AutoNumberOfTimeSteps` evaluates one extra RHS per step
      (`rk3.cpp:449-451`, +33%).
      Fixed 2026-10-27 (time-based I/O): the reaction stop criterion evaluates the RHS once per output, not per
      step.
- [x] **P-MBD-1** Static multibed blends recomputed in every RHS (`computeBulkSpeciesSink`,
      `computeChemisorptionTransportDerivatives`, `computeEnergyDerivatives`, the Ergun law in
      `updateVelocityAndPressure`): `fractionOfAdsorbent` is fixed after `initialize()`; precompute per-node tables.
      2026-10-27: the per-adsorbent coefficients (LDF rate, axial dispersion) are read into small tables outside the
      grid loops, and the chemisorption loops run only over the sites a component has. Against the single-bed code
      the RK3 runs are still 5–7% slower (BEA C7, Marx T15), mostly the non-uniform-grid WENO and dispersion
      stencils. A per-call WENO weight table did not change the BEA run time measurably but moved the monolith
      CVODE run by rounding alone (24,589 → 28,867 steps, 949 instead of 63 Newton failures), so it was dropped.
      2026-10-28: since the merge every run uses these kernels. In a BEA C7 RK3 profile (one adsorbent)
      `updateVelocityAndPressure` takes 8% and `computeBulkSpeciesSink` 3% of the samples, an upper bound for the
      blends there; multi-adsorbent runs are not profiled.
      Done 2026-10-28: `Column::computeAdsorbentBlends` fills per-node tables (bulk-sink prefactor, pore loading
      conversion, contact areas, pressure-drop law, LDF and dispersion coefficients, heats of adsorption and
      chemisorption) in `initialize()` and `readJSON()`; call it again after changing adsorbent properties (Python:
      `compute_adsorbent_blends`). Six kernels read them; the Ergun gradient of the other boundary modes still
      blends per adsorbent (state-dependent). The sums keep their order, but the simpler loops change the FMA
      fusion in unrelated code: RK3 outputs differ at rounding level (states ≤ 1.5e-13 of the column maximum; the
      monolith, CO2-N2 polishing, and MOR are byte-identical), CVODE within its tolerance (C6 layered 1500 → 1211
      steps, states ≤ 1e-4). RK3 run time: CO2-C3H8 gradient −6.0%, CO2-N2 polishing −5.0%, monolith −4.0%,
      CO2-N2 AB mixture −2.6%, BEA C7 −2.6%.
- [x] **P-IO-1** The `column.json` restart checkpoint (`ReadColumnFile`) is serialized with nlohmann `dump(4)` at every
      output step (`src/breakthrough.cpp`, `src/swing_adsorption.cpp`). For marx-T15 with CVODE (`WriteEvery` 1,
      750 outputs) this is 0.47 s of the 3.1 s wall time (15%). Decide whether every N outputs and the end
      suffice. The dump also includes every cache (about 3× the restart data).
      2026-10-28: marx2 T15 (CVODE, `WriteInterval` 1 s, 4001 outputs) takes 10.4 s, of which 5.5 s integration.
      `column.json` is 155 kB; serializing and writing it with `dump(4)` takes 0.68 ms (2.7 s per run, ~26%),
      compact `dump()` 0.54 ms (91 kB), the state vectors alone 0.15–0.18 ms (20–33 kB). Only `ReadColumnFile`
      reads it.
      Done 2026-10-28: `column.json` (and `failed_state.json`) hold only the state (concentrations, loadings,
      temperatures; 34 kB for T15). The `Breakthrough` constructor recomputes the derived fields after
      `ReadColumnFile` (precompute, rates, precompute, as at a swing phase start); older full files still read.
      A restart reproduces the original output except the rates, which lag one RK3 stage (dq/dt and dc/dt
      ≤ 1.4e-5 of their column maximum, equilibrium loadings 9e-11). Still written at every output. marx2 T15
      10.3 → 8.5 s (−17%), data files byte-identical.
- [x] **P-PY-1** Python: `plots.metadata()` runs `np.loadtxt` on up to 2000 frames to read one number; `export_mp4`
      renders every frame twice; `all_column_figures` re-parses all files per metric (~22×); `engine="python"` in
      `base_plotly.py:53` is unnecessary for `\s+`.
      2026-10-27: `metadata()` takes the frame times from the index (read while indexing), `export_mp4` builds each
      frame once, and `engine="python"` is gone. Open: `all_column_figures` (needs a block cache or a
      `column_animation` that takes parsed blocks).
      Done 2026-10-28: `column_animation` takes an optional keyword-only `block_cache`, which `all_column_figures`
      shares over the metrics: each data file is parsed once (T15: 53 → 5 reads), with identical figures.
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

- [x] **R-PHY-1** `computeReactionDerivatives` has four copies of the stoichiometry loop differing only in the target
      span.
      Done 2026-11-02: one loop over the two sides (reactants with sign −1, products +1) with the phase switch
      inside; outputs byte-identical on the reaction examples and Stadler.
- [x] **R-ITG-2** Template the duplicated integrator helpers (`updateStateRK`, `reactionAutoStopReached`,
      `breakthroughConverged` are identical for both column types), and share the step loop (init ramp,
      write/print cadence, integrator switch, failure dump) between `Breakthrough::run` and
      `SwingAdsorption::run`.
      2026-10-27: with one column type every helper exists once (R-GEN-1); the shared step loop is open.
      Dropped 2026-11-02: deprecated.
- [x] **R-GUI-1** `mainwindow.cpp` has ~950 lines of `addNumberField` calls that a `bind(&Config::member)` helper
      would collapse.
      Done 2026-11-02. By then the file had 87 field calls, 61 of them plain one-member setters. Member-pointer
      overloads of `addNumberField`, `addOptionalNumberField`, and `addComboField` with `selected<Config>()` (the
      selected column, simulation, component, or reaction) make 43 of them one-liners; the bed and swing-phase fields
      (indexed) and the setters with side effects stay. Test `LabInput.MemberBoundFieldsShowAndEditTheSelectedItem`.
- [x] **R-PY-1** Python dead code: `_resolve_grid_index`, `_rows_at_grid`, `dt_guess`, `_trace_kwargs`,
      `_safe_minmax`, `infer_components_from_files`, `parse_header_metadata`, the unreachable derived-mole-fraction
      path in `_y_values`, duplicate `Cdot`/`Dcdt` specs, `MANIFEST.in`.
      Done 2026-10-27, with `HeaderComponentInfo` (used only by `parse_header_metadata`), the `column_data` argument
      of `_y_values`, the `Dcdt` pop in `webapp/plots.py`, and the imports only they used.
- [x] **R-GEN-1 Drop the single-bed code path.** `column.cpp` is 82% identical to `column_multibed.cpp`, and `Column`
      is `MultibedColumn` with one adsorbent (`DebugForceMultibed` and `single_adsorbent_equivalence_tests.cpp`
      show it). Removes `column.{h,cpp}`, `compute.cpp`, `sorption.cpp`, `transport.cpp`, half of `rk3.cpp`,
      `cvode.cpp`, and `swing_adsorption.cpp`, and the class of "fixed in one twin only" bugs (SIRK3 layout
      B-ITG-4, pore-dot drift B-MBD-2, `chemisorptionDot` zeroing): 1,500–2,500 lines.
      Done 2026-10-27: `MultibedColumn` is now `Column` (`src/column.{h,cpp}`); the single-bed `Column`, its kernels
      (`compute.cpp` single-bed half, `sorption.cpp`, `transport.cpp`), `compute_multibed.*`, `DebugForceMultibed`
      (now an ignored input key), and `InputReader::isMultibed` are gone; the multibed kernels are split back into
      `compute.cpp`, `sorption.cpp`, and `transport.cpp`. `Breakthrough` and `SwingAdsorption` are no longer
      templates. Python: `Column`, `Breakthrough`, and `SwingAdsorption` are the former multibed classes, and
      `compute_weno` takes the node positions. Single-adsorbent results change at rounding level with RK3 and SIRK3
      and within the tolerance with CVODE (the multibed kernels sum in a different order). Run time against
      788fc61 (best of 3): BEA C7 RK3 4.67 → 4.99 s, Marx T15 RK3 2.38 → 2.54 s, T15 CVODE 2.99 → 3.10 s,
      monolith CVODE 9.29 → 5.59 s (P-ITG-5); see P-MBD-1. Fixed with it: B-ITG-4, P-ITG-6, the pore part of
      B-MBD-2.
- [x] **R-GEN-2** Replace the ~180 hand-written lines of copy constructor/assignment (`column.cpp`), which
      exist only to call `bindStateViews()`, by an aggregate data base class:
      `Column(const Column& o) : ColumnData(o) { bindStateViews(); }`.
      Done 2026-11-02: `ColumnData` holds every owned member, the enums, the static helpers, and the two
      constructors; `Column` adds the state views, the methods, a forwarding constructor, and the copy operations
      (194 lines removed). Outputs byte-identical on 11 cases (marx2 T15, BEA TSA, reactions, chemisorption, IAST).
- [x] **R-GEN-3** Dead code: `adsorbentGridPoints`; `Reaction::stoichiometricCoefficient`; the undefined
      `*BitString*` declarations in `special_functions.h`; `hypergeometric()`, `maxVectorDifference`, `maxTerms`; the
      never-set `*Flag` members and dead branches of `fit()` in `fitting.cpp`; the 80-line bisection duplicated in
      `isotherm.h:459-550` and `multi_site_isotherm.cpp:88-172`. Removed 2026-10-27 with R-GEN-1: the members
      `massFlux`, `facePressures`, `numberOfCalls`, `prefactorMassTransfer`, `adsorbentScaledVoidFraction`, and the
      two-argument `makeGeometry` overloads; ~~`computeWENO`~~ is used since F-ITG-1.
      Removed 2026-10-27: `Reaction::stoichiometricCoefficient`, the `*BitString*` declarations, `hypergeometric()`,
      `maxVectorDifference`, `maxTerms`, and the five `*Flag` members of `Fitting` with the `refittingFlag` branch of
      `fit()`. Open: `Column::adsorbentGridPoints` (a constructor argument, never read) and the duplicated bisection
      (a refactoring).
      Done 2026-11-02: `Column::adsorbentGridPoints` and its constructor argument removed (`InputReader` keeps its
      own, which it uses); the bisection is one function template, `inversePressureByBisection` in `isotherm.h`,
      called by `Isotherm` and `MultiSiteIsotherm`.

## Validation and reruns (V)

- [x] **V-PHY-6** Rerun Marx T6 and T18 and regenerate `marx_T1_T6.pdf` and `marx_T13_T18.pdf` (inputs already fixed,
      B-EX-7). Done in the 2026-10-30 rerun (both 4000 s with `InletPressureInletVelocity`); figures regenerated
      2026-11-02.
- [x] **V-ITG-1** Rerun `benchmark/breakthrough` (README results table, paper Figure `fig:benchmark` and the SI
      benchmark tables): the wall times predate commit `707e64c` (P-GEN-1), which speeds up RK3 and CVODE by
      different factors.
      Done 2026-10-28 with WENO3, the new CVODE defaults, a frozen binary, and a new CVODE Dense series in the
      precision study. The README and Figure `fig:benchmark` with its caption are updated; the text and SI tables
      are not (I-ITG-1).
- [x] **V-BC-1** 3.4.1 Boundary conditions: run `validation/3.4.1_boundary_conditions_marx_10bar` and analyze it,
      with the dilute control `validation/3.4.1_boundary_conditions_bea_dilute` (2026-10-19).
- [x] **V-CYC-2** 3.4.2 TSA: two-step cycles on the Marx T15 column, energy balance on
      (`validation/3.4.2_swing_adsorption`, 2026-10-19).
- [x] **V-GEN-2 Reruns after the 2026-10-26 fixes.** Results made before this date are out of date:
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
  - CVODE defaults (P-ITG-1, 2026-10-28): rtol 1e-5 and a 4 Δz/v step cap; every CVODE run that does not set
    `CVODERelativeTolerance`/`CVODEAbsoluteTolerance*`/`CVODEMaximumTimeStep` changes, including the
    default-tolerance runs of `benchmark/breakthrough` (the `components` and `methods` studies, and the speed-up
    panel e of `fig:benchmark`);
  - freeze a copy of the binary for the batch (other sessions rebuild `build/ruptura`).
  - 2026-10-28: `benchmark/breakthrough` rerun (V-ITG-1).
  - 2026-10-30/31: every simulation rerun with `run.py` (frozen binary). 3.4.1, 3.4.2, and 3.4.3 still use
    `ParticleDensity` 502 (not regenerated from marx2).
  - 2026-11-02: all data figures of the paper regenerated from the rerun, and the Lab screenshots recaptured
    (`simulations/20261102-102906-7e76c65a-Example-2-T15-CO-N-H`). The text and SI tables are not updated yet;
    statements that no longer match the figures:
    - Energy modeling: the simulated gas temperature now peaks far below the measured one (T17: 315 against
      327 K), not 1–6 K; t50 ratios (`tab:marx_t50_comparison`).
    - Benchmark (Results, "Numerical performance"): "These runs use first-order upwind convection" (now WENO3) and
      the e50 and wall-time values; convection paragraph values.
    - Sensitivity: S(t50, Q_ads) = +0.079 (text −0.16; the binary from before the temperature laws gives the
      same +0.079); h_fs now has |S| = 0.089 (ΔT_max) and is plotted (caption updated); t50 = 691 s; the D_ax and
      grid statements are based on upwind numerical dispersion.
    - Boundary conditions: `FixedPressureInletVelocity` now completes T15 and T21 (t50 −0.3% and −0.4%); the
      text says neither setting completes any Marx run.
    - Monoliths: the caption of `fig:stampi-bombelli` names 400 ppm and packed-bed panels that are not shown
      (the four panels are 4b, 4d, 5a, 5b, all 5.6%); figures now in s and K.
  - 2026-11-02 (later): 3.4.1 Marx, 3.4.2, and 3.4.3 regenerated from marx2 (`ParticleDensity` 456.6; the 3.4.3
    generator now scales the `Value` of the VantHoff affinity law) and rerun with the weekend binary; their figures
    regenerated. All the statements above, the SI tables (`tab:column_settings`, `tab:marx_t50_comparison`, the
    sensitivity, benchmark, and convection tables), and the SI boundary-condition and cycle sections are updated.
    Stadler: rerun in the weekend batch, and `output_time_offset` removed from the notebook (2026-11-02).
- [x] **V-GEN-3** 3.4.3 Sensitivity analysis: one-at-a-time study on Marx T15 (`validation/3.4.3_sensitivity`,
      2026-10-19).

## Paper writing (W)

- [x] **W-BC-1** Gas and liquid boundary-mode tables merged into one (Table S7, `tab:pressure_velocity_boundary_modes`):
      input setting × gas/liquid conditions, equations only; the legacy BC1–BC5 names are gone and the
      explanations are in the text.
- [x] **W-MBD-1** Heading "Construction of the multistacked axial geometry." is orphaned now that Algorithm S1
      floats. (2026-10-26: neither the heading nor the algorithm is in `main.tex` any more.)
- [x] **W-EX-1** Marx et al. used the multicomponent Sips equation for the ternary predictions (Sect. 4.2.1); marx2
      uses IAST. For the T15 feed at 298 K, multicomponent Sips gives 3.50 mol kg⁻¹ CO₂ against 3.70 (IAST). Marx
      et al. also fitted the feed velocity per experiment (values not published). Mention in the SI.
      Done 2026-11-02: SI `apx:marx`, paragraph "Mixture equilibrium and feed velocity".
- [x] **W-EX-2** Unused Stampi-Bombelli panels 4a, 4c, and 5c–5f: move to the SI or drop.
      Done 2026-11-02: 4c and 5c (monolith, 400 ppm) are SI Figure `fig:stampi_bombelli_400ppm` in the new
      subsection `apx:stampi_bombelli`; 5d–5f (dimensionless time) deleted from the paper folder. 4a is not used: it
      would show only the experiment, because `generate_runs.py` skips the 400 ppm packed bed
      (`packed_bed_400ppm_2mmol`, 2 × 10⁵ s); its PDF is still in `figures/stampi-bombelli/`.
- [x] **W-EX-3** Serna-Guerrero SI section removed.
- [x] **W-EX-4** `figures/serna-guerrero-fig7.pdf` is still in the paper folder, unused. It is tracked in the paper's
      git and has no identical copy elsewhere. Delete it if it is not needed. Deleted 2026-10-26 (not yet committed
      in the paper's git).
- [x] **W-EX-5** Unused Marx figures removed (`marx_B*.pdf`, `marx_figure_16.pdf`); identical copies remain in
      `marx2/output/pdf/`.
- [x] **W-GEN-2** SI Installation and Python API sections are still placeholders.
      Written 2026-11-02 (requirements, presets, conda, pip, CMake options, verification; module, classes, helper
      functions, step-wise runs, errors, and an example run against the current build). The main text now says
      nanobind instead of pybind11. Open: the repository URL and DOI of the release (TeX comment in
      `apx:installation`).
- [x] **W-GEN-3** Abstract and Conclusions are empty.
      Written 2026-11-02 from the Results as they stand; the empty Results subsections (liquid, MPD) are not
      claimed. Chemisorption and reactions (Stadler) added to both on 2026-11-02.
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
