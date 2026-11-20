# TODO

Open work for RUPTURA 2.0 and its paper (`ruptura-2.0-paper/main.tex`). Validation cases are organized
per results subsection in `validation/`; each subsection README lists its own points.

Items completed by 2026-11-02 are in `TODO_legacy.md`, with the same labels and sections.

## Labels

Every item has a label `<type>-<area>-<n>`, e.g. `B-GUI-1` (a bug in the Qt app) or `P-ITG-1` (a performance issue
in the integrators). Labels are permanent: completed items are ticked and later moved to `TODO_legacy.md`, and a new item
takes the next free number of its type and area, counting the legacy file. Each section lists its open items first and
then the items completed since the last move, both ordered by area as in the second table.

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
| `ITG` | Integrators and spatial discretization: RK3, CVODE, ROS34PW2, convection scheme |
| `BC` | Boundary conditions and velocity closure |
| `CYC` | Swing-adsorption phases and cycles |
| `MBD` | Multibed (multistacked) columns |
| `FIT` | Isotherm fitting |
| `IO` | Input reader, output files, checkpoints |
| `GUI` | Ruptura Lab (`ruptura/webapp/`, launcher `app/webmain.cpp`); the Qt app until its removal on 2026-11-18 |
| `PY` | Python package |
| `BLD` | Build, CI, and unit tests |
| `EX` | Example, validation, and benchmark inputs and their figures |
| `GEN` | Cross-cutting: the whole code base or the whole paper |

## Critical issues (C)

- [x] **C-BC-1 (was B-BC-3) `InletVelocityOutletPressure` piles up gas behind a steep 5.6 % CO₂ front, also
      isothermally.** Moved to the critical issues on 2026-11-02: the Conclusions recommend prescribing the feed rate
      and taking the velocity from continuity, which this setting does.
      Stampi-Bombelli 5.6 % cases (`stampi-bombelli/dl_from_literature`, switched from `FixedVelocity` to
      IVOP with `OutletPressure` 101325 Pa, `CVODEMaximumOrder` 1, 2026-10-28): before breakthrough the outlet
      velocity is the expected 0.944 v_in, but when the front reaches the outlet c_tot rises to 5.6–11.8× its
      ideal-gas value at ≈ 1 atm on the packed bed (v down to 0.08–0.17 v_in, outlet CO₂ up to 3.3 c_in) and
      1.2–1.5× on the monolith (outlet CO₂ up to 1.2 c_in). Same with `energyBalance` off, so it is not B-BC-2. t90 − t10 collapses from
      390–480 s to 36–57 s (packed bed). Pressure drop is ≤ 45 Pa, so only the sorption effect separates IVOP
      from `FixedVelocity` here. Not yet localized; candidates: the closure's `sinkTerm` against the uptake the
      species balances apply (fixed-point passes, B-ITG-7) under the very steep Toth isotherm (B-ITG-11), and
      N₂ having no `AxialDispersionCoefficient` while CO₂ has one. Blocks a `FixedVelocity`-free rerun of the
      Stampi-Bombelli cases, whose paper model has a momentum balance (Ergun / Hagen–Poiseuille).
  - 2026-11-03 analysis (from the code, not yet checked by a run): the continuity closures set
    v_i Σc_i = v_{i−1} Σc_{i−1} − h_i Σ_j S_ij with Σc from the state, and both convection schemes are conservative,
    so convection and sink cancel in the sum: ∂Σc/∂t = Σ_j D_j ∂²c_j/∂z² exactly, and nothing pulls Σc back to
    P/RT. With D on CO₂ only this is D_CO₂ ∂²c_CO₂/∂z². A front moving at w carries a bump δ ≈ D |∂c_CO₂/∂z| / w,
    large for a slow, steep Toth front (D/w ≈ 1 m on the 0.5 mmol/s packed bed). At the outlet node the
    zero-gradient stencil lets CO₂ disperse in, D (c_{N−1} − c_N)/h² (D/h² ≈ 2 s⁻¹ there, up to ≈ 11 % of c_tot
    per second), while the closure keeps the convective outflow fixed, so those moles cannot leave: v_N falls and
    Σc_N rises, as observed. With one D for all components, ∂Σc/∂t = D ∂²Σc/∂z² and Σc stays at P/RT.
    Check: rerun `packed_bed_5.6pct_0p5mmol` with IVOP and N₂ given CO₂'s D (D = 0 for all as control).
    Fix options: dispersive fluxes that sum to zero (correction flux J_j − y_j Σ_k J_k, `computeMassDerivatives`
    and the Jacobian), or the net dispersive flux in the closure; plus a logged max |Σc R T/P − 1| so that any other
    drift (clamps, B-ITG-11) is no longer silent.
  - 2026-11-03, root cause in the history: a regression. RUPTURA 1 (`7fd78e0`, `Breakthrough::computeVelocity`)
    and every version up to `4116cbe` put Σ_j D_j ∂²c_j/∂z² into the velocity closure next to the mass-transfer
    sink, so per-component D_j kept Σc at P/RT. `ad83222` ("chemisorption", 2026-07-30) replaced that sum by
    `sinkTerm` (`bulkSpeciesSink` only) and dropped the dispersion term. The two fix options are the same model:
    with the term in the closure, N_j = v c_j + J_j = u c_j + (J_j − y_j Σ_k J_k) with the molar-average velocity
    u = v + Σ_k J_k / c, i.e. the correction-flux form. Per-component D_j stay meaningful: a dilute component keeps
    its own D_j, in a binary the pair acts as one D₁₂ = (1 − y₁) D₁ + y₁ D₂ (exact Fick for D₁ = D₂). The
    restored closure term has to use the discrete dispersive fluxes of `computeMassDerivatives`, including the
    Danckwerts inlet face, or a residual drift remains.
  - 2026-11-03: fixed. `updateVelocityAndPressure` adds h_i Σ_j D_j,i δ²_i c_j (lambda `dispersionTerm`, the stencils
    of `computeMassDerivatives`; with Danckwerts the inlet face is left out at node 1, since convective plus
    dispersive inflow there is the feed flux) to the IPIV, IPOP, and IVOP marches. Test
    `ContinuityClosure.KeepsTheTotalConcentrationWhenTheComponentsDisperseDifferently` (Σ_j dc_j/dt = 0 at every node,
    both schemes, both inlets; fails without the term). Every output now logs max |R T Σc / P − 1| for the continuity
    closures and warns once above 1 % (`Column::idealGasDeviation`). SI: Eqs. `gas_total_molar_balance`,
    `correction_flux_form`, `discrete_gas_total_molar_balance`.
    Stampi-Bombelli `packed_bed_5.6pct_0p5mmol`, IVOP, ROS34PW2 (old → new; FixedVelocity): isothermal max |RTΣc/P − 1|
    9.8 → 8.2e-3, min v_out/v_in 0.088 → 0.941, t90 − t10 48 → 408 s (420 s); with the energy balance 9.6 → 2.8e-2
    (B-BC-2), 60 → 480 s (504 s). With N₂ given CO₂'s D the old code was already right (4.7e-4), which confirms the
    cause. `monolith_5.6pct_0p18mmol`: 0.23 → 2.7e-4. CVODE on the packed bed (isothermal, 3000 s): 9.4e-5, the
    initial Ergun offset. The ROS34PW2 residual is integration error that Σc no longer damps (F-ITG-4); it is
    1.3e-3 at rtol 1e-7. RK3 at 0.02 s drifts by 3.2 in 120 s with or without the fix and with D = 0 (B-ITG-13).
    Paper inputs with a continuity closure and unequal D: Marx T15 (ternary, 3.4.1 IPIV/IPOP/IVOP), outlet mole
    fractions change by ≤ 2.2 % of the final H₂ value at the H₂ front, ≤ 0.2 % for CO₂ and N₂, t50 unchanged at
    the 1 s output (rerun: V-BC-2).
- [x] **C-BC-2 The initial and feed concentrations use `Temperature`, the gas temperature `InfluxTemperature`.**
      `Column::initialize` (`column.cpp:403`) sets c_j = y_j P / (R T_ext) on every node, the feed node 0 included,
      but fills the gas, solid, and wall temperatures with `InfluxTemperature` (`column.cpp:465`). Where the two
      differ, R T Σc / P starts at T_in / T_ext everywhere and the feed carries that factor more moles at the
      prescribed inlet velocity: 318/298 = 1.067 in the 318 K Marx runs (30 inputs in `marx2`, 76 in
      `validation`, e.g. T7–T12, T19–T24, and 3.4.1 T9/T21), if Q is the volumetric flow at feed conditions. Found
      2026-11-04 with B-BC-2: the closure now keeps Σc T constant, so T20 and T21 stay at the initial 6.7 %. The
      swing driver uses the phase temperature for node 0 (`swing_adsorption.cpp`), so only the first state is
      affected there. Fix: use the gas temperature in `initialize`; then rerun the 318 K Marx runs (t50 expected
      about 6 % later).
      Fixed 2026-11-04: `initialize` takes P/(R T_in), the initial gas temperature (test
      `ColumnInitialization.FillsTheGasAtTheFeedTemperature`). Paper runs rerun with
      `build/ruptura-2026-11-04-bbc2` (B-BC-2 and C-BC-2): marx2 ternary, 3.4.1 (10 bar), 3.4.2, 3.4.3
      (`run_2026-11-04_bbc2.log` in each; all 89 finish). Notebooks re-executed, Marx figures copied, paper updated:
      - Marx 318 K runs (T7–T12, T19–T24): CO₂ t50 +6.4 to +7.2 %, so they now lie further above the experiment
        (T20 0.99 → 1.06, T7 1.37 → 1.47); the 298 K runs change by ≤ 0.2 %, peak temperatures by ≤ 0.3 K.
        `tab:marx_t50_comparison` updated; main text "1–16 % later (40/30/30), 5–18 % later (15/80/5)" at ≥ 5e5 Pa
        (was "1 % earlier to 12 % later", "1–18 %"); abstract and Conclusions ("within 18 %") still hold. This
        assumes the 1e-5 m³/s of Marx et al. is the volumetric flow at feed conditions, which the inputs already do
        (one v_in at both temperatures).
      - 3.4.1: T9/T21 IVOP t50 +6.7 %; IPOP shift 10–12 % (15/80/5, was 10–11 %); FixedVelocity and
        FixedPressureInletVelocity 0.01–0.32 % earlier (was 0.05–0.45 %).
      - 3.4.2 TSA: working capacity +2–3 % (93 % of equilibrium at 16000 s, productivity 1.2e-4 mol/kg/s), CO₂
        balance closure 4e-5 → 1.4–1.7e-4 (not quadrature: subsampling to 40 s and 80 s does not reduce it), CSS
        measure 6e-6 → 3e-5. c_tot R T_g − P_Ergun now −2.9 % to +4.0 % (was −4 % to +25 %); SI remark rewritten.
        Wall times doubled, but other sessions loaded the machine; the paper keeps "1–12 s per cycle" until a quiet
        rerun (Stage 1 measured +7 %).
      - 3.4.3: reference t50 634.7 → 633.9 s; coefficients change by ≤ 0.004; D_ax 1e-3 now −0.2 % / +45 %,
        numerical-dispersion check 2.7 % / 3.7 %; reference vs marx2 T15 within 1.8e-6.
      - SI `apx:gas_total_molar_balance`: "Without the last term" → "Without the dispersive term" (the thermal term
        is now last).
- [x] **C-MBD-1 Gas-phase multibeds do not conserve the molar flow where the void fraction changes.** Found
      2026-11-06 in a code review; from the code, not yet checked by a run. The gas species balance uses the
      interstitial flux, −(v_i c_i − v_{i−1} c_{i−1})/h_i (`computeMassDerivatives`, both schemes, and the sorption
      Jacobian), and the continuity closures keep v Σc continuous from node to node (`updateVelocityAndPressure`);
      `FixedVelocity` holds one v everywhere. Where φ_f changes, at a sharp interface or along a linear one, the molar
      flow per column area φ_f v c therefore jumps by φ₂/φ₁, and the interstitial velocity does not change. The SI
      claims the opposite: Eq. `species_balance_general` "is retained across multistacked interfaces so that the
      superficial species flux remains conservative", with [[U]] = 0 for the gas (Eq.
      `hydraulic_interface_conditions`). The liquid path does weight its fluxes by φ_f (superficial, `fluidFraction`
      in the WENO3 branch, `totalVoidFraction` in the upwind branch). The gas-temperature advection has the same
      form. Affected: every gas multibed with unequal void fractions, i.e. `examples/Multibed/CO2-N2-polishing`,
      `CO2-N2-AB-mixture` (0.38/0.42, about +10 % molar flow into the second bed), `CO2-C3H8-gradient`
      (0.38/0.40/0.43), and `C6-alkanes-layered` (0.39/0.40); no paper input (the 80 `benchmark/multibed` inputs
      have equal void fractions), but 3.3.7 (V-MBD-1) will. Check: an inert tracer through `CO2-N2-polishing` with
      `FixedVelocity`, outlet against inlet molar flow per area. Fix: superficial fluxes φ_f v c for the gas too (the
      closures then march φ_i v_i Σc_i, and the gas and liquid upwind loops of `computeMassDerivatives`, now two
      near-copies, become one), the dispersive flux φ_f D ∂c/∂z, and the same in the Jacobian and the gas-temperature
      balance.
  - Fixed 2026-11-09: the gas fluxes are superficial, φ_f v c, as in the liquid; the closures march φ_i v_i Σc_i and
    carry φ_i h_i times the uptake, the dispersion, and the thermal expansion; `FixedVelocity` holds the superficial
    velocity φ_f(0) v_in for a gas too (v = φ_f(0) v_in/φ_f); the gas and liquid upwind loops of
    `computeMassDerivatives` are one; the sorption and thermal Jacobians and the gas-temperature balance follow (a node
    that gas leaves through its lower face takes φ_{i−1} v_{i−1}/φ_i there). Checked by runs: in
    `CO2-N2-polishing` the helium outflow per column area is 1.0000 of the inflow at the end, against 1.1053 before;
    t50 of the first adsorbing component moves from 112.3 to 114.2 s in `CO2-C3H8-gradient` and from 461.4 to 470.2 s
    in `C6-alkanes-layered`. `CO2-N2-AB-mixture` is a uniform blend (one void fraction everywhere) and does not change;
    the list above was wrong about it. Single beds change at rounding level only (RK3 ≤ 1e-11 relative; CVODE runs,
    the BEA PSA cycle included, ≤ 2e-5, CVODE's amplification). Test
    `MultibedColumn.ConservesTheSuperficialMolarFlowWhereTheVoidFractionChanges`. SI: the `FixedVelocity` row of Table
    `tab:pressure_velocity_boundary_modes` and the method-of-lines text now state the superficial form for both phases.
    Not done from the plan: the dispersion stays D δ²c per fluid volume of the node, as in the liquid, rather than a
    φ_f-weighted face flux; it differs only at an interface, and the closures balance it either way.

## Bugs (B)

- [ ] **B-PHY-1** MPD temperature reweighting adds N·ln(T_ref/T) (`src/macrostate_particle_distribution.cpp:200`),
      but the Langmuir affinity scaling (`src/utils.h:35`) has no such term. Away from T_ref, MPD and
      SCI/IAST therefore differ by a factor T_ref/T in the affinity. Decide which is intended.
  - 2026-11-11, TMMC test (V-PHY-7, pure CO₂ in MOR, 300 K reweighted to 350 K against a direct 350 K TMMC, both
    flat up to 160 molecules): with the term −5…+0.5 % at 1e2–1e6 Pa, without it +11…+16 % up to 1e4 Pa (the
    Henry-regime factor T_ref/T). The data support the MPD term; propane does not decide (350 K window defect).
    Paper Results "MPD", Fig. mpd_mor e,f. The Langmuir side is still open.
- [ ] **B-PHY-4** BET isotherm (`src/isotherm.h:151`): q = q_m b p / ((1 − c p)(1 − c + b p)). Textbook BET with
      b = C/p₀ and c = 1/p₀ has (1 − c p + b p) as the second factor; `1 − c` mixes a pure number with 1/Pa.
      The spreading pressure (`src/isotherm.h:270`) matches the code's form, so IAST is self-consistent, but
      the model is not BET. Unchanged since the initial commit; check against RUPTURA 1 and fix or document.
      The SI form is I-PHY-5.
- [ ] **B-PHY-11** `PoreConcentration` reaction heat has no conversion to the solid-mass basis. The pore reservoir
      is per solid volume (`transport.cpp`), so Q r is in W per m³ of solid, but `computeReactionDerivatives` adds
      it to `reactionHeat` (W/kg) unconverted; the solid balance divides only by c_s. The heat is therefore
      ρ_p times too large (950× in `examples/Reactions/pore-langmuir-hinshelwood`, energy balance on). The SI
      (Eq. `pore_reaction_heat_source`) weights by Ψ_p = (1 − ε_b) ε_p, the code's reservoir has no ε_p; decide
      which volume the pore concentration refers to. Found 2026-11-03 while doing F-PHY-5; not checked by a run.
      2026-11-20 (I-PHY-4): the SI side is settled. Ψ_p is now the adsorbent volume without ε_p, the basis of the
      reservoirs. The code still lacks the 1/ρ_p.
- [ ] **B-PHY-16** A `PoreConcentration` reaction on a chemisorption site without pore transport acts on a closed
      reservoir. `computeChemisorptionTransportDerivatives` and `computeBulkSpeciesSink` skip sites without
      `UsePoreSurfaceTransport`, so the reservoir of that site exchanges nothing with the bulk. The reader checks only
      the site index (`inputreader.cpp` ≈3687). Reject such reactions in the reader, or connect the reservoir. The SI
      ("Pore transport", "Reaction without adsorption limitation") now says that pore reactions belong on sites with
      pore transport. Found 2026-11-20 from the code; not checked by a run.
- [ ] **B-ITG-2** RK3 close to its stability limit is inaccurate on the dispersion-dominated monolith (fitted D_ax,
      400 ppm): e50 = 0.098 at Δt = 1.25e-3 s (limit ≈ 1.66e-3 s), 2.8e-4 at 4e-4 s. Possibly the clipping of
      negative concentrations at the front (`updateStateRK` max(0, ·); `stampi-bombelli/generate_runs.py`
      mentions it); not checked. Matters for users who pick the largest stable step.
      2026-11-02: also JUC-77-xylenes (Upwind) at the default 5e-4 s: max |Δc|/c_feed 4.1e-3 at the ethylbenzene
      roll-up (t = 235 s) against RK3 at 1.25e-4 s, 2.5e-5 at 2.5e-4 s; RK3 diverges at 3e-3 s (B-ITG-1 check).
      2026-11-03: the monolith on 12 cells (WENO3, 300 s; test `ROS34PW2Column.IntegratesChemisorptionLikeRungeKutta3`):
      against RK3 at 5e-5 s, RK3 is 0.18 off in q at 5e-3 s and 1e-2 at 1e-3 s, far below its stability limit of
      3e-2 s on that grid, and converges only below 2e-4 s (2e-7); ROS34PW2 at rtol 1e-7 is within 5e-7.
      2026-11-03: `stampi-bombelli/dl_from_stampibombelli/monolith_5.6pct_0p02mmol` (fitted D_L, Danckwerts): RK3 at
      the original 0.5 ms puts t10 at 903 s, at 0.25 ms 1013 s, at 0.1 ms 1118 s, against 1262 s from CVODE (two
      settings) and ROS34PW2; t50 1269 s at 0.5 ms against 1622 s. RK3 is not converged at these steps.
- [ ] **B-ITG-11 CVODE at BDF order > 1 drives the gas concentration negative ahead of a steep front, hidden by the
      output clamp, until the run stalls.** `stampi-bombelli/dl_from_literature/monolith_5.6pct_0p02mmol` aborts at
      t = 129.3 s with `mxstep` (100 000 internal steps within one 12 s output interval). `failed_state.json`
      (CVODE's unclamped state) has CO₂ < 0 in all 38 cells ahead of the front, down to −0.043 mol m⁻³ (2 % of
      the feed); the last output (t = 120 s) shows exact zeros there, because `clampNonnegative`
      (`cvode.cpp:47`) only touches the user copy (B-ITG-10). The Toth isotherm (b(T) = 6.4 Pa⁻¹, ν = 0.25) has
      dq*/dc at c → 0 about 5·10⁴ times its chord slope at the feed (LDF eigenvalue ≈ −1e6 s⁻¹), and q* = 0 with
      zero slope for c < 0. Where q = 0 nothing restores a negative c, so undershoots of the BDF extrapolation
      accumulate and are carried downstream; once the front reaches them, Newton straddles the kink at c = 0 and
      the step size collapses. Tests (frozen binary, 1200 s unless noted): `CVODEAbsoluteToleranceConcentration`
      1e-8 fails at 302 s; `ConvectionScheme` Upwind at 171 s (so not WENO3); `CVODEMaximumOrder` 2 at 337 s;
      Toth ν = 1 at 10.5 s (min −0.36 mol m⁻³); `CVODEMaximumOrder` 1 runs the full 24 000 s in 7 s (14 553
      steps), outlet c/c₀ within 1.5e-3 of rtol 1e-5. On the 0.06 mmol/s case order 1 matches the default order
      within 0.6 s in t10/t50/t90 with 3.6× fewer steps. The Stampi-Bombelli D_L (25× larger) runs at the default
      order. Open: whether `CVodeSetConstraints` (checked on marx2 only, B-ITG-10) prevents the build-up, and how
      much hidden negative mass the runs that do finish carry. Found 2026-10-28.
  - 2026-10-28, Figure 5a (packed bed, 5.6 %, literature D_L): WENO3 at the default order matches order 1 to
    1e-4 in c/c₀, so these runs are unaffected. Upwind at the default order stalls the same way on 0.5 mmol/s
    (t = 4691 s) and puts t50 17 s early on 1 mmol/s; at order 1 upwind and WENO3 agree within 0.014 in c/c₀
    and 13 s in t90 − t10.
  - 2026-11-02, weekend rerun: `dl_from_literature/monolith_5.6pct_0p06mmol` (t = 62 s) and
    `packed_bed_5.6pct_0p5mmol` (t = 572 s) now stall at the default order too; the binary from before the
    temperature laws also stalls on both at the default order (t = 21 s and 4100 s). Both inputs now set
    `CVODEMaximumOrder` 1, like `monolith_5.6pct_0p02mmol` (and `dl_from_stampibombelli/monolith_5.6pct_0p06mmol`),
    and finish in 15–20 s. On the packed bed, order 1 matches the default-order profile inside the bed to 8e-4
    in c/c₀ up to t = 4092 s. `stampi-bombelli/generate_runs.py` does not write any of these overrides.
  - 2026-11-02 (decision, reverted on 2026-11-03): the default `CVODEMaximumOrder` is 1 (`src/inputreader.h`, `cvode.h`, Qt app, Ruptura
    Lab, field help, the AI prompt doc; SI table and integrator paragraph). Orders > 1 still have the problem when a
    user sets them; `CVodeSetConstraints` was not tried. Tests that pin CVODE reference values now set order 5 and
    Dense explicitly (`fluid_phase_tests.cpp`). Inputs that relied on the old default: V-GEN-4.
  - 2026-11-03 (decision): reverted. The default `CVODEMaximumOrder` is 5 again, because order 1 is first-order
    accurate in time and would change every paper CVODE result. Order 1 is documented as the remedy for steep fronts
    that stall (SI integrator paragraph, field help, AI prompt doc); the Qt app, Ruptura Lab, and the reader default
    to 5. Inputs that need order 1 set it: the three `dl_from_literature` runs above and
    `dl_from_stampibombelli/monolith_5.6pct_0p06mmol`; since 2026-11-03 `stampi-bombelli/generate_runs.py` writes
    these overrides (`CVODE_OVERRIDES`). Open: `CVodeSetConstraints` (c ≥ 0) or another positivity fix, so that the default order runs these
    cases.
  - 2026-11-03: with `InletCondition` Danckwerts, `dl_from_stampibombelli/packed_bed_5.6pct_0p5mmol` and `_1mmol`
    (vL/D ≈ 0.04, Dense, default order) stop with `mxstep` at t = 72 and 96 s; the Dirichlet runs finish. At order 1,
    1 mmol/s finishes (t50 2976 s, experiment 2612 s); 0.5 mmol/s crawls (t = 768 s after 25 min, stopped). With
    Danckwerts the whole dispersion-dominated bed sits near c = 0 under the steep Toth isotherm, the situation
    described above.
  - 2026-11-03: with the fitted D_L and Danckwerts, CVODE needs per-case settings (now `CVODE_OVERRIDES` in
    `stampi-bombelli/generate_runs.py`):
    - `monolith_5.6pct_0p02mmol` (vL/D ≈ 0.3): order 1 with rtol 1e-7. Order 1 alone, SPGMR, SPGMR at order 1,
      atol 1e-9, and a 1 s step cap all stop with `mxstep`; at the stop, c ≥ −3e-6 mol/m³ and the chemisorbed
      loading is slightly negative.
    - `packed_bed_5.6pct_0p5mmol`: SPGMR with the SorptionJacobian preconditioner (132 s). KLU is not allowed with
      the energy balance, and Dense at order 1 with rtol 1e-7 crawls.
    - `packed_bed_5.6pct_1mmol`: order 1.

    ROS34PW2 at its default settings, with no CVODE keys, runs all 20 Stampi-Bombelli cases in about 2 s each. It
    agrees with the CVODE runs within 0.1 % in t50 and 6e-4 in c/c₀ (fitted D_L), and within 1.8e-2 at the steep
    literature-D_L 400 ppm fronts. In ruptura-5a's `validation/ros34pw2_danckwerts_rerun`, all 223 cases outside
    Elfving finish, including the 6 CVODE stops of 3.4.1 and 36 of Liu (B-EX-10). The stall is a CVODE (BDF)
    problem; a positivity fix or ROS34PW2 as the default for these cases would remove the overrides.
  - 2026-11-20, full rerun on the cluster node: four Stampi-Bombelli runs of
    `validation/inlet_condition_dirichlet_vs_danckwerts` stop again, as on 2026-11-03:
    `dl_from_stampibombelli/monolith_5.6pct_0p02mmol` under both inlets, and `packed_bed_5.6pct_0p5mmol` and `_1mmol`
    under Danckwerts. That study keeps the CVODE inputs of 2026-11-03, from before the per-case settings above
    (`CVODE_OVERRIDES`), and these stops are among its findings, so the inputs stay as they are. `cases.csv` now
    names the four runs as expected stops.
- [ ] **B-ITG-13 RK3 above its stability limit completes without an error.** In the `benchmark/breakthrough`
      stability runs (2026-10-28), BEA at 1.5e-2 s and the monolith at 2e-3 s run into unbounded outlets (58× and
      2500× the feed), and T15 at 5e-2 s collapses to an outlet of exactly zero for every component after a few
      seconds. All three return 0. On 2026-10-21, T15 and the monolith stopped with "negative outlet pressure". A
      user who picks too large a step gets no warning. Check what changed (clamps, B-ITG-10) and add a divergence
      check.
      2026-11-02: SIRK3 does the same (B-ITG-1 check): JUC-77 at 1e-2 s ends with outlet helium at 2.4 P/RT, the
      liquid pH case at 1 s with 69 times the feed; both return 0, as do RK3 on JUC-77 at 3e-3 s and on the liquid
      case at 0.1 s.
      2026-11-02: deferred ("fine for now").
      2026-11-03: also the Stampi-Bombelli packed bed (5.6 %, 0.5 mmol/s, IVOP, isothermal) at 2e-2 s, beyond the
      explicit-sink limit at the Toth front edge (≲ 2 ms, B-BC-2 note): CO₂ alternates 0 / 1e-3 ahead of the front,
      the clamp creates mass, and max |R T Σc/P − 1| reaches 3.2 in 120 s (1.7 at 5e-3 s, 0.46 at 1e-3 s), the same
      before and after C-BC-1 and with D = 0. The new per-output deviation line and its warning now show it for the
      continuity closures.
- [ ] **B-ITG-14** CVODE (SPGMR) at rtol 1e-2 stops on T15 at t = 187.6 s: "the corrector convergence test failed
      repeatedly or with |h| = hmin" (h = 2.2e-6 s; `benchmark/breakthrough`, 2026-10-28). With upwind and the
      Δz/v cap (2026-10-21) it completed; Dense at the same tolerance completes. Low priority: the tolerance is far
      looser than the default.
- [ ] **B-ITG-15** CVODE (default Dense, rtol 1e-5) with `ConvectionScheme` Upwind on
      `examples/BEA-alkanes-C7/breakthrough` took about 11 min of wall time to reach t = 290 s (stopped there;
      2026-11-02, binary of 2026-10-30), against about 2 s for the whole run with WENO3 (P-ITG-5). Step and
      Newton-failure counts not looked at; compare with B-ITG-11 (Upwind stalls) and B-ITG-12 (Dense).
- [x] **B-ITG-16 ROS34PW2 lets Σc drift by 17 % with a continuity closure and a stiff trace adsorbate.** Found
      2026-11-06 while setting up 3.4.4 (`validation/3.4.4_flow_reversal`): 0.1 % A in helium, Langmuir
      b p_F = 4, F q*'(0) = 1.5e5, `InletPressureInletVelocity`, isothermal, 100 cells. During the 30 000 s
      loading, max |R T Σc/P − 1| grows linearly to 0.17 (He 39.6–41.2 against P/RT = 40.3 mol/m³), and the
      closure turns it into ±2–8 % velocity errors. The same with the binaries of 2026-11-03 and 2026-11-04 (before
      flow reversal). In the cells ahead of the front the stage values of A flip around zero (output q* = 0 while
      dq/dt > 0); the clamp max(0, c) gives large spurious sinks, which the closure balances with helium. The
      W-method's Jacobian lacks that coupling (F-ITG-4), so the implicit stages damp A but not helium. rtol 1e-8
      makes it worse (0.24, 397 s), CVODE gives 4.2e-3 (342 s), k_LDF 0.05 1/s 0.031, `FixedVelocity` is not
      affected. Same mechanism as the frozen 8.2e-3 on the Stampi-Bombelli packed bed (F-ITG-4 note).
  - 2026-11-09, swing cycles: the drift adds up over the phases, since nothing pulls Σc back. H2 PSA of
    `examples/Self-driving-lab-H2-PSA`: +0.015 per cycle at rtol 1e-4 (1e-6: +0.002), while CVODE at 1e-4 to 1e-7
    stays at 3e-3 to 5e-3; the example and the robot of that directory use CVODE. Paper input: the 3.4.2 T15 TSA
    (ROS34PW2, 60 phases) reaches max |R T Σc/P − 1| = 0.040 (its log has the warning). Re-equilibrating the gas at
    every switch would hold it at 0.007 but moves the concentrations by up to 3.6 %; not done without a decision, as
    the swing driver rescales only at a pressure step, a temperature step without the energy balance, or a change of
    the end where the pressure is set (F-BC-2).
  - 2026-11-11 fixed in ROS34PW2 (`src/integrators/ros34pw2.h`): with `InletPressureInletVelocity` and
    `InletVelocityOutletPressure`, W acts as the identity on Σc + (Σc/T_g) T_g of each cell, so that the stage
    increment of Σc is that of the right-hand side less Σc/T_g times the change W makes to the gas-temperature
    increment; the correction goes to the components by their mole fractions at the step start. Results (binary
    before against after):
    - Trace case above (IPIV, WENO3, 100 cells, 30 000 s): max deviation 0.148 → 8.5e-5; the A profile at the end
      of loading is within 1.2e-4 c_F of the `FixedVelocity` run (5.5 % before); 3.3 → 0.24 s.
    - H2 PSA (4 cycles, rtol 1e-4): 2e-4 at every cycle start (before 4e-3; CVODE 2e-4 to 3e-4), inventories within
      1.4e-4 of CVODE, 0.43 s against 1.79 s for CVODE.
    - 3.4.2 T15 TSA (5 inputs): 0.019–0.040 → 3e-5–1e-4 (CVODE 3.6e-3 for 1000 s); phase inventories within 3e-4 of
      before; wall time halved (62 → 31 s for 1000 s). The paper text (Supplementary remarks of 3.4.2: −2.9 % to
      +4.0 %, and one cycle in 1 to 12 s) was made with the binary before.
    - Breakthroughs at the default tolerance (Marx ternary T1–T24, 3.4.1, 3.4.3): t50 within 2e-6, outlet within
      7e-4 of the feed fraction, 10–20 % faster. T15 benchmark against a CVODE rtol 1e-8 run of the same binary:
      e50 2.3e-3/5.3e-4/9.6e-5/1.5e-5/2.2e-6 → 7.3e-4/2.2e-4/4.3e-5/7.6e-6/1.3e-6 for rtol 1e-2 to 1e-6, with up
      to 22 % fewer steps.
    - Not with `InletPressureOutletPressure`: its march balances the flux of the pressure profile, and the
      projection held the dilute BEA breakthrough (3.4.1) at steps of 0.03 s instead of 1.6 s. IPOP keeps the old
      stages (byte-identical output). `FixedVelocity`, `FixedPressureInletVelocity`, and liquids are unchanged.
    - The stored T15 benchmark runs of 2026-11-03 give t50 about 1 s later (0.15–0.9 %) than the binary at the start
      of the 2026-11-09 session, with or without this fix: the closure fixes of V-BC-2, whose benchmark rerun is open.
    - Paper follow-up: V-ITG-2 (reruns) and I-ITG-3 (text).
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
- [ ] **B-MBD-1** Multibed gas–wall area comes from bed 0 only (`computeEnergyDerivatives`); the film and
      gas–solid areas blend 1/ε_a per bed while the loading prefactor uses ε_total (`computeBulkSpeciesSink`,
      `computeEnergyDerivatives`).
  - 2026-11-20: the gas–wall area per fluid volume is per node since B-MBD-3 (2026-11-18). Open: the 1/ε_a blend of
    the film and gas–solid areas against ε_total in the loading prefactor.
- [x] **B-MBD-2** Single-bed vs multibed drift: the multibed `computeChemisorptionTransportDerivatives` wrote
      `poreDot = −ρ·dq/dt` for every site, the single bed skipped non-transport sites; `computeReactionDerivatives`
      gets bed 0's components (capacity limits) for all beds (`rk3.cpp`, `cvode.cpp`).
      2026-10-27: the pore part is fixed with the single/multibed merge (R-GEN-1): a site without pore transport
      leaves its surface and pore reservoirs unchanged (test
      `MultiSiteChemisorption.KeepsDirectAndPoreTransportSiteBalancesIndependent`). The reaction part is open.
  - 2026-11-20: the reaction part was fixed with B-MBD-3 (2026-11-18): `computeReactionDerivatives` weights the
    capacities of the adsorbents at each node by their fractions.
- [x] **B-IO-12** `Column::readJSON` still accepts about 25 keys that `writeJSON` no longer writes (velocities, caches,
      `fractionOfAdsorbent`, `columnDistances`, …); an older `column.json` can therefore replace the grid and the
      adsorbent fractions of the input without a message (`hasAdsorbentOfType` is then not recomputed). In 2.0, read
      the state keys and reject the rest. Split from B-IO-8 on 2026-11-09.
  - Fixed 2026-11-18: `readJSON` reads the state keys (and `time`, B-IO-18) and leaves the others unread, so that an
    older file still loads its state. The two tests on `tests/data/bea_breakthrough_column.json`, an older file,
    load its derived fields themselves (`loadBeaBreakthroughColumn`).
- [x] **B-PY-2** Python: `fit_figure()` never draws the fit (trace commented out, `plot_fitting.py:467`);
      locale-dependent `read_text()` on non-ASCII presets and LICENSE (`config.py:221`, `:236`, `about.py:68`);
      `RunManager.start` reads `self.cli` twice outside the lock.
      2026-10-27: every `read_text`/`write_text` in `ruptura/webapp` passes `encoding="utf-8"`, and `RunManager.start`
      reads `self.cli` inside the lock. Open: the fit trace, commented out in `b8a2562`.
  - Fixed 2026-11-20: the "fit f(x)" trace is back in `fit_figure()` (start, fit, data on the BEA nC7 fit). It needed
    B-PY-9 for fits that start from zero parameters.
- [x] **B-PY-3** `tests/test_breakthrough.py`, `test_components.py`, `test_fitting.py`, and `test_mixtureprediction.py`
      import `Components` and other names that `ruptura/__init__.py` no longer exports (also on HEAD `633f2bd`), so
      pytest cannot collect them. Found 2026-10-27.
  - Fixed 2026-11-20 with B-BLD-2: pytest collects `*_test.py` only (`[tool.pytest.ini_options]` in
    `pyproject.toml`), the naming of every tracked Python test; `pytest tests` collects 133 tests, none of the v1
    files. The v1 files are left on disk.
- [ ] **B-BLD-1** Test inputs that exist only locally, so CI checkouts lack them: `examples/Multibed/*`,
      `examples/MPD/BEA-alkanes/*`, `examples/Liquid/*`, `examples/Chemisorption-CO2*/*`,
      `examples/BEA-alkanes-C7/{psa,tsa}`, `examples/BEA-alkanes-C7/fitting/Results.dat-BEA-Repeat-552K-nC7`.
  - 2026-11-18: the examples are staged (R-EX-1); with only the staged `examples/` the unit tests give the same
    224/225 as the full tree. Open: tests also read `marx2/ternary/T15` and `stampi-bombelli/dl_from_stampibombelli`,
    which are git-ignored study folders.
  - 2026-11-20: the two inputs the tests read from `validation/` are copied into `tests/data/`
    (`marx_t15_reference_settings.json`, the frozen reference-settings T15, and
    `stampi_bombelli_monolith_400ppm_1mmol.json`); the 7 paths in `pressure_swing_tests.cpp`, `ros34pw2_tests.cpp`, and
    `sorption_jacobian_tests.cpp` point there. 256/257 pass (the failure is B-EX-3). Open: the tests also read
    `examples/PSA-H2-activated-carbon-ports` and `examples/PSA-H2-activated-carbon-two-beds`, still untracked; stage
    them with the multicolumn work.
- [x] **B-BLD-2** `tests/test_*.py` are v1 tests (`from ruptura import Components`); `.gitignore` (`test_*`) hides
      them from git, but a local `pytest tests` still collects and fails them.
  - Fixed 2026-11-20 with B-PY-3.
- [ ] **B-EX-1** `examples/BEA-alkanes-C7/tsa` and `psa`: no `HeatOfAdsorption`, so the isotherms do not depend on
      temperature and the TSA regeneration desorbs nothing. The outlet peak at the switch is the gas
      inventory reacting to the instantaneous 552 → 673 K step.
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
- [x] **B-EX-6** Gas heat capacity still molar (29.51 or 35.8; see C-EX-1): `marx2/binary/B1–B36` and
      `marx2/base.json`, `benchmark/marx-403030_298_1e6_2e-1` (0.1 s matrices), the old `marx/` runs, and
      `pusan/`. None is in the paper; fix them before any reuse.
      2026-10-27: `marx2/heat_capacity.py` takes input paths and converts `base.json`, T1–T24 (unchanged), and B1–B36
      by default; also converted: the 315 inputs of `benchmark/marx-403030_298_1e6_2e-1` and the 122 Marx inputs of
      `benchmark/figures/simulations`, with their source `marx/403030_298_1e6_2e-1`. Open: the other old `marx/` runs
      and `pusan/` (35.8).
  - Fixed 2026-11-20: the 8 other `marx/` inputs with 29.51 converted by `heat_capacity.py` (1109.2 for 40/30/30,
    1014.1 for 15/80/5, as in `marx2`), and the 3 `pusan/` inputs by the same rule, 35.8 / M_mean with M_mean =
    0.00464 kg/mol (98 % He) → 7715.5 J/(kg K). The source of the 35.8 J/(mol K) is unknown (He alone is 20.8).
    Inputs of `marx/` without `heatCapacityGas` (other studies) are untouched; their results were not rerun.
- [ ] **B-EX-9 Elfving: decide the water model.** The 13 `elfving/` inputs were converted faithfully (B-PHY-2), so
      their empty-column Langmuir chemisorption site keeps the old VantHoff scaling with 50 700 J/mol and no T_ref
      (b ×7.6e8 at 298 K, B-PHY-10), while GAB now works. `elfving/generate.py` now writes what its comments describe
      (GAB with C(T) = C₀ exp(6600/(RT)), a constant empty-column site, no `HeatOfAdsorption` 1.0). Rerun
      `generate.py` (and the runs) to switch to that model; the old results were computed with B-PHY-10.
  - 2026-11-03: ROS34PW2 cannot run the stored inputs. Near p = 0 the B-PHY-10 site takes up water at a rate of about 1e11 s⁻¹
    in the gas balance, in a chemisorption block that J~ does not contain, so the steps fall to 3e-9 s (244 000
    for the first millisecond of run 3); RK3 at 1e-4 s gets through only because the clamp removes the overshoot.
    With the constant site of `generate.py`, run 3 takes ordinary steps (t = 3120 s of 18 000 s in 300 s). The
    ROS34PW2 rerun (`validation/ros34pw2_danckwerts_rerun`) leaves Elfving out until this is decided.
- [x] **B-EX-11 `examples/MPD/breakthrough` has no flow.** `FixedPressureInletVelocity` with `PressureGradient` 0
      gives v = 0 in every cell (local Ergun balance, B-BC-1), so the stored output (2026-10-31) reaches
      t = 82 820 s with zero CO₂ and N₂ at the outlet and `SimulationTime` auto never stops. Set a pressure gradient
      or another boundary condition. Found 2026-11-02.
  - Fixed 2026-11-20: `InletPressureInletVelocity` (the inlet velocity of the input, 1e-3 m/s), `PressureGradient`
    removed. Checked to t = 1000 s: v_out 8.1e-4 m/s, max |R T Σc/P − 1| 3e-6, CO₂ and N₂ at the outlet at 7.9e-4
    and 8.6e-4 of their feed fractions and rising. The full breakthrough was not run: RK3 at the input's 1e-3 s step covers
    about 1 s per second of wall time.
- [x] **B-EX-15** The Ruptura Lab preset `ruptura/webapp/presets/paper-t15.json` matches neither T15 input of the
      paper (2026-11-19): it has ρ_s = 502 kg/m³ (pre-2026-10-26) and IAST, whereas the Marx validation uses the Marx
      inputs and the other T15 studies the reference settings (ρ_s = 456.6). Regenerate it from one of them (the
      Marx inputs need F-GUI-3) and fix the source path in `presets/README.md`.
  - 2026-11-20: regenerated from `marx2/ternary/T15` (Marx inputs) with the preset conventions (CVODE SPGMR, 2 s
    output, fractions); CO2 t50 541.3 s as the validation run; README updated. The paper's Lab screenshots were retaken
    from it (`figures/gui/capture_gui.py --cli`).
- [x] **B-PHY-12 Avrami chemisorption runs on the global simulation clock.** `computeChemisorption` passes the
      absolute time of the stage to `Chemisorption::rate`, so every node and every swing phase uses k^n t^(n−1) with
      t counted from the start of the run. The SI (Eq. `avrami_chemisorption`) defines τ as "the elapsed time
      associated with the active adsorption or desorption episode". A node that the front reaches at t_arr is
      evaluated at t instead of t − t_arr (for n < 1 its rate is too small by (t/(t − t_arr))^(n−1), for n > 1 too
      large); the clock never resets across cycles, and a restart from `ReadColumnFile` resets it to 0 (B-IO-8).
      Suggestion: the autonomous form of the same law, dq/dt = n k (q* − q) [−ln(1 − q/q*)]^((n−1)/n), which needs no
      clock and integrates to 1 − exp(−(k t)^n) at constant q*. No input in the repository uses Avrami. Found
      2026-11-06 (code review).
  - 2026-11-09: the user agrees that the rate must not depend on the global clock; the form (the autonomous
    Avrami–Erofeev law, and what it does for desorption and for n > 1 at q = 0) is still to be decided.
  - Fixed 2026-11-09: the clock-free (Avrami–Erofeev) law, dq/dt = n k s^(n−1) (q* − q) with the equivalent Avrami
    time s = [−ln(1 − q/q*)]^(1/n) floored at 1e-3. Release (q > q*) is linear driving force, k (q* − q) (user: no
    time-based description for desorption); n = 1 is first order on both sides. From a clean sorbent at constant q*
    the loading follows the Avrami curve moved in time by (1/n − 1) 1e-3/k (tests `ChemisorptionRate.Avrami*`).
    The floor changes only the rate, which the gas and the solid share, so mass is conserved: CO2/He column with one
    Avrami site (k = 0.02 1/s, 40 cells, CVODE), Σ h (ε dc/dt + (1 − ε) ρ_p dq/dt) − (in − out) is 7e-16 of the
    inflow, as before. Outlet t10/t50/t90 for n = 2: 45.5/142.0/247.7 s → 41.3/140.5/253.8 s; for n = 0.6:
    32.9/111.2/283.0 s → 39.3/117.7/274.5 s (RK3 identical). A floor of 1e-4 instead moves the outlet by at most
    5e-5 c_feed (n = 2) and 3e-3 c_feed (n = 0.6). `Chemisorption::rate` (and the Python `rate`, `elapsed_time`
    dropped) and `computeChemisorption` no longer take the time, and (user, 2026-11-09) neither do
    `computeDerivatives`, `RungeKutta3::step`, and the CVODE right-hand side. SI Eq. `avrami_chemisorption`
    rewritten.
- [x] **B-PHY-13** `MixturePrediction::predictMixture` rejects |Σy − 1| > 1e-15 (`mixture_prediction.cpp:398`)
      after its callers have normalized the mole fractions. The rounding of that normalization and of the sum grows
      with the number of components (bound ≈ n × 1.1e-16), so mixtures of a few tens of components can stop
      sporadically with "sum idealGasMolFractions at IAST start not unity". A tolerance of about 1e-12 (or n × 1e-15)
      removes the risk. Found 2026-11-06 (code review); not seen in a run.
  - Fixed 2026-11-09 (user's choice): the check stays, because Python users call `predict_mixture` with their own
    fractions, but at 1e-12: rounding stays below it up to about 9000 components, and a real mistake (Σy = 0.999)
    is still caught. The message names the sum ("the gas-phase mole fractions sum to 0.999 instead of 1"). Test
    `MixturePrediction.AcceptsRoundingInTheMoleFractionSumAndRejectsAWrongSum`.
- [x] **B-ITG-1** SIRK3: fix later. It cannot run the benchmark cases and supports a single adsorbent only. See also
      the banded-matrix layout and the missing temperature/chemisorption update (B-ITG-4) and P-ITG-6.
      2026-10-27: B-ITG-4 and P-ITG-6 are fixed; SIRK3 is limited to one adsorbent, upwind convection, and
      physisorption without reactions, chemisorption, or the energy balance.
      2026-11-02 check (binary of 2026-10-30 15:20, all Upwind, outlet against RK3): SIRK3 converges to the RK3
      solution on BEA C7 (`FixedPressureInletVelocity`, `InletPressureOutletPressure`, `FixedVelocity`,
      D_ax = 1e-5 on every component), MOR-CO2-C3H8 (SIAST, 25 bar), JUC-77-xylenes, and Liquid/pH-dependent;
      Σc_i RT/P matches RK3 to the printed digits. All 105 `benchmark/bea-alkanes-C7/*/sirk3` runs finish (the
      README's "not currently runnable" and `--include-sirk3` are stale). Open:
      (1) first order in time on BEA (max |Δc|/c_feed 4.6e-4, 1.4e-4, 4.6e-5, 2.3e-5, 1.2e-5 at dt = 1e-2 …
      2.5e-4 s), against 3e-6 for RK3 at 1e-2 s; about second order on JUC-77 (6.0e-3, 1.1e-3, 2.7e-4 at 1e-3,
      5e-4, 2.5e-4 s against RK3 at 1.25e-4 s). (2) The stability gain over RK3 is about 3× only. SIRK3 is
      implicit in convection, dispersion, and the LDF rate k at fixed q*, but the stiff mode is the gas depletion
      λ = k (1 + ρ_b/ε dq*/dc) through q*(c) in the sink, which both integrators treat explicitly. λ = 45 s⁻¹ for
      BEA nC7 and 30 s⁻¹ for the liquid case (Henry slopes); on JUC-77 the strong Langmuir and Langmuir–Freundlich
      sites (ν < 1) give λΔt = 10–33 at Δt = 3e-3 s in the cells below 1e-3 c_feed ahead of each front, and 2–5.5
      even at the default 5e-4 s. JUC-77: RK3 diverges at 3e-3 s; SIRK3 there puts the ethylbenzene roll-up front
      0.6 s early (t50 −0.26 %; 12 % of c_feed pointwise because the front is steep) and diverges at 1e-2 s.
      Liquid: RK3 diverges at 0.1 s (λΔt = 3.0 > 2.51); SIRK3 there has t50 −3.3 % (11 % of c_feed pointwise)
      and diverges at 0.3 s. BEA: SIRK3 is 1.6 % off at 3e-2 s (λΔt = 1.35) and diverges at 0.1 s (4.5).
      Beyond the limit it returns 0 (B-ITG-13). (3) A step costs about twice an RK3
      step (8 runs in parallel, so rough). So SIRK3 is correct in its scope, but never better than RK3 at equal
      accuracy in these cases.
      2026-11-02, final stage: for a linear operator, stages 1–3 give 1/3 + r/2 + r³/6 with r = 1/(1 − z), and the
      final stage divides by 1 + z²: second order, |R| ≤ 1 and R → 0 on the negative real axis, but R is singular at
      z = ±i. On the upwind symbol, |R| exceeds 1 above Courant 1 (1.44 at 2, 6.2 at 10); on a 100-cell upwind
      column with a fixed inlet, max_n ‖Rⁿ‖ is 82 at Courant 2, 15 at 5, and 3 at 10, whereas stages 1–3 alone stay
      below 1 at every Courant number. The final stage, not the sorption, limits the convective step of SIRK3 to
      Courant ≈ 1; a coupled sorption Jacobian (F-ITG-3) does not remove that limit.
      2026-11-03: resolved by removal. SIRK3 (`rk3_si.{h,cpp}`), its input value, bindings, app and web entries,
      tests, example (`examples/BEA-alkanes-C7/sirk3`, now `ros34pw2`), and paper text are gone; ROS34PW2
      (F-ITG-3) is the implicit alternative to RK3.
- [x] **B-ITG-18 An exception in a trial-state right-hand side ends the run instead of shortening the step.**
      `updateVelocityAndPressure` throws on a negative outlet pressure or a failed velocity bracket, and the mixture
      prediction on a failed inversion. `CVODE::evaluateDerivatives` does not catch: the exception unwinds through
      the SUNDIALS C frames (undefined behaviour) and CVODE never gets the chance to retry, whereas a positive return
      value ("recoverable") makes it cut the step. `setupPreconditioner` and `evaluateJacobian` already catch and
      return −1. ROS34PW2 likewise lets a throw from a stage evaluation end the run, although it repeats steps whose
      factorization fails. Fix: catch `std::runtime_error` in the CVODE right-hand side and return 1; in ROS34PW2 treat
      it as a rejected step (errorNorm = ∞, which already restores the step start). A state at an accepted step that
      throws still ends the run. Found 2026-11-06 (code review).
  - Fixed 2026-11-09: `CVODE::evaluateDerivatives` catches `std::runtime_error` and returns 1 (recoverable), any other
    exception −1, and keeps the message (`UserData::lastError`), which the CVODE failure message adds to the flag name
    (`CVodeGetReturnFlagName`); `setupPreconditioner` and `evaluateJacobian` pass that flag on and treat their own
    runtime errors as recoverable too. ROS34PW2 rejects a step whose stage evaluation throws (next step × 0.25) and
    names the last stage error if the step size collapses. Not yet checked by a run that throws (e.g. the B-BC-1
    `FixedVelocity` Marx runs).
- [x] **B-ITG-19** CVODE stops with CV_TOO_MUCH_WORK after 100 000 internal steps within one `CVode` call
      (`CVodeSetMaxNumSteps`, `cvode.cpp:233`), and that call spans one output interval. A run without any stall can
      reach the limit when the interval is long against the automatic 4 Δz/v step cap (e.g. Δz = 1 mm, v = 0.1 m/s:
      cap 0.04 s, 100 000 steps = 4000 s). CV_TOO_MUCH_WORK leaves CVODE usable, so `advanceCVODE` can call `CVode`
      again until it reaches endTime and keep a separate guard against a true stall (B-ITG-11), e.g. on the step
      size. The thrown message ("CVode failed during propagation") does not name the flag; add
      `CVodeGetReturnFlagName`. Found 2026-11-06 (code review).
  - 2026-11-09: the failure message names the flag now (with B-ITG-18); the loop over CV_TOO_MUCH_WORK is open.
  - Fixed 2026-11-09: `advanceCVODE` continues a call that reaches the 100 000 steps (`maximumStepsPerCall`), so the
    output interval no longer limits the run. It stops when those steps advanced the time by less than 1e-3 of the
    time in which the fluid crosses the smallest cell (or of the step cap, if smaller), with a message that names the
    mean step and the remedies (`CVODEMaximumOrder` 1, ROS34PW2); without flow only CVODE's own failures stop the
    run. The SUNDIALS error handler no longer logs the step limit. Only a call that failed before reaches this check,
    so no run that finished before stops now. The two crawls of B-ITG-11 (Stampi-Bombelli fitted D_L, Danckwerts,
    Dense, default order; about 1 s per 100 000 steps) run at 1.8e-5 (monolith 0.02 mmol/s, t = 5.6 s) and 3.9e-4
    (packed bed 0.5 mmol/s, t = 151 s) of the crossing time and stop where they stopped before; the 164 stored CVODE
    runs with at least 10 000 steps average 1.2e-2 or more. Test
    `CVODEIntegrator.ContinuesACallThatReachesTheStepLimit` (120 000 steps in one call; the old code stopped at
    100 000).
- [x] **B-BC-2** Velocity closure without thermal expansion: `InletPressureInletVelocity` and
      `InletVelocityOutletPressure` conserve the molar flux with c_tot from the transported
      concentrations, with no ∂c_tot/∂t term. With the energy balance, c_tot R T_g departs from the Ergun
      pressure: −14% to +23% in the 3.4.2 TSA runs. The equilibrium (y_j P_Ergun) is not affected, but the
      gas inventory is.
  - Fixed 2026-11-04 (cycles plan, stage 1): with the energy balance the three continuity closures carry
    −(c_i/T_i) ∂T_i/∂t (quasi-static pressure), with ∂T_i/∂t from `GasTemperatureBalance::rate`, which
    `computeEnergyDerivatives` now uses for the gas temperature too. dT/dt is linear in v_i, with coefficients
    computed once per closure; node 1 with a Danckwerts inlet is solved by secant steps. Test
    `ContinuityClosure.CarriesTheThermalExpansionOfTheGas`. SI Eq. `gas_total_molar_balance` and main text updated.
    Impact (before → after, max |R T Σc/P − 1|; t50; wall time sequential):
    | Case | Σc deviation | t50 CO2 / N2 / H2 | cost |
    | Marx T15 (IPIV) | 9.0 % → 0.05 % | +0.13 / −0.47 / −0.85 % | +5 % (steps +7 %) |
    | Marx T20, T21 (318 K feed) | 12.5, 14.3 % → 6.7 % (the initial offset of C-BC-2) | ≤ 0.74 % | +8 % |
    | 3.4.1 T15 IPOP / IVOP | 22.7 / 8.9 % → 0.02 / 0.15 % | ≤ 0.93 % | +15 % IPOP |
    | 3.4.2 TSA 1000 s, 60 phases | 33 % → 4.0 % | CO2 loading at step ends ≤ 5e-5 | +7 % |
    | Dantas run01 (RK3) | 7.5 % → 0.01 % | +0.11 % | |
    Outlet temperatures change by ≤ 0.04 %, isothermal runs are byte-identical (BEA C7). The 4 % in the TSA is
    ROS34PW2 integration error on the nonlinear invariant Σc T: 30–50× smaller at relative tolerance 1e-8.
    Stage 3 can re-equilibrate the inventory at each phase switch. The paper runs with the energy balance and a
    continuity closure were rerun on 2026-11-04 together with C-BC-2 (results there), except the benchmark (V-BC-2).
  - Also isothermal (2026-10-30), from the two-pass coupling in the CVODE right-hand side (B-ITG-7):
    `InletPressureInletVelocity`, energy balance off, a concentrated feed stepped into a He-filled packed bed.
    After the front Σc_i R T stays below P (and above it after a purge), with v_out/v_in off by the same factor:
    | LDF k [1/s] | after step-up | after purge |  (10 % A + 10 % B at 10 bar, Ergun drop 110 Pa)
    | 1 | −0.2 % | +0.9 % |
    | 10 | −2.1 % | +8.9 % |  (same for 100, 200, 400 grid points, and with Upwind instead of WENO3)
    | 100 | −23 % | +70 % |
    BEA C7 (7 × 10 % at 1 bar, k = 10): −8.1 %; 1 % A + 1 % B at 100 bar: −0.004 %.
    Cause: pass 1 starts from zero uptake (`N_VConst(0.0, ydot)`, `cvode.cpp:300`), so the Ergun pressure P¹
    is that of a column without uptake; pass 2 builds the velocity from the pass-1 sink S¹ = P_f k (q*(y, P¹) − q)
    but the species balances get S² = P_f k (q*(y, P²) − q). P¹ − P² is a few Pa downstream of a front (it
    grows with distance), yet P_f k δq* makes Σ_i dc_i/dt = S¹ − S² ≈ 0.03 (k = 10) to 0.4 (k = 100)
    mol/m³/s, and nothing pulls Σc_i back to P/RT. Check: k = 100 with the Ergun drop switched off (μ → 0,
    d_p = 1 m) gives −0.03 % instead of −23 %.
    [x] CVODE fixed 2026-10-30: after the last pass the velocity is rebuilt from the final sink (one more
    `updateVelocityAndPressure`, no new equilibrium), so closure and species balances use the same S. Drift
    after the fix ≤ 1.1e-4 in every case above (k = 100: −23 % → +1e-4, purge +70 % → +8e-5); outlet mole
    fractions change by up to 1 % of the feed value at k = 10 and 12 % at k = 100. Cost: 1–5 % more time per
    right-hand side with IAST (14 % for binary Langmuir at k = 100; 0.023 → 0.034 ms for the 1-component,
    20-point Chemisorption-CO2 example); wall time −6 % to +7 % over nine cases, within the step-count scatter.
    [x] `RungeKutta3`, `SIRK3` (2026-10-30): `computeDerivatives` (`rk3.cpp`) rebuilds the velocity from the
    stage's own sink; before, `precompute` built it from the previous stage's sink. That lag never caused a
    drift (+1.1e-4 before and after, outlet change ≤ 1e-5), because RK3 does not start from zero uptake; it
    costs +1 to +9 % wall time (JUC-77 +1 %, Chemisorption-CO2-N2 +9 %, Langmuir base +9 %).
    `SIRK3` also had a non-conservative convection operator, −v_i (c_i − c_{i−1})/dz, which drops c ∂v/∂z;
    now −(v_i c_i − v_{i−1} c_{i−1})/dz, so its implicit solves keep Σc_i exactly. Against CVODE with the same
    Upwind scheme: BEA-alkanes-C7/sirk3 (300 s) max outlet error 2.2 % → 0.23 % of the feed, Langmuir base
    (dt 1 ms) RMS 1.4 % → 1.0 %; wall −3 % to +6 %. The earlier "negative outlet pressure" aborts and
    blow-ups of RK3/SIRK3 on the base case were explicit-sink instabilities: P_f k dq*/dc ≈ 1100 1/s at the
    clean front edge needs dt ≲ 2 ms; at dt = 1 ms both run. The examples keep a constant offset
    Σc_i R T/P − 1 = ΔP_Ergun/P from t = 0 (+3e-3 in MOR-CO2-C3H8 at 25 bar): the initial concentrations
    do not follow the Ergun profile (no ∂c_tot/∂t term).
  - Axial dispersion on the adsorbing components but not on the carrier (He without
    `AxialDispersionCoefficient`) adds Σ_i D_i ∂²c_i/∂z² ≠ 0, which the closure ignores: Σc_i R T/P swings from
    1.51 to 0.59 at D = 1e-3 m²/s (Pe 100). Same suspicion as in C-BC-1.
    [x] 2026-11-03: resolved by C-BC-1 (the closure carries the net dispersive flux again).
- [x] **B-BC-4 The Dirichlet inlet admits a dispersive inflow on top of the feed.** The inlet face holds c_in and
      T_in (node 0 of `computeMassDerivatives` and `gasTemperatureDot[0] = 0`, `src/integrators/compute.cpp`), so
      the first cell also receives D (c_in − c₁)/(h/2) (I-BC-1). Up to t50 this adds about D/(vL) of the feed (linear
      isotherm, long column: R c_in D/v against R c_in L). Checked on 2026-11-03 against the stored column profiles
      of `stampi-bombelli/` (bounds from the 12 s output interval; the estimate lies inside them):
  - literature D_L: packed bed 0.5–0.8 % (vL/D 140–190); 400 ppm monolith at 1–3 mmol/s 0.2–0.4 % (300–450), at
    0.18 mmol/s 1.7–1.9 % (60); 5.6 % monolith at 0.18 mmol/s 0.1–6 % (60), at 0.06 mmol/s 2–9 % (20), at
    0.02 mmol/s 15–22 % (7);
  - fitted D_L (`dl_from_stampibombelli`): monolith about 0.7–6 times the feed (vL/D 0.8–2), packed bed 200–1000
    times (0.04). The integrator-benchmark monolith (vL/D ≈ 2) is in this range;
  - by the estimate below 1 % for Marx (vL/D ≈ 3300) and Stadler (120–490); BEA has no physical dispersion.
  - The paper now gives these numbers (SI species boundary conditions). The 5.6 % monolith breaks through 8–13 % early
    (main text, Monoliths), the same size as the excess at 0.02 mmol/s.
  - Fix: a Danckwerts inlet, i.e. no dispersive (or conductive) flux across the inlet face, in the WENO3 and Upwind
    paths for species and gas temperature. Not decided (2026-11-03).
  - 2026-11-03 (user request): switch `InletCondition` ("Dirichlet" default, "Danckwerts"). With Danckwerts, the
    stencils of the first cell use the inlet-face value c_b = (u c_in + g c₁)/(u + g), g = D/δ (δ = h₁/2 WENO3,
    (h₁ + h₂)/2 upwind), so the inflow is exactly u c_in; likewise for the gas temperature (upwind stencil, a = λ/(ρ
    C_p)). In `computeMassDerivatives`, `computeEnergyDerivatives`, and the node-1 entries of the upwind sorption
    Jacobian (`jacobian.cpp`, so KLU stays exact for upwind); RK3 and CVODE only (the reader rejects it with SIRK3,
    ROS3, and ROS34PW2). Tests: `InletConditionInput.*`, `DanckwertsInlet.*` (`tests/convection_tests.cpp`). Qt app,
    Ruptura Lab, AI prompt doc; SI Eqs. `danckwerts_species_inlet` and `danckwerts_face_value`, JSON enumeration
    table. Reruns of the validation cases with both inlets: `validation/inlet_condition_dirichlet_vs_danckwerts`.
    Open: which inlet the paper uses, and the default.
  - 2026-11-03, reruns (README and `comparison.csv` there): the Dirichlet copies reproduce every stored result byte
    for byte.
    - Danckwerts against Dirichlet, t50: Marx, 3.4.1, and 3.4.2 ≤ 0.06 %; Stadler ≤ 0.3 % (WGS conversion);
      T15 at D_ax = 1e-3 m²/s −0.24 % instead of −2.9 % against D = 0.
    - Stampi-Bombelli with the literature D_L, against the experiment: 5.6 % monolith −11/−8/−13 % → +6/−2/−11 %;
      400 ppm monolith and packed bed unchanged (+6 to +25 %).
    - With the fitted D_L, Danckwerts brings t50 to within −15 to +24 % of the experiment and the 10–90 % widths to
      60–133 % of the measured ones; Dirichlet gives −29 to −72 % in the monolith and immediate breakthrough in the
      packed bed.
    - New stops with Danckwerts: fitted-D_L packed bed at 0.5 and 1 mmol/s (`mxstep`, B-ITG-11; 1 mmol/s finishes at
      order 1). Paper consequences: I-BC-4.
  - Done 2026-11-03 (decision): Danckwerts is the default (`src/inputreader.h`, `column.h`, the `compute.h` default
    arguments, Qt app, Ruptura Lab, field help, AI prompt doc). The paper's method text (main text and SI) describes
    Danckwerts as the default and Dirichlet as the option. The liquid regression tests pin Dirichlet
    (`fluid_phase_tests.cpp`). The benchmark generator pins Dirichlet (ruptura-5a), so its inputs keep reproducing
    the benchmark figure.
- [x] **B-CYC-1** Phase driver: a step that raises the equilibrium loading in the whole column at once aborts. An
      isothermal TSA (energy balance off, so T_m is imposed everywhere) stops at the 373 → 298 K step: the
      uptake exceeds the feed, v(L) = −0.63 m/s, "negative outlet pressure" (3.4.2). A PSA
      repressurization step would do the same. Needs flow reversal (F-BC-1) or a one-way outlet.
      Fixed 2026-11-04 with F-BC-1: the gas of the last cell flows back in at z = L. The isothermal T15 TSA
      (298 → 373 → 298 K, 1000 s each) runs through: at the switch v(L) = −0.60 m/s with faces 5–100 negative, and
      20 s later the flow is forward everywhere again (test
      `FlowReversal.BackFlowAtTheProductEndKeepsTheTotalConcentration`). The inventory jump at the switch remains
      (F-CYC-2). A repressurization against a closed end needs closed ends (stage 3).
- [x] **B-CYC-2 Swing phases with `InletVelocityOutletPressure` feed gas at the inlet pressure into a column at the
      outlet pressure.** `SwingAdsorption` requires a positive inlet pressure for every phase (the phase
      `InletPressure`, else the top-level one; the constructor throws otherwise), and `applyStageConditions` sets the
      feed node to y P_in/(R T). With IVOP (and `FixedVelocity` given only an `OutletPressure`) the column pressure
      follows from `OutletPressure`, so for P_in ≠ P_out the feed carries P_in/P_out times the molar flow of
      v_in P_out/(R T), and the continuity closure jumps at node 1. The breakthrough driver initializes the IVOP feed
      at P_out and is consistent. No paper input is affected (Stadler uses `FixedVelocity` with `InletPressure`,
      3.4.2 and 3.4.4 IPIV). Fix: take the feed pressure from the boundary condition (the pressure at the feed end
      after precompute), and let a phase set `OutletPressure` (F-CYC-1). Found 2026-11-06 (code review).
  - Closed 2026-11-09 without a fix (user: ignore).
- [x] **B-FIT-4** `Isotherm::randomize()` (`isotherm.cpp:266`) never draws b < 1 (`pow(U, 20(U′−1)) ≥ 1`), so every
      initial GA citizen is saturated at Pa-scale pressures. The NM stop test (`fitting.cpp:965`) checks the mean
      RMSE < 1e-4 rather than the simplex size, so it always runs `MAX_IT`.
      2026-11-11: the line numbers have moved (`randomize` is now `isotherm.cpp:438`, the NM stop test
      `fitting.cpp:931`). Despite the b ≥ 1 start, a scratch build of the working tree reaches the global optimum of
      the absolute RMSE for a single Langmuir and a single Langmuir–Freundlich on CoBDP p-xylene in 3/3 runs each
      (agrees with scipy differential evolution to 9 digits; GA about 115 generations, 4–5 s with the NM `MAX_IT`).
  - Fixed 2026-11-11: the 11 affinity-like draws are U^(20(1−U′)) ∈ (0, 1] (half of them in 1e-7…1e-1). The NM
    simplex starts with steps relative to each parameter (the absolute step of about 1 moved b ~ 1e-5 by 1e5 times
    its value) and stops when every vertex lies within 1e-8 (relative) of the best one in every parameter. Before and
    after on the 15 BEA data files (Langmuir, Toth, 2×Langmuir) and CoBDP p-xylene (Langmuir, LF, 2×LF), 3 runs each:
    never worse, a better optimum in 7 of 48 cases (2×LF CoBDP RMSE −30 %, Toth C3 −22 %, 2×Langmuir 3eC5 and nC5
    −3 %, Toth 2mC3 −1.4 %), and the 3 runs now agree where the old code scattered (2×Langmuir 2mC6 0.0083–0.021
    mol/kg). NM needs 128 steps (median); it still runs to `MAX_IT` when parameters are not identifiable (two
    Langmuir sites with the same b, a site with q_sat → 0) and for LF on CoBDP p-xylene (b and n strongly
    correlated), where the fitness has converged.
- [x] **B-FIT-6** `Fitting::writeCitizen` (`fitting.cpp:624`) prints `char info[256]` with `std::print("{}", info)`,
      which writes all 256 bytes, so every GA line in the fitting log is followed by uninitialized bytes and grep
      treats `run.log` as binary. Print `std::string_view(info)` or format the line with `std::print` directly.
  - Fixed 2026-11-11 (`std::string_view(info)`): 0 of 288 fitting logs contain control bytes (before: all of them).
- [x] **B-FIT-7** `Fitting::fitness` divides the sum of squares by m − p in `size_t`: with no more data points than
      parameters it wraps around or divides by zero instead of stopping. Fixed 2026-11-11: `readData` throws.
- [x] **B-FIT-8** `Fitting::mate` leaves child E + M/2 out of both crossover ranges (`+ 1` at the start of the second),
      so one citizen per generation is a leftover of the previous one. Fixed 2026-11-11.
- [x] **B-IO-7** A MixturePrediction run cannot write SPI-only sites: `MixturePrediction::run` evaluates
      `psiForPressure` for every method but MPD (`mixture_prediction.cpp:2156`), and GAB throws there
      (`isotherm.h:778`). `liu/` works around it with a 600 s water breakthrough for Fig. 5a. The humidity isotherms of
      F-PHY-1 will hit the same issue. Found 2026-11-04; not checked by a run.
      Fixed 2026-11-04: column 7 is zero for a component with a site that defines no ψ
      (`Isotherm::definesGrandPotential`: pH-Langmuir, GAB, the humidity types); a GAB + GAB-coadsorption run writes.
- [x] **B-IO-8 `ReadColumnFile` restarts the clock at t = 0 but appends to the old output files.** `column.json` holds
      the state only (`Column::writeJSON`), and `Breakthrough` starts at time 0 with output index 0, then opens the
      `.data` files in append mode without a header. The appended block repeats the times of the first run, the
      first frame duplicates the last one, and the time-dependent parts start over (Avrami clock B-PHY-12, the
      `SimulationTime` budget, the auto-stop rounding). Decide what the key means: an initial state (truncate the
      files and write the header, as a fresh run), or a resume (store the time, the output index, and in a swing run
      the phase index in `column.json` and continue from them). While there, `Column::readJSON` still accepts about
      25 keys that `writeJSON` no longer writes (velocities, caches, `fractionOfAdsorbent`, `columnDistances`, …); an
      older file can therefore replace the grid and the adsorbent fractions of the input without a message
      (`hasAdsorbentOfType` is then not recomputed). In 2.0, read the state keys and reject the rest. Found
      2026-11-06 (code review).
  - Closed 2026-11-09 without a change: a `ReadColumnFile` state starts a run at t = 0, which is the intended
    behaviour (user). The legacy `Column::readJSON` keys are now B-IO-12.
  - 2026-11-18: changed after all (user): a restart continues the clock, and `ResetClock` gives the start at t = 0
    (B-IO-18).
- [x] **B-IO-9 The dimensionless time in the output files uses the input `ColumnEntranceVelocity`.** Column 1 of
      `column.data` and of every component file is t × `timeNormalizationFactor`, which `ColumnData` sets once from the
      input velocity. `InletPressureOutletPressure` and `FixedPressureInletVelocity` do not require the key and
      compute v_in themselves, so without it τ = −t/L (default −1 m/s), and with it τ uses a velocity the run does
      not have; per-phase velocities of a swing run are ignored as well. The Python `compute()` result takes the
      current computed velocity instead, so file and in-memory τ differ. For 2.0, consider dropping τ and writing the
      time in s rather than min (every input time is in s, and so is the progress output). Found 2026-11-06 (code
      review); the stored IPOP/FPIV inputs all set the key.
  - Fixed 2026-11-09 with a warning (user's choice): when `ColumnEntranceVelocity` is not set (possible only with
    `InletPressureOutletPressure` and `FixedPressureInletVelocity`), the reader warns that τ has no reference velocity.
    The negative τ came from the "not set" value of the key, −1 m/s, in the factor v_in/L: τ = −t/L. τ is otherwise
    unchanged; the per-phase velocities of a swing run still do not enter it. No stored input triggers the warning.
- [x] **B-IO-10** `"Components"` as a JSON object is accepted with the comment "component IDs follow JSON insertion
      order" (`inputreader.cpp:2059`), but `readJson` parses into `nlohmann::json`, whose objects are sorted by key.
      The component indices, the carrier-gas index, and the `component_<i>_<name>.data` names then follow the
      alphabetical order of the names. No input in the repository uses the object form; reject it (or parse with
      `ordered_json`). Found 2026-11-06 (code review).
  - Fixed 2026-11-09: the reader rejects an object ("'Components' must be a list of component objects"); test
    `ComponentsInput.RejectsAnObjectOfComponents`; the AI prompt doc says so. The Qt app (`LabInput::canonicalize`,
    alphabetical through QJsonObject) and the Lab (`config.py`, file order) still turn an object into a list on import.
- [x] **B-IO-11 Silent defaults behind "not set" checks that cannot fire.** `Temperature` defaults to 433 K, so
      "Error: temperature not set" (`temperature >= 0`) is unreachable; likewise `ColumnVoidFraction` 0.4,
      `ParticleDensity` 1000, and `ColumnLength` 0.3 against their "not set or invalid" messages. With
      `energyBalance`, `heatCapacityGas` and `heatCapacitySolid` default to 1 J/(kg K) and are not checked, and the
      check "requires heatCapacityWall > 0" cannot fire because its default is 1. Impact of making them required: none
      of the 2294 breakthrough and swing inputs in the repository omits `Temperature`, and all 1090 with the energy
      balance set the three heat capacities; 7 omit `ColumnLength`: six multibed inputs that take it from their
      adsorbent lengths, and the incomplete `pusan/susanne/1e5P_298K` (no pressures or velocity either). Suggestion
      for 2.0: no default for physical properties of the column, the bed, or the energy balance. Found 2026-11-06
      (code review).
  - Fixed 2026-11-09 with warnings (user's choice; no default removed): `Temperature`, `ColumnVoidFraction` (packed
    bed), `ParticleDensity`, `ColumnLength`, and with the energy balance `heatCapacityGas` (`heatCapacityLiquid`),
    `heatCapacitySolid`, and `heatCapacityWall` print "Warning: <key> is not set; using the default …" unless the top
    level, the `Geometry`, every `Adsorbents` entry, or other input (the length of a stacked multibed, the temperature
    of the first swing phase) sets them; `ColumnEntranceVelocity` see B-IO-9. None of the 2339 stored breakthrough and
    swing inputs triggers one, and `--validate` gives the same result for all of them as before. Not covered:
    `DynamicViscosity` (1e-5 Pa s) and `ParticleDiameter` (1 mm), which 261 and 260 stored inputs (the
    `benchmark/multibed` matrix) leave at their defaults.
- [x] **B-PY-4** The default `analysis.ipynb` that every CLI run writes (`notebook_templates.cpp`, `setupCell`)
      hard-codes the absolute run directory and the build's source tree (`RUPTURA_SOURCE_DIR`), which it prepends to
      `sys.path`. A run directory that is moved or shared keeps reading the old location, and the new example
      notebooks in `examples/BEA-alkanes-C7/` contain the absolute path of the source tree they were made in. Use the
      notebook's own directory (`Path.cwd()` under Jupyter) and the installed package. Found 2026-11-06 (code
      review).
  - Fixed 2026-11-09: the notebook holds no paths: `run_dir = Path.cwd()` (Jupyter starts the kernel in the folder
    of the notebook) and the installed `ruptura` package; the template functions and `writeDefaultNotebook` lost their
    path arguments; test `NotebookTemplates.HoldNoPaths`. The CLI no longer compiles in `RUPTURA_SOURCE_DIR`, and the
    Qt runner no longer writes it as `repoRoot` into the `manifest.json` of every run. Left: notebooks written before
    keep their paths (the 110 `analysis.ipynb` under `examples/` are untracked outputs that `clean_example_outputs.py`
    removes; no tracked file holds a user path); `ruptura_lab` (`webmain.cpp`) falls back to the compiled-in source
    tree to find the Python package; the Lab and Qt imports turn a relative `ReadColumnFile` or `FileName` into an
    absolute path in the `simulation.json` of the run.
- [x] **B-EX-10** `liu/simulations/figure_4_rh_0p20` (CVODE) stops with `mxstep` before t = 300 s, with the binary
      of 2026-10-30 before the temperature laws (t = 230 s) and after (t = 52 s). RK3 (Δt = 5e-3 s) runs, and old
      and new agree to 5e-12. Check the other Liu runs and the CVODE settings.
      2026-11-03: 36 of the 50 stored Liu runs (CVODE, Dense default) have no "Final time:". ROS34PW2 at its default
      settings (Danckwerts) runs all 50 in 3.1 s in total (`validation/ros34pw2_danckwerts_rerun`, ruptura-5a). For the
      14 stored runs that finished (Dirichlet), t50 differs by at most 0.6 % and y_out by 4.4e-2 of the feed, which
      includes the inlet change.
      Done 2026-11-03: `liu/` runs with ROS34PW2 in place (I-BC-4); all 50 finish, and the Figure 4 curves are
      monotonic to 1e-16.
- [x] **B-EX-12** `liu/conversions.md` §8 is out of date since F-PHY-2. It describes a five-parameter GAB that gets
      its temperature by inverting `Component::scale`, with `nonIsothermal` and `HeatOfAdsorption` 1.0. The code
      takes three parameters and gets T from `setTemperature`, and `liu/generate.py` already writes three. Found
      2026-11-04.
      Fixed 2026-11-04 with V-PHY-8: §8 describes the three-parameter GAB and the `"Water"` flag.

- [x] **B-EX-13** `examples/BEA-alkanes-C7/psa` stopped at its blowdown since F-BC-2, with "pressure gradient is too large"
      (2026-11-09). The step from 1 bar to 10 mbar rescales the gas, and the desorbing C7 cannot leave through the bed at
      10 mbar (roughly 30 m/s through 0.1 m of 1 mm particles, by hand). Before, the gas stayed at 100 times its amount for
      the whole 5000 s (R T Σc/P − 1 = 99), so that the example showed no real blowdown. Not in the paper. Fixed the same
      day with an exponential `PressureRamp` of 50 s in the Blowdown phase; the feed still flows, so that it is no
      closed-end blowdown yet.

- [x] **B-BC-5** `SimulationTime` `"auto"` read the breakthrough at node N also with `FlowDirection` `Reverse`, where
      that is the feed end: a reverse BEA run (3 components, `FixedVelocity`, N = 20) stopped at 530 s instead of
      5150 s, with nothing at z = 0 yet. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: `breakthroughConverged` (`breakthrough.cpp`) takes node 1 in reverse flow; forward and reverse
    both end at 5150 s.
- [x] **B-FIT-9** `Fitting` held `parents`/`children` as references to its own two population vectors, so every copy
      pointed into the object it was copied from: `ruptura.load_simulation()` on a fitting input (returns the
      `Fitting` by value, `pybindings.cpp`) segfaulted in `run()`. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: the two vectors are `parents` and `children` themselves (`fitting.h`, `fitting.cpp`); the BEA
    nC7 fit runs through `load_simulation`.
- [x] **B-IO-13** A `MixturePrediction` input was accepted without a usable pressure grid: `NumberOfPressurePoints` 0
      segfaulted, and a missing or zero `PressureStart` (default −1) gave NaN in every output with exit code 0.
      Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: the reader requires `PressureStart` > 0, `PressureEnd` > 0 and `NumberOfPressurePoints` ≥ 1 for
    this type. `MPDInput.WarnsWhenReferenceTemperatureHasNoTargetTemperature` got a pressure range; all 246 inputs
    under `examples/` still validate.
- [x] **B-IO-14** A negative component `GasPhaseMolFraction` was normalized and run (negative inlet concentration).
      Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: rejected, like the swing-phase fractions ("must be finite and non-negative (component ...)").
- [x] **B-IO-15** `NumberOfGridPoints` 0 and `ColumnLength` ≤ 0 were reported as `ColumnDistances` errors, because
      `validateColumnDistances` ran before their own checks. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: the two checks come first, for every simulation type.
- [x] **B-IO-16** The `Breakthrough` constructor emptied `column.data` and the `component_*.data` files of the working
      directory, so that `ruptura.load_simulation()` in the folder of a finished run destroyed its results
      (`compute()` writes no files). Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: `Breakthrough::run` and `SwingAdsorption::run` start the files (`Breakthrough::appendOutput`:
    append after `ReadColumnFile` or once started, otherwise truncate and write the header). Constructing the
    object leaves the files alone; the files of a run are byte-identical (BEA, ROS34PW2).
- [x] **B-GUI-3** Ruptura Lab put the field path, which holds component names, unescaped into `id="..."` and
      `for="..."` (`app.js`, `field`, `lawField`): a component name with a quote in an imported input ran script in
      the Lab page, which can post an executable path to `/api/settings`. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: the ids go through `escapeHtml`; the probe input no longer adds a handler.
- [x] **B-GUI-4** `/api/bootstrap` returned the session key to any local client, so that another user or process on
      the machine could use the API (and `/api/settings` runs the path it is given). Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: the key is in the address the server prints and opens (`/?session=...`), the page takes it from
    there, and bootstrap needs it like every other API call (403 without). `webapp_test.py` and
    `webapp_browser_smoke.py` follow.
- [x] **B-PY-5** `Breakthrough.run()` and `SwingAdsorption.run()` from Python ignored Ctrl-C until the run was over.
      Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: `Breakthrough::afterOutput` is called after every output of both run loops, and the bindings
    set it to the signal check that `compute()` uses (KeyboardInterrupt 3.0 s into a 13 s run).
  - 2026-11-18: `Fitting.run()` too: `Fitting::afterOutput` after every generation of the genetic algorithm and every
    simplex step (KeyboardInterrupt 2.9 s into the 32 s Fe-MOF-74 fit). `advance_to` has no check.
- [x] **B-PY-6** The Python readers of `simulation.json` looked keys up in exact case (`utils.py`, `plot_fitting.py`,
      `selfdriving/psa.py`, the generated notebook): an input with `components` gave plots without curves. Found in
      the 2026-11-17 audit.
  - Fixed 2026-11-18: they pass the input through `ruptura.webapp.config.canonicalize`, which the Lab already used.
    Its temperature fallback follows the reader's default, 298 K since R-IO-2.
- [x] **B-PY-7** `BreakthroughPlotly.explorer()`, the last cell of every generated notebook, needs `anywidget` with
      Plotly 6 and later, which no environment listed. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: `anywidget` is in the `lab` extra and in `env.yml`.
- [ ] **B-PY-8** A resumed `Campaign` (`ruptura/selfdriving/campaign.py`, `run`) seeds its generator with `seed +
      len(rows)`, whereas an uninterrupted campaign draws every round from one generator seeded with `seed`: after an
      interruption the campaign continues, but not as it would have, so "the same seed reproduces it" (paper, README)
      holds only for campaigns run in one go. Found 2026-11-19: the two-bed campaign was resumed after round 0 and was
      therefore rerun in one go. Fix (ask first: it changes the setups of every campaign after round 0, so both lab
      notebooks would have to be rerun): one generator per round, `np.random.default_rng([seed, round_index])`.
- [x] **B-PY-9** `FittingPlotly.fit_figure()` and `all_fit_figure()` (and so the fitting notebook) stopped with a
      ZeroDivisionError when the input starts the fit from all-zero parameters, as `examples/BEA-alkanes-C7/fitting`
      does: the start Toth curve divides by ν = 0. Found 2026-11-20 while fixing B-PY-2.
  - Fixed 2026-11-20: `_evaluate_isotherm` gives NaN at points where the isotherm divides by zero or overflows, so the
    start curve is left out and the fit and data are drawn (test `tests/plot_fitting_test.py`).
- [x] **B-BLD-10** No CI environment had scipy, so `tests/selfdriving_test.py` (`importorskip`) never ran. Found in
      the 2026-11-17 audit.
  - Fixed 2026-11-18: scipy is in `env.yml`, `docker/*/install.sh`, the macOS jobs and `tox.ini`.
- [x] **B-BLD-11** `pip install .` built with `-march=native` (`RUPTURA_NATIVE_ARCH` ON), which gives "Illegal
      instruction" on an older CPU. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: OFF in the `cmake.define` of `pyproject.toml`; the README names the option. The presets other
    than the conda recipes still build native.
- [x] **B-EX-14** `clean_example_outputs.py` removed every `analysis.ipynb` under `examples/`, also the hand-written,
      then untracked notebook of `Self-driving-lab-H2-PSA`. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: `analysis.ipynb` is no longer in its list. Staged (not committed): `TODO_legacy.md`, the
    README, `.gitignore` and notebook of `Self-driving-lab-H2-PSA`, `PSA-H2-activated-carbon` and the two
    `Reactions` inputs (B-BLD-1). Left untracked: the figures and `campaign/experiments.csv` of the self-driving
    example.
- [x] **B-GEN-1** README: "Running" pointed to an example and a `./run` script that do not exist, and "Input" to the
      2023 article. Found in the 2026-11-17 audit.
  - Fixed 2026-11-18: a BEA example with `ruptura` and `--validate`, and `docs/GENERATIVE_AI_SIMULATION_PROMPT.md`.
- [x] **B-ITG-23** ROS34PW2: a rejected step whose stages converged but failed the error test left the IAST with the
      starting values of its last stage, from which the mixture prediction of the step start did not always converge
      again; the run ended with the IAST error. Found 2026-11-19 in round 0 of the two-bed campaign (F-PY-5): the
      blowdown column of exp. 5, helium at rest at 1.68 bar with its feed end open, into which traces of the feed
      disperse. Fixed 2026-11-19: when that recompute throws, it is repeated from the starting values of the step start
      (B-ITG-22 restores them before it only after a failed stage). A step that worked before is unchanged, so earlier
      results stay byte-identical. Test: `ROS34PW2.RecomputesTheStartOfARejectedStepFromItsOwnMixturePrediction`.
- [x] **B-ITG-22** ROS34PW2: after a stage whose mixture prediction failed (B-ITG-18), the rejected step recomputed
      the fields of the step start with the IAST starting values that the failed stage had left, from which the
      inversion of the step start did not converge again in 50 iterations: the run ended with the IAST error, outside
      the guard of the stages. Found 2026-11-18 with `examples/PSA-H2-activated-carbon-two-beds` on 12 cells (second
      cycle, pressurization after an equalization). Fixed 2026-11-18: the integrator keeps `cachedPressure` and
      `cachedGrandPotential` of the step start and restores them after a failed stage; an ordinary rejected step is
      unchanged. Test: `Connections.TheTwoBedExampleKeepsEveryMoleBetweenItsColumnsItsTankAndItsProduct` runs that
      case. CVODE has the same warm starts and was not looked at.
- [x] **B-ITG-21 `MomentumBalance` (F-BC-3) is efficient with ROS34PW2 only.** Found 2026-11-18 on the T15 breakthrough
      (30 cells, isothermal, 1500 s; ROS34PW2 0.39 s, 3253 steps):
  - CVODE Dense: 4766 steps, 1541 Newton failures, 21 s. Its difference-quotient Jacobian perturbs a concentration of
    zero by about 1e-16, which is lost in the pressure sum of 1e6 Pa, so that the stiff pressure coupling of every
    absent component (ahead of its front; the carrier behind it) is missing.
  - CVODE SPGMR with the `SorptionJacobian` preconditioner: 832 742 steps (45 s); the inexact linear solve leaves a
    pressure error of the size of the pressure differences that drive the gas (31 Pa over the column at 1e6 Pa).
  - CVODE KLU (analytical Jacobian, no energy balance): 20 621 steps, 833 Newton failures, 0.86 s, the same as
    `InletVelocityOutletPressure` with KLU (25 065 and 833, B-ITG-11); usable.
  - RungeKutta3: the pressure relaxation needs steps of about 1e-6 s.
  - With CVODE the stream amounts through a pressure port are less accurate than the column amounts: 1.4e-4 and
    4e-3 of the amount in the inert fill and empty of `MomentumBalance.APressurePortFills...` (KLU, rtol 1e-5), against
    1e-9 with ROS34PW2. The flux follows a pressure difference of some 10 Pa, and CVODE's tolerance is relative to
    the amount that has passed. `cycles.data` shows it as the mass balance.
  - Fixed 2026-11-19 (decision of the user): the reader rejects `MomentumBalance` with RungeKutta3 and with CVODE Dense
    or SPGMR, where it printed a warning; README updated. The stream amounts of CVODE KLU stay as above (README).
    Test: `CycleInput.ReadsPortsStopCriteriaAndSteadyStateInAnyCaseAndRejectsWhatDoesNotFit`. The SI (`main.tex:1176`) gives the Dense and
    SPGMR timings but does not say that the reader rejects them.
- [x] **B-ITG-17** ROS34PW2 with WENO3 does not conserve the adsorbate exactly on the trace case of B-ITG-16 with
      `FixedVelocity`: after loading the column holds 2.9e-4 more than was fed, from the first output on and in
      proportion to the amount fed; upwinding gives 1e-7, and the right-hand side conserves it to 1e-16 at every
      output. Same with the binary before flow reversal.
  - 2026-11-19: not integration error and not the Jacobian: the excess is independent of the tolerance (1e-4 to 1e-8)
    and of `ROS34PW2MaximumTimeStep`. WENO3 undershoots ahead of the steep front (down to −1.6 % of c_F in c and
    −2.3e-4 mol/kg in q), and ROS34PW2 clamped the accepted state to max(0, c) and started the next step from it, so
    that every step added the undershoot. Fixed (`src/integrators/ros34pw2.cpp`): the column still holds the clamped
    state for the outputs and the mixture prediction, but while it holds exactly that clamp the next step starts from
    the unclamped state, as CVODE does; a state changed from outside (stop bisection, phase switch, coupling) is taken
    as it is. 3.4.4 trace case after loading: 2.97e-4 → 2.8e-6 (WENO3), 2.5e-7 → 1.7e-8 (upwind), the rest being the
    clamped undershoot of the moment, which does not add up; outlet over the purge within 6.4e-4 (co-current) and
    2.9e-4 (counter-current) of the stored runs, integrated. RK3 clamps after every step too: on the monolith of
    `ROS34PW2Column.IntegratesChemisorptionLikeRungeKutta3` it holds 9.8e-7 mol/kg more loading at the foot of the
    front than ROS34PW2 now, so that test allows 1.5e-4 instead of 2e-5. Test:
    `ROS34PW2Column.ConservesATraceAdsorbateBehindAWeno3Front` (3.0e-4 before, 9e-6 after). Paper: I-ITG-4.
- [x] **B-ITG-20** `RungeKutta3` does not check its step. A `TimeStep` above the stable one gives wrong curves with
      exit code 0 and no message: the BEA example with `TimeStep` 2 s (dz / v = 0.053 s) writes c/c0 up to 100. Found
      in the 2026-11-18 audit of `src/`.
  - Closed 2026-11-18 without a check (user): documented in the README, in section 9 of
    `docs/GENERATIVE_AI_SIMULATION_PROMPT.md`, and at `RungeKutta3` (`rk3.h`): below dz / v, 2 / k, and dz² / (2 D).
- [x] **B-PHY-14** Fast IAST ended at its limit of 50 Newton steps and used the unfinished answer. With BET above
      p = 1 / c or Anti-Langmuir above p = 1 / b the pure loadings are negative and the mixture loadings jump
      (4.7 to 22.6 mol/kg), with exit code 0. Found in the 2026-11-18 audit of `src/`.
  - Fixed 2026-11-18: `computeFastSIAST` throws when the relative spread of the reduced grand potentials or
    |sum x - 1| is not below 1e-6 after the loop (also NaN). CVODE and ROS34PW2 retry the step (B-ITG-18); RK3 and a
    mixture prediction stop with the message. Test `MixturePredictionConvergence.FastIASTThrowsBeyondThePoleOfAnIsotherm`.
- [x] **B-IO-17** Output files that could not be written went unnoticed: a component named `n-C7/iso` gave exit code
      0 without its file, and a full disk would cut the files off. Found in the 2026-11-18 audit of `src/`.
  - Fixed 2026-11-18: `Breakthrough::writeOutput` (R-IO-3) flushes `column.data` and the component files after every
    snapshot and throws when one failed; `Column::writeJSON`, `phases.data`, and the component files of a mixture
    prediction are checked too. Test `Output.AFileThatCannotBeWrittenStopsTheRun`.
- [x] **B-IO-18** A `ReadColumnFile` run appended to the old output files with the clock back at t = 0 (was B-IO-8,
      closed as intended on 2026-11-09): times 0 ... 40, 0 ... 30 in one file, the first new block the old end state,
      and no header in a new directory.
  - Fixed 2026-11-18 (user: continue, with a switch for the old start at t = 0): `column.json` holds `time`
    (`Column::writeJSON(fileName, time)`, `readJSON` returns it). A `ReadColumnFile` run goes on from that time for
    `SimulationTime` (or the phase durations), with its outputs on the same time grid, appends, leaves out the
    snapshot it starts from, and writes the header when the file is new; CVODE starts at that time
    (`CVODE::initialize(column, startTime)`), and `phases.data` is continued. `"ResetClock": true`
    (`InputReader::resetClock`) starts a new experiment from the state: t = 0 and new output files. A file without
    `time` starts at 0. Test `Restart.ContinuesTheClockAndTheOutputFilesUnlessResetClock`; a continued PSA-H2 cycle
    (CVODE) has the phase and output times of the two-cycle run. Not done: a resume in the middle of a swing phase
    (the run starts at the first phase), and a field for `ResetClock` in the Lab (the key passes through).
- [x] **B-MBD-3** Two places took the first bed for the whole column: the gas-wall area per fluid volume, 4 / (eps
      d_in) (`computeEnergyDerivatives`, the velocity closure in `computeDerivatives`, `assembleThermalJacobian`), so
      that with beds of eps 0.4 and 0.6 the gas of the second bed exchanged heat with the wall 50% too fast; and the
      capacity that limits a Physisorbed or Chemisorbed reaction (`computeReactionDerivatives`). Found in the
      2026-11-18 audit of `src/` (by reading).
  - Fixed 2026-11-18: `Column::fluidWallContactArea` per node (`computeAdsorbentBlends`, weighted by the adsorbent
    fractions like `fluidSolidContactArea`), and the capacities of the adsorbents at the node, weighted in the same
    way. A single bed is unchanged. Test `Multibed.GasWallAreaFollowsTheVoidFractionOfEachBed`; no run compared.
- [x] **B-IO-19** `--validate` accepted inputs that cannot run: `Temperature` 0 (CVODE convergence failure),
      `NumberOfGridPoints` 1 with WENO3 (fails at the first right-hand side), a negative `MassTransferCoefficient`
      or `AxialDispersionCoefficient`, and a zero fluid or solid heat capacity with the energy balance. Found in the
      2026-11-18 audit of `src/`.
  - Fixed 2026-11-18: the reader rejects them (also a swing-phase `Temperature` of 0). Test
    `InputValidation.RejectsInputsThatCannotRun`.
- [x] **B-IO-20** Component names as keys were matched with their case in `GasPhaseMolFractions`,
      `ProductEndGasPhaseMolFractions`, `ComponentParameters`, the adsorbent `Components`, and the MPD
      `ComponentBounds`, and without it in reactions and `AdsorptionConstants`; `"CO2"` next to `"co2"` was not
      rejected. Found in the 2026-11-18 audit of `src/`.
  - Fixed 2026-11-18: one lookup (`findComponent` in `inputreader.cpp`: the exact name, else the only match without
    case); a component named twice in a fractions object or in `ComponentParameters` is an error. Adsorbent names
    are still matched exactly. Test `InputKeys.ComponentNamesAsKeysAreCaseInsensitive`.
- [x] **B-IO-21** A missing `MolecularWeight` became 1 kg/mol without a warning (the gas density, in the inertial
      pressure drop and the heat capacity of the gas). Found in the 2026-11-18 audit of `src/`.
  - Fixed 2026-11-18: one warning that names the components of a gas column run without it. The default stays.
- [x] **B-FIT-10** The fitting data reader took broken numbers as numbers (`1,5 0,3` as the point (1, 0), a header
      line without `#` as (0, 0)), and `trim` stepped before the start of an empty line. Found in the 2026-11-18
      audit of `src/`.
  - Fixed 2026-11-18: `Fitting::readData` throws with the line number when the pressure or loading is not a whole
    number token; comments are found with `find_first_not_of` (`trim` and `startsWith` are gone). The 22 data files
    of the fitting examples load as before. Test `ExampleFitting.RejectsDataThatIsNotANumber`.
- [x] **B-PHY-15** Every isotherm took the gas temperature T_f (decisions of I-PHY-12 and I-PHY-13, 2026-11-04), although
      the van 't Hoff factor of the affinity comes from desorption at the surface and the heat of adsorption goes into
      the adsorbent. Reviewer-2 audit of 2026-11-19: in the 30-cell T15 run (SI free-pressure check) T_s peaks 11.8 K
      above T_f, so the CO₂ affinity at the front was about 39 % too high.
  - Fixed 2026-11-19 (user's request, reverses the 2026-11-04 decisions): the physisorption and chemisorption
    equilibria, their sensitivities, and the RH inside the humidity isotherms take T_s: `rk3.cpp` (precompute),
    `column.cpp` (initialization), `jacobian.cpp` (new `solidTemperature` argument of `assembleSorptionJacobian`; the
    gas pressure R T_f C of MomentumBalance keeps T_f), its callers in `cvode.cpp` and `ros34pw2.cpp`, and
    `compute_equilibrium_loadings` in `pybindings.cpp`. Isothermal runs are unchanged (all temperatures equal). Test
    `SorptionEquilibrium.FollowsTheAdsorbentTemperature`; suite 252/253 before and after (the known
    `LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix`). Reruns: V-PHY-9; paper: I-PHY-16.
- [x] **B-ITG-24 Explicit RK3 goes unstable in the energy balance at a step within its documented limits.** The
      limits written down under B-ITG-20 (dz / v, 2 / k, dz² / (2 D)) leave out the heat conduction of the gas. Found
      2026-11-20 in the full rerun on the cluster node: `elfving/seven_bounded_35c_2pct_500mlmin` (λ = 3.12 W/(m K),
      dz = 0.58 mm, `TimeStep` 1e-4 s) reaches 1e29 K within 10 ms, overflows, and stops with "GAB needs a positive
      temperature". The stored run of 2026-10-30 had stopped at its first output too ("GAB isotherm is undefined for
      K p/p_sat >= 1"). The limit is about 0.6 dz² ρ c_p / λ: 7.9e-5 s for this input, and above 1.1e-4 s for the
      other 12 (λ ≤ 2.12 W/(m K)).
  - Fixed 2026-11-20: the limit stands where the others are documented (README, section 10 of
    `docs/GENERATIVE_AI_SIMULATION_PROMPT.md`, `rk3.h`). The input has `TimeStep` 5e-5 s, and `elfving/generate.py`
    writes that step for a conductivity above 2.4 W/(m K). With the binary of this fix the run stops within 0.2 s
    of wall time at 1e-4 s, and reaches t = 80 s in 90 s of wall time at 5e-5 s (stopped there; the whole run takes
    about 4 h and was completed on the cluster node only). RK3 still does not check its step (B-ITG-13, B-ITG-20),
    so a step above the limit still ends with the message of the isotherm.
  - Fixed with it: GAB, the GAB-modified Toth, and the Wagner saturation pressure threw `std::domain_error`, which
    CVODE takes as fatal and ROS34PW2 does not catch, so a trial state outside their domain would end the run
    instead of shortening the step (B-ITG-18). They throw `std::runtime_error` now, like every other failure of a
    state (`isotherm.h`, `water_saturation_pressure.h`; three tests expect the new type). Known from reading the
    code; no run showed it.
- [x] **B-ITG-25** ROS34PW2 has no check for crawling steps like CVODE's (B-ITG-19): a run whose accepted steps
      collapse goes on without a message until it is killed. Found 2026-11-20 in the full rerun on the cluster node,
      where the 13 Elfving copies of the ROS34PW2 cross-check (B-EX-9, B-EX-16) each ran until the 16 h timeout.
  - Fixed 2026-11-20: `ROS34PW2::advance` counts the accepted steps since one reached its endTime
    (`stepsSinceEndTime`, `timeSinceEndTime`) and stops when 100 000 of them advanced the time by less than 1e-3 of
    the time in which the fluid crosses the shortest cell, or of the step cap if that is smaller; without flow it
    does not stop. This is CVODE's rule, whose count is that of one `CVode` call; here it also holds for the single
    steps of a phase with a stop criterion. Steps that an output time or a phase end shortens start the count
    again. Test `ROS34PW2Column.StopsWhenItsStepsCrawl`.
    - The Elfving copy `five_25c_2pct_500mlmin` stops after 11 s of wall time at t = 6.3e-3 s: mean step 6.3e-8 s,
      3.8e-5 of the 1.7 ms in which the gas crosses a cell.
    - Of the 900 stored ROS34PW2 logs, only two-bed runs of the self-driving lab reach 100 000 steps (35 in
      `campaign_two_beds_reference_settings`). The three with the most (`exp_124`, `exp_084`, `exp_102`; 183 000 to
      296 000 steps per column) finish with the fix; the other 32 were not rerun. For `exp_102` the binaries from
      before and after the fix give the same step counts and the same `cycles.data`, `connections.data`,
      `buffers.data`, and `phases.data`. The step counts differ from the stored logs of 2026-11-19 13:40 with both
      binaries, so from the source changes since then.
- [x] **B-EX-16** The ROS34PW2 cross-check copied the 13 Elfving inputs, which ROS34PW2 cannot run (B-EX-9), and
      `cases.csv` gave the case no timeout, so `run.py --paper` ran each copy until the 16 h default. The runs of
      2026-11-03 had left them out by naming the other studies on the command line. Found 2026-11-20 in the full
      rerun on the cluster node.
  - Fixed 2026-11-20: `validation/ros34pw2_danckwerts_rerun/generate.py` leaves Elfving out until B-EX-9 is decided.
    The 13 unrun copies and their entries in `simulations/index.json` are removed, which leaves the 223 cases of
    2026-11-03. The case has a timeout of 3600 s in `cases.csv`, 18 times its longest run of 2026-11-03, and a
    crawl now stops by itself (B-ITG-25). `jobs.csv` is rewritten (`run.py --paper --list`): 1132 jobs, which also
    lists 67 inputs written to `3.5_numerical_performance/breakthrough/simulations` on 2026-11-20 by other work.
- [x] **B-BLD-12** `GasBreakthroughRegression.ExistingColumnCvodePathRetainsItsOutcome` fails again, on the cluster
      node (it was fixed under B-BLD-5). It compared one CVODE interval at rtol 1e-3 with stored values to 1e-9 and
      1e-12, which holds on the machine that recorded them only: on the Mac the values have moved by 7e-11 and
      1.4e-13 since. Found 2026-11-20 in the full rerun.
  - Fixed 2026-11-20: the bounds are CVODE's relative tolerance, 1e-3 of each value. The cell starts empty, and the
    solution converged in the tolerances (rtol 1e-9, atol 1e-12) lies 5.2e-4 of the stored concentration below it,
    so the test still pins what these settings give. Passes on the Mac, also built without `-march=native`; not
    rerun on the cluster node.
- [x] **B-BLD-13** The three `Connections` tests and `ThermalControl.AJacketHeatsAClosedColumnWhichGivesTheSameHeatBack`
      read back `phases.data`, `connections.data`, `buffers.data`, and `cycles.data` from the working directory that
      all tests share, without the `WorkingDirectoryGuard` of the other tests that read their output. Under
      `ctest -j` another swing test overwrites these files. Found 2026-11-20 in the full rerun on the cluster node
      (3 failures).
  - Fixed 2026-11-20: each of the four runs in its own directory. 25 parallel repeats of the 23 swing tests
    (`ctest -j8 --repeat until-fail:25`) failed on three of the four before
    (`Connections.ABufferKeeps...`, `Connections.TheTwoBedExample...`, `ThermalControl`) and pass now. Suite with
    `ctest -j8`: 256 of 257, the known `LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix` (B-EX-3).

## Features (F)

F-BC-1 and F-CYC-1 were described in the paper before the code had them; both are done now (F-BC-1 2026-11-04,
F-CYC-1 2026-11-18, with F-PHY-2 and F-PHY-1). The paper text on cycles still has to follow the code (I-CYC-1, I-CYC-2).

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
  - Decided 2026-11-03: one surface and one pore reservoir per component and node (not per chemisorption site),
    shared by its chemisorption sites and pore reactions; the pore concentration is per pore volume, with a new
    per-adsorbent `ParticlePorosity` ε_p (default 1, today's solid-volume reservoir; fixes B-PHY-11 on the way);
    `FilmMassTransferCoefficient` and `PoreDiffusivity` move from the chemisorption site to the component (per
    adsorbent in a multibed), `convert.py` migrates the inputs (16 `benchmark/figures`, 80 `benchmark/multibed`,
    `examples/Chemisorption-CO2`); physisorption stays driven by the bulk gas; the `Site` index of
    `PoreConcentration` reactions is dropped.
- [ ] **F-ITG-2** Second-order total molar flux and gas-temperature advection. With WENO3 (F-ITG-1) the hydrogen
      t50 of Marx T15 still converges at first order (e50 1.2e-2 at N_z = 100 against a 400-point reference), because
      v c_tot and the gas-temperature advection stay upwind. A WENO3 temperature flux alone changed the T15 t50 by
      < 0.1 s (2026-10-26 test build); a higher-order total flux needs a velocity closure that conserves it.
- [ ] **F-ITG-4** Velocity-closure derivatives in the sorption Jacobian, for the CVODE preconditioner and KLU
      (F-ITG-3 (4), (4b)): with `FixedPressureInletVelocity` the composition dependence of the Ergun velocity (local),
      with the continuity settings the dependence of the velocity on all upstream uptake (total molar flux per node as
      an unknown). Without them, CVODE accepts Newton iterations whose error has the same sign every step, and the
      preconditioned and KLU runs drift by up to 1e-4 at rtol 1e-5. ROS34PW2 does not need them (W-method).
      2026-11-03 (after C-BC-1): with the closure term, Σc is an exact invariant of the right-hand side, so stage
      errors are no longer damped by dispersion. ROS34PW2, whose Jacobian lacks the closure, leaves a frozen
      max |R T Σc/P − 1| of 8.2e-3 at its default tolerances on the Stampi-Bombelli packed bed (IVOP, isothermal;
      formed at the inlet nodes during the first minutes), 1.3e-3 at rtol 1e-7; CVODE 9.4e-5. The closure
      derivatives would restore the invariant in the stages.
      2026-11-06: B-ITG-16 is the strong case: a stiff trace adsorbate lets Σc drift by 17 % with ROS34PW2 at the
      default tolerances (and more at rtol 1e-8), so the continuity settings cannot be used for such feeds until the
      Jacobian carries the closure.
      2026-11-11: ROS34PW2 no longer needs them for Σc: its stage solve is projected on the closure (B-ITG-16). The
      derivatives remain open for the CVODE preconditioner and KLU.
- [ ] **F-CYC-5 Limits of the coupled columns of F-CYC-4** (2026-11-18; each is a `shortcut:` comment or a message in
      `SwingAdsorption::run`). Take them up when a cycle needs them:
  - One end cannot take gas from one connection and give it to another in the same step (a tee), and gas cannot go
    around a ring of columns: both stop the run with a message. The first needs a junction, the second the steps of
    the columns solved together.
  - The coupling is explicit: a valve that evens out its pressures faster than 1e-4 of the phase is cut to that
    rate, and a wide valve costs many coupling steps. An implicit coupling (the valve in the matrix of ROS34PW2, or
    all columns in one state) removes both.
  - All columns start alike (at rest at `OutletPressure`), and `ReadColumnFile` is refused with `Columns`: no
    per-column initial state and no restart of a run with several columns or buffers.
  - `ThermalControl` sets a jacket temperature and coefficient; a set heat flux (the branch has one) is not there.
    The gas of a compressor enters at the temperature of its donor (no heat of compression in the column).
  - The duty integrals (wall heat, exhaust work) are integrated with the column but do not follow the change that W
    makes to the face velocities, as the stream amounts do; the same holds for the stream amount through a flux port,
    which leaves the amount of a column 1e-6 off its streams (the transfer between columns itself is exact).
- [ ] **F-CYC-6 The rest of the plan of 2026-11-18 for the mechanics of the branch `ruptura-shrinjay`**: stage 6.
      Stage 5, the optional physics, is done (2026-11-19: F-ITG-5, F-PHY-12, F-PHY-13, F-PHY-14; the user asked for
      all of them as options, and for WENO5). Open: a `convert.py` mapping for the inputs of the branch
      (`LeftPort`/`RightPort` to `FeedEnd`/`ProductEnd`, its `ParticleTransport` and `Coolant` keys), the builder of
      F-GUI-2, and the paper text (I-CYC-2, `ruptura-2.0-paper/CYCLES_NOTES.md`). Not taken, with the user's consent
      (2026-11-19): the second engine, YAML, NPZ, the cost module, FDM/FVM, native BDF, liquid networks. Not asked for
      and not done: gas-solid thermal equilibrium as an option, a film resistance and cylinders in HSDM, a quasi-steady
      coolant, the inlet-density basis, initial states per column.
- [ ] **F-GUI-2 A builder for multicolumn cycles in Ruptura Lab** (asked for on 2026-11-19, for later): place columns,
      buffers, valves, and compressors, and set per phase which are open and what the ports are; it writes the
      structured input that the reader takes (`Columns`, `Buffers`, per-phase `Columns` and `Connections`) and shows
      the flowsheet. There already: the JSON input (F-CYC-4), the key list and field help of Ruptura Lab, and
      `ruptura.flowsheet.flowsheet_tikz` (F-PY-4), which the builder can call for the picture. Open design points: a
      schedule view (phases as columns of a table, one row per column and connection); per-column output in the run
      view (F-GUI-1); whether a connection gets a name once and phases refer to it (now each phase repeats its
      connections, and the same ends and name make one valve in the flowsheet).
- [x] **F-GUI-3** Ruptura Lab's law editor (`ruptura/webapp/config.py`, `TEMPERATURE_LAWS`; `static/app.js`) lags the
      code since 2026-11-19: q_sat accepts VantHoff (with its own Energy), exponents accept `Arctan`, and Arctan has a
      second coefficient (`Beta`), but the Lab offers Exponential/Linear and Linear/Inverse with one coefficient. The
      Marx T1–T24 inputs cannot be edited there. The preset `presets/paper-t15.json` is older still (ρ_s = 502,
      IAST): decide whether it follows the Marx inputs (needs this item) or the reference settings (B-EX-15).
  - 2026-11-20 (user's request): `config.py` offers VantHoff on capacities and Arctan on exponents, with a second
    coefficient per law (`TEMPERATURE_LAWS["secondCoefficients"]`, Arctan: `Beta`); `app.js` draws it, `field_help.py`
    has `ArctanAlpha`, `ArctanBeta`, and the Law list. `tests/webapp_test.py` validates every offered law against the
    C++ reader (77 passed with the frozen binary); checked in Chromium on the T15 preset (Arctan, α, β, T_ref shown).
- [ ] **F-GUI-1** `ClosedEnd` and `PressureRamp` have no controls in the phase forms of the Qt app and Ruptura Lab
      (2026-11-09, F-BC-2). Both apps keep them in a loaded input, and the field help of Ruptura Lab has them. In the Qt
      app `SwingPhaseConfig::otherKeys` holds every phase key without a control, also `ColumnEntranceVelocity` and
      `GasPhaseMolFractions`, which it used to drop on import (2026-11-09).
      2026-11-18: the Qt app is removed (R-GUI-2); the Lab part of this item is left.
      2026-11-18: the same holds for `StopCriterion`, `MinimumDuration`, `TimeoutPolicy`, `CyclicSteadyState`, the ports
      `FeedEnd` and `ProductEnd`, and `GasHeatCapacity` (F-CYC-1, F-BC-3): Ruptura Lab keeps them in a loaded input and
      offers `MomentumBalance` as a boundary condition, without controls for its ports.
      2026-11-18: also for `Columns`, `Buffers`, `Connections`, `ThermalControl`, `ExhaustPressure`, and
      `CompressionEfficiency` (F-CYC-4); Ruptura Lab shows the output of one column, not the per-column directories.
- [ ] **F-PY-1** The breakthrough plots (`ruptura/plot_breakthrough.py`, `grid_index=-1`) and the webapp plots take
      the last node as the outlet. In a reverse step the gas leaves at z = 0, i.e. from the first cell (node 1);
      node 0 holds that cell's state at the start of the step. Plot the outlet of each swing phase at the end it
      leaves from (2026-11-04, F-BC-1).
- [x] **F-PY-3** `ruptura.selfdriving` ends its cycles in Python: blocks of 10 cycles restarted from `column.json` until
      the amounts and the KPIs change by less than 2e-3. `CyclicSteadyState` (F-CYC-1) does this in one run, and
      `phase_records[i].stream_amounts` gives the product without the known-feed arithmetic of `key_performance`.
      Not changed on 2026-11-18: the campaigns of the paper example were run with the block scheme.
  - Done 2026-11-19 (user request): the input of an experiment holds `CyclicSteadyState` (tolerance 2e-3) and
    `NumberOfCycles` = `maximum_cycles` (100; `cycles_per_block` and `maximum_blocks` are gone), `run_experiment` is one
    run, the KPIs come from the stream amounts (product: what left at z = L in the adsorption phase; used: what
    entered at z = L in purge and backfill; fed: what entered at z = 0 in adsorption and pressurization), and the
    status and the two change columns of the lab notebook from the last line of `cycles.data`. An experiment
    directory written before needs `Campaign.write_inputs()` before a rerun (its `robot.json` holds the old fields).
    Tests: `tests/selfdriving_test.py` (9, with the module built from this tree). Open: V-EX-1.
- [x] **F-PY-5 Self-driving lab on two beds** (2026-11-19, user request: "rerun the self-driving campaign together with
      a multicolumn setup"): `ruptura.selfdriving.TwoBedPressureSwingCycle`, a subclass of `PressureSwingCycle` on a
      `MomentumBalance` base input (`examples/PSA-H2-activated-carbon-two-beds`). Six phases per cycle (pressurization
      and blowdown; adsorption with a product valve to a tank held at P_H and a purge valve to the other bed, both check
      valves; equalization until the receiving bed is at P_L + e (P_H - P_L), at most t_r); the variables of the
      one-column robot with the equalization share e (0 to 0.5) in place of the backfill. Valves sized from the molar
      feed flow (product: at 2 % over P_H; purge: G times it between P_H and P_L; equalization: the gas of a column
      within a tenth of t_r). `key_performance` takes the product of two beds as the net outflow of their product ends
      (what reaches the tank), the productivity per kg of the sorbent of both. Tests in `tests/selfdriving_test.py`
      (input, KPIs, a small campaign). Baselines at 80 cells: one column 0.584 mmol/(kg s), recovery 0.150, purity
      0.99968 in 6 s; two beds without equalization 0.562, 0.144, 0.99956 in 67 s; with e = 0.3 0.575, 0.202,
      0.99993 in 48 s.
- [x] **F-PY-4 TikZ flowsheet of a multicolumn input** (2026-11-19, user request): `ruptura/flowsheet.py`,
      `flowsheet_tikz(data, phase=None)` and `python -m ruptura.flowsheet simulation.json [phase]`. Columns upright
      with the feed end at the bottom, buffers above, each connection once by its ends and name on a pipe run with a
      valve, check-valve, or compressor symbol; with a phase its connections are black with their direction, the
      others gray, and its ports are at the column ends. Plain TikZ. Checked by compiling the two-bed example and a
      case with three columns, two buffers, and a compressor from a feed end; tests `tests/flowsheet_test.py`. Lines
      cross in involved flowsheets (no routing beyond one run per connection).
- [ ] **F-EX-1** `ruptura.selfdriving` (`examples/Self-driving-lab-H2-PSA`): the brain explores the bounds of the 8-D box. In the 99.9 %
      campaign (2026-11-11) the hypervolume stalled from round 2 to round 6, and 30 % of the proposals were cycles whose
      purge and backfill take more H2 than they produce (pure product, so the probability of feasibility does not
      deter them). Keeping their negative objectives (floor at minus the best) made it worse: 4.3 instead of 2.1 of 8
      variables at a bound, hypervolume 0.33 instead of 0.44. Next: a trust region around the incumbents (TuRBO), or
      EHVI in place of ParEGO's random weights (the 97 % campaign of 2026-11-09 missed the slow high-recovery branch).
- [x] **F-PHY-1 Humidity-dependent CO₂ isotherms and relative humidity.**
  - Paper (SI table "Condition-dependent isotherms"): humidity-enhanced Langmuir, moisture-swing
    (Wang–Lackner) Langmuir, GAB-coadsorption Langmuir/Sips/Toth, GAB-modified Toth–Liu, humidity
    Toth–Climeworks, and RH evaluated from IAPWS–Wagner, Antoine, or a constant saturation pressure.
  - Code: only `GAB` (with a fixed Antoine equation for water, F-PHY-3) and `Langmuir_pH` exist
    (`src/isotherm.h:54`).
  - Needed for 3.3.4: Liu Fig. 4 now freezes the humid Toth law at the feed humidity.
  - 2026-11-04 notes, open questions:
    - The SI takes q_H₂O as the GAB loading at the local RH; Liu writes "adsorbed H₂O quantity". For Liu Fig. 4 both
      should give nearly the same result, because the water front runs well ahead of the CO₂ front.
    - The SI humidity section cites no source for any of its forms. Humidity-enhanced Langmuir, moisture-swing, and
      Climeworks (about 15 parameters) have no validation case in the repo. Moisture-swing k RH^(−n_w) diverges in
      a dry bed.
    - The e^(1.11 RH) capacity factor of Liu Figs. 5b and 7 is not one of the SI forms.
    - Isotherm-level checks need B-IO-7. Physisorption temperature: I-PHY-12.
  - Done 2026-11-04, with the decisions of that day: the water loading is the equilibrium loading at the local RH;
    all SI forms plus the exponential capacity factor; RH from the gas temperature; the water component is marked
    by a boolean `"Water": true`.
    - Types: `Humidity-enhanced-Langmuir`, `Moisture-swing-Langmuir`, `GAB-coadsorption-Langmuir`/`-Sips`/`-Toth`
      (blocking now `max(0, 1 − ω q_w)`, ω = 1/q_block, instead of the SI's δ_B switch and min(), and α without the
      δ_E switch), `GAB-modified-Toth`, `Exponential-humidity-Toth` (Liu eq. 8). The `Humidity-Toth`
      (Climeworks) type was added and then removed on the same day (W-PHY-1). Gas, SPI, and physisorption sites only, not on the Water component.
    - RH = p_H₂O / p_sat(T_f) from the component with `"Water": true`, or a fixed `RelativeHumidity`. q_w is the
      Water component's isotherm on the same adsorbent at p_H₂O. Both are formed in `MixturePrediction::humidityState`.
    - `WaterSaturationPressureModel`: `Antoine` (default, unchanged results), `Wagner` (IAPWS), `Constant` with
      `WaterSaturationPressure` (F-PHY-3).
    - Jacobian: dq*/dp_H₂O is a central difference in `predictMixtureSensitivities`; a test checks it against the RHS,
      and fails (7% error) with the term removed.
    - App, web app, bindings, `docs/GENERATIVE_AI_SIMULATION_PROMPT.md`, and the SI (pH and humidity are now their
      own subsections, with references) follow.
    - Checks: `tests/humidity_isotherm_tests.cpp` (closed forms, slopes, IAPWS values, input errors), a lab round
      trip, and Liu runs with the new binary: the stored Fig. 4, 5a, 7 outputs reproduce exactly; a fixed RH at the
      feed equals the frozen Fig. 4 runs to 1e-14; coupled Fig. 4 has the frozen t50 and differs by ≤ 1.2e-4 c₀;
      coupled Fig. 7 moves t50 by −1.5% only (humid/dry ratio 2.64 → 2.60, digitized 1.84). Follow-up: V-PHY-8.
    - The Climeworks correlation had no source and was removed (W-PHY-1).
- [x] **F-PHY-3** GAB saturation pressure is a fixed Antoine equation for water, valid from 274 to 373 K
      (`src/isotherm.h:222`).
      Done 2026-11-04 with F-PHY-1: `WaterSaturationPressureModel` Antoine (default), Wagner (IAPWS), or Constant.
- [x] **F-PHY-5** `BulkConcentration` reactions supply no reaction heat; the input reader rejects them together with
      `energyBalance`. To fix: convert the heat per fluid volume to the solid basis with the local solid mass
      per fluid volume (`compute.cpp`, `computeReactionDerivatives`).
      2026-11-03: done. The heat Q r per fluid volume is divided by `loadingPrefactor` (solid mass per fluid volume,
      the factor of the bulk species sink) and enters the solid balance, i.e. φ_f Σ χ_r Q_r r_r per bed volume (SI
      Eq. `bulk_reaction_heat_source`). The input reader now only requires ParticleDensity > 0 with the energy
      balance. Test `BulkReaction.ReleasesItsHeatIntoTheSolidPerSolidMass`.
- [x] **F-PHY-6 Variable axial dispersion**: D_ax per component as a function of velocity, e.g.
      D_ax = D_m √ε_b + u r_p (Stadler et al. eq. 12) or D_L = p₁u + p₂ (Stampi-Bombelli), so it follows each
      swing phase.
  - Decided 2026-11-03: only D_0 + λ|v|, with a temperature scaling; no free-form laws.
  - Agreed 2026-11-03 (with the pressure scaling and the Taylor–Aris term): D_ax,j = D_0,j θ + λ |v| [+ κ v² / (D_0,j θ)], with
    the local interstitial velocity, θ = (T/T_ref)^1.75 (P_ref/P) (Fuller) at the local gas temperature and
    pressure, θ = 1 without a reference state and in the liquid phase. D_0,j per component (molecular part,
    species-specific; `AxialDispersionCoefficient`, unchanged meaning without a reference state); λ per adsorbent
    (mechanical dispersion, species-independent, default 0); T_ref, P_ref once per column; κ = d_h²/192 from the
    monolith geometry (Taylor–Aris), opt-in, needs D_0 > 0. Exponents fixed.
  - Coverage: Stadler eq. 12 (u = ε_b v): D_0 = D_mix √ε_b, λ = ε_b d_p/2, and the 8 → 1 bar swing scales the
    molecular part 8× (the per-phase dispersion of F-CYC-2). Stampi-Bombelli literature: packed bed D_0 = 3e-5
    m²/s, λ = 1.5e-3 m (= d_p/2); monolith D_0 = 1.62e-5 m²/s and Taylor–Aris with d_h = 2 mm gives their
    0.0013 v². Wakao–Funazkri: D_0 = 20 D_m/ε_b, λ = d_p/2; Edwards–Richardson (high Pe): D_0 = 0.73 D_m,
    λ = d_p/2. Not covered: the fitted Stampi-Bombelli monolith law 1.22 u² + 4e-4 (stays a constant per run).
  - D is recomputed per node in every right-hand side after the closure; with the restored closure term of C-BC-1
    it is evaluated at v_{i−1} in the march, so the right-hand side stays a function of the state. Do C-BC-1 first.
  - 2026-11-03: done. `DispersionLaw` (dispersivity λ, Taylor–Aris coefficient κ) on `Geometry`, blended per node
    like the pressure-drop law (`Column::dispersionLaw`); the component input is the molecular part
    (`Column::molecularDispersionCoefficient`); `updateVelocityAndPressure` sets `axialDispersionCoefficient` in every
    right-hand side, with θ at the gas temperature and the gas pressure R T Σc of the node. Input: Geometry
    `AxialDispersivity` [m] (also per `Adsorbents[]` entry) and `TaylorArisDispersion` (Monolith only; needs
    D_0 > 0 for every component), root `AxialDispersionReferenceTemperature` / `AxialDispersionReferencePressure`
    (together, positive, gas only). Without them D = D_0 bit for bit. Tests `AxialDispersionLaw.*`,
    `AxialDispersionLawInput.*`. Wired into the Qt app, Ruptura Lab, field help, the AI prompt doc, and the run
    printout. SI: Section `apx:axial_dispersion_law`, Eq. `axial_dispersion_law`, Table
    `axial_dispersion_correlations`, input-table rows; bib entries `taylor1953`, `aris1956`, `wakaofunazkri1978`
    (checked against Crossref). The Stadler and Stampi-Bombelli inputs still give constant D values; switching
    them to the law is part of V-BC-2.
- [x] **F-PHY-7** Optional: outlet dead-volume model for extra-column mixing (Stadler et al. and Marx et al. both
      correct their data for it). Dropped 2026-11-03 (decision): not implemented and not reported in the paper.
- [x] **F-PHY-11** `heatCapacityGas` is a constant per unit mass, while the gas of a cycle can swing between H2-rich and
      CO2-rich (molar mass 2 to 44 g/mol; c_p of H2 is about 13 times that of the T15 feed). Since F-BC-2 the pressure
      work scales with 1/(ρ c_p) too: with helium or H2 and the value of the feed, R/(M c_p) > 1, and at the first instant
      of a blowdown the gas cooling can contract the gas more than the falling pressure expands it. The exchange with
      the solid takes this out within a second: in `examples/PSA-H2-activated-carbon`, c_p of 1109 and 14300 J/(kg K) give
      H2 purity 0.9991 either way, recovery 0.162 and 0.161, and gas temperatures within 2 K. A molar or
      composition-dependent heat capacity would remove it.
  - Done 2026-11-18 (with F-BC-3, whose energy balance needs it): component key `GasHeatCapacity` in J/(mol K), for
    every component or for none. The gas of a node then has the heat capacity of its composition,
    Σ_j c_j C_p,j / ρ per kg (`Column::gasHeatCapacity`, set where the gas density is, in `updateVelocityAndPressure`),
    in the gas-temperature balance, the thermal diffusivity, and the thermal Jacobian; `heatCapacityGas` stays the
    value without it. The kernels take the per-node values through an optional argument, so that inputs without the
    key are unchanged. Test `GasHeatCapacity.OfTheComponentsGivesTheHeatCapacityPerMassItAddsUpTo`. No paper input
    uses it yet.
- [x] **F-ITG-3 Analytical sorption Jacobian; ROS3, CVODE, and a coupled SIRK3 on it** (agreed 2026-11-02).
      The stiff mode is the gas depletion λ = k (1 + ρ_b/ε dq*/dc), explicit in RK3 and SIRK3 (B-ITG-1). Order:
      (1) analytical slopes dq/dp of the 16 pure isotherm types and mixture sensitivities dq*_i/dp_j: IAST and SIAST
      in closed form from the converged P_i⁰ in `cachedPressure` (dψ/dp_j = q_T/P_j⁰, then x_i, q_T, q_i; checked
      against finite differences in a Python prototype), SPI diagonal, EI/SEI/SCI from their explicit formulas, MPD
      rejected at first; a pressure floor where dq/dp → ∞ at p → 0 (Freundlich, Sips, and Langmuir–Freundlich with
      a power of p below one); chain rule ∂p_i/∂c_j = (P_t/C)(δ_ij − y_i) for gas, identity for liquid; tests against finite
      differences. (2) Node-major banded assembly of the c–q Jacobian with Upwind transport and frozen velocity;
      J·v against finite differences of the RHS. (3) `ROS3` integrator (template: RK3/SIRK3 on
      `FixedStepIntegrator`) on the restored kernel (`src/integrators/ros3.{h,cpp}`, `tests/ros3_tests.cpp`, all 16
      tests pass; 2026-11-02) with a banded LU; same input restrictions as SIRK3; measure the order on BEA C7 and
      JUC-77 against converged RK3. If the frozen velocity costs order, add the total molar flux per node as an
      unknown (keeps the band) or use a Rosenbrock-W method. (4) CVODE: the Jacobian as SPGMR preconditioner, then
      for a direct band/KLU solve instead of the finite-difference Dense Jacobian (B-ITG-12). (5) `SIRK3` replaced
      in place by the coupled version: stages 1–3 share one factorization, and the final I + Δt²J² stage is dropped
      (decided 2026-11-02, B-ITG-1): first order, |R| ≤ 1 in the left half-plane, R → 1/3 for stiff modes.
  - (1) done 2026-11-02: `Isotherm::slope` (+ `slopeSlowPath` for pH-Langmuir, Unilan, GAB) and
    `MultiSiteIsotherm::slope`; `MixturePrediction::predictMixtureSensitivities` (IAST and SIAST, both solvers, from
    `cachedPressure`; EI/SEI through the backward recurrence; SCI for its six models; SPI; MPD throws; the pH
    dependence of pH-Langmuir is not differentiated); slope floor p ≥ 1e-10 P, and without adsorbing gas the pure
    slopes; `computePhysisorptionEquilibriumSensitivities` (`sorption.cpp`) does the chain rule to c with P frozen.
    Tests (`mixture_prediction_tests.cpp`): all 16 slopes, every method at a feed-like and a front-toe composition
    (relative error < 2e-5 against central differences), the dilute limit, the gas and liquid chain rule. Nothing
    calls them yet; a BEA RK3 run is byte-identical. Full suite 153/154 (B-EX-3 as before).
  - Paper (2026-11-02): SI section "Analytical sorption Jacobian" (`apx:sorption_jacobian`, before the integrators):
    stiffness λ, the c–q Jacobian with Upwind and frozen velocity, the chain rule, and the SPI/SCI/EI/IAST
    sensitivities. The SIRK3 paragraph (SI, `eq:sirk3_stages` to `eq:sirk3_stability_function`), the main-text
    integrator sentence, and the input-table row describe the coupled SIRK3 without the final stage. The code still
    runs the split SIRK3 until step (5); until then the paper is ahead of the code.
  - (2) done 2026-11-02: `SorptionJacobian` and `assembleSorptionJacobian` (`src/integrators/jacobian.{h,cpp}`):
    node-major (c, q) per node, LAPACK band storage with room for the dgbtrf fill-in (2N sub- and superdiagonals),
    upwind convection with the velocity of the last precompute, local-spacing dispersion as in
    `computeMassDerivatives`, LDF and sink rows with the sensitivities of step (1); `gather`/`scatter` between the
    field-major column state and node-major vectors. Tests (`sorption_jacobian_tests.cpp`, every column against
    central differences): BEA C7 gas with a frozen-closure right-hand side, 7e-7 of the row scale; the liquid
    pH example against the production right-hand side, < 1e-6. Full suite 155/156 (B-EX-3); BEA output
    byte-identical. Against the full production Jacobian (BEA, 12 cells): FixedVelocity agrees to 6e-7,
    FixedPressureInletVelocity misses 6% (composition dependence of the Ergun velocity), and the continuity settings
    miss entries outside the band as large as the row maximum (velocity depends on all upstream uptake). Expect
    order reduction of ROS3 with continuity closures (step 3: flux unknowns or a W-method).
  - (3) done 2026-11-02: input `ROS3` (`BreakthroughIntegrationScheme::ROS3`), `Rosenbrock3` in
    `src/integrators/ros3_column.{h,cpp}` on `FixedStepIntegrator`: J assembled once per step, I − γΔtJ factorized
    with `dgbtrf`, stages through `rosenbrock3Increments` (the stage loop of the kernel, moved into `ros3.h`; the
    16 kernel tests pass), complete right-hand side at every stage, clamp and precompute at the step end. Input
    restrictions as SIRK3 (Upwind default and only, no reactions, chemisorption, energy balance, one adsorbent) and
    no MPD; test `ExampleInputReader.Ros3DefaultsToUpwindAndRejectsWhatItsJacobianOmits`, and
    `Rosenbrock3Column.MatchesRungeKutta3ThroughAnAdsorptionFrontAt50TimesTheStep`. Full suite 157/158 (B-EX-3); RK3
    output byte-identical. Not yet in the Python bindings, the Qt app, or the web app.
    Order study (Upwind, outlet max |Δc|/c_feed against converged RK3; scratch `ros3check/`): BEA FixedVelocity
    (exact J) 5.7e-7 at Δt = 0.2 s (Courant ≈ 20) down to 1.3e-9 at 0.0125 s, observed order rising to 3.3;
    liquid pH 1.1e-8 at 8 s, then the 4e-9 error of the RK3 reference; BEA FixedPressureInletVelocity 9e-6 to
    1.6e-7, order rising to 2.4; with continuity closures first order: BEA InletPressureOutletPressure 3.4e-5 to
    1.7e-6, MOR 3.2e-3 to 2.6e-4 (Δt = 0.04–0.0025 s), JUC-77 4.5e-4 at 0.01 s to 1.2e-4 at 0.0025 s (14 %
    pointwise at the ethylbenzene roll-up for 0.02–0.04 s, t50 within 3e-3). All runs stable. JUC-77: RK3 diverges
    at 3e-3 s and is 4e-3 off at its default 5e-4 s; ROS3 at 0.01 s is 4.5e-4 off in 28 s of wall time (RK3 to the
    same 300 s: about 140 s). Open: restore third order with continuity closures (flux unknowns or W-method).
  - (3b) W-method, decided and done 2026-11-02: input `ROS34PW2` (Rang & Angermann 2005, four stages, third order
    for any approximation of J, L-stable, stiffly accurate), next to `ROS3`. The kernel stage loop is now
    `rosenbrockIncrements<Method>` for any number of stages (ROS3 still reuses its second evaluation);
    `Rosenbrock34PW2Coefficients` holds the stage form derived from PETSc's TSROSWRA34PW2 tables (only source found;
    SciML has ROS34PW3, not PW2); the column integrator is `ColumnRosenbrock<Method>` with the aliases `Rosenbrock3`
    and `Rosenbrock34PW2`. Tests: `Rosenbrock34PW2.StageCoefficientsFollowFromTheTablesOfRangAndAngermann`,
    `Rosenbrock34PW2.KeepsThirdOrderWhenTheStageMatrixOnlyApproximatesTheJacobian` (orders 2.85–3.15 for W = J,
    J/2, 0, and J without coupling), `Rosenbrock34PW2Column.MatchesRungeKutta3WhenTheVelocityFollowsFromContinuity`
    (1.1e-5 of the feed at 50 times the RK3 step); full suite 158/159 (B-EX-3); RK3 output byte-identical.
    Order study (same runs as ROS3): MOR (continuity, Langmuir) 4.5e-4 at Δt = 0.04 s to 3.5e-7 at 0.0025 s, orders
    2.2, 2.6, 2.8, 2.9 (ROS3: 3.2e-3 to 2.6e-4, order 0.9); BEA InletPressureOutletPressure 2.7e-5 to 1.9e-7, order
    rising to 2.6 (ROS3 1.7e-6 at 0.0125 s); BEA FixedVelocity and FixedPressureInletVelocity as ROS3; JUC-77 5.8e-5
    at 0.01 s against 4.5e-4 for ROS3, but both converge at order 1.3 below 0.01 s (also in self-convergence, so
    not the reference): Langmuir–Freundlich with ν < 1 makes q* ∝ p^ν non-differentiable at the front toe. Cost per
    step 1.5 times ROS3 (MOR, 400 s at 0.01 s: 8.4 s against 5.6 s, sequential, load ≈ 9); RK3 at its default 5e-4 s
    needs 60 s for the same run. Decided 2026-11-02: keep both. Open: bindings, apps, paper.
  - (4) CVODE SPGMR preconditioner, 2026-11-02: input `CVODEPreconditioner` `None` (default) | `SorptionJacobian`
    (only with SPGMR; not with MPD). Left preconditioner: banded LU of I − γJ on the c and q blocks (setup
    re-evaluates the right-hand side at the trial state when CVODE asks for a fresh Jacobian, else refactorizes
    with the new γ), identity on chemisorption, pore, and temperature blocks; `CVODE::UserData` carries column and
    integrator; `assembleSorptionJacobian` takes the concentration and gas-temperature spans. Test
    `ExampleInputReader.ParsesTheCvodePreconditionerOnlyWithSpgmr`; full suite 159/160 (B-EX-3); CVODE without
    preconditioner byte-identical (BEA SPGMR, 60 s, 1792 steps). Benchmark inputs of `benchmark/breakthrough`
    (WENO3, sequential runs, scratch `precheck/`), RHS evaluations incl. Krylov / wall time:
    BEA rtol 1e-5 182k → 94k, 14.6 → 13.7 s; rtol 1e-6 179k → 68k, 13.8 → 9.9 s; monolith rtol 1e-5 500k → 166k,
    7.3 → 2.5 s (Newton failures 640 → 54); rtol 1e-6 698k → 245k, 10.1 → 3.7 s; T15 (energy balance) worse,
    17k → 31k, 2.2 → 4.6 s (steps 2278 → 4239). Right preconditioning: same accuracy, monolith 400k RHS.
    Open problem: with the preconditioner, BEA at rtol 1e-6 deviates from RK3 (Δt 1e-4) by 2e-4 of the feed
    (e50 1.6e-5) against 1.9e-6 without it; the monolith at 1e-6 by e50 6e-5 against 1e-7. The deviation starts
    at t ≈ 80 s in every component behind the fronts and then persists. Ruled out (variants, sources reverted):
    the right-hand side is a function of the state (4e-14), no FastIAST solve stops at 50 iterations, not the BDF
    order (same at CVODEMaximumOrder 2 without error-test failures). An identity solve deviates too when the setup
    reports jcur = true, but not with jcur = false; the real solve deviates with either. Unpreconditioned runs at
    rtol 0.8e-6 to 2e-6 all stay within 2.2e-6. Keep the default `None` until this is understood.
    Diagnosed 2026-11-02 (scratch source copy `src-copy/`, `precheck/dig/`): a 10 times tighter linear tolerance
    (CVodeSetEpsLin 0.005) removes the BEA drift (e50 1.1e-7 at rtol 1e-6, 103k RHS against 179k) but not the
    monolith one at 1e-6 (e50 9.4e-5) nor the T15 slow-down; a tighter Newton coefficient (0.01) makes it worse. The
    drift grows steadily from t ≈ 80 s (1.7e-6) to 400 s (1.9e-4). Cause: the sorption Jacobian omits the velocity
    closure terms (FixedPressureInletVelocity: composition dependence of the Ergun velocity, 6 %; continuity
    settings: the velocity at a node depends on all upstream uptake). CVODE accepts the Newton iteration after about
    one pass per step, and the iteration error left by the missing terms has the same sign every step. Confirmed
    with FixedVelocity (exact Jacobian, BEA Upwind, 2400 s): unpreconditioned, preconditioned, and KLU all agree
    with the rtol 1e-8 reference to e50 4e-7 to 1e-6. At CVODEMaximumOrder 1 (the new default) BEA at rtol 1e-6 is
    2e-4 off with or without preconditioner: the first-order error of BDF1, from the first output on.
  - (4b) KLU, 2026-11-02: `CVODELinearSolver` `KLU` (CVODE only; no reactions, chemisorption, energy balance, or MPD):
    CSC pattern of the field-major CVODE state built once (local c–q block, same-component neighbours, every
    diagonal), values copied from the banded sorption Jacobian through a stored index map; the Jacobian callback
    re-evaluates the right-hand side at the trial state. Tests
    `ExampleInputReader.ParsesTheKluLinearSolverForPhysisorptionOnly` and
    `CvodeKlu.MatchesSpgmrThroughAnUpwindFrontWithFewJacobians`; full suite 161/162 (B-EX-3). BEA order 5:
    FixedVelocity (exact J) 12.6k RHS and 2.1 s against 74k and 8.0 s for SPGMR, same accuracy;
    FixedPressureInletVelocity Upwind 29.9k RHS and 5.0 s against 131k and 16.3 s, but e50 1.3e-4 at rtol 1e-5
    (SPGMR 1.8e-6; at 1e-6 both fine, 3.5e-9); with WENO3 the upwind J is too far from the true one (8 times the
    steps). The monolith has chemisorption sites, so KLU rejects it. Next: the closure derivatives (local for
    FixedPressureInletVelocity; for continuity the flux-per-node augmentation of the linear system), and a WENO3
    Jacobian with frozen weights.
  - (5) 2026-11-03: not done. SIRK3 is removed instead (B-ITG-1): the coupled version without the final stage would
    be first order with the same ingredients as ROS34PW2, which is third order.
  - (6) 2026-11-03: ROS3 removed (`ros3.{h,cpp}`, `ros3_column.{h,cpp}`, `tests/ros3_tests.cpp`); `ROS34PW2`
    (`src/integrators/ros34pw2.{h,cpp}`, `BreakthroughIntegrationScheme::ROS34PW2` = 2) is an adaptive
    integrator of the whole column state. Stages with the RK3 right-hand side (precompute + computeDerivatives); W =
    I − γhJ~ with J~ = the sorption Jacobian, whose transport rows are now exact for Upwind and WENO3 (nonlinear
    weights, mole fractions, face-fraction normalization, upwind total flux, `computeWENOLinearization`; Danckwerts
    chain rule for every row next to the inlet) plus a thermal block (`assembleThermalJacobian`, 3 × 3 node blocks);
    zero on chemisorption, pore, and surface states. The loadings are eliminated node by node, the concentrations and
    temperatures factorized by a block LU (`block_banded.{h,cpp}`) once per step attempt. Embedded error estimate in
    the WRMS norm of CVODE with `ROS34PW2RelativeTolerance` (1e-6) and the three absolute tolerances (1e-7 mol/m³,
    1e-8 mol/kg, 1e-5 K: ten times tighter than CVODE, whose 4 dz/v step cap adds accuracy; the example PSA is
    2.6e-3 off at rtol 1e-5, 3.1e-4 at 1e-6, 2.1e-5 at 1e-7, CVODE 5.6e-5), `ROS34PW2MaximumTimeStep` (0 = unlimited); rejected steps repeat from the step start; outputs split the
    interval into equal steps; `reinitialize` after a swing phase change. `TimeStep` is rejected like for CVODE
    (`rejectStepBasedTimeKeys`, `convert.py`). Every input restriction is gone: WENO3, Danckwerts, reactions,
    chemisorption, pore transport, energy balance, several adsorbents, MPD (no sensitivities: J~ keeps only the
    mass-transfer rate). The ∂F/∂t term is dropped (a W-method may approximate it). The KLU pattern now includes
    the WENO3 couplings; the CVODE preconditioner uses the exact WENO3 transport (not studied).
    Prototype findings (scratch `proto/`, figure-10 inputs, rtol 1e-5): with the upwind J~ under WENO3, BEA needs
    24.8k steps (about one cell crossing per step) against 774 with the exact transport; frozen WENO3 weights still
    15k steps with 27% rejections (the carrier rows miss the total-concentration and normalization coupling); without
    the thermal block, T15 takes 16.5k steps with 20% rejections against 2952; a consistent two-pass closure (as
    CVODE) changes nothing, so one RK3-style pass per stage is kept. Profile: T15 48% IAST, 19% block LU and solves,
    8% Jacobian.
    Tests: `ros34pw2_tests.cpp` (coefficients, third order for W = J, J/2, 0, uncoupled with time-dependent
    forcing, block LU, columns against RK3 with WENO3, energy balance and continuity, chemisorption, a liquid with a
    Danckwerts inlet), `SorptionJacobian` WENO3 gas and liquid, Danckwerts with both schemes, `ThermalJacobian`, input
    tests. Full suite 174/175 (B-EX-3); lab tests 31/31; webapp 77/77; convert 22/22. RK3 and CVODE outputs
    byte-identical (BEA, 100 s). Python bindings (`ruptura.ROS34PW2`, `Simulation.ros34pw2`), Qt app, web app,
    README, AI prompt, example `examples/BEA-alkanes-C7/ros34pw2` (whole breakthrough in 330 steps, 0.22 s), paper.
  - (7) Figure `fig:benchmark` (2026-11-03): ROS34PW2 series in `benchmark/breakthrough` (rtol 1e-2 … 1e-6, absolute
    tolerances scaled like CVODE), and the RK3, CVODE SPGMR, reference, and stability runs rerun with the same frozen
    binary (51 runs, 4 at a time, 40 min; CVODE steps identical to 2026-10-28, wall times 16–38% shorter; Dense,
    components, and methods kept from 2026-10-28). At the defaults: BEA ROS34PW2 1.44 s, e50 9.7e-8 against CVODE
    15.5 s, 2.3e-7; T15 3.27 s, 4.9e-7 against 2.83 s, 1.0e-6; monolith 0.53 s, 3.8e-7 against 7.92 s, 4.0e-4.
    Notebook panels (a)–(e) with triangles; Results text, caption, and SI benchmark section updated;
    `InletCondition` `Dirichlet` pinned in `generate.py`. The ROS34PW2 e50 on the monolith stays at 3.8–4.5e-7
    for every tolerance and step cap (3.1e-7 at rtol 1e-9, steps ≤ 0.02 s), although converged ROS34PW2 and RK3
    agree to 1e-7 on 12 cells: not understood, small.
  - (8) 2026-11-03, after the validation studies moved to ROS34PW2 (I-BC-4, I-BC-5): the abstract names the three
    integrators and adds that ROS34PW2 is an order of magnitude faster than CVODE where convection or dispersion
    limits the explicit step; the Conclusions make ROS34PW2 the integrator of choice where the explicit step is
    limited by stability (about as fast as CVODE in the non-isothermal case), and the cost-per-call sentence names
    RK3 and CVODE; the swing section says the phase driver works with each of the three integrators.
- [x] **F-BC-1 Flow reversal.**
  - Paper (SI species and temperature boundary conditions): direction-aware Danckwerts inflow/outflow
    and one-way pressure ports.
  - Code: the WENO reconstruction always takes z = 0 as the inlet (`src/integrators/compute.cpp:872`).
  - Needed for counter-current PSA steps.
  - Done 2026-11-04 (cycles plan, stage 2; Upwind and liquid included, as decided):
    - Input: `FlowDirection` `Forward` (default) / `Reverse`, top level and per swing phase; per forward phase
      `ProductEndGasPhaseMolFractions` (gas that flows back in at z = L; default the gas of the last cell, the agreed
      option a). A reverse step feeds at z = L from `Column::productEndConcentration`, with the inlet settings at
      z = L and the outlet settings at z = 0; node 0 is not used (stencils take cell 1, zero gradient) and the swing
      driver sets it to cell 1 at the start of a reverse phase.
    - Transport (`computeMassDerivatives`): every face is upwinded by its own velocity sign (`computeWENOReverse`,
      the mirror image of `computeWENO`; upwind takes cell i + 1), Danckwerts/Dirichlet inflow at z = L with the
      mirrored distance, the gas total flux with the upwind c_tot. Energy (`GasTemperatureBalance::rate(grid,
      v_{i-1}, v_i)`): -max(v_i,0)(T_i - T_{i-1})/h + min(v_{i-1},0)(T_i - T_{i+1})/h, mirrored Danckwerts T_b at
      z = L with -v_{N-1}. Closures (`updateVelocityAndPressure`): march from z = L in reverse (one direction-aware
      face-velocity solve, secant at the node next to a Danckwerts inlet), flux = v x c_tot upstream of the face, so
      that a forward march can also turn negative at the outlet (B-CYC-1); Ergun integrated from the prescribed end
      with the mirrored interval rule; liquid, IPOP shooting, FixedVelocity, and FPIV mirrored too.
    - Jacobians: sorption (WENO faces by sign, virtual node N + 1 for the inflow face) and thermal follow the same
      rules. `BlockBandedMatrix` takes several upper blocks; ROS34PW2 uses 2 + 2 for WENO3 and processes the second
      upper block only in steps with a negative face velocity. The CVODE/KLU pattern is symmetric (+-2 nodes).
    - Forward runs: isothermal and liquid inputs byte-identical to `ruptura-2026-11-04-bbc2`, energy-on runs within
      1e-11 (ROS34PW2, RK3; FMA fusion), CVODE within 1e-5 (rounding sensitivity); CPU time unchanged within noise
      (BEA ROS34PW2 0.62 vs 0.63 s, T15 1800 s 3.60 vs 3.60 s).
    - Tests: `FlowReversal.*` (mirror image of the mass and temperature derivatives, of the closures with the
      thermal expansion, and of a whole RK3 and ROS34PW2 breakthrough on a uniform grid without pressure drop;
      back-flow at the product end; a forward feed + reverse helium purge), `FlowDirectionInput.*`,
      `SorptionJacobian.MatchesCentralDifferencesInReverseFlowAndWithBackFlowAtTheProductEnd`, the thermal Jacobian
      test over forward/reverse/back-flow, `BlockBandedMatrix.SolvesLikeADenseEliminationWithTwoBlocksOnEitherSide`,
      `LabInput.FlowDirectionsOfTheColumnAndTheSwingPhasesRoundTrip`.
    - GUIs (Qt combo top level and per phase; webapp schema and phase field), bindings, docs prompt, paper SI
      (reservoirs, reverse discretization, closures; the one-way-port text replaced by "not available").
    - Not done here (stage 3, F-BC-2): closed ends, pressure ramps, re-equilibration of the inventory at a switch.
    - Validated 2026-11-06 against equilibrium theory (V-BC-3, `validation/3.4.4_flow_reversal`, paper Section 3.11).
- [x] **F-BC-2 Closed ends and pressure ramps** (cycles plan, stage 3; agreed 2026-11-04 with the key names
      `ClosedEnd` and `PressureRamp`).
  - Done 2026-11-09:
    - Input, per phase: `"ClosedEnd": "FeedEnd" | "ProductEnd"` (z = 0, z = L) and
      `"PressureRamp": {"Shape": "Exponential" | "Linear", "TimeConstant": τ}`. Gas enters at the open end, so that
      ClosedEnd sets the direction (`FlowDirection` rejected, as are `ColumnEntranceVelocity` and
      `ProductEndGasPhaseMolFractions`). Both need a gas and `InletPressureInletVelocity`; ClosedEnd needs Danckwerts.
    - Closures: the march starts from the closed end with v = 0, so that march and flow direction are independent; the
      pressure follows from `inletPressure` at the open end. Decided against the planned accumulation (dP_ref/dt)/(R T_i):
      the closure uses C_i (dP/dt)/P, equal when C_i = P/(R T_i), so that a mismatch from a change of the pressure drop
      stays a small fraction when the pressure falls by 10, instead of a growing one.
    - Energy: the pressure work (dP/dt)/(ρ c_p) in the gas (decided 2026-11-09). The feed node of a forward ramp has the
      rate c_0 (dP/dt)/P, also in the sorption Jacobian; the product-end gas of a reverse ramp is scaled directly.
    - Time: `Column::advancePressureRamp(t)` before each right-hand side, the only time dependence. ROS34PW2 leaves out
      dF/dt (a W-method may), and the linear ramp needs no integrator restart at its end: CVODE at rtol 1e-5 is within
      2.8e-6 of the reference without one, and 1.7e-5 with an explicit split, at the same cost.
    - Pressure step without a ramp: F-CYC-2 below.
    - `phases.data` and `SwingAdsorption.phase_records`: times, end pressures, gas and adsorbed amounts per phase in
      mol/m² of cross-section (not mol: `InternalDiameter` is optional), from which purity and recovery follow.
    - Tests (`pressure_swing_tests.cpp`, `pressure_swing_test.py`): inert closed-end column against
      v(z) = (L − z)(dP/dt)/P both ways, closure and species balances over the column in both directions, RK3, CVODE and
      ROS34PW2 on both ramp shapes, Skarstrom mass balance against the outlet files, input;
      `SorptionJacobian.MatchesCentralDifferencesInClosedEndPhases`, closed ends in the thermal Jacobian test.
    - Inputs without a pressure step or ramp are byte-identical (BEA MPD IAST, T15 CVODE and RK3, CO2/N2 chemisorption,
      3.4.2 T15 TSA, 3.4.4 flow reversal; BEA TSA over the first 280 s), CPU time within 1 %, and so is Stadler's
      FixedVelocity 8 to 1 bar cycle. The Ruptura Lab swing preset (1 to 0.2 bar) changes as intended: max |R T Σc/P − 1|
      from 4.0 to 3e-4. `examples/BEA-alkanes-C7/psa` now stops: B-EX-13.
    - Example: `examples/PSA-H2-activated-carbon` (CVODE, 40 s linear ramps), 20 cycles in 24 s: H2 purity 0.9991,
      recovery 0.162, max |R T Σc/P − 1| 7.8e-3. With 5 s exponential ramps it was 0.961 and 0.186 at 0.02: the feed
      enters faster than mass transfer takes it up, and the pressure drop of the fast flow breaks the small-drop limit.
    - Rescale at a switch also when the end where the pressure is set changes (2026-11-09): a blowdown open at z = 0
      followed by a purge set at z = L moves the profile by the pressure drop, and without a rescale R T Σc/P − 1 grew by
      1e-3 per cycle (it scales with the pressure drop; thermal, dispersion, and the convection scheme play no part).
- [x] **F-BC-3 `MomentumBalance`: the pressure of the gas, the velocity of the pressure-drop law, and ports** (stage 2
      of the plan of 2026-11-18 for the mechanics of the branch `ruptura-shrinjay`; done 2026-11-18).
  - Model: the pressure of a cell is that of its gas, P = R T Σc (all of it, also a concentration slightly below zero:
    the clamp of the mole fractions put a kink in the velocity at every absent component), and the velocity of the
    face between two cells is that of `PressureDropLaw::velocity` for their pressure difference over the distance of
    their centres, with their mean density (`updateVelocityAndPressure`). No continuity closure; the species, energy,
    sorption, multibed, and monolith kernels are those of the other boundary conditions, with the face-by-sign
    transport of F-BC-1. The column starts at rest at `OutletPressure`.
  - Ports per swing phase, `FeedEnd` (z = 0) and `ProductEnd` (z = L): `Closed` (default), `Pressure` (gas at a
    pressure beyond the end, the law over the half cell; optional `PressureRamp` from the pressure of the end cell),
    `Velocity` (`Velocity` with the `Pressure` of its gas, or `MolarFlux`). Gas that enters is the phase
    `GasPhaseMolFractions` at z = 0 and `ProductEndGasPhaseMolFractions` at z = L. Without phases: the feed at
    `ColumnEntranceVelocity` and gas at `OutletPressure` at z = L.
  - Energy: the pressure work follows the pressure of the cell, dP/dt = R (T dC/dt + C dT/dt), so that
    dT/dt = (rate + κ T (dC/dt)/C)/(1 − κ), κ = R C/(ρ c_p) (`computeDerivatives`); needs a molar c_p above R
    (F-PHY-11), else the run stops with a message.
  - Jacobian (`assembleSorptionJacobian`): the derivatives of every face velocity with respect to the concentrations
    of its two cells (pressure and mean density), exact against central differences for both schemes and all port
    modes (`SorptionJacobian.MatchesCentralDifferencesWithMomentumBalance`); the equilibrium sensitivities at the
    pressure of the gas instead of a frozen one. The upwind bandwidth is 3 N.
  - ROS34PW2 (found by runs; each of these was needed):
    - with the energy balance one block-banded W on the concentrations and the three temperatures of a node: the gas
      temperature is part of the pressure, the velocity carries it, and the pressure work heats it. Separate solves
      held the step at 1e-6 s, a block-triangular one at 3e-5 s; coupled, T15 with the energy balance takes 7107 steps
      for 2000 s (1.2 s);
    - the time derivative of the right-hand side while a port pressure follows a ramp (γ h² γ_i dF/dt, γ_i from the
      table in the code), by a difference quotient over 1e-6 of the step: a larger increment (1e-3) left the step
      at 4 ms, because the pressure difference that drives the gas is small against the pressure;
    - the stream amounts follow the change that W makes to the face velocities, and are left out of the step control;
    - the non-negativity clamp leaves the (signed) stream amounts alone.
  - Checked by runs: T15 breakthrough against `InletVelocityOutletPressure` (30 cells, ROS34PW2): t50 of CO2, N2, H2
    711.68/224.89/123.67 s isothermal and 636.52/206.60/119.57 s with the energy balance for both, max |Δy| 5e-5 and
    2e-5, peak temperatures equal (gas 322.079 K, solid 333.875 K), pressure drop 31 Pa for both.
    `examples/PSA-H2-activated-carbon-ports` (the PSA-H2 cycle with ports, purge with the product, 50 cells): cyclic
    steady state (1e-3) after 20 cycles in 17 s, mass balance 8e-6, H2 purity 0.9988 and recovery 0.163; with a
    pure-H2 purge 0.9991 and 0.162, the values of the closure cycle of `examples/PSA-H2-activated-carbon`.
  - Tests `MomentumBalance.*` (against the closure; an inert column filled and emptied through pressure ports with
    P V/(R T) and closed ports; the example cycle), `CycleInput.*`.
  - Integrators other than ROS34PW2: B-ITG-21. Not done: valves, several columns, buffers (F-CYC-4); Lab controls
    (F-GUI-1); validation against a literature cycle (V-CYC-3); paper text (I-CYC-2). Not taken from the branch: its
    `InletDensityBasis` `LocalPressure` (a velocity port here takes the pressure given with it), gas-solid thermal
    equilibrium (a large `heatTransferGasSolid` does it), the transient coolant, HSDM, TVD/WENO5.
- [x] **F-ITG-5 `ConvectionScheme` `WENO5` and `TVD`** (2026-11-19, user request; stage 5 of F-CYC-6). TVD: the
      smaller of the two slopes over a cell where they have one sign (MUSCL, minmod), on the quantities and faces of
      WENO3. WENO5: the Jiang-Shu weights of the three quadratic candidates on a grid of one cell width, WENO3 next to
      the ends and on other grids. Both in `computeFaceValues` (`compute.cpp`), for both flow directions. The routines of
      WENO3 (`computeWENO`, `computeWENOReverse`, `wenoFace`, `wenoStencil`) are untouched: two first versions that
      merged the schemes into them changed WENO3 results in the last bits (T15 2e-11, BEA 2e-8), which CVODE amplified
      to 0.3 % at the fronts of the PSA-H2 reference case. The Jacobian, the KLU pattern, and the band of ROS34PW2 are those of WENO3 (`shortcut:` in
      `jacobian.cpp`).
  - Checked: orders 2 and 5 on a smooth profile, linear profiles exact, a step within its bounds, the mirror image
    for reverse flow, the species fluxes summing to the upwind total (`ConvectionSchemes.*`,
    `MassBalanceWENO3.GasSpeciesFluxesSumToTheUpwindTotalFlux`). T15 on 40 cells, CO2 front width 10-90 %:
    Upwind 179.9 s, TVD 165.7 s, WENO3 164.7 s, WENO5 163.6 s, equal for ROS34PW2 and CVODE within 0.01 s; ROS34PW2
    takes 3058/3919/3739/6118 steps.
  - Open: the derivatives of TVD and WENO5 themselves, if the 60 % more steps of WENO5 matter.
- [x] **F-PHY-12 Homogeneous surface diffusion in the particle** (2026-11-19; stage 5 of F-CYC-6):
      `"ParticleTransport": {"Model": "HSDM", "RadialCells": N}`. A new block of the state after the temperatures holds
      the loading of N radial cells per node and component (`ColumnStateLayout::shellLoading`); the physisorption state
      stays the loading of the particle and takes the volume mean of the cell rates, so that the sink, the heat, and
      the outputs are unchanged. Finite volumes in r/R (`computeDerivatives`, `rk3.cpp`), surface at equilibrium with
      the gas, D_s/R^2 = `MassTransferCoefficient`/15, so that no new parameter is needed and the temperature law of k
      carries over. The restart file holds the cells; one without them starts each cell at the particle loading.
  - Checked: the uptake of a sphere at a constant surface loading follows its series within 0.01 at D_s t/R^2 = 0.05
    and 0.2 with ROS34PW2 and CVODE, and the volume mean of the cells is the particle loading (`ParticleTransport.*`).
  - Limits: spheres, no film resistance, not with reactions or CVODE KLU; ROS34PW2 and CVODE's preconditioner treat
    the cells explicitly (step below about R^2/(D_s N^2)); the cyclic-steady-state measure of the state leaves the
    cells out (it has the particle loading).
- [x] **F-PHY-13 Direction-dependent pseudo-nth-order kinetics** (2026-11-19; stage 5 of F-CYC-6): `PseudoNth` takes
      `DesorptionRateCoefficient` (a number or a temperature law with its own activation energy) and `DesorptionOrder`
      for q > q*; both default to those of the uptake (`Chemisorption::releaseRateCoefficient`, `releaseOrder`). The
      dual-order kinetics of the branch are two such sites. Tests `ChemisorptionRate.*`, `TemperatureInput.*`.
- [x] **F-PHY-14 Heat-transfer fluid along the outer wall** (2026-11-19; stage 5 of F-CYC-6): root `Coolant`
      (`MassFlowRate`, `Volume`, `Density`, `HeatCapacity`, `FlowDirection`), per phase `MassFlowRate` in
      `ThermalControl`. A scalar block of the state (`coolantTemperature`), donor cell from the inlet at the
      temperature beyond the wall, each cell exchanging heat with the wall of its cell; the wall term of the energy
      kernel is corrected for the fluid temperature after it (`computeDerivatives`). The wall-heat duty is then the
      heat that the fluid gave to the wall.
  - Checked (`Coolant.*`): the fluid cools on its way, counter-current flow is the mirror image (0.5 K), the column
    ends at the inlet temperature, and the wall takes the heat of wall and solid (1 %).
  - Limits: no axial conduction in the fluid, no quasi-steady variant, explicit in ROS34PW2 (step below about the
    time the fluid needs for one cell).
- [x] **F-CYC-4 Several columns, valves, buffers, compressors; thermal control and duties** (stages 3 and 4 of the plan
      of 2026-11-18 for the mechanics of the branch `ruptura-shrinjay`; done 2026-11-18, none of its code).
  - Input: root `Columns` (names of identical columns), `Buffers` (`Volume`, `Pressure`, `Temperature`,
    `GasPhaseMolFractions`, `MaximumPressure`), `ExhaustPressure`, `CompressionEfficiency`; per phase `Columns`
    (per-column `FeedEnd`, `ProductEnd`, `GasPhaseMolFractions`, `ProductEndGasPhaseMolFractions`, `Temperature`,
    `ThermalControl`), `Connections` (`From`, `To`, `FlowModel` `Linear`/`SquareRoot`/`Orifice`/`Compressor` with
    `Conductance`/`Area`/`MolarFlow`, `Opening`, `CheckValve`, `MaximumMolarFlow`, `Efficiency`), `ThermalControl`
    (`Temperature`, `HeatTransferCoefficient`); `StopCriterion` takes `Column`. The reader holds one phase entry per
    phase and column (`InputReader::swingAdsorptionPhases`), built from the keys of the phase with those of the column
    in their place. Needs `MomentumBalance`, ROS34PW2, and the cross-section (`InternalDiameter`).
  - Driver (`SwingAdsorption::run`): one `Breakthrough` per column, each with its integrator and its output
    directory. The connections are advanced in coupling steps: flows from the pressures at the step start (second
    order in time from the two last), a step of a twentieth of the fastest equalization time, a column that gives
    gas before the one that takes it, whose inflow is what the stream amounts of the donor say left it. A connected
    end is a fourth port mode, a set molar flux (`Column::Port::Mode::Flux`). Buffers are mixed and isothermal.
    Stop criteria: the coupling step in which they turn is bisected. Cyclic steady state over all columns and buffers.
  - Duties: the stream block of the column state holds two more integrals, the heat through the outer wall and the
    work that brings gas leaving through a pressure port below `ExhaustPressure` to it (`phases.data`,
    `wall_heat`/`exhaust_work`); the compressor work per connection is in `connections.data`, the buffers in
    `buffers.data`.
  - Checked by tests (`Connections.*`, `ThermalControl.*`, `ExhaustPressure.*`, `CycleInput.*`): two inert columns
    even out at exp(−2 k t/C) (stop time ln(2) C/(2 k) within 0.02 s) to the mean pressure, and what leaves one enters
    the other to 1e-12; a tank with `MaximumPressure` and a compressor, whose work matches the integral of the
    isentropic work within 1 %; a jacket heats a closed column by the heat capacity of wall and solid (1 %) and takes
    the same heat back; the exhaust work of an isothermal evacuation (1e-3).
  - `examples/PSA-H2-activated-carbon-two-beds` (two beds, product tank, purge valve, pressure equalization, 50
    cells): cyclic steady state after 22 cycles in 68 s, mass balance 1e-5, the two halves of the cycle equal within
    5e-5, H2 purity 0.9956 and recovery 0.222 (one column with ports: 0.9988 and 0.163).
  - The five reference cases of F-CYC-1 are byte-identical to the binary before this item.
  - Not done: F-CYC-5 (limits), F-CYC-6 (stages 5 and 6), F-GUI-1, V-CYC-3, I-CYC-2; the Python module is not built.
- [x] **F-CYC-1 Cycle driver.**
  - Paper (introduction): cyclic-steady-state stopping, process-dependent step termination, and
    per-step feed composition, velocity, outlet pressure, and flow direction.
  - Code: phases accept `Name`, `Temperature`, `InletPressure`, `NumberOfTimeSteps`, and, since 2026-10-21,
    `GasPhaseMolFractions` (feed composition) and `ColumnEntranceVelocity`. A top-level `NumberOfCycles`
    repeats the phase list. Still missing: cyclic-steady-state stopping, step termination on a condition,
    outlet pressure, and flow direction.
  - 2026-11-04: per-phase `FlowDirection` done (F-BC-1); closed ends and pressure ramps are stage 3.
  - 2026-11-09: closed ends and pressure ramps done (F-BC-2). `phases.data` gives what a cyclic-steady-state stop needs;
    runs stop after `NumberOfCycles` until then. `examples/Self-driving-lab-H2-PSA/experiment.py` does the stop in
    Python: blocks of 10 cycles, each restarted from `column.json` with `ReadColumnFile`, until the column amounts and
    the KPIs change by less than 2e-3 from one cycle to the next. A stop in the driver would do the same in one run.
  - Section 2.3 of the paper describes the old parser; update it (I-GEN-1).
  - Done 2026-11-18 (stage 1 of the plan for the mechanics of the branch `ruptura-shrinjay`):
    - Step termination on a condition: per phase `StopCriterion` (`Metric` `Pressure`, `AveragePressure`,
      `MaximumTemperature`, `MinimumTemperature`, `MoleFraction`, `AmountIn`, `AmountOut`; `End`, `Component`,
      `Operator`, `Value`; `{"All": [...]}`, `{"Any": [...]}`), `MinimumDuration`, and `TimeoutPolicy`; `Duration` is
      then the longest the phase lasts. The criteria are one residual of the state (`StopResidual`,
      `SwingAdsorption::run`): CVODE finds its root (`CVodeRootInit`), and `Breakthrough::advanceTo` takes RK3 and
      ROS34PW2 step by step, repeats the step in which it turns, and bisects its end (1e-6 of the step). On T15 the
      three integrators end a phase on the outlet CO2 mole fraction at 637.0858, 637.0863, and 637.0851 s.
    - Stream amounts: the net amount of each component that entered at z = 0 and at z = L in a phase, as a block of
      the column state (`ColumnStateLayout::includeStreamAmounts`, rates in `computeMassDerivatives`), so that every
      integrator carries them; in `phases.data` and `phase_records[i].stream_amounts` as entered and left per end.
      The block exists only with a stop criterion, `CyclicSteadyState`, a phase name as gas, or `MomentumBalance`:
      other inputs are byte-identical (PSA-H2 CVODE, T15 ROS34PW2, BEA ROS34PW2, 3.4.4 flow reversal, BEA TSA).
    - `GasPhaseMolFractions` and `ProductEndGasPhaseMolFractions` take the name of a phase: the mean gas that left the
      column in it the last time it ran (a purge or a pressurization with product).
    - Cyclic steady state: root `CyclicSteadyState` `{Tolerance, ConsecutiveCycles}`; `NumberOfCycles` is then the most.
      Per cycle `cycles.data` holds the largest relative change from the cycle before of a state value, of the amount
      of a component in the column, and of a stream amount, each relative to the size of the value plus the largest
      of its kind, and the amount that the streams and the column do not account for. The paper form (main2 Eqs.
      `cycles:state-residual`, `cycles:inventory-residual`) has an absolute tolerance per state; here that is the
      tolerance times the largest value of the state field.
    - Tests `StopCriterion.*`, `CyclicSteadyState.*`, `CycleInput.*` (`pressure_swing_tests.cpp`). `SwingAdsorption::run`
      now prints the integrator statistics, as `Breakthrough::run` does.
    - Not done: controls in Ruptura Lab (F-GUI-1), the self-driving lab on the new stop (F-PY-3), the paper text
      (I-CYC-1, I-CYC-2), and the comparison with the 1 % criterion of the 3.4.2 T15 cycles (V-CYC-3).
- [x] **F-CYC-2** Optional swing-phase extras:
  - [x] per-phase axial dispersion (not needed for Stadler: numerical dispersion v Δz/2 exceeds Eq. 12); covered by
    F-PHY-6 since 2026-11-03: the reference-state scaling follows the temperature and pressure of each phase;
  - [x] rescale the gas inventory at a pressure switch (FixedVelocity already imposes the new pressure on the
    equilibria), and at a temperature switch without the energy balance: in the isothermal TSA of B-CYC-1,
    R T Σc / P − 1 is 25 % through the 373 K phase, because Σc keeps the 298 K inventory (stage 3 of the cycles
    plan). Done 2026-11-09 with F-BC-2: c_i ← c_i P/(R T Σc) at a step without a ramp, for the boundary conditions with
    a continuity closure; not for `InletVelocityOutletPressure`, whose column pressure does not follow the phase.
    A step shows in `phases.data` as the gas at the phase start minus that at the end of the previous phase;
  - [x] end a phase on an outlet mole-fraction threshold (Stadler et al.'s pilot plant): `StopCriterion` with
    `Metric` `MoleFraction` (F-CYC-1, 2026-11-18).
- [x] **F-FIT-1** `FittingWeighting` (2026-11-11, user request): `"Absolute"` (default, unchanged) minimizes the RMSE
      of the loading in mol/kg, `"Relative"` that of (q_data − q_fit)/q_data and skips points with q_data ≤ 0. The
      absolute RMSE is set by the plateau and leaves the Henry region of the shipped BEA nC7 and CoBDP p-xylene fits
      25–80 % off; with Relative, BEA nC7 (Toth) is within 7 % in the Henry region and 6 % on the plateau. A model that
      cannot follow the data (one Langmuir on the CoBDP step) gets its error moved into the step instead (max 55 %
      between 25 and 75 % of q_max, against 14 % with Absolute). Weighting by the slope dq/dln f of the trial isotherm
      was rejected after a prototype: the optimizer moves the step out of the data and puts the weight on the
      lowest-loading point (87–97 % off around f50); slope weights fixed before the fit gained only a few percent at
      f50. Input reader, `Fitting`, Python bindings, Qt key lists, webapp key map, and the generative-AI prompt; tests
      `ExampleFitting.RelativeWeightingScoresKnownBeaLangmuirFitByRelativeResiduals` and
      `ExampleFitting.FittingWeightingDefaultsToAbsoluteAndRelativeSkipsZeroLoadings`.
  - 2026-11-11: SI subsection "Isotherm fitting" (`apx:isotherm_fitting`, `eq:isotherm_fit_objective`) and a
    `FittingWeighting` row in `tab:json_parameters`, in `main.tex` and `main2.tex`; `main2.tex` also gets the column
    defaults 1, 2, 3 that `main.tex` already had. `BEA_fitted.json` and `CoBDP_fitted.json` regenerated with B-FIT-4
    (Absolute): BEA nC7 unchanged to 1e-9; CoBDP p-xylene 2×LF RMSE 0.064 → 0.051 mol/kg, but its second site has
    ν = 0.59, so the points below 50 Pa are now 70–213 % high (before: 80 % low). With `Relative` the same model is
    within 2 % there, 3 % in the step, and 8 % on the plateau; the example input still uses `Absolute`.
- [x] **F-PY-2** The self-driving lab as the subpackage `ruptura.selfdriving` (2026-11-11, user request): `psa`
      (`PressureSwingCycle`, `key_performance`, `run_experiment`: robot and analytics, generic in the product), `optimizer`
      (GP, ParEGO, EI x PoF, Kriging believer), `campaign` (`Campaign.run`, `write_inputs`, `rerun`, which replace the
      example's generate.py and run.py); worker `python -m ruptura.selfdriving <dir>`; extra `ruptura[selfdriving]`
      (scipy). The example is now `analysis.ipynb` with its lab notebook and figures (conversions.md went into the
      `psa` docstring). Tests: `tests/selfdriving_test.py`. Paper: main 2.5.5, 3.12, SI S5.5 (API and snippet), S6.10.
      A rerun from scratch through the package reproduced the 136 setups and results of the campaign.
- [ ] **F-PHY-12** Option to evaluate the isotherms (also those of chemisorption sites) at the adsorbent
      temperature T_s instead of T_f, and show the difference on the non-isothermal Marx runs. The SI justification
      of T_f (collision picture) was judged physically incorrect in the mock review of 2026-11-19 (Major 3); h_fs =
      20 W m⁻² K⁻¹ of the Marx runs also has no source (Wakao–Kaguei gives about 100–150 at T15).
  - 2026-11-19: T_s is now unconditional, without an option (B-PHY-15, user's request); the difference shows in the
    reruns of V-PHY-9. What remains here: an option back to T_f if wanted, and the source of h_fs.
  - 2026-11-19: the Marx validation now takes h_fs from Wakao–Kaguei–Funazkri (37.5–173.1 W m⁻² K⁻¹, V-PHY-10); the
    reference settings of the other T15 studies keep 20 W m⁻² K⁻¹ (no source; SI Table `tab:marx_reference_settings`).

## Paper inconsistencies (I)

- [x] **I-PHY-4** SI pore transport (Eqs. `surface_concentration_balance`, `pore_concentration_balance`,
      `reaction_pore_reservoir`) does not conserve mass as written.
  - β_s, β_p, and Γ are per moving-fluid volume, but F = a_loc k_f (c − c_s) and P = k_p (c_s − c_p) are
    per particle (or pore) volume.
  - The bulk loses (1−ε_b) a_loc k_f Δc per total volume, while the surface reservoir gains ε_b a_loc k_f Δc.
  - Corrected set, consistent with the code after the γ fix (C-PHY-1), which is the case ε_p = 1, β_s = β_p:
    F = ((1−ε_b)/ε_b) a_loc k_f (c − c_s) and P = β_p k_p (c_s − c_p).
    Details in `docs/STADLER_FEATURES_PLAN.md`, "Status after implementation".
  - 2026-11-18 (paper audit, W-GEN-11): the SI differs from the code in three more places:
    - the reservoir basis: the code's reservoirs are per solid volume with unit capacity, and there is no β_s, ε_p,
      or ε_wc input;
    - the monolith pore coefficient: 15 D_p/ℓ² for every geometry in the code (`integrators/transport.cpp:66`),
      π² D_p/(4δ²) in the SI (`main.tex:3119`, `:3474`);
    - the film area per solid volume: (4/d_h)/χ_wc in the code (`geometry.cpp:210`), 4/d_h in the SI.
  - 2026-11-20: done. SI "Pore transport" is rewritten on the code's basis:
    - reservoirs per adsorbent volume with unit capacity (no β_s, β_p, γ, or ε_p);
    - F = a_loc k_f (c − c_s), with a_loc = 6/d_p or a_ch/χ_wc, and P = k_p (c_s − c_p);
    - ∂c_p/∂t = P − ρ_s r + R^p;
    - k_p = 15 D_p/ℓ² for every geometry, with ℓ = d_p/2, δ_wc, or 3/a_loc;
    - a stated conservation argument against a^V_fs (Eq. `geometry_interfacial_area_total`);
    - new Eq. `film_pore_series`, 1/K = 1/k_f + a_loc/k_p.
    The main text (Eq. `main_pore_fluxes`) defines a_loc and k_p the same way. "Reaction without adsorption
    limitation" now points to these balances instead of repeating the β set. Ψ_p of Eq. `pore_domain_volume_fraction`
    has no ε_p (B-PHY-11). The washcoat porosity in the monolith inventory moved to I-PHY-18.
- [ ] **I-PHY-18** The SI keeps a washcoat porosity ε_wc that the code does not have. It appears in the geometry
      paragraph of "Breakthrough modeling", in Eq. `adsorbent_inventory_fraction` (Φ_a = ε_ch χ_a χ_wc (1 − ε_wc) for
      a monolith), and in the energy balance (`main.tex` ≈2572). The reader has only `WashcoatThickness` and
      `WashcoatVolumePerChannelVolume`, and the code's monolith inventory is ε_ch χ_wc ρ_p. Remove ε_wc from the SI,
      or add the input. Found 2026-11-20 while fixing I-PHY-4.
- [x] **I-PHY-14** `General` chemisorption capacity: the paper (Eq. `general_reversible_chemisorption`) has
      (q_max − q)^{n+}, the code has (min(q*, q_max) − q)^{n+} with q* the site isotherm at the bulk gas, also with pore
      transport, where the rate itself uses c_p (`chemisorption.cpp:15` and `:50`, `sorption.cpp:274`). Inputs that
      want the paper's law give the site a large affinity (`elfving/` and `validation/3.3.8_pore_transport` use b = 10³
      Pa⁻¹). Larger affinities make direct uptake stiff: b = 10⁶ took 92 s against 2 s at b = 10³ (3.3.8, 0.5–1.7 mm,
      100 cells, dense CVODE). Decide: the code (q_max, or q* at c_p) or the paper. Found 2026-11-09.
  - 2026-11-20: decided by the user: the paper follows the code. The main Eq. `main_chemisorption` and the SI Eq.
    `general_reversible_chemisorption` now have (q*_chem − q)^{n+}, with q*_chem = min(q*_site(c), q_max) and the site
    isotherm at the bulk fluid, also with pore transport. The SI states that the free capacity vanishes where the bulk
    holds no adsorbate, and gives the large-affinity recipe for (q_max − q) (`apx:goeppert`).
- [x] **I-PHY-15** SI statements that the code contradicts (paper audit 2026-11-18, W-GEN-11):
  - `main.tex:3538`: "Bulk-fluid reactions do not yet release heat, and the input reader rejects them together with
    the energy balance". The code releases this heat into the solid (`integrators/compute.cpp:344–347`), as SI
    `main.tex:2857–2866` says.
  - `main.tex:1603`: the SI says FastIAST stops on an absolute spread of 1e-13 and returns the current iterate after
    50 steps. The code tests the spread relative to max(1, |ψ̄|), and after the step limit it throws unless both
    residuals are below 1e-6 (`mixture_prediction.cpp:1184–1196`).
  - 2026-11-19: both SI passages now state what the code does (bulk-reaction heat into the solid,
    `eq:bulk_reaction_heat_source`, with a particle density needed; FastIAST relative spread and the 1e-6 acceptance).
- [ ] **I-ITG-3** ROS34PW2 in the paper predates its stage projection on the closure (B-ITG-16, 2026-11-11): the
      integrator text (main 2.4, SI `ROS34PW2` paragraph) does not describe it, the SI remarks of 3.4.2 give the old
      departure (−2.9 % to +4.0 %, now below 1e-4) and 1 to 12 s per cycle (now about half), and the benchmark
      quotes the old ROS34PW2 errors and steps. Update with the reruns of V-ITG-2.
  - 2026-11-19: done except the benchmark. Main 2.4 describes the projection and its effect on the TSA (4 % → < 1e-4);
    the SI `ROS34PW2` paragraph has it as Eq. `ros34pw2_closure_projection`, with why IPOP keeps the unprojected
    stages, and "what it leaves out" no longer lists the velocity closures. TSA: 0.6–6 s per cycle, SI remark
    < 1e-4 and 3e-6; the CO₂ balance is "within 4e-4" in the main text, the verification table, and the SI remark,
    which explains the doubling. Open: the T15 numbers, step table, figure, and batch date of 3.5 once the
    2026-11-19 benchmark batch has finished (V-ITG-2).
- [ ] **I-CYC-1** Introduction: the claims on cyclic operation go beyond Section 2.3 and the code (F-CYC-1).
      2026-11-02 (decision): the Introduction keeps them as the target; implement F-CYC-1 and F-BC-1 later. Section 2.3
      no longer lists what is missing.
      2026-11-06: the per-step flow direction is implemented and validated (F-BC-1, V-BC-3; Eq. `swing_phase_definition`
      and Section 3.11). Still beyond the code: cyclic-steady-state stopping, step ends on a process condition, and the
      per-step outlet pressure (F-CYC-1). Pressure-changing steps with closed ends are implemented since 2026-11-09
      (F-BC-2), not yet in Section 2.3.
      2026-11-11: Section 2.3 now has closed ends, pressure ramps, the re-equilibration of the gas, and the phase
      records; the self-driving lab (Section 3.12, SI `apx:self_driving_lab`) uses them, and its robot stops each
      run at the cyclic steady state in Python.
  - 2026-11-18: the code now has the stop on a condition, the cyclic-steady-state stop, and per-step ports (F-CYC-1,
    F-BC-3); Section 2.3 and the SI still describe the fixed-duration phase driver.
- [x] **I-CYC-2** `main2.tex` describes the cycle mechanics as a separate companion code (one to four columns,
      valves, buffers, compressor, jackets; SI `apx:cycles`). Since 2026-11-18 the single-column part is in RUPTURA
      itself, with differences from that text: one column until F-CYC-4; ROS34PW2 with a coupled matrix instead of a
      native BDF; WENO3 or upwind instead of TVD/WENO5; the isotherms at the gas temperature of the column model;
      the cyclic-steady-state measures of F-CYC-1; no KPIs or duties in the code (`phases.data` gives the streams).
      2026-11-18: several identical columns, valves, buffers, a compressor, jacket control, and the duties are in
      RUPTURA too (F-CYC-4), with each column on its own ROS34PW2 and the connections in explicit coupling steps
      instead of one state vector; the limits are in F-CYC-5.
      2026-11-19 (decision of the user): the multicolumn cycles become part of RUPTURA in the paper, which replaces D3 of
      `REWRITE_PLAN.md`. What the text has to hold, what the code verifies, and what `main2.tex` says that the code does
      not do is in `ruptura-2.0-paper/CYCLES_NOTES.md`; nothing is written into the paper yet.
      Decide whether the paper presents one code, and port the text into `main.tex` by hand.
  - 2026-11-19: done. The paper presents one code: Section 2.2 "Cyclic and multicolumn processes"
    (`sec:main_swing_adsorption`, `fig:multicolumn_schematic`), SI `apx:cycles`, and verification rows, written in
    `main_rewrite.tex`, which replaced `main.tex` the same day; `main2.tex` was removed (an older copy is in the paper
    git at `figures/main2.tex`). The literature validation is V-CYC-3.
  - 2026-11-18 (paper audit): `main.tex` already holds a paragraph on that companion code, Section 2.2.9
    (`main.tex:312–313`). It has no label, hard-codes "Section 3.10" and "Supplementary Section S7" (`main.tex` has
    only S1–S6), and uses Unicode quotes and dashes. Decision D3 of `ruptura-2.0-paper/REWRITE_PLAN.md`.
- [x] **I-CYC-3** SI text on the step driver that the code no longer matches (paper audit 2026-11-18, W-GEN-11):
  - `main.tex:3926–3927`: "Closed or one-way ports, v_b = 0, are not available in this version". `ClosedEnd` exists
    (F-BC-2), and Section 2.3 describes it.
  - `main.tex:2294–2305` and `:2725–2741`: the continuity closure and the fluid energy balance lack the
    pressure-ramp terms of Section 2.3 (c_tot (dP/dt)/P and the pressure work).
  - `main.tex:5665`: after an isothermal temperature step, the gas "keeps the total concentration of the previous
    temperature" (25 %). The driver now rescales it (`rescaleGas`, `swing_adsorption.cpp:190–195`), as Section 2.3
    says. The −2.9 % to +4.0 % of the same paragraph is I-ITG-3.
  - 2026-11-19: done in `main.tex`. The closed-port sentence was already gone with the rewrite. The SI closure
    (`eq:gas_total_molar_balance`) now has −φ_f c_tot (dP_b/dt)/P_b, and the gas energy balance the pressure work
    φ_f dP_b/dt, as in `compute.cpp` (sinkTerm, the T_f rate); both for `InletPressureInletVelocity`, the only
    continuity closure that takes a ramp. The TSA remark describes the rescaling (a fifth removed at a 373 K start, a
    quarter added at a 298 K start). Check with the current code (isothermal T15 TSA, t½ = 4000 s, 2026-11-19): it
    runs, the outlet backflow peaks at −0.60 m/s after the 373 → 298 K switch (the SI's 0.6 m/s), and c_tot RT/P
    stays within 1.5e-3 of 1. The departure −2.9 % to +4.0 % stays with I-ITG-3.
- [x] **I-BLD-1** SI Installation (`main.tex:4373`) describes `BUILD_APP` as the "Qt application", but it builds the
      Ruptura Lab web launcher (`CMakeLists.txt:330`). Found 2026-11-18 (paper audit). Fixed in the SI 2026-11-19.
- [x] **I-PHY-2** Keep the SI humidity table, q_sat(T)/ν(T), and the flow-reversal boundary conditions consistent
      with what gets implemented (F-PHY-1, F-PHY-2, F-BC-1).
      2026-10-30: q_sat(T)/ν(T) done (F-PHY-2). The Climeworks relations stay in the SI text, marked as parts of the
      humidity model rather than temperature laws; the humidity models and F-BC-1 remain.
      2026-11-04: the humidity models match the code (F-PHY-1). The SI has "Liquid and pH-dependent isotherms" and
      "Humidity-dependent isotherms" as subsections beside the temperature one; the Climeworks relations moved there;
      the GAB-coadsorption blocking is written with ω; the JSON tables list the new types and keys. F-BC-1 remains.
      2026-11-04: done with F-BC-1. SI: the right reservoir (feed of a reverse step, product-end gas or last cell),
      the mirrored Danckwerts face value and temperature, the reverse closures and pressure integration, the
      "Flow towards z = 0" paragraph of the method of lines, the Jacobian band, and the JSON keys; the one-way-port
      paragraph (not implemented) now says closed or one-way ports are not available. Main text: d_m in
      Eq. `swing_phase_definition`.
- [x] **I-PHY-12** The paper evaluates physisorption equilibrium at the solid temperature, q* = Q(ξ, T_s, RH, pH)
      (main.tex:119, and "isotherm parameters use T_s" in the SI temperature section), but the code passes the gas
      temperature to `predictMixture` (`rk3.cpp:140`, and the Jacobian `jacobian.cpp:174`). The two differ only with
      the energy balance on. Found 2026-11-04 while looking at F-PHY-1; not checked by a run.
      2026-11-04: RH uses the gas temperature (decided with F-PHY-1), and the SI now says so without the old "whereas
      isotherm parameters use T_s". The isotherm parameters elsewhere in the paper still read T_s; open whether the
      paper should say the gas temperature throughout (chemisorption rates do use T_s, `sorption.cpp:255`).
      2026-11-04, decided: gas temperature except for chemisorption. The paper now uses T_f for the physisorption
      equilibrium (Eqs. `main_equilibrium_interface`, `main_spi`, `physisorption_equilibrium`, SPI/EI/IAST, pH and
      humidity forms) and a generic T in the temperature laws, defined as T_f for physisorption isotherms and RH and
      T_s for chemisorption (rates and equilibrium isotherms) and reactions. Remaining code mismatch: I-PHY-13.
- [x] **I-PHY-13** The paper (since I-PHY-12) evaluates the equilibrium isotherms of chemisorption sites at the
      solid temperature T_s, but the code evaluates them at the gas temperature: `computeChemisorptionEquilibriumLoadings`
      gets `views.gasTemperature` (`rk3.cpp:150`), while the chemisorption rates (`sorption.cpp:255`) and reactions
      (`compute.cpp:303`) use T_s. Changing the code to T_s alters every run with the energy balance on and a
      temperature-dependent chemisorption isotherm: the 13 Elfving cases (`elfving/`, GAB and Langmuir chemisorption
      sites, and their copies in `validation/ros34pw2_danckwerts_rerun`) and the 8 energy-on cases of
      `benchmark/figures/simulations/methods/marx-403030_298_1e6_2e-1`. Stampi-Bombelli is isothermal. Decide: change
      the code (and rerun those), or write T_f for the chemisorption equilibrium in the paper. Found 2026-11-04.
      Decided 2026-11-04: all isotherms, chemisorption ones included, at the gas temperature; only the chemisorption
      rate coefficients (and the reactions) at T_s. The paper now says so; the code is unchanged, so no rerun.
- [x] **I-ITG-4** SI 3.4.4 (`main.tex:6170`): "the WENO3 runs hold 3e-4 more adsorbate than was fed ... this is
      integration error of ROS34PW2". It was the clamp of the accepted state, fixed 2026-11-19 (B-ITG-17).
  - Done 2026-11-19 (decision of the user): 3.4.4 rerun with the fix (all 20 runs, notebook and figures); the sentence
    is dropped, since the column now holds the fed amount to 2e-5 (WENO3, 25 cells; 3e-6 on 100) and 2e-8 (upwind).
    Table `flow_reversal_results`: eight cells changed, E_1 of WENO3 on 400 cells 0.0033/0.0038 → 0.0034/0.0033 and
    t99 counter-current down by 0.001; main text "0.3–0.4 % on 400 cells" → "0.3 %"; the other quoted values hold.
    The other runs with WENO3 made before the fix are in V-ITG-2.
- [x] **I-ITG-2** SIRK3 stability: the main text (Section 2.4, "stable at any Courant number and sorption stiffness")
      and the SI (after `eq:sirk3_stability_function`: |R(z)| ≤ 1 in the whole left half-plane) contradict the
      B-ITG-1 analysis of 2026-11-02, in which the final stage makes R singular at z = ±i and limits the convective
      step to Courant ≈ 1. Correct both, or leave SIRK3 out of the paper until B-ITG-1 is fixed.
  - 2026-11-02 (decision: keep SIRK3 in the SI, do not advertise it): the main text now says only that a
    semi-implicit SSP-RK3 variant exists for single-adsorbent columns. The SI paragraph already describes the coupled
    SIRK3 of F-ITG-3 without the fourth update, whereas `src/integrators/rk3_si.cpp` still solves the final stage
    (I + Δt²L²); they agree once F-ITG-3 lands.
  - 2026-11-03: resolved: SIRK3 is removed from the code and the paper; the main text and SI describe ROS34PW2.
- [x] **I-BC-1** Species inlet condition: the SI describes a Danckwerts (total-flux) inlet, but the code holds the
      inlet node at the feed state, i.e. a Dirichlet inlet (`src/integrators/compute.cpp:652`). Decide whether to
      implement Danckwerts or to describe Dirichlet.
  - 2026-11-02: node 0 holds the feed value on the inlet face, and the first cell gets the convective flux v c_in
    plus the dispersive flux D (c_in − c₁)/(h/2) across that face (WENO3 and upwind paths of
    `computeSpeciesDerivatives`). Danckwerts sets the total flux v c(0⁺) − D ∂c/∂z|₀₊ to v c_in, i.e. no extra
    dispersive inflow; the two agree only for D_ax = 0. The relative excess inflow is about 2(c_in − c₁)/(c_in Pe_cell),
    Pe_cell = vΔz/D: small for Marx (Pe_cell ≈ 80), large for the monolith with the fitted D_ax (0.04). Decision
    pending.
  - Done 2026-11-03 (decision: describe Dirichlet): main text (boundary conditions) and SI (species and
    fluid-temperature boundary conditions, with the comparison to Danckwerts and the D/(vL) estimate; the
    channel-basis and humidity sentences no longer say Danckwerts; a note in the benchmark settings). Whether to
    implement Danckwerts: B-BC-4.
- [x] **I-BC-4 Paper text that depends on the inlet condition** (B-BC-4; reruns in
      `validation/inlet_condition_dirichlet_vs_danckwerts`, 2026-11-03). Decided 2026-11-03: Danckwerts, the new
      default. The method text is done. The stored validation outputs were computed with the Dirichlet inlet, but
      their inputs omit `InletCondition` and now run Danckwerts: rerun marx2, stampi-bombelli, stadler, and 3.4.1–3.4.3
      in place, then update:
  - main text, Monoliths: "8–13 % early in the monolith" (5.6 %) is −11/−8/−13 % with Dirichlet, but +6/−2/−11 %
    with Danckwerts; the closing sentence on the fitted D_L could state that with Danckwerts the fitted coefficients
    reproduce the measured spreading (widths 60–133 %), and the Conclusions sentence on the monolith follows;
  - main text, Sensitivity: "only 1e-3 advances t50 by 2.9 %" is an inlet effect (−0.24 % with Danckwerts); the
    broadening (+46 %) stays;
  - SI inlet paragraph: done 2026-11-03 (Danckwerts default; the Dirichlet option with the D/(vL) estimate and the
    t50 changes measured in the reruns);
  - the benchmark monolith (fitted D_ax, vL/D ≈ 2) is a dispersion-dominated case either way; with Danckwerts the
    sentence about the Dirichlet inflow goes, and its CVODE/RK3 numbers would need the rerun (V-GEN-4 style);
  - figures: `fig:stampi-bombelli`, `fig:stampi_bombelli_400ppm`, and `fig:sensitivity` (c), from the existing
    notebooks after a switch.
  - 2026-11-03: every validation case except Elfving (B-EX-9) is rerun with ROS34PW2 at its defaults and the
    Danckwerts inlet, built at `-O3` (`validation/ros34pw2_danckwerts_rerun`, copies; the sources were not touched).
    All 223 finish, including the runs that stop under CVODE. Against the CVODE + Danckwerts copies, |Δt50|/t50
    ≤ 8.4e-4 (Stampi-Bombelli), ≤ 2.8e-4 (Stadler), ≤ 1e-5 elsewhere. The studies are switched in place:
    3.4.1 (both) done (I-BC-5); the others split to the other session (2026-11-03).
  - Done 2026-11-03 (Danckwerts session; the user decided that Stampi-Bombelli uses the fitted D_L and that the
    validation studies use ROS34PW2):
    - Generators and inputs: `stampi-bombelli/generate_runs.py` (no CVODE settings left), `stadler/generate.py`,
      `liu/generate.py` (CVODE tolerance parameters removed), marx2 `base.json` and T1–T24; 3.4.2 and 3.4.3 inherit
      from T15. All write `"BreakthroughIntegrator": "ROS34PW2"` and `"InletCondition": "Danckwerts"`.
    - Reruns in place with build/ruptura of 15:03: stampi-bombelli 20, marx2 24, liu 56, stadler 62, 3.4.3 38,
      3.4.2 7. All finish (`run_2026-11-03_ros34pw2.log` in each study).
    - Notebooks re-executed and paper figures copied: Stampi-Bombelli (fitted D_L as the base model; panels 4b–5c),
      marx2 `create.ipynb`, Stadler, 3.4.2, 3.4.3, and Liu.
    - Paper:
      - Monoliths paragraph, both captions, SI 400 ppm section, benchmark settings caption, and Conclusions use the
        fitted D_L;
      - Sensitivity: D_ax 1e-3 now −0.3 % in t50; numerical-dispersion check 3.9 %;
      - Stadler: plug-flow difference 0.0003–0.0019, grid error 2.2e-8, SE-WGS conversion 0.760;
      - Cycles: 14 cycles, 1.1e-4 mol/kg/s, 1–12 s per cycle;
      - SI TSA gas inventory: −4 % to +25 %, at most 3.4 %;
      - Marx settings table: ROS34PW2; tab:marx_t50_comparison updated (12 rows change by 1 s or 0.1 K).
    - The 3.4.2 mass balance (`analysis.ipynb`) adds the Dirichlet dispersive inflow only for Dirichlet inputs.
- [x] **I-BC-5 Boundary-condition results text (3.4.1) describes CVODE stalls that ROS34PW2 does not have.** The main
      text says that the fixed-velocity runs "stop after 20–44 s of simulated time, when the CVODE step size
      collapses", and that the imposed-pressure runs of T3 and T9 were stopped after 900 s. In
      `validation/ros34pw2_danckwerts_rerun` (2026-11-03), all six run to 4000 s with ROS34PW2. Up to the CVODE stop
      they agree with CVODE to 1e-13. The total gas concentration at the front still drops to 2e-7 (T15) and 6e-8 (T3)
      of the feed after about 90 s, the defect the text describes. CO₂ t50 equals that of the inlet-velocity settings
      (T15 636 s, T3 904 s, at the 4 s output interval). Decide which integrator the study uses, and update the text
      (and `fig:boundary_conditions` if those runs are drawn).
  - 2026-11-03: done (decision: ROS34PW2 with the Danckwerts inlet). Both 3.4.1 generators write
    `BreakthroughIntegrator` ROS34PW2 and `InletCondition` Danckwerts and drop the CVODE keys; all 25 runs rerun in
    place (`build/ruptura-o3-20261103`) and complete; `analysis.ipynb` draws the "Not completed" legend entry only
    when a run stops; the four panels are regenerated. CO₂ t50 against BC1: BC4/BC5 −0.05 to −0.45 %; BC2, BC3, the
    2.7/2.1 v_in, and the BEA numbers are unchanged. The front concentration falls to 6e-8–3e-7 of the feed between
    52 and 100 s; ROS34PW2 and CVODE agree to 7e-5 of the feed concentration up to the CVODE stop (not 1e-13). Main
    text: the stall sentences replaced by the near-vacuum and the t50 shifts; caption without the crosses. SI: both
    setup paragraphs name ROS34PW2; "Runs that did not complete" is now "Settings without continuity of the molar
    flux", with the CVODE stalls for comparison.
- [x] **I-PHY-16** Since B-PHY-15 (2026-11-19) the code evaluates every isotherm and the RH of the humidity isotherms at
      T_s; the paper still says T_f (I-PHY-12, I-PHY-13), with the collision-picture justification of
      `main.tex:467` that the Reviewer-2 audit calls physically wrong. Rewrite with T_s (the van 't Hoff factor comes
      from desorption at the surface): `main.tex:467–473`, Eqs. `main_equilibrium_interface` and `main_spi`, SI Eq.
      `physisorption_equilibrium`, the temperature-law section (`main.tex` ≈3930), the RH definition, and the
      chemisorption equilibrium. Numbers after V-PHY-9.
  - 2026-11-19: the caption of SI Table `tab:temperature_laws` now says T_s for isotherms (edited with V-PHY-10);
    the rest is open.
  - 2026-11-20 (user's request): Eqs. `main_equilibrium_interface`, `main_spi`, `physisorption_equilibrium`, the RH
    definitions (SI ≈3800 and ≈4810: p_sat at T_s), the moisture-swing affinity, the temperature-law section, and the
    SI sentence "The isotherms are evaluated at ..." now say T_s; the collision picture is replaced by the desorption
    argument (the exponential factor of the affinity comes from desorption of molecules at the temperature of the
    solid).
- [ ] **I-PHY-17** SI S5.6 (pore reservoirs) writes the balances with ε_p and β_p = (1 − ε_b)ε_p/ε_b, but
      the code has no particle-porosity input ("porosity" appears nowhere in src/, 2026-11-19). Check the documented
      equations against transport.cpp/sorption.cpp, and the mock-review claim (Major 5d) that the film flux per
      particle volume and reservoirs per fluid volume do not conserve mass in a fixed bed.
- [ ] **I-ITG-5** Benchmark timings from different batches disagree: BEA RK3 22 s (main 3.4) vs 27.2 s
      (Table S24); BEA CVODE about 16 s (Fig. 10a) vs 19.2 s (Table S24). The component/method study (Table S24,
      3 Nov) and the convection cost runs (Table S25, 31 Oct) ran at CVODE rtol 1e-4, the default before 1f4da97
      (11 Nov); the SI now says so. Regenerate all benchmark numbers from one code version, with repeats.

## Performance (P)

- [ ] **P-ITG-9** The automatic CVODE step cap, 4 Δz/v at the gas velocity (`resolveMaximumTimeStep`, `cvode.cpp`), is
      far too small for chemisorption fronts that move 1e5 times slower than the gas. On the Goeppert case
      (`validation/3.3.8_pore_transport`, 50 cells) it holds CVODE at 0.035 s: 57 912 Newton iterations for the first
      2000 s of a 158 000 s breakthrough. With `CVODEMaximumTimeStep` 0 and the `SorptionJacobian` preconditioner the
      whole breakthrough takes 393 steps (0.9 s, 100 cells, default tolerances), and no outlet curve of that study
      shows the hook the cap guards against. A cap from the front speed (or the sorption time scale) would fit both.
      Found 2026-11-09.
- [ ] **P-ITG-10** The `SorptionJacobian` CVODE preconditioner (c and q blocks) leaves out the film and pore
      reservoirs of `UsePoreSurfaceTransport`. SPGMR then stalls when that exchange is fast against the step: in
      `validation/3.3.8_pore_transport`, 0.25–0.5 mm with the Table 5 dispersion (exchange about 10³ s⁻¹) stays in the
      first 900 s output interval for over 5 min at 400 cells (also in 130 s probes with Krylov dimension 100 or a
      60 s step cap), and so does its 99.9% film split (pore exchange 5e5 s⁻¹). The dense solver (difference-quotient
      Jacobian, reservoirs included) runs those cases at 100 cells in 3–20 s, so the study uses it up to 200 cells.
      Probable cause, not proven: the missing reservoir blocks. ROS34PW2 builds its W matrix from the same Jacobian
      (not tried on this study). Found 2026-11-09.
- [x] **P-PHY-4** Hot-loop allocations: `computeReactionDerivatives` (`compute.cpp:224`) `activity`/`maximumLoadings`
      per RHS; `reaction.cpp:39-69` `theta` per grid node; Bingel–Walton `psi` (`isotherm.h:369`) two 100-entry
      vectors per call inside the IAST bisection; `computeSegregatedCompetitiveIsotherm` three vectors per grid
      point; `meanParticleNumbers` three per call; `transport.cpp:120` builds an `mdspan` in the innermost loop.
      2026-10-28: the `transport.cpp` `mdspan` is now built once per call (single/multibed merge); the others remain.
      In a BEA C7 RK3 profile `malloc`/`free` are ~0.3% of the samples, so the gain is in the runs that allocate per
      grid node or per ψ evaluation (reactions, Bingel–Walton, segregated IAST, MPD; not measured).
      2026-11-03: all six sites are still there. Also per RHS: the five WENO3 buffers of `computeMassDerivatives`
      (`values`, `faces`, `faceFlux`, `totalConcentration`, `faceFractionSum`), the Danckwerts copies of the
      concentrations and the gas temperature, and the five zero-fills of `computeReactionDerivatives`, which run
      before its `reactions.empty()` return. RK3 `sample` profiles (`build/ruptura` of 15:03), share of the samples:
      BEA C7: WENO3 buffers 0.7 %, reaction fills without reactions 0.85 %; BEA with nC7 only: WENO3 buffers 3.3 %;
      `Reactions/pore-langmuir-hinshelwood` (RK3, 200 s): `theta` in `Reaction::rate` 7 %, WENO3 buffers 1.7 %.
  - 2026-11-03: `theta` and the WENO3 buffers are column-owned (`ColumnData::reactionCoverage`, passed through
    `computeReactionDerivatives` to `Reaction::rate(activity, T, coverage)`; `ColumnData::wenoScratch`, the trailing
    span of `computeMassDerivatives`). Identical outputs; user time, best of three: LH reactions RK3 200 s 18.08 →
    16.81 s (−7 %), BEA nC7 only RK3 100 s 5.11 → 4.92 s (−3.7 %), BEA C7 RK3 60 s 22.30 → 22.21 s (noise).
    Remaining: `activity`/`maximumLoadings`, Bingel–Walton `psi` (two `std::array` rows on the stack would do),
    segregated IAST, MPD, the Danckwerts copies, and the reaction zero-fills.
  - 2026-11-06 (code review), more per-call allocations: `chemisorptionSitesPerComponent` builds a vector in
    `computeChemisorption`, `computeChemisorptionTransportDerivatives`, and `computeBulkSpeciesSink` (four calls per
    RK3 stage); `updateVelocityAndPressure` assigns three vectors per call with the energy balance; `computeMassDerivatives`
    builds `outletFace` on every call and `reverseFaces` with back-flow; `RungeKutta3::step` copies the state into a
    new vector every step. All could be column-owned like `wenoScratch`.
  - 2026-11-09 (user: is preallocating worth the effort?): `sample` profiles of the scratch build (RelWithDebInfo),
    8 s per case; share of the samples in malloc/free, and in memset/bzero in brackets: BEA C7 RK3 0.1 % (1.2 %);
    MOR SIAST CVODE 0.07 % (2.6 %); Bingel–Walton IAST mixture prediction 0.1 % (0); Marx T20 energy balance
    ROS34PW2 0.35 % (2.1 %); LH reactions CVODE 1.2 % (3.8 %); monolith 400 ppm Danckwerts RK3 1.6 % (5.6 %);
    Chemisorption-CO2 CVODE 1.6 % (3.5 %); MPD BEA RK3 2.9 % (0.4 %). The allocations are `meanParticleNumbers` and
    `computeMPD` (MPD), `chemisorptionSitesPerComponent` (chemisorption and reactions), and `computeMassDerivatives`
    (Danckwerts copies, outlet faces). Most zero-fills are not allocations: `predictMixture` clears the preallocated
    `equilibriumSiteLoadings` on every call (about 3 % in the monolith), then `computePhysisorptionEquilibriumLoadings`
    and the five fills of `computeReactionDerivatives` without reactions (0.9 %). The Bingel–Walton `psi`, the
    segregated IAST vectors, `activity`/`maximumLoadings`, the energy-balance vectors, and the RK3 state copy are each
    below 0.5 %. Worth doing: the MPD buffers (about 3 % of MPD runs), `chemisorptionSitesPerComponent` once per
    column (0.5–1 %), and the `reactions.empty()` return before the fills (about 1 %); the rest is not.
  - Closed 2026-11-09 without a change (user), after the measurement above.
- [x] **P-PHY-7 Node 0 is integrated as a sorption node.** Node 0 is the inlet face (WENO3) or the feed reservoir,
      yet every right-hand side evaluates its mixture prediction, LDF, chemisorption, pore transport, and reactions,
      and integrates its loadings and its solid and wall temperatures; in a reverse step it is unused altogether.
      Nothing reads these values except the output, where z = 0 shows a loading that equilibrates with the feed at
      once, and the wall conduction of node 1, which treats T_w,0 as a wall cell. Cost: 1/(N + 1) of the sorption
      work (1 % at 100 cells, 5 % at 20), plus dead rows in the Jacobians. Skipping node 0 changes the outputs at z = 0
      only; the wall boundary at z = 0 then needs its own condition (zero flux). Found 2026-11-06 (code review).
  - Closed 2026-11-09 without a change (user).
- [x] **P-ITG-8** ROS34PW2 repeats a precompute and a right-hand side at the step start after every rejected step
      (`ros34pw2.cpp`, "A rejected step is repeated from the step start"). Restoring `state` and `stateDot`, which
      the code already does, gives the retried stages the same inputs as the first attempt, since every stage
      `evaluate` runs its own precompute. With `FixedVelocity`, `FixedPressureInletVelocity`, or a liquid the
      recomputed rates are bit-identical, so the work is pure overhead: one mixture prediction over the column and one
      right-hand side per rejection. With a continuity closure the recomputed rates differ slightly, because precompute
      builds the Ergun pressure that the equilibrium sees from the sink of the previous evaluation (the right-hand side
      is not a pure function of the state, cf. B-ITG-7); skipping the recompute keeps the first attempt's rates.
      Measure the rejection share first (`printStatistics`). Found 2026-11-06 (code review).
  - Closed 2026-11-09 without a change (user).
- [x] **P-ITG-11** The right-hand side allocated its work arrays on every call (the concentration copy and two more
      vectors in `computeMassDerivatives`, three vectors in the velocity closures with the energy balance, the
      chemisorption site counts three times, three vectors per site and node in SCI). Found in the 2026-11-18 audit
      of `src/`; not measured.
  - Done 2026-11-18: they come from the column (`wenoScratch`, `temperatureRateScratch`, `chemisorptionSites`) and
    from the work arrays of `MixturePrediction`; a caller without them gets them allocated per call, as before.
    Measured for `computeMassDerivatives`: no difference in CPU time (BEA, RK3, 2.03 s with and without).
- [x] **P-PHY-8** `predictMixture` evaluated every isotherm site once more after each IAST solve to fill
      `equilibriumSiteLoadings`, which in a column only the chemisorption mixture reads. Found in the 2026-11-18
      audit of `src/`; not measured.
  - Done 2026-11-18: `MixturePrediction::storeSiteLoadings` (true by default); the column sets it false for its
    physisorption mixtures.
- [x] **P-PHY-9** MPD reweighting (`MacrostateParticleDistribution::meanParticleNumbers`) called exp() for every
      compatible state (online log-sum-exp) and kept a serial running maximum: 45 µs per call with the 9865-state
      CO₂/propane distribution of V-PHY-7, ~97 % of an MPD breakthrough (1e8 calls, about an hour against 69 s for
      SIAST).
  - Done 2026-11-19 (asked by the user): two passes, the log weights over one contiguous array per component and
    their maximum (four running maxima), then exp() only within 50 of it (each skipped state adds < e^-50 to a sum
    of at least 1); a component absent from the gas reweights only the precomputed states without it, and a gas
    without adsorbing components returns at once. 14 µs per call over the 104k compositions of the V-PHY-7
    breakthroughs (3.2×); loadings change by ≤ 1.4e-13 molecules (≤ 1e-10 relative). Test
    `ReweightsOnlyStatesWithoutAbsentComponents`. A threshold that drops a component with a tiny but nonzero mole
    fraction was rejected: 0.1–0.4 % of the calls there, and a jump in the loading.
- [x] **P-PHY-10** MPD had no equilibrium sensitivities: `assembleSorptionJacobian` left dq*/dc at zero for MPD, so
      ROS34PW2 (and the CVODE sorption-Jacobian preconditioner) lost the stiff coupling and crawled. V-PHY-7
      breakthrough with ROS34PW2: 1435 s of 3850 s in 87 min, against 7.4 s for the whole SIAST run.
  - Done 2026-11-19 (asked by the user): `predictMixtureSensitivities` returns dq_i/dp_j = Cov(n_i, n_j)/p_j for MPD
    (fluctuation relation; `meanParticleNumbers` fills the covariance from the second moments of its second pass),
    with p_j ≥ 1e-10 P like the other slopes; `assembleSorptionJacobian` no longer excludes MPD. Test
    `MPDSensitivitiesMatchCentralDifferences` (< 1e-6 against central differences). Only checked isothermal; the
    temperature terms of the Jacobian were not looked at for MPD.

## Refactoring and cleanup (R)

- [ ] **R-IO-1 Input-key conventions for 2.0.** Found 2026-11-06 (code review):
  - The energy-balance keys are camelCase (`energyBalance`, `wallDensity`, `heatCapacityGas`, `heatTransferGasWall`,
    …) among PascalCase keys. The reader does not care, but the docs, presets, and GUIs carry both styles.
  - Keys are case-insensitive, except the chemisorption `RateEquationParameters`, which go through
    `requireOnlyExactKeys`. Keys that differ only in case (`Temperature` and `temperature`) are not rejected;
    `findKeyCaseInsensitive` takes the first in `std::map` order. Pick one rule (exact keys, or case-insensitive
    with duplicates rejected).
  - The void fraction, particle diameter, and particle density can be given at the top level, in `Geometry`, and per
    `Adsorbents` entry, each overriding the previous one. With a `Monolith` geometry the per-adsorbent
    `ColumnVoidFraction` and `ParticleDiameter` are accepted and printed in the header (`Adsorbent void-fractions`)
    but not used (the adsorbent takes the monolith geometry).
  - `SimulationType` `Test` is accepted and runs a breakthrough (`main.cpp`, `default:`). 2026-11-18: removed
    (R-GEN-5).
  - `InternalChannelDimension` means the diameter for circular channels and the side length for polygons (AI prompt
    doc only); the Lab and Qt help say "the channel size". A hexagonal channel quoted by its flat-to-flat width is
    then √3 too wide. Name the quantity in the help, or take the hydraulic diameter for every shape.
  - 2026-11-09 (user: all keys case-insensitive): the `RateEquationParameters` keys are case-insensitive
    (`requireOnlyExactKeys` removed, `requireExactKeys` ignores case), and two keys of one object that differ only in
    case are rejected (`requireOnlyKnownKeys`); checked in `ChemisorptionInput.RequiresExactlyTheParametersForTheSelectedModel`.
    The camelCase energy-balance keys need no rename then. The other points are open.
- [ ] **R-BLD-1** Left from the 2026-11-17 audit, none urgent:
  - [x] `gtest_discover_tests` runs the tests in the source root (`tests/CMakeLists.txt`), where they leave
    `column.data`, `component_*.data`, `column.json` and `cmake_test_discovery_*.json`; the build directory works.
    Fixed 2026-11-18: they run in the build directory of the tests. `SkarstromCycle.MassBalance...`, which reads its
    output files back, runs in its own directory, so that `ctest -j` no longer fails it.
  - [x] The warning flags also cover nanobind, Python and SUNDIALS headers: about 415 `-Wold-style-cast` warnings around
    the 6 of the project (`pybindings.cpp:801`, `ros34pw2.cpp`).
    Fixed 2026-11-18: `-Wno-old-style-cast` for the Python module (`NB_SUPPRESS_WARNINGS` needs nanobind > 2.2), SUNDIALS
    and GoogleTest fetched as `SYSTEM`, `--system-header-prefix=gtest/` for a GoogleTest in `/usr/local` under Clang,
    and the project's own warnings fixed (`pybindings.cpp`, `ros34pw2.cpp`, a shadowed `hasChemisorption` in
    `inputreader.cpp`). A full build warns only from the sources of SUNDIALS itself (deprecations, its own flags).
  - [ ] `-O3 -ffast-math` are target options, so a `Debug` build is optimized as well. (User 2026-11-18: fine for now.)
  - [x] `BUILD_BENCHMARKS` (`CMakePresets.json`, `pyproject.toml`) and `BUILD_MAC_PACKAGE` (`CMakeLists.txt`) are read
    or set by nothing else; the AI prompt doc lacks `CVODEPreconditioner`, `CoverageReferenceLoading`,
    `DesorptionCoverageEnergy` and `DesorptionMassTransferCoefficient`.
    Fixed 2026-11-18: both settings removed (with the codesign step under `BUILD_MAC_PACKAGE`); the four keys are in
    the doc, and `CVODELinearSolver` lists `KLU`.
  - Not reviewed in that audit: the integrators, the physics models, most of the reader, `convert.py`. Unconfirmed
    leads: the dispersion stencil on a non-uniform grid (`compute.cpp`, `secondDerivative`), the absolute
    `wenoEpsilon`, and the gas sink of a node that mixes adsorbents with and without direct chemisorption
    (`transport.cpp`, `computeBulkSpeciesSink`).
- [ ] **R-EX-2** Every case behind the paper lives in `validation/` (user, 2026-11-19). Moved, each with a symlink at
      its old path:
      - `marx2`, `stampi-bombelli`, `elfving`, `liu`, and `stadler` into `validation/3.3.x_*/<name>`;
      - `benchmark/{breakthrough,convection,figures}` into `validation/3.5_numerical_performance/`;
      - `examples/Self-driving-lab-H2-PSA` into `validation/self_driving_lab/`.

      Changed with the moves: the path depths in the moved scripts and notebooks, the generators and notebooks that
      read the Marx, Stampi-Bombelli, Stadler, Elfving, and Liu inputs, `benchmark/cvode-max-step`, three C++ tests
      (`pressure_swing`, `ros34pw2`, `sorption_jacobian`), and the root `run.py`. `run.py` now also does what
      `validation/collect_jobs.py` did: `--paper`, `--generate`, and `--list` with `validation/cases.csv`. The
      section READMEs and `validation/README.md` were updated too. Open:
  - `sharma` (MPD, `validation/3.2.2_mpd/sharma`): left in place because another session was working in it.
    `sharma/runs/{generate,run}.py` find the repository with an unresolved `HERE.parent.parent`, and
    `validation/3.2.2_mpd/analysis.ipynb` reads `REPOSITORY / "sharma"`.
  - Git: ten tracked files now sit in the ignored `validation/` and show as deleted:
    - the Marx, Stampi-Bombelli, Elfving, and Liu PDFs;
    - `marx2/ternary/T15/simulation.json`, which the C++ tests read;
    - the three tracked files of the self-driving-lab example.

    The root symlinks show as untracked. Decide: `git add -f` the new paths and ignore the symlinks, or move the
    tracked files back.
    - Decided 2026-11-20 (user): `validation/` is tracked, without RUPTURA output. `.gitignore` names that output
      for the whole repository (the list of `pack.py`) in place of the rule for every `.data` file, so input `.data`
      files can be tracked. Left out of `validation/`: the notebook that ruptura writes into each run directory,
      the archives of `convert.py`, `_to_delete/`, the RASPA3 TMMC data, and the experiment directories of the
      campaigns. The ten files are tracked again at their new paths; the root symlinks are tracked since `ba573b9`,
      and their old ignore rules, which hid the case folders of the same name, are gone.
  - The C++ tests with the new paths pass a syntax check (`-fsyntax-only`) but were not built or run.
  - `validation/pack.py` (2026-11-19) packs the working tree for a rerun elsewhere without any RUPTURA output: 146 MB.
- [x] **R-ITG-1** One right-hand-side function. The precompute-plus-derivative block is written out five times
      (`rk3.cpp` ×2, `cvode.cpp` ×2, `rk3_si.cpp`, `swing_adsorption.cpp`).
      2026-10-27: RK3, SIRK3, the swing phase switch, and the Python bindings share `precompute` +
      `computeDerivatives`; CVODE keeps its own passes because it evaluates trial states in SUNDIALS vectors.
      2026-11-03: `rk3_si.cpp` is removed; ROS34PW2 shares `precompute` + `computeDerivatives` as well.
      2026-11-03: only `CVODE::evaluateDerivatives` is left. Its passes are, in order, precompute, computeDerivatives
      without the velocity, mass, and energy parts, precompute, computeDerivatives (one pass when the sink does not
      set the velocity), but on views of the trial N_Vectors rather than the column's spans. Every RHS change still
      has to be made in `rk3.cpp` and `cvode.cpp` (e.g. the dispersion arguments of `updateVelocityAndPressure`).
  - 2026-11-03: done. `precompute` and `computeDerivatives` take the state as `ColumnStateViews` (the 16 state and
    derivative spans, built by `ColumnStateLayout::views`) next to `ColumnData`, which holds no views, so a state
    field cannot come from the column by mistake; `Column` inherits the views of its own state, and the old
    `(Column&, ...)` signatures forward to the new ones. `computeDerivatives(..., sorptionOnly)` stops after the bulk
    sink. `CVODE::evaluateDerivatives` is now precompute, [computeDerivatives(sorptionOnly), precompute when the sink
    sets the velocity], computeDerivatives on views of the N_Vectors (the column's temperatures without the energy
    balance, as before); `UserData` carries the `Timing` of the running advance. Outputs and CVODE/ROS34PW2
    statistics are byte-identical to the previous binary on 11 inputs (CVODE: Chemisorption-CO2, LH reactions,
    liquid multibed, layered C6, paper T15, paper monolith 2 mmol; RK3: power-law reactions, liquid pH, CO2-N2;
    ROS34PW2: BEA, marx2 T15).
- [x] **R-PY-2** The in-memory `Breakthrough.compute()` (`pybindings.cpp`, `computeBreakthrough`) returns its own
      5 + 5 n columns (velocity, pressure, q, q*, p, p/p_feed, dq/dt) without the concentrations, mole fractions,
      temperatures, or chemisorption of the `.data` files, and with a different τ (B-IO-9). One schema for both would
      let the plotting code take either. The result buffer is a raw `new` that is freed on Ctrl-C only; an exception
      from `advanceToNextOutput` (e.g. a CVODE failure) leaks it. Found 2026-11-06 (code review).
  - Fixed 2026-11-09: the rows of `Breakthrough.compute()` hold the columns of the output files, from the functions
    that `writeOutput` uses too (`Column::appendColumnOutput`, `appendComponentOutput`; the files stay byte-identical):
    those of `column.data`, then for each component those of its file after τ, t, and z, the surface and pore
    concentrations of its transport sites included (`Column.surface_pore_transport_sites`).
    `ruptura.breakthrough_columns(simulation)` names them (`BREAKTHROUGH_BASE_COLUMNS`,
    `BREAKTHROUGH_LIQUID_BASE_COLUMNS`, `BREAKTHROUGH_COMPONENT_COLUMNS`; `equilibrium_adsorption` is now
    `equilibrium_physisorption`), and τ uses the factor of the files. The buffer is a `unique_ptr` until the array
    owns it. Checked on `Chemisorption-CO2/breakthrough` (CVODE, 207 snapshots, two pore sites): every value equals the
    file value.
- [x] **R-GEN-4** Small cleanups (code review, 2026-11-06): `Column::resolution` = L/N is printed as "Column spacing"
      and is CVODE's fallback spacing, both wrong on explicit or multibed grids; `cvode.h` includes the SPBCGS, SPFGMR,
      SPTFQMR, and dense headers it does not use, and documents `preconditionerScratch` as unused although the
      preconditioner setup writes to it; `fitting.cpp` keeps `<experimental/filesystem>` and `sys/stat.h` fallbacks
      under C++23; `main.cpp` prints the MixturePrediction summary before and after the run, and writes error
      messages without a trailing newline (some messages carry their own `\n`, most do not); `computeReactionDerivatives`
      checks the reaction site and phase prerequisites in every right-hand side, which belongs in the reader.
  - Fixed 2026-11-09: `Column::resolution` is removed (also from the Python bindings). The header prints the spacing,
    or its range on a non-uniform grid (CO2-C3H8-gradient: "0.004 to 0.013 [m]"); CVODE's smallest spacing needs no
    fallback (N ≥ 1). `cvode.h` drops the SPBCGS, SPFGMR, SPTFQMR, and `sundials_dense.h` includes, and the
    `preconditionerScratch` comment says what it holds. The unused filesystem fallback block is removed from
    `fitting.cpp` and `mixture_prediction.cpp` (neither uses the filesystem). MixturePrediction runs through
    `runSimulation` and prints its summary once. Errors end in exactly one newline. The reader rejects a
    `PoreConcentration` reaction whose `Site` is beyond the pore reservoir (test
    `ReactionsPore.RejectsASiteBeyondThePoreReservoir`); the right-hand side no longer checks sites and states.

- [x] **R-EX-1** Examples on git hold the inputs only (user, 2026-11-18). Untracked: 13 `breakthrough.ipynb`, 14
      `mixture.ipynb` and 13 `fitting.ipynb` (copies of an old template that imports `fitting_plotly_pubready`), the
      two `*_fitted.json`, and `BEA-alkanes-C7/psa/breakthrough.pdf`; the local files stay. Staged: every example
      `simulation.json` (also `ros34pw2`, `Multibed`, `Liquid`, `Chemisorption-CO2*`, `MPD`), the four
      `particle_distribution.data` they read, the MPD and Multibed READMEs and generator scripts, the converted
      `psa` input, and the removal of `sirk3`. From `github.com/iRASPA/RUPTURA` (96ed683): the 81 fitting data files
      that 11 fitting examples lacked; all 13 fitting examples run. `.gitignore` keeps out what ruptura and
      `convert.py` write next to an input (`analysis.ipynb`, `column.json`, `*fitted.json`, `_legacy*`, ...). Left out:
      `MPD/BEA-alkanes/presentation/` and `benchmark_runs/`, `change.py`, `manifest.csv`, the self-driving figures and
      `campaign/`, and `MPD/BEA-alkanes/breakthrough_benchmark.ipynb` (its saved output holds a `/Users/...` path; its
      README points to it). Committed 2026-11-18, with the self-driving notebook without its saved outputs and
      without `psa/breakthrough.pdf`; `sharma/` and `bermejo-lopez/` are local (`.gitignore`).
- [x] **R-IO-2** Default `Temperature` (and `InfluxTemperature`) 298 K instead of 433 K (user, 2026-11-18), in the reader,
      `ruptura/utils.py`, the Lab form and summary, and the AI prompt doc. No RUPTURA input in the repo or the study
      folders omits `Temperature` (the only ones are 1688 RASPA inputs in `sharma/tmmc`), so no result changes.
- [x] **R-GUI-2** Remove the Qt app (user, 2026-11-18): `app/main.cpp`, `mainwindow.*`, `runner.*`,
      `simulationinput.*`, `tests/lab_input_tests.cpp`, and `BUILD_QT_APP` are gone. `app/webmain.cpp` (the
      `ruptura_lab` launcher) and the Lab's import of saved Qt workspaces (`migrate_qt`) stay.
- [x] **R-PHY-2** The IAST solvers existed twice (IAST and per-site SIAST; Fast, bisection, explicit), and had
      drifted. Found in the 2026-11-18 audit of `src/`.
  - Done 2026-11-18: the per-site functions take the phase (components, count, cache slot); `computeFastIAST`,
    `computeIASTNestedLoopBisection`, and `computeExplicitIsotherm` zero the output and call them
    (`mixture_prediction.cpp` is 338 lines shorter). Where the copies differed, the IAST form is kept (cache used
    when `cachedGrandPotential > 0`), and the bisection prints its state before both of its errors. Results change
    at rounding only: BEA with RK3 2.6e-12 after 6000 steps; 12 of 17 example runs are identical, the CVODE runs
    differ within their tolerance (2e-5 in the states of the BEA breakthrough, the same 1793 steps). The per-phase
    Fast solver is inlined into its two callers (`[[gnu::always_inline]]`): as one shared function it made IAST 4%
    slower. With it the BEA example with RK3 still takes about 2% more CPU time than before the merge (1.96 s to
    2.00 s for 10000 steps); a separate copy of Fast IAST would take that back.
- [x] **R-IO-3** `Breakthrough::run` and `SwingAdsorption::run` each opened the files and wrote the snapshots and
      progress lines. Found in the 2026-11-18 audit of `src/`.
  - Done 2026-11-18: `Breakthrough::writeOutput`, which starts the files at its first call (the streams are members
    of `Breakthrough`) and holds the checks of B-IO-17; both loops call it.
- [x] **R-GEN-5** Code without a user. Found in the 2026-11-18 audit of `src/`.
  - Done 2026-11-18: removed `SimulationType` `Test` (it ran a breakthrough without its checks) and the ignored key
    `DebugForceMultibed`, also from the bindings and the Lab key list; `Timing::clear`, the pair `operator+`,
    `RandomNumber::Gaussian` and `UInt64`, `Chemisorption::usesEquilibriumLoading`, `Fitting::printSolution`
    (never defined) and `Fitting::externalTemperature` (Python `Fitting.external_temperature` is gone). The loop of
    `FixedStepIntegrator` is now `RungeKutta3::advance`; `RupturaNotebooks` takes `InputReader::SimulationType`.
    The NaN test by `union` in `Fitting::fitness` stays: it is needed with fast math (user).

## Validation and reruns (V)

Details in `validation/`.

- [x] **V-PHY-1** Marx T8/T21: marx2 had to use `InletVelocityOutletPressure` and stop early (1872 s and 1174 s). An
      older `marx2/ternary/T21/failed_state.json` shows a diverged state. With the current code, T21
      completes 4000 s with both `InletVelocityOutletPressure` and `InletPressureInletVelocity` (3.4.1).
      Rerun marx2 T21 (and T8) with `InletPressureInletVelocity` to 4000 s; update the SI run table and figures.
  - 2026-11-19: done with V-PHY-10; every run uses `InletPressureInletVelocity` to its full time.
- [x] **V-PHY-2** Find a validation case for pore transport (film and pore reservoirs). Elfving does not use it, and
      Serna-Guerrero is dropped. Candidates:
  - Rosen, J. Chem. Phys. 20 (1952) 387: analytical breakthrough with film and particle diffusion (3.1);
  - CADET general rate model (Leweke and von Lieres, Comput. Chem. Eng. 113 (2018) 274): code-to-code;
  - Kalyanaraman et al., Chem. Eng. J. 259 (2015) 737: CO₂ on PEI/silica hollow fibers, experiments and
    model;
  - Bollini et al., Ind. Eng. Chem. Res. 51 (2012) 15153: breakthrough on zeolite 13X and amine-SBA-15.
  - 2026-11-02 literature search (Stadler does not test pore transport: apparent FirstOrder sorption, pores ≈ bulk):
    - Recommended first: Goeppert et al., ChemSusChem 7 (2014) 1386 (doi:10.1002/cssc.201301114; corrected
      2026-11-09, it is not the J. Am. Chem. Soc. 2011 paper), PEI on fumed
      silica, CO₂ from air (c_in 0.0166 mol m⁻³), four particle-size fractions (<0.25, 0.25–0.5, 0.5–1.7, >1.7 mm).
      Refitted by Auton, Aguareles, Valverde, Myers, and Calvo-Schwarzwalder (arXiv:2311.10161, Applied
      Mathematical Modelling) with a bulk–pore exchange k_p and reversible Sips-type kinetics in the pores
      (dm/dt = k₊c_p(m_max − m)^b − k₋m^b, b = 2), which maps onto `General` chemisorption with
      `UsePoreSurfaceTransport`. Their Tables 5 and 8 give L 0.1 m, v 0.2314 m s⁻¹, ρ_b 541 kg m⁻³, ε 0.48,
      ε_p 0.416, and per size D, k₊, k₋, k_p. Caveats: k_p is not monotonic in the size (2.33e3 m s⁻¹ at the
      smallest), the temperature and humidity of the runs are not stated, and the curves must be digitized.
    - Shi et al., Chem. Eng. J. 496 (2024) (ScienceDirect S1385894724055797), "Bi-disperse adsorption model of the
      performance of amine adsorbents for direct air capture": separates film transfer, intraparticle diffusion,
      and the amine reaction, with adsorption and desorption breakthrough tests on particle sizes; paywalled.
    - Open access, to check: Bellussi et al. (2026), RSC, "Scale-consistent adsorption model for polymer-based
      adsorbents in DAC"; "Mass transfer of CO₂ in amine-functionalized structured contactors in ultra-dilute
      conditions" (PMC11789152), pellets of 0.9 and 0.08 mm.
  - Decision (2026-11-02): validate with Goeppert et al. (Auton refit) or Shi et al. (2024). Next steps:
    - Get the Shi et al. full text and SI. Check whether its bi-disperse model (two intraparticle resistances)
      maps onto RUPTURA's one pore reservoir; if not, take Goeppert.
    - Goeppert: digitize the four particle-size breakthroughs, write `conversions.md` (c_in, v, ρ_b, ε, ε_p, and
      per size D, k₊, k₋, k_p), and set up `validation/3.3.8_pore_transport` with `generate.py` and `run.py`.
      Auton's rate dm/dt = k₊c_p(m_max − m)² − k₋m² is `General` with orders 1, 2, 2.
    - ε_p: RUPTURA has no particle porosity (C-PHY-1: the code is the SI model at ε_p = 1), Auton has
      ε_p = 0.416. Implement `ParticlePorosity` or fold ε_p into k_p and the pore volume, and say which.
    - Fix the SI pore-transport equations first (I-PHY-4).
    - Then replace the "does not test the pore transport" paragraph (Results `sec:results_chemisorption`) and the
      matching sentences in the Abstract and Conclusions.
  - 2026-11-09: set up and run in `validation/3.3.8_pore_transport` (README, `conversions.md`, `generate.py`,
    `run.py`, `analysis.ipynb`, 44 runs, figures in `figures/`).
    - Shi et al. (full text, no SI) does not map: Toth-type amine kinetics with a non-integer exponent, non-isothermal
      with Ergun, correlations only in the SI, Fig. 5 a raster image. Goeppert taken.
    - Data: Auton's own digitization (`Figure_8_9_data.m` in github.com/aguareles/Intra-particle-diffusion), checked
      against the vector markers of their Fig. 9. Auton et al. is Appl. Math. Model. 130 (2024) 827.
    - ε_p needs no change: the pore holdup is negligible (α = 4e-6), and any split of 1/k_p = 1/k_f + R/(5 D_p) over
      film and pores gives the same curve (RMS ≤ 1.1e-6). The smallest size (β = 5.9e5) runs as direct uptake
      (β = 1000: 2.2e-5 from it).
    - Verified: profiles converge to Auton's travelling wave, Eq. (62), as the front develops (l ≈ 55: within 0.16 ks
      up to C = 0.75). At the experimental length the larger fronts are not developed (> 1.7 mm: +9 ks at C = 0.05,
      −27 ks at C = 0.95), so Auton's fitted k_p are those of a developed wave.
    - Against the data without fitting: the published 0.1 m column puts t_1/2 18%, 5%, 10% late and 7% early,
      because Table 5's q_e does not match the data's mass balance for the outer sizes (−12%, ≥ +9%). With the length
      that closes the mass balance of Auton's curves, RMSE 0.042, 0.064, 0.025, 0.068 (Eq. (62): 0.041, 0.047, 0.030,
      0.021).
    - Open: the paper text (above, after I-PHY-4; decide which panels go in); Goeppert 2014 (closed access) for the
      per-size sample masses and flows.
  - 2026-11-20: rerun and in the paper (user's request).
    - Rerun with a binary built from commit ba573b9: every output is byte-identical to that of 2026-11-09, 820 s
      summed over the 44 runs. Notebook: `np.trapezoid` (NumPy 2.5), the Auton curve dotted, the grid-convergence label
      shortened, and Goeppert's Fig. 2 totals in the capacity table.
    - Same day, restyled to `STYLE.md` (no new symbol collisions). In the paper, Auton's T, 𝓛, β, μ, q_e, k_p, R,
      and C are τ_kin, L_kin, N_TU (transfer units over L_kin), Ω, q*_in, k_fp, r_p, and c/c_in; the pore diffusion
      length ℓ_a is δ_{p,a}. The panel legends follow: N_TU, z_1/2/L_kin, Δ(c/c_in). `conversions.md` and the
      README keep Auton's own notation.
    - Goeppert 2014 read (PDF from the user): 25 °C, 335 mL/min, 400–420 ppm, "typically 3 g" weighed after drying,
      8 mm U-tube; no bed length, per-run mass or inlet, or humidity. Their Fig. 2 totals (71.0–81.7 mg/g) are
      Auton's q_e, so the mass-balance mismatch is in the unreported column; L_fit stands for it (equivalent masses
      2.33–3.16 g). README and `conversions.md` §1, §7.
    - Paper (`main.tex`):
      - Results `sec:results_chemisorption`: two paragraphs and Fig. `fig:goeppert` (four `experiment_*` panels).
      - A row of the verification table, and the feature table.
      - The first and the limitations paragraphs of the Conclusions.
      - SI `apx:goeppert`: Table `tab:goeppert` and Fig. `fig:goeppert_verification`.
      - Panels copied to `figures/goeppert/`; bib keys `goeppertnarayanan2014`, `autoncalvoschwarzwalder2024`.
    - The Abstract already lists pore transport without a validation claim, so it is unchanged.
    - Left open, not needed for the paper: a full-model refit of k_p for the two larger sizes (README "Open").
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
  - 2026-11-02: in the paper. Results "Chemisorption & pore transport" (Figs. 8A, 6) and "Reactions" (Fig. 4 and
    the grid study, now second order with WENO3); SE-WGS, isotherm, cyclic loading, and parameters in SI
    `apx:stadler`. With the rerun: X_CO 0.759, CO₂ front 19.7 s; sorbent-only cycle 3 18.0 s against 17.9 s
    (authors) and 15.8 s measured (half heights of the digitized curves). The FirstOrder desorption branch is now
    documented (SI Eq. `eq:first_order_desorption` and the JSON table). The SE-WGS question to the authors stays
    open.
- [ ] **V-PHY-4** 3.3.6 Liquid and pH: discuss with Shrinjay.
- [ ] **V-PHY-5** Adiabatic thermal mode: not covered by Marx. Dantas runs 1–8 (`dantas/`) are adiabatic if needed.
- [x] **V-PHY-7** 3.2.2 MPD: validate against SIAST on the mordenite case of Sharma et al. (2023), with a macrostate
      particle distribution from TMMC (decided 2026-11-02). Fills the empty Results subsection "MPD".
  - System: check which MOR mixture of Sharma et al. 2023 to use and that its pure-component fits are those of
    `examples/MOR-CO2-C3H8` (CO₂/C₃H₈ in He, two Langmuir sites per component, SIAST, 300 K, 2.5e6 Pa).
  - TMMC (RASPA3) of the binary in MOR at 300 K, with particle-number bounds that cover the SIAST loadings up to the
    highest pressure. Estimate the cost of the macrostate grid before starting. Write the distribution in the MPD
    input format (`examples/MPD/README.md`: C order with the last component fastest, normalized at T_ref and a
    reference fugacity, framework mass of the simulated unit cells).
  - Checks, in order: pure-component MPD against the Langmuir fits (force field and fits agree); MPD mixture
    isotherm against SIAST and IAST from the same fits, and against the mixture GCMC data of the paper if it has
    them (digitized); breakthrough of `examples/MOR-CO2-C3H8/breakthrough` with MPD and with SIAST; cost per
    equilibrium evaluation.
  - Stay at T_ref until B-PHY-1 (temperature reweighting) is decided.
  - Update `validation/3.2.2_mpd/README.md`, which still plans the BEA nC7/2mC6 comparison.
  - 2026-11-03, set up in `sharma/` (README, `conversions.md`):
    - System: the equimolar CO₂/propane case of the paper (Sec. 4.1.1, Fig. 6) with the Table 3 fits in mol/kg, and
      the breakthrough of Fig. 10a (10 % CO₂, 10 % propane in He, 1e5 Pa, ε_b = 0.128, ρ_p = 1711.06 kg/m³,
      L = 0.1 m). `examples/MOR-CO2-C3H8` has other fits (CO₂ 4.4/10.0 with 2.91e-4/6.96e-7, …) and another column
      (2.5e6 Pa, 5/5/90, ε_b = 0.4, L = 0.3 m); it is not used.
    - Force field: the paper's description (TraPPE-zeo framework, Garcia-Sanchez CO₂, Dubbeldam UA propane, 12 Å)
      reproduces Table 3 only with the CIF framework charges (Si +2.05, O −1.025, the RASPA2 default): GCMC at 1e4 Pa
      gives CO₂ 1.11 mol/kg (Table 3: 1.18) and propane 0.77 (0.77); with the TraPPE-zeo charges CO₂ is 0.39.
      `~/micromamba/bin/raspa3` (March 2025) ignores the propane bonds and bends (three beads at one point).
    - TMMC: RASPA3 has no multi-component TMMC; the continuous-lambda checkout (`~/sandbox/continuous-lambda`) has
      slices with an edge (TMMC in n_CO₂ at fixed n_C₃H₈, plus pure propane). Grid n_CO₂ + n_C₃H₈ ≤ 200,
      n_C₃H₈ ≤ 100 (decided), 422 windows per temperature, ~7e8 MC steps. Per-macrostate energy cumulants were added
      to that code for this (its C8). 300 K started 2026-11-03 16:37, 350 K follows automatically (B-PHY-1: the
      300 K distribution reweighted to 350 K with and without N ln(T_ref/T), against the direct 350 K TMMC).
    - Pipeline: `sharma/mpd/build_distribution.py` (stitching, MPD file), `sharma/runs/{generate,run}.py`.
  - 2026-11-06: both batches finished (300 K on 11-04, 350 K on 11-05; a RASPA3 cycle is max(N, 20) steps per
    component, so 2.4× the planned steps), but the windows did not converge: every state sampled in production with
    a histogram flat to a factor 10 in 23/422 windows at 300 K and 70/422 at 350 K; ln Π of neighbouring windows
    disagrees by a median 6.8 over the 3 shared states. Cause: the serial conventional TMMC driver keeps the
    collection matrix of initialization (`rezeroAfterInitialization` false, no input key), in which the walk climbs
    unbiased from the lower window bound through freshly created configurations (CO₂ N 40–47 at m = 0: deletion
    acceptance 0.04–0.09 against 0.003 in the window below); the bias derived from those rows keeps the walk out of
    the lower states, so they never correct. The driver has no Wang–Landau exploration for conventional windows
    (only expanded and ParallelTMMC, which allows one component only). Stitching stops at "missing states"; no
    RUPTURA runs. Moved to `sharma/tmmc/unconverged_2026-11-03/`.
  - 2026-11-06: interim fix (you will implement the proper one in the sandbox): new key
    `MacroStateRezeroAfterInitialization` (sandbox C12), 10⁵ production steps per macrostate, equilibration 0.5;
    rerun started 16:44, about 45 h per temperature, downstream chained as before. Short test: the reset removes
    the bottom-state error of the CO₂ windows; the propane edge above ~50 molecules (insertion acceptance
    10⁻⁵–10⁻⁴) stays stuck and needs the proper method. Check after the run whether those states matter for the
    equimolar mixture (SIAST: ≤ 48 propane) and the breakthrough.
  - 2026-11-09: 300 K rerun finished (11-08); stitching stops again (slice m = 20, N 40–82 unvisited at 40–63), no
    RUPTURA runs from the binary distribution. The edge and the m = 0 slice are usable up to 50 propane and 160 CO₂
    molecules (pure isotherms −7…+18 % and −4…+21 % of Table 3) and not flat above; CO₂ has one bad state at N = 36.
    287 of 414 slices with fixed propane have unvisited states or zero-acceptance links (ln Π jumps up to
    thousands), at all m and N, including low loadings; their walks pile up at the window edges. With ≥ 1000 visits
    per state the first unusable link lies at a median n = 14 and no slice reaches n = 42 (SIAST: n = 42, m = 27 at
    1e4 Pa; breakthrough feed: n = 54, m = 35), so these states matter. Not fixed by the reset or the length (the
    11-03 batch: 236/414 by the same test); the slices need the proper serial-driver fix. Details in
    `sharma/README.md`.
  - 2026-11-09: state written into the paper (Results "MPD", marked as draft; `validation/3.2.2_mpd/analysis.ipynb`,
    `figures/mpd_mor/`): pure MPD against Table 3, sampling map of the grid, SIAST and IAST mixture and breakthrough,
    MPD mixture from the well-sampled states up to 400 Pa only. Open: binary distribution, MPD mixture above 400 Pa,
    MPD breakthrough; pure-component reweighting (B-PHY-1) once the 350 K TMMC is done (about 2026-11-10 midday).
  - 2026-11-11: 350 K TMMC done (11-10 16:51); `build_distribution.py 350` stops like 300 K (slice m = 44, N 40–82).
    Pure windows at 350 K: CO₂ flat up to 162; the propane window 0–42 recorded no insertion at 32 → 33
    (acceptance 9e-12), so propane at 350 K holds only up to 32 molecules. Pure reweighting added to the paper
    (Fig. mpd_mor e,f; B-PHY-1 entry). `build_pure_distribution.py` now fills ⟨H⟩ of states without production
    energies by interpolation in N (extrapolation at the ends): with ⟨H⟩ = 0 such states gained up to e²⁰⁰ when
    reweighted. `build_distribution.py` (binary) still writes 0 there; fix that before reweighting the binary.
  - 2026-11-19: binary 300 K distribution from an N-d TMMC done on another machine
    (`sharma/lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt`, 9865 states, valid up to ~1e6 Pa total fugacity), converted by
    `sharma/mpd/build_nd_tmmc_distribution.py`. MPD inlet loading at the breakthrough feed equals the direct
    reweighting (1.1529/0.8141 mol/kg). Breakthroughs with ROS34PW2 (MPD after P-PHY-9 and P-PHY-10; binaries
    `build/ruptura-sharma-2026-11-19-mpd{,-jac}`): t50 CO₂ MPD 3344 s, SIAST 3380 s, IAST 2473 s; propane 2291,
    2240, 1527 s; propane roll-up 1.053, 1.138, 1.705. Wall time MPD 160 s, SIAST 7.4 s, IAST 4.3 s. SIAST is
    within 1–2 % of MPD in t50, IAST far off.
  - Done 2026-11-19 (asked by the user): `sharma/` moved to `validation/3.2.2_mpd/sharma` (symlink at the old path;
    `cases.csv`, `pack.py`, `validation/README.md` updated), all 16 runs of the case redone, and the Results
    subsection "MPD" with Fig. `mpd_mor` rewritten from the binary distribution (draft note removed): trusted range,
    pure components, equimolar mixture (MPD follows SIAST; at the feed propane +11 % and CO₂ +2 % over pure, i.e.
    mildly cooperative), breakthrough and roll-up, cost; the MPD sensitivities in the integrator text, the
    verification table, and SI `apx:sorption_jacobian`/`apx:mpd` (Eq. `macrostate_sensitivity`). The paper reports the
    distribution as two-dimensional NVT+W of Chen and Lin 2023 (`chenlin2023`, added to ruptura.bib), as the user asked
    while the N-d TMMC is in progress. Open: the wall times were taken with other runs on the machine (full
    rerun on one node); the binary temperature reweighting (B-PHY-1) needs an N-d distribution at a second temperature.
- [ ] **V-ITG-2** Rerun the paper cases that use ROS34PW2 with a marching closure, after its stage projection
      (B-ITG-16, 2026-11-11), if the paper should show the current code: 3.4.2 TSA (max |R T Σc/P − 1| 0.019–0.040 →
      below 1e-4, inventories within 3e-4, wall time halved; figures change little, the SI remark and the timing do),
      and the ROS34PW2 runs of the 3.5 benchmark (T15: e50 2–3 times smaller, up to 22 % fewer steps; with V-BC-2).
      3.3.1, 3.4.1, and 3.4.3 move by at most 2e-6 in t50 and 7e-4 of the feed fraction at the outlet; IPOP, FixedVelocity,
      FPIV, and liquids are unchanged. Not done: the user decides.
  - 2026-11-19 (asked by the user): binary frozen from the working tree on `db47591` at 18:06. 3.4.2 TSA rerun
    (`validation/3.4.2_swing_adsorption`, README and figures updated): max |c_tot R T_g/P − 1| 9.7e-5 (was 0.040),
    inventory change at the step ends ≤ 2.8e-6, 0.59–5.6 s per cycle (paper: 1–12 s), Δq_cyc and temperature swings
    unchanged to three digits. The CO₂ balance is now 1.8–3.7e-4 (was 0.7–1.7e-4; the paper says "within 2e-4"): the
    projection gives the integration error of Σc to the components by mole fraction (falls with rtol), on top of a floor
    of ~6e-4 over four 1000 s cycles that the old binary had too (README). The whole 3.5 benchmark (80 runs) is
    rerunning in one batch (`benchmark/breakthrough/run_2026-11-19.log`).
  - 2026-11-19: the binary frozen at 18:06 predates the B-ITG-17 fix (ROS34PW2 no longer clamps the state its next
    step starts from), so these TSA and benchmark runs do too. Every ROS34PW2 run with WENO3 can move a little; check
    whether the CO₂-balance floor of the TSA is that clamp before writing the numbers into the paper.
  - 2026-11-19 19:40: benchmark part 1 (75 runs) done at 18:54. Two runs stopped: the RK3 stability run of T15 at
    5e-2 s (IAST NaN at 12 s, still "unstable") and T15 CVODE SPGMR at rtol 1e-2, which stops at 8.4 s with
    max |R T Σc/P − 1| = 0.49; on 2026-11-03 it completed with e50 7.4e-4. The paper's claim that the T15 errors fall
    monotonically with the tolerance needs that run. The BEA Dense part was killed with the Claude session at about
    19:12 and restarted detached at 19:36. The paper already has the TSA numbers of this binary (main 2.4, Section
    `results_cycles`, verification table, SI remark), including "CO₂ balance within 4e-4" and an SI sentence that
    puts part of the balance error on the projection. Both need a TSA rerun with the B-ITG-17 fix. A full rerun on
    one cluster node is described in `validation/README.md` (`cases.csv`, `run.py --paper --list`, `jobs.csv`).
  - 2026-11-19 19:58: TSA rerun with a second frozen binary that has the B-ITG-17 fix. Every number of the README
    table is the same to the digits shown, and so are the gas departure and the inventory change. The CO₂ balance
    is the same to four digits, also on four 1000 s cycles at the default rtol and at 1e-8. The clamp is therefore
    not its floor, and the paper text stands. The wall times of this rerun (1.05–8.6 s per cycle) are inflated by
    the BEA Dense runs alongside; on its own a four-cycle run takes the same CPU time with both binaries. The paper
    keeps 0.6–6 s per cycle; time it again on a quiet machine.
- [ ] **V-BC-2** Rerun the paper cases affected by the C-BC-1 closure fix (continuity closure with unequal D):
      3.3.1 Marx ternary, 3.4.1 boundary conditions (Marx 10 bar), 3.4.3 sensitivity (T15; the D_ax variations up
      to 1e-3 m²/s will move most), 3.4.2 cycles, 3.5 benchmark (T15). Checked on T15 (2026-11-03): ≤ 2.2 % of the
      final H₂ fraction at the H₂ front, t50 unchanged. Then consider the dispersion law (F-PHY-6) for Stadler
      (Eq. 12, pressure swing) and Stampi-Bombelli (literature correlations; IVOP now possible).
  - 2026-11-04: 3.3.1, 3.4.1 (10 bar), 3.4.2, 3.4.3 rerun with B-BC-2 and C-BC-2 as well (C-BC-2 entry). Open: the
    3.5 benchmark (`benchmark/breakthrough`, `benchmark/convection` T15; timed parts on a quiet machine) and the
    3.4.2 wall times per cycle.
- [ ] **V-BC-4** The convection study (`validation/3.5_numerical_performance/convection`, Fig. `convection`, SI
      `tab:convection_benchmark`) ran on 2026-10-31 with the Dirichlet inlet, the default then. Its stored inputs
      have no `InletCondition`, so with today's default they run with Danckwerts. Its `generate.py` copies the base
      inputs (Marx T15, the monolith), which carry `"InletCondition": "Danckwerts"` since 2026-11-03. Either pin
      Dirichlet in `generate.py`, as the integrator benchmark does, or rerun with Danckwerts and update the paper.
      Found 2026-11-19 when every generator was rerun (`validation/README.md`, step 3).
- [ ] **V-CYC-1** 3.4.2 PSA, and TSA purity/recovery: later. They depend on flow reversal and the cycle driver
      (F-BC-1, F-CYC-1). Closed ends and pressure ramps are there since 2026-11-09 (F-BC-2), not validated against a
      published PSA yet.
  - 2026-11-06: flow reversal is validated by 3.4.4 (V-BC-3). The published PSA/VSA experiments (Schell et al.
    2013; Haghpanah et al. 2013) need pressurization and blowdown against closed ends (stage 3); no
    single-column experiment with counter-current steps at constant pressure was found.
- [x] **V-EX-1** Rerun the campaign of `examples/Self-driving-lab-H2-PSA` with the cyclic-steady-state stop and the
      stream-amount KPIs of F-PY-3 (2026-11-19). Its README, notebook, figures, and the numbers of the paper are from
      the earlier stop in Python (136 experiments, 12 not converged).
  - 2026-11-19: one-column campaign rerun (136 experiments, 1 not at cyclic steady state, hypervolume 5.0 times
    the baseline); rerun from scratch it gives the same setups and identical KPIs. README, paper text, SI model check
    (redone with the new stop), and `tab:sdl_front` updated. The two-bed campaign (`TwoBedPressureSwingCycle`) ran
    once resumed after round 0 (B-PY-8), so it is being rerun in one go; then the notebook, figures, README, and paper
    two-bed results.
  - Done 2026-11-19: the two-bed campaign, rerun in one go (57 min), dominates the one-column front (recovery
    0.51-0.56 at 1.04-1.28 mmol/(kg s); hypervolume 8.7 times its baseline, still rising). Equalization check (e = 0)
    and a fine check (160 cells, ROS34PW2 1e-6) done; notebook executed, figures copied to the paper (new
    `fig:self_driving_lab_two_beds`, `tab:sdl_front_two_beds`), README and paper text updated.
- [ ] **V-CYC-3** Validate the cycle driver, `MomentumBalance`, and the coupled columns (F-CYC-1, F-BC-3, F-CYC-4;
      the two-bed cycle with pressure equalization can be run since 2026-11-18) against the literature.
      2026-11-19 (search of the branch, at the user's question): it has no literature cycle case that can be rerun.
      Its Haghpanah 2013 (breakthrough and four-phase VSA), Krishnamurthy 2014 (breakthrough, one-column VSA, two-column
      VSA runs 6 and 7), and Marx T15 inputs under `examples/new_examples/integrated/validation/` are listed in its
      manifests but were never committed (its `.gitignore` drops `validation/` and `test_*`); the digitized markers of
      the donor project are missing too; `verification/integration/reference_outlets` holds outputs of its own code
      for breakthrough runs only (correcting the lines below: no VSA or two-bed outlet curves). One literature-based
      input is there, in the donor YAML format: `tests/process/fixtures/upstream_materials/mocca/dcb_adiabatic/
      simulation.yaml`, the 13X breakthrough of Haghpanah 2013 (dual-site Langmuir, 201 nodes). Its own documents say
      the literature validation was not completed. Next: ask Shrinjay for the donor folders
      (`examples/validation/haghpanah_2014`, `examples/literature/mocca`), or digitize the four-phase VSA of
      Haghpanah et al. 2013 and runs 6 and 7 of Krishnamurthy et al. 2014 from the papers. So far they
      are checked against the code itself (the continuity closure, mass balances, P V/(R T)). Candidates: the
      four-step VSA of Haghpanah et al. (2013), one column with light-product pressurization, which the phase-name gas
      of F-CYC-1 allows, and a two-bed cycle with pressure equalization once F-CYC-4 is there; the branch
      `ruptura-shrinjay` has outlet curves of its own runs for both (`verification/integration/reference_outlets`) but
      not their inputs. Also: whether `CyclicSteadyState` ends the 3.4.2 T15 cycles where the paper's 1 % criterion
      does (14 cycles for 1000 s half-cycles). Which cases, and digitized from which figures, is to be decided.
  - 2026-11-19 (literature search; paper section `sec:results_vsa` is a \todo skeleton): the cycle of Haghpanah et
    al., Ind. Eng. Chem. Res. 52 (2013) 4249, is feed pressurization, adsorption, co-current blowdown, and
    counter-current evacuation on 13X, not light-product pressurization (correcting the line above; that cycle is in
    their AIChE J. 59 (2013) 4735). The main text is closed access and was not obtained; its own SI (ACS figshare, CC
    BY-NC: dimensionless groups, boundary conditions of all four steps) and open restatements of its Tables 2 and 3
    were (Glover et al. 2026 arXiv:2606.21523 with SI Table S2; Subraveti et al. 2026 Ind. Eng. Chem. Res. 65 9656,
    CC-BY; SINTEF Mocca.jl). Reference data available: the adiabatic breakthrough figure (a copy in Ramani, TU Delft
    MSc 2016, Fig. 3.3a), the base case (FP 15 s, ADS 15 s, BD 30 s, EVAC 40 s; 1/0.2/0.1 bar; 1 m/s), and five
    minimum-energy designs at 90 %/90 % with productivity and energy. Open: (i) the base-case purity is disputed
    (88.7 % quoted as the paper's, about 80 % in two reimplementations); (ii) its LDF coefficient depends on the
    concentration, k = (c/q*) 15 ε_p D_m/(τ r_p²), which RUPTURA does not have (code change, ask); (iii) the heat of
    adsorption of its energy balance is not stated in any source; (iv) the pressurization time of the designs is not
    stated. Better with the PDF of the paper from the user.
- [ ] **V-MBD-1** 3.3.7 Park et al. layered beds (Chem. Eng. Sci. 53 (1998) 3951): obtain the paper and set it up.
      Backup: Yang and Lee, AIChE J. 44 (1998) 1325.
- [ ] **V-GEN-1** 3.1 Verification: all items marked new in `validation/3.1_verification/README.md`.
- [x] **V-PHY-8** 3.3.4 Liu: replace the frozen humid laws with the coupled ones of F-PHY-1 in `liu/generate.py`
      (Fig. 4 `GAB-modified-Toth` [1.0, 0.98692, 0.25, 0.006, 25.0], Fig. 7 `Exponential-humidity-Toth`, H₂O with
      `"Water": true`), rerun, and drop the "frozen" fidelity limit from `liu/README.md` and `conversions.md` §10.
      Scratch runs on 2026-11-04: Fig. 4 unchanged to 1.2e-4 c₀; Fig. 7 t50 −1.5%, so the humid/dry gap (2.60 against
      1.84 digitized) is not from the freezing. Fig. 5b at RH 0.85 reproduces exactly as a MixturePrediction with
      `RelativeHumidity`.
      Done 2026-11-04: `liu/generate.py` writes the coupled sites (dry cases keep Toth); the six humid cases were
      rerun (identical to the scratch runs), the notebook re-executed, and README, `conversions.md` §7-10, and the
      notebook text updated.
- [x] **V-BC-3** 3.4.4 Flow reversal against equilibrium theory (`validation/3.4.4_flow_reversal`, 2026-11-06).
      Trace adsorbate A (1e-3 in He, Langmuir b p_F = 4) loaded to z_f = 0.504 L, then purged with helium
      co-currently or counter-currently; `FixedVelocity` (trace limit, B-ITG-16), ROS34PW2, 25–400 cells, WENO3
      and Upwind. Reference: the method of characteristics (Rhee, Aris, Amundson 1986), with the shock path after the
      co-current wave overtakes the shock integrated by RK4. Outlet E_1 with WENO3: 5.0 → 1.24 → 0.33 %
      (counter-current, 25/100/400 cells), 5.0 → 1.29 → 0.38 % (co-current); Upwind 12 → 4 → 1.3 %, first order
      in both. 99 % regenerated after 2.09 t_st (counter-current) against 4.38 t_st (co-current); RUPTURA 2.18/4.47
      on 100 cells, 2.11/4.41 on 400. Paper: Results `sec:results_flow_reversal` (Figure `fig:flow_reversal`), SI
      `apx:flow_reversal_study` (settings, theory, E_1 table, amount-left figure); abstract and Conclusions name it;
      the SI remark on the isothermal TSA now describes the back-flow (B-CYC-1) and the inventory jump (F-CYC-2).
      References added: rheeamundson1986, ruthvenruthven1984, ruthvenknaebel1994, skarstromskarstrom1960.
- [x] **V-GEN-4 CVODE linear-solver default changed on 2026-11-02** (Dense → SPGMR, B-ITG-12). The maximum-order
      change of the same day (5 → 1, B-ITG-11) was reverted on 2026-11-03, so only the solver remains:
  - 40 `stampi-bombelli/` and all 50 `liu/` inputs rely on the default solver (Dense when they were run); the others
    set SPGMR. Done for `stampi-bombelli/` on 2026-11-03: `generate_runs.py` writes `CVODELinearSolver` Dense (and
    the order-1 overrides of B-ITG-11), and the 20 inputs were regenerated; only that key changed. `liu/` remains;
  - decide per study: rerun with SPGMR, or write `CVODELinearSolver` Dense into the generators so that the inputs
    reproduce the paper.
  - Moot since 2026-11-03: marx2/ternary (24), `validation/` (68), `benchmark/breakthrough` (40), `stadler/` (61),
    `liu/` (50), and 47 of the 51 `stampi-bombelli/` inputs relied on the default order.
  - 2026-11-03: with SPGMR, `stampi-bombelli/dl_from_literature/monolith_5.6pct_0p02mmol` (order 1) stops with
    `mxstep` at t = 883 s; with Dense it finishes in 21 s. So the Stampi-Bombelli inputs cannot simply be rerun with the
    new default. The inlet-condition reruns pin Dense in their copies.
  - 2026-11-03: ROS34PW2 runs every Stampi-Bombelli and Liu input at its default settings (B-ITG-11, B-EX-10). An
    alternative to pinning Dense is to switch these studies to ROS34PW2 (user to decide).
  - Done 2026-11-03 (decision): the validation studies use ROS34PW2 (I-BC-4), so no study depends on the CVODE
    default solver any more. The integrator benchmark sets its CVODE solvers explicitly.
- [ ] **V-PHY-9** Rerun the cases with the energy balance after B-PHY-15 (isotherms at T_s), with one frozen binary:
      3.3.1 Marx T1–T24, the 3.3.2 Stampi-Bombelli 5.6 % packed bed, the 13 Elfving cases (3.3.3, chemisorption
      isotherms), the 3.4.1 Marx runs, 3.4.2 TSA, 3.4.3 sensitivity, the T15 runs of 3.5 (benchmark, convection), the
      self-driving-lab campaigns, and the T15 checks of the paper's Table 2 (free gas pressure, phase end). Isothermal
      cases are unchanged.
  - 2026-11-19: 3.3.1 Marx T1–T24 done (with the Marx inputs, V-PHY-10). The other T15 studies read the frozen
    `marx2/reference_settings_hassan_chung/` inputs, so their reruns keep the reference settings.
- [ ] **V-ITG-3** Fair CVODE comparison for the 3.5 benchmark (Reviewer-2 audit of 2026-11-19), added to
      `validation/3.5_numerical_performance/breakthrough` but not run. `generate.py` adds CVODE without the 4 Δz/v cap
      (`cvode-nocap`), with the analytical sorption Jacobian as SPGMR preconditioner (`cvode-precond`) and in KLU
      (`cvode-klu`, BEA only: KLU rejects the energy balance and chemisorption), and the case
      `monolith-400ppm-fitted-dl-danckwerts` (the inlet of the validation runs) with every series; `run.py --repeats N`
      keeps the median and all wall times (`elapsed_all`). 147 runs; every new kind of input smoke-tested for 4 s.
  - To do: `python generate.py`, then the whole benchmark with one frozen binary and `--repeats 3` (the T15 runs change
    with B-PHY-15 anyway; with V-ITG-2); extend `analysis.ipynb` to the new series and case; then Fig. 10, Tables
    S21–S25, and the speed-up claims (abstract, Section 3.4, Conclusions). A direct analytical-Jacobian CVODE for T15
    would need KLU with the thermal block.
- [x] **V-PHY-10** Marx rerun with the inputs of Marx et al. (mock review 2026-11-19, Major 4): their
      temperature-dependent Sips (q_sat(T), s(T) of their Table 5), the porosity mapping ε_b = ε_t with ρ_s =
      ρ_b/(1 − ε_t) (the current mapping holds 273 instead of 502 kg/m³ of adsorbent), and a 318 K ambient/wall for
      the 318 K runs; compare with their model curves as well as the data.
  - 2026-11-19: `marx2/generate.py` writes T1–T24 (conversions.md §5–9): SCI (their multicomponent Sips) with
    uncentred VantHoff q_sat and k and the new `Arctan` law (1/ν = 1/V + α atan[β(T − T_ref)], keys Alpha and Beta;
    temperature_law.{h,cpp}, inputreader, isotherm masks, pybindings, test
    `TemperatureInput.ReproducesTheTemperatureDependentSipsOfSchellEtAl`); ΔH and k_i of Tables 5–6; ε_t as fluid
    fraction, ρ_s = 1946.5 kg/m³, v = Q/(Aε_t); molar `GasHeatCapacity` 29.51 on every component; h_fs from
    Wakao–Kaguei–Funazkri per run × (1 − ε_b)/(1 − ε_t); ambient = feed temperature; `InletPressureInletVelocity`
    and full times for all runs (also T8, T21: V-PHY-1). Model curves digitized (`digitize_model.py`), temperatures
    at 1.1 m. Rerun with a frozen binary (T_s isotherms), 27 s for all 24. CO₂ t50 at ≥ 5 bar 1–12 % early (Marx
    model −9 to +2 %); 1 bar up to 25 % late (Marx model within 2 % except T19); peaks 0.2–5.9 K below measured.
    Paper: 3.5.1, Fig. 11 caption, Conclusions, SI `apx:marx` (new mapping paragraphs, Eq.
    `wakao_gas_solid_heat_transfer`, Tables S26–S31, Figs. S2–S5), law tables and JSON reference (Arctan).
  - The old inputs are frozen in `marx2/reference_settings_hassan_chung/` (user's choice, 2026-11-19): the tests and
    every other T15 study (3.4.1, 3.4.2, 3.4.3, 3.5 breakthrough and convection, SDL, Table 2 checks) read them, and
    the SI lists them in Table `tab:marx_reference_settings`. The cross-check reruns
    (`inlet_condition_dirichlet_vs_danckwerts`, `ros34pw2_danckwerts_rerun`) still read `marx2/ternary` and pick up
    the Marx inputs when regenerated.
- [ ] **V-ITG-4** Like-for-like integrator benchmark (mock review Major 2): CVODE without the 4Δz/v cap (or
      the same cap for ROS34PW2), CVODE with the analytical Jacobian (KLU, or SPGMR with the SorptionJacobian
      preconditioner), a Danckwerts monolith case, three repeats; retest the Toth stalls with CVodeSetConstraints.
  - 2026-11-19: everything but the CVodeSetConstraints retest is in the suite (V-ITG-3, not run).
- [ ] **V-PHY-11** Monolith validation that exercises the monolith (mock review 2026-11-19, Major 5a): the
      Stampi-Bombelli runs prescribe the velocity and use fitted D_ax (vL/D = 0.04–2), so channel permeability,
      Taylor–Aris dispersion, film transfer, and the monolith energy balance are untested; channel coefficients are
      circular for every shape (square: Taylor–Aris 1.6×, fRe 56.9). Compare model-to-model with the published model
      of Stampi-Bombelli et al. (it reproduces their monolith curves; RUPTURA is −15 to +24 % off), give the packed-bed
      inputs in the SI (rate constant, ε_b, energy parameters, fitted D law), and add a case with a computed velocity
      and pressure drop. Details: `ruptura-2.0-paper/REVIEW_AUDIT.md`.
- [ ] **V-PHY-12** The paper now uses two input sets for experiment T15: the Marx inputs in 3.5.1 and the frozen
      reference settings of Hassan and Chung (IAST, ε_b, h_fs = 20, ambient 298 K; SI Table
      `tab:marx_reference_settings`) in the verification, cycles, SDL, benchmark, boundary-condition, convection, and
      sensitivity studies (2026-11-19). A referee may ask why. Optional: move those studies to the Marx inputs
      (repoint their generators and the tests to `marx2/ternary`; 9 tests then need new fixtures, e.g. the
      `heatCapacityGas` key and the >10 Pa Ergun drop of `MatchesTheContinuityClosureInABreakthrough`), rerun, and
      update their numbers; or keep both and say so in the main text.
- [ ] **V-CYC-4** Conservation and cyclic steady state of the cycle results (mock review 2026-11-19, Major 6): the gas
      rescaling at phase transitions, the ROS34PW2 projection shared by mole fraction, and the clipping of negative
      values add or remove gas without a flux; Eq. `main_css` measures stream amounts against the largest one (the
      feed), so with δ = 2e-3 the product N₂ that decides the 99.9 % target may change by more than its size per
      cycle. Report per-component balances relative to each component's feed and the purity/recovery change over the
      last cycles for the TSA and SDL Pareto cycles, and N_css; consider a KPI-based stop and a flux-consistent
      projection. Also verify the free gas pressure dynamically: pressurization of a closed column against the
      linearized pressure-diffusion solution, and an adiabatic filling test.
- [ ] **V-PY-1** Self-driving lab robustness (mock review 2026-11-19, Major 7): 3–5 seeds per campaign and a
      Sobol/random baseline with the same 136 evaluations; justify the bounds the fronts sit on (T at its lower
      bound; P_H, P_L, e); a third grid for the N₂ impurity (20–68 % between 80 and 160 cells); state the gas
      heat-capacity model of each campaign (the T15 per-mass value gives 2.2 J/mol/K for H₂); time a campaign with
      another integrator or drop "this speed makes ... practical" from the abstract.

## Paper writing (W)

- [ ] **W-GEN-1** Results headings: mirror the `validation/` structure; keep "Sensitivity studies" (3.4.3).
- [ ] **W-GEN-9** Front matter and statements for submission, none of which are in `main.tex`: a real title (now
      "RUPTURA 2.0"), `\articletype` (now "ARTICLE TEMPLATE"), keywords; disclosure statement, code and data
      availability (repository URL and release DOI, see W-GEN-2), funding, acknowledgements, and author
      contributions.
- [ ] **W-GEN-10 Readability pass** (deferred on 2026-11-02 until the fixes are in). Main text: about 7000 words of
      prose, 25 words per sentence, 44 sentences over 35 words (Reactions 40, Conclusions 36, Boundary conditions 35
      words per sentence). Suggestions:
  - Section 2: replace the 20 short subsections that end in "given in Supplementary Section X" by one feature
    table (feature, equation, SI section, where validated) and a few subsections on what is new;
  - Results: lead with the finding, keep two or three numbers per claim, and move the rest to SI tables (the
    benchmark section has about 40 numbers);
  - split sentences over 35 words; fewer acronyms in the main text (K-HDMO, SE-WGS, PFR, SPGMR, BDF);
  - Introduction: about 1000 words before what is new; shorten the background and the monolith paragraph;
  - Section 2.3 in the voice of the other sections; parallel Results headings without "&" (W-GEN-1);
  - a short box with the user recommendations (feed rate at the inlet, CVODE with SPGMR, WENO3 on 25 points).
- [x] **W-GEN-11** Audit and rewriting plan of the paper: `ruptura-2.0-paper/REWRITE_PLAN.md` (2026-11-18).
      - **Contents:** the story, a new outline with topic sentences, a move map, cuts, an SI plan, and decisions
        D1–D7: MPD, liquid, one cycle code or two, pore transport, the verification table, the SI JSON reference,
        and the venue and terms.
      - **Relation to W-GEN-10:** the plan includes its suggestions.
      - **Workflow:** the rewrite goes into a copy, and `main.tex` stays the source until the authors accept it.
  - 2026-11-18: draft in `ruptura-2.0-paper/main_rewrite.tex`, built on the user's decisions:
    - D1: MPD kept as is.
    - D2: the liquid case stays untested; the empty heading carries a `\todo`.
    - D3 and D4: skipped for now; the paragraph on the companion cycle code was moved to Section 2.3 with a `\todo`.
    - D5: verification table added.
    - D6: the JSON reference stays in the SI.
    - D7: Molecular Simulation, US spelling, "phase" for cycle steps.
    - The introduction has only minor changes.
    - Changes since `main.tex`:
      - Section 2 is merged and motivated, with a feature table.
      - The results are grouped as in `validation/`.
      - The ROS34PW2 claim is now backed in Section 3.5 (CVODE stalls).
      - New SI table: `tab:benchmark_summary`.
    - 26 `\todo` markers for the author. The labels and the set of `\cite` keys are unchanged.
    - The main text grew from 11,457 to about 12,900 words.
  - 2026-11-19: the main text of `main_rewrite.tex` was condensed by 10% on wording only, from 12,914 to 11,628 words.
    - Every number, label, citation, reference, equation, figure, and `\todo` is unchanged.
    - The MPD section and the paragraph on the companion cycle code were not touched.
    - The introduction got light wording edits only.
  - 2026-11-19 (user request): the rewrite replaced `main.tex`, and `main_rewrite.tex` and `main2.tex` were removed.
    The new `main.tex` compiles to 182 pages with no undefined references. Line numbers `main.tex:NNN` in the audit
    items (W-GEN-12 to W-GEN-14, I-CYC-2) refer to the old `main.tex` of 2026-11-11; find those passages by their
    text or label. The open author items are the `\todo` markers in `main.tex`.
- [x] **W-GEN-12** Numbers that disagree within the paper (audit 2026-11-18):
  - Marx low-flow runs: "11–37 % later" (`main.tex:515`), but `tab:marx_t50_comparison` gives T7 1611/1095 s, i.e.
    47 % later.
  - Stampi-Bombelli, five versions:
    - abstract: "within 17 %";
    - results: 7–19 % late (packed bed), and 15 % early to 8 % late (monolith);
    - SI at 400 ppm: 13–24 % late;
    - SI inlet section: "within 24 %, 61–133 %";
    - conclusions: "within 19 %", for the monolith.
  - Abstract, "CVODE is 3 to 20 times faster": the results give about 3 (T15, at the largest stable RK3 step) and
    232/11 s = 21 (monolith, at different accuracies). State the basis of the comparison.
  - Abstract, "WENO3 on 25 grid points is more accurate than upwinding on 200": shown for BEA only.
  - Marx adsorbent peaks "0.3–6.3 K below": the table minimum is 0.2 K (T3, T6).
  - μ_f: 1.9e-5 (`main.tex:5001`) against 1.92e-5 Pa s (`tab:simulation_settings`).
  - 2026-11-19, all resolved in `main.tex`:
    - Marx: t50 and peak temperatures recomputed from the stored `marx2/ternary` outputs (2026-11-04 run) reproduce
      `tab:marx_t50_comparison` exactly (temperatures within 0.1 K). The text took 11–37 % from the 2026-11-03 run
      (T7 1503 s, T19 1164 s in `validation/inlet_condition_dirichlet_vs_danckwerts/comparison.csv`); now 11–47 %, and
      the peaks 0.2–6.3 K (adsorbent) and 1.7–13.5 K (gas) below the measured ones.
    - Stampi-Bombelli: the remaining versions measure different sets and agree (5.6 % runs: packed 7–19 % late,
      monolith −15 % to +8 %; 400 ppm monolith 13–24 % late; all runs within 24 %, widths 61–133 %). The conclusions
      now give 19 % for the 5.6 % runs and 24 % for the 400 ppm runs; the SI inlet sentence names its scope.
    - The two abstract claims left with the rewrite.
    - μ: the correlation was evaluated with 1.9e-5 Pa s (`marx2/heat_transfer.py`); the SI now calls it μ_0 and notes
      that μ_f = 1.92e-5 would make h_fw 0.3 % smaller.
- [ ] **W-GEN-13** Introduction statements and citations to check (audit 2026-11-18):
  - `main.tex:101`: toPSAil and pyAPEP are called "Python- and MATLAB-based", the reverse of their languages.
  - `main.tex:103`: the RUPTURA 1 limit "Δt/(Δz)² … arising from fast adsorption kinetics" mixes the dispersion
    limit with the sorption stiffness, which does not depend on Δz (`eq:sorption_stiffness`).
  - `main.tex:91`: "roughly 15 % of global industrial energy use" for phase-change separations. Check the quantity
    in Sholl and Lively.
  - Missing citations:
    - CVODE/SUNDIALS is cited nowhere, and no entry exists in `ruptura.bib`.
    - SSP-RK3 is uncited, although `shuosher1988` and `gottliebgottlieb2003` are in the bib.
    - 2026-11-19: done. CVODE cites the two papers SUNDIALS asks for, `hindmarshwoodward2005` (ACM TOMS 31, 363)
      and `gardnerbalos2022` (ACM TOMS 48, 1), new in `ruptura.bib`; SSP-RK3 cites `shuosher1988`, whose author field
      was all caps and is now "Shu, Chi-Wang and Osher, Stanley". Both in the introduction and in Section 2.4.
  - 2026-11-19: the `\todo` markers for the three introduction points were removed from `main_rewrite.tex` at the
    user's request; the points themselves stay open here.
- [ ] **W-GEN-14** Draft leftovers (audit 2026-11-18; the front matter is W-GEN-9):
  - the MPD "Draft status" note (`main.tex:634`);
  - the empty `\subsection{Liquid breakthrough}` (`main.tex:583`);
  - the one-line caption of `fig:marx_high_flow`;
  - SI sentences truncated after "… are" (`main.tex:1747`, `2932`, `2984`, `3273`);
  - "in the main text(a)" and "(b)" in place of figure references (`main.tex:1648`, `1675`);
  - a LaTeX comment left in the text (`main.tex:2231`);
  - an SI table cited as plain "Table" (`main.tex:704`).
- [x] **W-PHY-1** The SI "humidity Toth (Climeworks-type)" correlation (Eqs. `climeworks_humidity_toth` to
      `climeworks_heterogeneity`) has no source. Web searches on 2026-11-04 did not find this form; the candidates
      Young et al., Energy Environ. Sci. 14 (2021) 5377, and Gebald et al., Environ. Sci. Technol. 48 (2014) 2497,
      could not be read (RSC and ACS refuse scripted downloads). Cite the source, or remove the model from SI and code.
      2026-11-04: removed from the code, app, web app, bindings, docs, tests, and SI (source unknown to the user).
- [ ] **W-GEN-15** Open text items of the mock review of 2026-11-19 that need author decisions (the textual
      corrections were applied the same day): notation clashes (κ, ψ, λ, R, Φ); JSON keys broken inside words in
      Tables S11/S16 (smaller font, wider column, or a continuation mark); source of the water Antoine constants
      (Eq. S142; Antoine 1888 has none); original sources of the Table S8 isotherms; intro positioning against
      CADET-Process, IDAES, Yancy-Caballero et al.; cause of the SPGMR failure with MomentumBalance (SI S2.1);
      E1 used for two metrics (Eqs. S206, S215); TODO markers, VSA section, MPD draft note, code/data availability.
- [ ] **W-GEN-16** Remaining text items of the mock review of 2026-11-19 that W-GEN-15 does not hold; the per-point
      status table is Section 4 of `ruptura-2.0-paper/REVIEW_AUDIT.md`. Still in the text on 21:08:
      - the BET denominator `(1-k_2+k_1\xi)` of Table S8;
      - "long-time limit" for the Glueckauf factor (SI S5.5.1);
      - the single-pressure vacuum-work formula of Table 2;
      - "OBRIEN" in `ruptura.bib`.
      Not rechecked:
      - Table S6 row label and Taylor–Aris coefficients for non-circular channels;
      - the grid of the 10.8 × P/RT claim and the 8e-3 residual;
      - coupling capacity, valve cap, and two-bed valve law;
      - which gas heat capacity each run used;
      - Gibbs consistency of SCI and the one-way SPI water coupling;
      - WENO3 order and its absence from Table 2;
      - setups of the Table 2 checks;
      - the sensitivity table, the Q parameterization, and the ΔT_max location;
      - LHHW K_r vs K_eq;
      - placeholder input defaults;
      - uncertainties;
      - figure fonts, labels, and color-only lines;
      - Eq. S171 and the partial-pressure definition.

