# What RUPTURA lacks for Stadler et al. (2022), and how to add it

Status: gaps 1–5 are implemented (2026-10-21); see "Status after implementation" at the end. The sections
below are the original analysis. Case setup: `stadler/`. Validation notes:
`validation/3.3.5_reactions_stadler/README.md`.

Two experiments run today: the Fig. 5 isotherm, and the sorbent-only adsorption step on a clean bed
(Figs. 6 and 8A). Everything else is blocked by four gaps. The table below also lists a mass-balance
error in the existing pore model, found during this analysis, and items that looked blocking but are not.

| # | Gap | Needed for | Blocking | Proposed fix | Size |
|---|---|---|---|---|---|
| 1 | A reaction among non-adsorbing gases cannot exchange with the bulk gas | Figs. 4, 7, 8B | yes | Catalyst-weighted reaction on the bulk concentration (§1, option A) | small |
| 2 | No catalyst amount and no catalyst placement along z | Figs. 7, 8B (Fig. 4: only the amount) | yes | `CatalystDensity` per `Adsorbents` entry, and at top level (§2) | small |
| 3 | No coverage-dependent desorption (Eqs. 7–10) | cycles, Figs. 6–8 | yes | Optional desorption branch in the chemisorption `FirstOrder` law (§3) | small |
| 4 | Swing phases cannot change feed composition or flow | cycles, Figs. 6–8 | yes | `GasPhaseMolFractions` and `ColumnEntranceVelocity` per phase, and CVODE re-setup at each switch (§4) | small–medium |
| 5 | The pore/chemisorption coupling does not conserve mass | not used by Stadler with option A | no, but it is a bug | Make γ consistent with the reservoir volumes (§5) | small (fix), large (paper model) |
| 6 | Repeating N cycles | cycles | no | Already possible: list the phases N times, as `validation/3.4.2_swing_adsorption/generate.py` does | none |
| 7 | Axial dispersion per phase | cycles | no | Numerical dispersion is larger than Eq. 12 (§6) | optional |
| 8 | Outlet dead volume | all breakthroughs | no | The notebook convolves the outlet with the measured CO tracer | optional |
| 9 | Run time: 248 s for one 900 s step with the dense CVODE solver | cycles | no | `"CVODELinearSolver": "SPGMR"`: 9.6 s, outlet within 3e-4 | input change |

## 1. Reaction coupling (gap 1)

### What the code does now

- `Reaction::Phase` is `Physisorbed`, `Chemisorbed`, or `PoreConcentration` (`src/reaction.h`).
- The WGS species CO, H₂O, and H₂ do not adsorb, so adsorbed-phase reactions cannot carry them.
- Pore and surface concentrations are state blocks per chemisorption site: the layout has 2 + 3S
  component blocks (`ColumnStateLayout`, `src/column.h:29`).
  - A component exchanges with its site's pores only if it has a `General` or `Elovich` site with
    `UsePoreSurfaceTransport` (`Chemisorption::usesSurfacePoreTransport`, `src/chemisorption.h:47`).
  - For any other component, a `PoreConcentration` reaction changes a pore state that never reaches the
    bulk gas.
- The paper's SI describes something else: one pore phase per adsorbent and component, with particle
  porosity ε_p and factors β_p and γ (`main.tex`, Supplementary "Pore transport", Eq.
  `chemisorption_pore_capacity_factor`). The code has no ε_p.

### Stadler's model is pseudo-homogeneous

Their Eqs. 15–16 carry a pore phase, but the fitted model does not depend on it. An ideal plug-flow
reactor needs no pore phase and no dispersion:

$$\frac{dX}{dW} = \frac{k_\infty e^{-E_a/RT}\,P^2 y_{\mathrm{CO},0}^2\left[(1-X)(S/G-X) - X^2/K_{\mathrm{eq}}\right]}{F_{\mathrm{CO},0}}$$

Integrated over W = 0.5 g of catalyst, it reproduces the authors' digitised Fig. 4 model lines to
|ΔX| ≤ 1e-4 at 225, 250, and 275 °C. Against the measured points, the largest deviations are 0.038,
0.030, and 0.025.

The time scales agree. Bulk↔pore exchange is (1−ε_b)ε_p k_LDF/ε_b = 225–1000 s⁻¹ for
D_eff = 2e-7 to 9e-7 m² s⁻¹. The fastest WGS pseudo-first-order constant is 26 s⁻¹ (Fig. 4, 275 °C,
S/G = 4.5), and CO₂ uptake is about 2 s⁻¹.

So a reaction acting directly on the bulk gas, weighted by the local catalyst mass, is the same model
the authors fitted.

### Option A (recommended): catalyst-weighted reaction on the bulk concentration

- New phase `"BulkConcentration"`. Activities are the bulk concentrations c_j, in mol m⁻³; the rate
  laws (`GeneralPowerLaw` and the others) are unchanged.
- The rate r is per kg of catalyst. The source in the bulk balance, per fluid volume, is
  $S_j = \nu_j\,\rho_{\mathrm{cat}}(z)\,r/\varepsilon(z)$, where ρ_cat is the catalyst mass per bed volume (gap 2).
- Code:
  - `computeReactionDerivatives` (`src/integrators/compute.cpp:182`): a new case writes into a new array
    `reactionBulkSource`, of size nodes × components. It applies an inventory limiter like
    `limitedPoreRate`, on c_j.
  - `computeBulkSpeciesSink` (`src/integrators/transport.cpp`, and `compute_multibed.cpp:629`): subtract
    `reactionBulkSource`. It must be included before `updateVelocityAndPressure`, so that non-equimolar
    reactions change the velocity. The WGS is equimolar.
  - The RHS is assembled in several places, each for the single and the multibed column:
    - `cvode.cpp:372` and `:500`;
    - `rk3.cpp:143` and `:197` (`computeDerivatives`), plus `precompute` at `:52` and `:95`;
    - `swing_adsorption.cpp:24` (`refreshColumnState`);
    - the initialisation in `column.cpp:445` and `column_multibed.cpp:669`.

    `rk3_si.cpp` has no reaction or transport terms at all. If phase 2 of
    `docs/SIRK3_ANALYTIC_JACOBIAN_PLAN.md` lands first, this becomes one place.
  - Heat: `reactionHeat` is per kg of solid. Convert with ρ_cat/(ρ_p(1−ε)). It is irrelevant for Stadler
    (isothermal).
- Input:

  ```json
  {"Phase": "BulkConcentration", "Style": "GeneralPowerLaw",
   "Reactants": ["CO", "H2O"], "Products": ["CO2", "H2"], "ForwardOrders": [1, 1], "BackwardOrders": [1, 1],
   "Kinetics": {"forwardRateCoefficient": k0, "forwardActivationEnergy": 51845.0,
                "equilibriumConstant": 0.0131675, "gibbsFreeEnergy": -38061.9}}
  ```

  - With concentrations as activities, k0 = 2.01e5 · (RT/1e5)² · 1e3/3600, in m⁶ mol⁻¹ kg⁻¹ s⁻¹; the
    catalyst mass then enters through ρ_cat.
  - The (RT)² factor makes k0 specific to the temperature. That is fine for isothermal runs.
  - Optional: an `"ActivityBasis": "PartialPressure"` key (in Pa), so that pressure-based rate laws need
    no conversion.
- Paper: the SI says "Homogeneous reactions in the bulk gas or liquid are not included." It needs a short
  paragraph on catalyst-weighted pseudo-homogeneous coupling, as the third coupling next to "with" and
  "without adsorption limitation".
- Fig. 4 then runs on a 25 mm column (the catalyst quarter) without sorbent, at ρ_cat = 100 kg m⁻³ and
  FixedVelocity v = 0.101 m s⁻¹. The PFR integration above is its exact verification counterpart: a
  Danckwerts or dispersion-free limit is not needed at Pe ≈ 1000.

### Option B: a particle pore phase for every gas (the SI model)

- Replace the per-site surface and pore blocks with c_s and c_p per adsorbent and component, for every
  component, when enabled.
- Input:
  - `ParticlePorosity` ε_p;
  - film and pore coefficients per component, moved out of the chemisorption `RateEquationParameters`.
- Chemisorption and `PoreConcentration` reactions then act on this one pore phase. Its volume is defined,
  so a catalyst density per pore volume is meaningful.
- Touches:
  - both state layouts and `transport.cpp` (single column and multibed);
  - `computeChemisorption`, reactions, `writeOutput`, and `writeJSON`/`readJSON`;
  - the input reader;
  - the examples `Chemisorption-CO2/breakthrough` and `Reactions/pore-langmuir-hinshelwood`;
  - every RHS assembly listed under option A, and the SI.
- Recommended only when a validation case actually resolves pore diffusion. The TODO search for one is
  open: Rosen 1952, the CADET general rate model, Kalyanaraman 2015. Stadler does not need it.

## 2. Catalyst amount and placement (gap 2)

- Reactions are global. In the multibed column, `computeReactionDerivatives` receives
  `physisorptionMixtures.front().components` and acts at every node, with no section weighting
  (`cvode.cpp:561`, `rk3.cpp:217`).
- Add `CatalystDensity`, in kg of catalyst per m³ of bed, default 0:
  - per `Adsorbents` entry (`inputreader.cpp:1797` key list), giving
    ρ_cat(z) = Σ_a χ_a(z) ρ_cat,a with the existing `fractionOfAdsorbent`. That handles interfaces
    automatically.
  - at top level for a single column.
- Only `BulkConcentration` reactions use it. The meaning of adsorbed-phase reactions is unchanged.
- Figs. 7 and 8B:
  - Two sections, 25 mm and 75 mm long. The first carries the catalyst: 0.51 g, i.e.
    ρ_cat = 102 kg m⁻³.
  - Its sorbent share of the 9.35 g is not stated. Assume uniform packing: 2.34 g in the first quarter
    and 7.01 g downstream.

## 3. Coverage-dependent desorption (gap 3)

Stadler, Eqs. 6–10: dq/dt = k (q* − q), with

- k = k_ads for adsorption;
- k = k¹_des exp[β_des (q/q_max)/RT] for desorption.

Parameters: k_ads = 0.1125 s⁻¹, k¹_des = 5.38e-4 s⁻¹, β_des = 74,000 J mol⁻¹, and
q_max = q*(20 bar) = 0.3421 mol kg⁻¹.

- Put it in `Chemisorption::rate`, in the `FirstOrder` case (`src/chemisorption.cpp:23`). The rate
  already receives the solid temperature. One change there serves the single column, the multibed column,
  and every integrator.
- New optional `RateEquationParameters` keys for `FirstOrder`, parsed at `inputreader.cpp:367`, where
  `requireOnlyExactKeys` becomes known-keys with `MassTransferCoefficient` required:
  - `DesorptionMassTransferCoefficient` (k¹_des, 1/s);
  - `DesorptionCoverageEnergy` (β_des, J/mol);
  - `CoverageReferenceLoading` (q_max, mol/kg).

  If they are absent, the law stays symmetric (current behaviour).
- Do not reuse the site's `MaximumLoading` for q_max: it caps q* (`chemisorption.cpp:15`).
- Stadler input:
  - CO₂ carries a `ChemisorptionSites` entry with Freundlich [0.0140596, 4.54545] and `FirstOrder`
    {0.1125, 5.38e-4, 74000, 0.3421};
  - the `PhysisorptionSites` entry becomes the zero-Henry placeholder.
  - Check that the chemisorption mixture accepts a Freundlich site.
- The adsorption-only runs can stay on physisorption, where k_ads is the LDF coefficient. A one-step
  regression test should show that the two give identical outlets.
- Numerics: Freundlich has an infinite slope as p → 0, which is everywhere during regeneration. Expect
  stiffness; use SPGMR (gap 9).

## 4. Swing phases with their own feed and flow (gap 4)

- Phase keys are `Name`, `Temperature`, `InletPressure`, and `NumberOfTimeSteps` (`inputreader.cpp:1347`).
- `applyStageConditions` (`swing_adsorption.cpp`) does three things:
  - sets T and P_in;
  - rewrites the inlet node from each component's single `initialGasMoleFraction`;
  - calls `refreshColumnState`.
- Add per phase:
  - `GasPhaseMolFractions`: an object of component name to y, carrier included, summing to 1. Default:
    the `Components` values.
  - `ColumnEntranceVelocity`: sets `column.columnEntranceVelocity` and `timeNormalizationFactor`.
  - Optional: `AxialDispersionCoefficients` per component (§6).
- These go into `SwingAdsorptionSubStage` (`swing_adsorption.h:17`). `applyStageConditions` writes the
  inlet node from the stage fractions; the inlet node is Dirichlet (`dc/dt = 0` at node 0,
  `compute.cpp:652`), so that is all a feed change needs.
- CVODE:
  - `CVODE::reinitialize()` only calls `CVodeReInit` (`cvode.cpp:273`).
  - The step cap Δz/v and the absolute tolerances are set once, at setup (`cvode.cpp:205-215`).
  - Regeneration flows 4× faster (v = 0.404 against 0.101 m s⁻¹), so the cap must be recomputed at each
    switch, and the tolerances refreshed for the new feed.
  - Give `reinitialize` the column and redo both.
- Pressure switch: no change needed for Stadler.
  - With `FixedVelocity`, the pressure profile is imposed from `InletPressure` (`compute.cpp:596`), and
    partial pressures are the renormalised mole fractions times that pressure (`compute.cpp:437`).
    Equilibria therefore switch to 1 bar instantly, as Stadler assumes.
  - The 8 bar gas inventory is flushed out within one residence time (0.25 s). That affects molar flows,
    not the dry fractions that are compared.
  - Optional: rescale interior concentrations by P_new/P_old at the switch, so that c and p stay
    consistent.
- Flow direction: Stadler does not state it; assume co-current. Counter-current regeneration needs the
  flow-reversal feature (TODO, "Flow reversal").
- One Stadler cycle:
  - adsorption: 900 s, 8 bar, feed of Fig. 6 or 7, v = 0.101 m s⁻¹;
  - regeneration: 2400 s, 1 bar, H₂O/N₂ = 0.4/0.6, v = 0.404 m s⁻¹.

  Five cycles span 16,500 s; cycle 3 is compared with Fig. 8.

## 5. The existing pore model loses mass

In `transport.cpp`:

- The bulk loses `fluidSolidPerFluidVolume · k_f (c − c_s)`.
- The surface gains `solidFluidPerSolidVolume · k_f (c − c_s)`.
- The pores exchange `P = k_p (c_s − c_p)` with the surface.

These three conserve mass only if both reservoirs have the solid volume (1−ε) per bed volume. The pore
then loses γ·q̇ with γ = `loadingPrefactor` = ρ_p(1−ε)/ε, while the solid gains ρ_p(1−ε)q̇ per bed volume.
The mismatch is a factor (1−ε)/ε: for ε = 0.38, 1.63 times more leaves the pores than arrives on the solid.

Numerical check: a copy of `examples/Chemisorption-CO2/breakthrough` with FixedVelocity, one General site,
and no physisorption. On the discrete measure that RUPTURA conserves (nodes 1..N, weight Δz):

| run | closure, ρ_p(1−ε)/ε on the solid | closure, ρ_p((1−ε)/ε)² on the solid |
|---|---|---|
| slow flow, fast uptake (t = 59.9 s) | 0.750 | 1.0001 |
| example settings (t = 299 s) | 0.915 | 1.0012 |

The model is therefore consistent only with the wrong solid weight: a quarter of the fed CO₂ disappears in
the first run.

- Minimal fix: γ = ρ_p, consistent with reservoirs per solid volume. This applies in
  `computeChemisorptionTransportDerivatives`, its multibed version, and the inventory limiter in
  `computeChemisorption` (`sorption.cpp:175`).
- Proper fix: option B, with explicit ε_p, β_p, and γ, so that code and SI agree.
- Either way, the SI's β_s/β_p equations must also conserve mass. As written, the film flux in the bulk
  and in the surface balance use different volume bases.
- A related inconsistency: for pore-transport sites, q* is evaluated from the bulk mole fractions, while
  the rate uses c_p (`computeChemisorptionEquilibriumLoadings` receives the bulk `moleFraction`).
- Any balance taken with trapezoid weights also counts the inlet node. Its solid and reservoirs take up
  material from the pinned feed concentration, which is an O(Δz) artefact of the Dirichlet inlet and
  applies to every sink.

## 6. Not blocking, but worth knowing

- Numerical dispersion: convection is first-order upwind (`compute.cpp:659`), which adds v Δz/2 =
  5e-5 m² s⁻¹ on 100 cells. That is 7 times Eq. 12's 7.5e-6 m² s⁻¹, so per-phase physical dispersion
  barely matters. The breakthrough fronts are set by k_ads and by the tubing; refining the grid is the
  lever, if any.
- Solver: SPGMR is 26 times faster than the dense default on the clean-bed step (9.6 s against 248 s),
  with the outlet within 3e-4 in mole fraction. Use it in `stadler/generate.py` for the cycles.
- Cycle 1 capacity loss (irreversible carbonate) is outside the model, as in the paper.

## 7. Order of work, with a check at each step

1. Gap 3, the desorption law.
   - Unit test of `Chemisorption::rate` in both branches.
   - Regression: an adsorption-only step through physisorption and through chemisorption must give
     identical outlets.
2. Gap 4, the swing phases.
   - Stage 3: Figs. 6 and 8A, five sorbent-only cycles.
   - Checks: CO₂ balance over each cycle; cycle-to-cycle convergence (Stadler reports a periodic state
     from cycle 4 in the pilot plant); cycle 3 against the measured curve and against the authors' curve
     in Fig. 8A, both through the measured CO response.
3. Gaps 1A and 2: Stage 2, Fig. 4.
   - The steady outlet X_CO against the PFR integral, which must agree to discretisation error, and
     against the digitised model lines (1e-4).
   - The equilibrium limit at long residence time against Moe.
4. Stage 4: Figs. 7 and 8B, with the multibed column (catalyst quarter plus sorbent), using gaps 1–4.
5. Gap 5, independent of Stadler: fix γ, then decide on option B.

## 8. Decisions needed

1. Reaction coupling: option A (catalyst-weighted bulk reaction; recommended) or option B (pore phase for
   every gas).
2. Where the catalyst amount lives: `CatalystDensity` on `Adsorbents` entries and at top level
   (recommended), or a z-window on each reaction.
3. Desorption law: in chemisorption `FirstOrder` (recommended; it is what TODO.md planned) or in the
   physisorption LDF.
4. Pore mass balance: the minimal γ = ρ_p fix now, or straight to option B.

## Status after implementation (2026-10-21)

### What was implemented

- **Gap 5, pore mass balance.**
  - Fix: γ = ρ_p (`transport.cpp`). In the multibed column (`compute_multibed.cpp`), γ is the
    solid-volume-weighted ρ_p. Both inventory limiters were changed the same way.
  - Closure on the copy of `examples/Chemisorption-CO2/breakthrough`: 0.750 → 1.0001.
  - Test: `ChemisorptionTransport.FilmPoreAndChemisorptionExchangeConservesMass`.
- **Gap 3, desorption kinetics.** Optional `FirstOrder` keys `DesorptionMassTransferCoefficient`,
  `DesorptionCoverageEnergy`, and `CoverageReferenceLoading`, in `Chemisorption::rate`.
- **Gap 1, bulk-gas reaction.**
  - `Reaction::Phase::BulkConcentration`. Its source per fluid volume (`reactionBulkSource`) enters the bulk
    balance through `computeBulkSpeciesSink`, in every RHS assembly.
  - In the single-column CVODE RHS, the reaction is now evaluated before the final sink, as RK3 and the
    multibed path already did. The result is bit-identical for the existing reaction phases: the
    `GasBreakthroughRegression` values are unchanged.
  - The input reader rejects the phase together with `energyBalance` (no heat coupling yet) and with SIRK3.
- **Gap 2, catalyst placement.** A reaction key `"Adsorbents": [names]` restricts a reaction to multibed
  sections, weighted by their local fraction. The catalyst amount goes into `forwardRateCoefficient`, so no
  `CatalystDensity` is needed.
- **Gap 4, swing phases.**
  - Per-phase `GasPhaseMolFractions` and `ColumnEntranceVelocity`.
  - Top-level `NumberOfCycles`.
  - `CVODE::reinitialize(column)` re-derives the Δz/v step cap at each switch.

All of these have unit tests. The suite passes 141 of 159; the 18 failures are the same tests that fail on
the clean checkout.

### Results (`stadler/stadler_validation.ipynb`)

- **Fig. 4.** RUPTURA is within 0.0007–0.0035 of the plug-flow conversion over all 51 runs.
- **Grid study (275 °C, S/G = 4.5, no physical dispersion).** The conversion error halves with every
  doubling of the grid: 3.5e-3 at 25 cells, 5e-5 at 1600. The upwind convection is first order, as
  expected.
- **Sorbent-only cycles.** Periodic from cycle 2 on: the bed-average loading goes 0.1446 → 0.0177 mol/kg.
- **Sorbent-only, Fig. 8A, cycle 3.** The CO₂ delay behind CO through the measured response is 17.9 s,
  the same as the authors' model; 16.3 s measured.
- **SE-WGS cycles.** The WGS acts only in the catalyst quarter, and the conversion follows the Table 3
  kinetics: X = 0.753, against a plug-flow 0.761. The paper reports 0.896, and its CO₂ front is earlier
  (10.6 s against RUPTURA's 19.8 s). A lower flow cannot explain both, because it would delay the front
  further. The stated SE-WGS conditions do not reproduce Fig. 8B; confirm them with the authors.

### Numerical dispersion

The dispersion is not an unavoidable property of the model. It comes from the first-order upwind
difference of the convective term (`computeMassDerivatives`, `compute.cpp:659`), which adds v Δz/2. The
grid study shows it vanishing as 1/N. It can be reduced in two ways:

1. **Refine the grid.** The cost grows with N. With CVODE's automatic cap Δz/v, the step count also grows
   with N, because the step is capped by the smallest cell.
2. **Use a higher-order convection scheme.** RUPTURA already has a WENO reconstruction for the energy
   balance (`computeWENO`), and a TVD/WENO flux for the species balance would give second order or better.
   This is a design change; not done.

For Stadler it does not matter: the fronts are set by k_ads and by the tubing response, which is about
30 s wide. Only the catalyst quarter needs a fine grid, because the conversion responds to back-mixing.

### Are the paper's pore-transport equations right?

No, not as written. The SI equations `surface_concentration_balance`, `pore_concentration_balance`, and
`reaction_pore_reservoir` mix volume bases:

- β_s, β_p, and Γ are defined per moving-fluid volume: β_p = (1−ε_b)ε_p/ε_b and
  Γ = γ ρ_s r with γ = (1−ε_b)(1−ε_p)/ε_b.
- F = a_loc k_f (c − c_s), with a_loc = 6/d_p, is a flux per particle volume.
- P = k_p (c_s − c_p) is an LDF rate per pore volume.

Per unit total volume, the bulk loses a^V_fs k_f Δc = (1−ε_b) a_loc k_f Δc (Eq. `film_controlled_species_sink`),
while the surface reservoir gains ε_b β_s ∂c_s/∂t = ε_b a_loc k_f Δc − … . The two agree only at ε_b = 0.5.
The pore equation conserves mass between surface and pores, but ∂c_p/∂t = P/β_p relaxes with k_p/β_p
instead of k_p.

A consistent set, written per unit moving-fluid volume as the β's are:

$$\mathcal{F} = \frac{1-\varepsilon_b}{\varepsilon_b}\,a_{\mathrm{loc}}k_f(c-c_s),\qquad
\mathcal{P} = \beta_p\,k_p(c_s-c_p),$$
$$\beta_s\frac{\partial c_s}{\partial t} = \mathcal{F}-\mathcal{P},\qquad
\beta_p\frac{\partial c_p}{\partial t} = \mathcal{P}-\Gamma+\beta_p\mathcal{R}^p,\qquad
\Gamma = \frac{1-\varepsilon_b}{\varepsilon_b}\,\rho_p\,r^{\mathrm{chem}},$$

The bulk sink is S = 𝓕 per fluid volume, which is (1−ε_b) a_loc k_f (c − c_s) per total volume. Dividing by
the β's:

$$\frac{\partial c_s}{\partial t} = \frac{1-\varepsilon_b}{\varepsilon_b\beta_s}a_{\mathrm{loc}}k_f(c-c_s) - \frac{\beta_p}{\beta_s}k_p(c_s-c_p),\qquad
\frac{\partial c_p}{\partial t} = k_p(c_s-c_p) - \frac{\rho_p}{\varepsilon_p}r^{\mathrm{chem}} + \mathcal{R}^p.$$

The code after the γ fix is exactly this set with β_s = β_p = (1−ε_b)/ε_b, i.e. both reservoirs hold the
solid volume and ε_p = 1:

- ∂c_s/∂t = a_loc k_f (c − c_s) − k_p (c_s − c_p);
- ∂c_p/∂t = k_p (c_s − c_p) − ρ_p r + R^p.

Supporting ε_p < 1 needs a `ParticlePorosity` input, entering as ε_p in the surface exchange and as ρ_p/ε_p
in the pore sink.
