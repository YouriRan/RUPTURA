# 3.3.5 Heterogeneous reactions: Stadler et al. (2022), sorption-enhanced water-gas shift

Reference: T. J. Stadler, L. J. Bender, P. Pfeifer, Dynamic simulation of a compact sorption-enhanced
water-gas shift reactor, Front. Chem. Eng. 4 (2022) 1000064, doi:10.3389/fceng.2022.1000064.
PDF: `stadler/fceng-04-1000064.pdf`.
Related: Stadler et al., Energies 14 (2021) 355 (laboratory setup); Stadler et al., Processes 10 (2022)
1160 (numerical approach).
Source: `stadler/`, with `stadler_validation.ipynb`, `generate.py`, `run.py`, `digitize_reference.py`,
`reference_digitized/`, `conversions.md`, and `figures/`.
Status (2026-10-21): every experiment of the paper runs.
In the paper (2026-11-02): Results "Chemisorption & pore transport" (Figs. 8A and 6, `fig:stadler_sorption`) and
"Reactions" (Fig. 4 and the grid study, `fig:stadler_reactions`); the SE-WGS cycles (Figs. 7, 8B), the isotherm,
and the cyclic loading in SI `apx:stadler` (`fig:stadler_si`, `tab:stadler_parameters`). The notebook writes the
panels at print size into `ruptura-2.0-paper/figures/stadler/`. This case does not test pore transport (V-PHY-2).
- Figs. 4–8 are digitized; Fig. 4 has 51 steady-state runs plus a grid study, and Fig. 5 runs.
- Figs. 6, 8A: five sorbent-only cycles. Figs. 7, 8B: five SE-WGS cycles.
- The features this needed are listed below and in `../../docs/STADLER_FEATURES_PLAN.md`.

## Reference model (Table 1, Eqs. 4–13)

- Isothermal at 523 K. Pressure drop is negligible and the superficial velocity constant. Dirichlet
  inlet, zero-gradient outlet.
- Bulk gas:
  - axial dispersion $D_{\mathrm{ax},i} = D_{\mathrm{mix},i}\sqrt{\varepsilon_b} + u\,r_p$, with $u$ superficial;
  - exchange with the particle pores by LDF, $k_{\mathrm{LDF},i} = 15 D_{\mathrm{eff},i}/r_p^2$;
  - no film resistance.
- Particle pores ($\varepsilon_p = 0.5$):
  - WGS source $w_{\mathrm{cat}}\,\nu_i\,r_{\mathrm{WGS}}$ and CO₂ sorption sink $(1-w_{\mathrm{cat}})\,a_{\mathrm{CO_2}}$;
  - catalyst and sorbent particles are lumped into one phase with catalyst mass fraction $w_{\mathrm{cat}}$.
- WGS rate: $r = k_\infty e^{-E_a/RT}\,(p_{\mathrm{CO}}\,p_{\mathrm{H_2O}} - p_{\mathrm{CO_2}}\,p_{\mathrm{H_2}}/K_{\mathrm{eq}})$.
  - $k_\infty = 2.01\times10^{5}$ mol bar⁻² g⁻¹ h⁻¹, $E_a$ = 51.845 kJ mol⁻¹.
  - $K_{\mathrm{eq}} = \exp(4577.8/T - 4.33)$ (Moe).
- CO₂ isotherm (Freundlich, 523 K, dry): $q_{\mathrm{eq}} = k_{\mathrm{Fr}}\,p^{1/n}$, with
  $k_{\mathrm{Fr}}$ = 0.177 mmol g⁻¹ bar⁻¹ᐟⁿ and $1/n$ = 0.220.
- CO₂ sorption kinetics:
  - adsorbing: $a = k_{\mathrm{ads}}(q_{\mathrm{eq}} - q)$ with $k_{\mathrm{ads}}$ = 0.1125 s⁻¹;
  - desorbing: $a = k_{\mathrm{des}}(q_{\mathrm{eq}} - q)$ with $k_{\mathrm{des}} = k_{1,\mathrm{des}}\exp[\alpha_{\mathrm{des}}(q/q_{\max})/RT]$,
    $k_{1,\mathrm{des}} = 5.38\times10^{-4}$ s⁻¹, $\alpha_{\mathrm{des}}$ = 74,000 J mol⁻¹, $q_{\max} = q_{\mathrm{eq}}$(20 bar).
- Only CO₂ adsorbs; the interaction of H₂O with the sorbent is neglected.

## Experiments usable for validation

Laboratory reactor with two slits, each 2 mm × 50 mm × 100 mm.

| Fig. | Content | Conditions |
|---|---|---|
| 4 | Steady-state CO conversion vs steam-to-gas ratio (0.5–4.5), with equilibrium lines | 0.50 g catalyst diluted in SiC in the first quarter; 8 bar; 225, 250, 275 °C; 2000 mL min⁻¹; 10% CO in N₂ |
| 5 | CO₂ isotherm | 523 K, dry |
| 6, 8A | Cyclic breakthrough on sorbent only; CO acts as a non-adsorbing tracer | 9.70 g sorbent. Adsorption: 15 min at 8 bar, CO/CO₂/H₂/H₂O/N₂ = 0.05/0.05/0.05/0.2/0.65. Regeneration: 40 min at 1 bar, H₂O/N₂ = 0.4/0.6 |
| 7, 8B | Cyclic SE-WGS breakthrough | 0.51 g catalyst mixed with sorbent in the first quarter, pure sorbent downstream (9.35 g sorbent). Adsorption feed CO/H₂O/N₂ = 0.05/0.2/0.75; regeneration as above |

- Parameters were fitted from cycle 2 on, because cycle 1 loses capacity to irreversible carbonate
  formation. Fig. 8 compares the third cycle.
- Outlet data are dry-basis volume fractions, and they include strong back-mixing in the downstream
  tubing.
- This paper does not give the flow rates of the breakthrough experiments. Stadler et al. (2021),
  Table 2, uses 2000 mL min⁻¹ (STP) for adsorption and 1000 mL min⁻¹ for regeneration, with the same
  15 min / 8 bar and 40 min / 1 bar / 40 % H₂O protocol. These values are assumed here.
- The measured and simulated plateaus in Figs. 6 and 8 are about 0.066. For the stated feed, the dry
  value is 0.05/0.8 = 0.0625.
- Figs. 9–12 (six-chamber pilot plant with threshold-triggered switching) are simulations only and are
  out of scope.

## How the reference model maps onto RUPTURA

- Freundlich isotherm, isothermal operation, and `FixedVelocity` with a Dirichlet inlet, as in the
  reference.
- WGS (Eq. 4): a `BulkConcentration` reaction with `GeneralPowerLaw`.
  - Forward orders 1/1 on CO and H₂O; backward orders 1/1 on CO₂ and H₂; $k_b = k_f/K$.
  - The catalyst amount per fluid volume goes into `forwardRateCoefficient` (`conversions.md` §9).
  - The authors' Fig. 4 model is an ideal plug-flow reactor without pores; it is reproduced to 1e-4.
  - In SE-WGS, `"Adsorbents": ["catalyst-zone"]` restricts the reaction to the first quarter of a
    two-section multibed column.
- Moe's $K_{\mathrm{eq}}$ maps exactly: a VantHoff `EquilibriumConstant` with `Value` exp(−4.33) = 0.01317 and
  `HeatOfReaction` = 4577.8·R = 38,062 J mol⁻¹.
- CO₂ sorption (Eqs. 6–10): a Freundlich chemisorption site with the `FirstOrder` law and its desorption
  branch (`DesorptionMassTransferCoefficient`, `DesorptionCoverageEnergy`, `CoverageReferenceLoading`).
  The pore LDF of Eq. 13 is about 10⁴ times faster than k_ads, so the pores follow the bulk.
- Cycles: `SwingAdsorption` with per-phase `GasPhaseMolFractions`, `ColumnEntranceVelocity`, and
  `InletPressure`, repeated with `NumberOfCycles`.
- Not modelled:
  - per-phase D_ax (RUPTURA holds one D_ax per component);
  - the outlet dead volume. The notebook passes the simulated dry outlet through the measured CO
    response instead, with no fitted model.

## Staged plan

- [x] Stage 1: Fig. 5 isotherm with MixturePrediction (`figure5_isotherm_CO2_523K`); matches Table 3 to
      machine precision.
- [x] Stage 1b: the sorbent-only adsorption step on a clean bed (`figure6_sorbent_only_fresh_bed_8bar_523K`).
  - Compared with cycle 1 of Fig. 6 and cycle 3 of Fig. 8A; the mass balance closes to 1e-4.
  - The CO₂ delay behind CO is 15.5 s at the bed outlet and 20.3 s through the measured tracer response.
    The measurement gives 62 s (cycle 1, with carbonate formation), then 18.8, 15.8, 12.4, and 12.3 s
    (cycles 2–5). The authors' model gives 17.9 s (cycle 3).
- [x] Stage 2: Fig. 4 steady-state conversion (`figure4_wgs_conversion/`, 51 runs over the catalyst
      section only).
  - RUPTURA is within 0.0007–0.0035 of the plug-flow conversion, so it lies on the authors' model lines.
  - Grid study (`figure4_grid_study_548K_SG4p5/`, physical dispersion off): the error falls as 1/N, from
    3.5e-3 at 25 cells to 5e-5 at 1600. This is the first-order numerical dispersion of the upwind scheme.
  - 2026-11-02, WENO3 (default since 2026-10-27): the error falls with the square of the grid spacing, from
    1.8e-4 at 25 cells to 3.8e-8 at 1600; the 51 runs are within 0.0005–0.0028 of plug flow.
- [x] Stage 3: Figs. 6 and 8A, five sorbent-only cycles (`figure6_8a_sorbent_only_cycles_8bar-1bar_523K`).
  - Periodic from cycle 2 on.
  - CO₂ delay behind CO through the measured tracer: 19.7 s in cycle 1, then 17.7, 18.0, 17.7, and
    17.2 s. The measurement gives 62.1 s (cycle 1, carbonate formation), then 18.8, 15.8, 12.4, and
    12.3 s.
  - Fig. 8A, cycle 3: RUPTURA 17.9 s, the authors' model 17.9 s, the measurement 16.3 s.
  - 2026-11-02 rerun (WENO3, half-height times of the digitized curves): cycles 1–5 give 20.2, 17.7, 18.0, 17.7,
    and 17.2 s; measured 62.1, 18.8, 15.8, 12.4, and 12.3 s; the authors' model 17.9 s in cycle 3.
- [x] Stage 4: Figs. 7 and 8B, five SE-WGS cycles (`figure7_8b_sewgs_cycles_8bar-1bar_523K`). See
      "SE-WGS result" below.
- [x] `conversions.md`: flow at STP to velocity, bed and sorbent mass, $k_{\mathrm{Fr}}$, Eq. 12 with Fuller
      diffusivities, dry basis, and the WGS rate to a `BulkConcentration` rate.
- [x] Digitize Figs. 4, 5, 6, 7, and 8 into `reference_digitized/`.
  - Read from the vector paths in the PDF, not from pixels.
  - The Fig. 5 Freundlich fit matches Table 3 to 0.02 %; the Fig. 4 equilibrium lines match Moe to 1e-4.
- [x] `generate.py` and `run.py` for the runnable cases; extend them as stages 2–4 become possible.
- [ ] Verification counterpart: Wehner–Wilhelm conversion (3.1).

## SE-WGS result (Figs. 7, 8B)

- The WGS acts only in the catalyst quarter, which has 40 cells; downstream, H₂ is constant. The run is
  periodic from cycle 2 on.
- With the stated conditions (2000 mL min⁻¹, 0.51 g catalyst, Table 3 kinetics), the steady CO conversion
  is 0.753 (carbon basis; 0.759 in the 2026-11-02 rerun). A plug-flow reactor with the same kinetics gives 0.761.
- The paper reports 0.896, and the authors' Fig. 8B model plateaus agree with that.
- The measured CO₂ front is 10.6 s behind the Fig. 8A CO response; the authors' model gives 10.7 s and
  RUPTURA 19.8 s (19.7 s in the 2026-11-02 rerun). The measured front is even earlier than in the sorbent-only cycle 3 (16.3 s), although
  less CO₂ is produced.
- The measured unconverted CO rises after the CO₂ and slowly: t₅₀ = 59 s and t₉₀ = 118 s, against 42 s and
  90 s for the inert CO of Fig. 8A. The conversion is therefore higher at the start of the step and falls as
  the sorbent saturates.
  - Eq. 4 cannot produce this: far from equilibrium (X_eq = 0.996), CO₂ removal does not change the rate.
    That needs CO₂-inhibited kinetics.
  - The authors' model shows the slow CO rise anyway, so their Fig. 8B simulation is not Eq. 4 with
    Table 3 as printed.
- The measured CO₂ front moves earlier every cycle (t₅₀ = 66.7, 52.4, 45.3, 37.9, 36.5 s for cycles 1, 3, 13,
  43, 63). The SE-WGS sorbent loses capacity faster than the sorbent-only bed.
- Sensitivity runs (three cycles, 10 cells over the catalyst quarter, in the session scratchpad; not kept):
  - 1400 mL min⁻¹: X = 0.83, but the front moves later (68.6 s). A lower flow lowers the CO₂ production X·ṅ.
  - A 1.7× more active catalyst: X = 0.87, front at 59.9 s.
  - The same catalyst with 40 % less sorbent capacity: front at 52.9 s, against 52.5 s measured.
- Fig. 8B therefore needs a faster reaction (more active catalyst, or a CO order below 1 at 5 % CO, which
  Fig. 4 at fixed 10 % CO barely constrains) and a smaller cyclic capacity. The flow alone cannot explain it.
- Conclusion: the stated SE-WGS conditions and parameters do not reproduce Fig. 8B, and RUPTURA follows
  what is stated. Confirm the conditions of the 2022 SE-WGS runs with the authors before using Figs. 7
  and 8B quantitatively. The sorbent-only Figs. 6 and 8A, and Fig. 4, are consistent.
