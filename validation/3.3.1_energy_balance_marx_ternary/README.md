# 3.3.1 Energy balance: Marx et al. (2014), ternary experiments T1–T24

Reference: D. Marx, L. Joss, N. Casas, J. Schell, M. Mazzotti, Adsorption 20 (2014) 493–510.
Source: `marx2/` (moved here from the repository root on 2026-11-19), with `generate.py` (writes the inputs of the model
of Marx et al., `conversions.md` §5–9), `run.py`, `ternary/T1`–`T24` (inputs and output), `reference_digitized/` (experiments
extracted from the vector markers in the PDFs, `digitize_paper.py`, and the model curves of Marx et al.,
`digitize_model.py`), `create.ipynb` (figures and `t50_comparison.csv`), and `output/pdf/`.
`marx2/reference_settings_hassan_chung/` holds the frozen inputs of the reference settings, which only the C++ unit
tests read; since 2026-11-20 every paper study of T15 uses `marx2/ternary` (TODO V-PHY-12). Do not regenerate them. The paper figures in
`figures/marx/` are identical to the files in `output/pdf/`.
Status: text and figures are in the paper (Results "Energy modeling" and SI "Marx benchmark").

## Points

- [x] SI tables updated to the marx2 inputs, with all 24 runs in the experiment table. S1 (selected
      low-flow runs) was removed, so the complete results are now Figs. S1–S4.
- [x] T6 and T18 used $k_{\mathrm{LDF,CO_2}}$ = 0.07 and 0.15, but Marx et al. (Table 6) give 0.04 and 0.11
      for these mixtures at 298 K. Corrected in the inputs and the SI table; all 24 runs now follow
      Table 6.
- [x] Rerun T6 and T18 and regenerate `marx_T1_T6.pdf` and `marx_T13_T18.pdf` via `create.ipynb`
      (done with the full rerun below, 2026-10-20).
- [x] Gas–wall heat-transfer coefficient (2026-10-20). Every run used 25.6 W m⁻² K⁻¹, the Hassan & Chung value
      for T21 only. All 24 runs now use the correlation of Marx et al. (Eq. 3), h_L = 41.13 Re₀^0.32
      exp(−6d_p/d_i) λ_ax/d_i, with Re₀ from the initial helium density and the superficial feed velocity, and
      λ_ax = 0.0400 W m⁻¹ K⁻¹ recovered from their Fig. 6 (`marx2/heat_transfer.py`,
      `marx2/conversions.md`): 25.4–77.1 W m⁻² K⁻¹. All runs rerun with `marx2/run.py --rerun`.
- [x] `ParticleDensity` = 502 kg m⁻³ is the bulk density in Marx et al.; their particle density is
      850 kg m⁻³. Combined with $v = Q/(A\varepsilon_t)$, the solid capacity per fluid volume is about 10%
      above the reference model (`../3.4.1_boundary_conditions_marx_10bar/conversions.md`). Check this
      and state the mapping.
      2026-10-26: mapping stated in `marx2/conversions.md` §6, $\rho_p = \rho_b\varepsilon_b/((1-\varepsilon_b)
      \varepsilon_t) = 456.6$ kg m⁻³, written into T1–T24 by `marx2/particle_density.py`. Not rerun yet.
- [x] Rerun T1–T24 for the 2026-10-26 fixes (particle density, energy-conserving gas–wall exchange; see
      `../../TODO.md`, "Reruns after the 2026-10-26 fixes") and update the SI tables and figures. Superseded by the
      rerun with the inputs of Marx et al. (2026-11-19, below).
- [x] $t_{50}$ table (simulation/experiment, CO₂), now SI Table `tab:marx_t50_comparison`, with the Energy-modeling
      text quantified accordingly (2026-10-20). With the correlation h_fw and the gas heat capacity per unit
      mass (1014.1 / 1109.2 J kg⁻¹ K⁻¹, `marx2/conversions.md` §5; rerun 2026-10-20):
      - 1 bar (T1, T7, T13, T19): 1.34, 1.46, 1.15, 1.22 (late; h_fw barely changed at 1 bar);
      - 15/80/5 at 5–20 bar: 1.04–1.17 (late);
      - 40/30/30 at 5–20 bar: 0.98–1.06.

      With the molar C_p,g (29.51): 15/80/5 1.04–1.15, 40/30/30 0.97–1.05 (t50 shifts of 0–9 s).
      With the old h_fw = 25.6: 1 bar 1.33, 1.46, 1.14, 1.22; 15/80/5 0.94–1.03; 40/30/30 0.82–0.92.
      Largest temperature at 1.1 m: adsorbent within −2.9 to +2.2 K of the thermocouple, gas 1.1–5.6 K below.
- [x] Optional: rerun the 318 K experiments with an ambient temperature of 318 K (`Temperature`) to
      quantify the effect described in the SI note. Done with the 2026-11-19 rerun: the ambient is the feed temperature. It may explain part of the larger 1 bar deviation at
      318 K (T7 1.46 vs T1 1.33; T19 1.22 vs T13 1.14).
- [ ] Optional: feed-velocity sensitivity for T1, T7, T13, and T19. Marx et al. fitted the feed velocity
      because the flow-controller uncertainty at low mass flow can be as large as the set point.
- [x] Temperatures are plotted at $z = L$, but thermocouple T5 is at 1.1 m. Sample the simulation at 1.1 m.
      2026-11-19: `create.ipynb` interpolates the column profiles at 1.1 m.
- [ ] Binary experiments B1–B36: the figures exist (`figures/marx/marx_B*.pdf`). Move them to the SI or drop.
- Adiabatic operation is not covered: all Marx runs exchange heat with the surroundings. Dantas runs 1–8
  (`../../dantas/`) are adiabatic if the adiabatic mode needs a case.
- [x] 2026-11-19, rerun with the inputs of the model of Marx et al. (TODO V-PHY-10; `marx2/generate.py`,
      `marx2/conversions.md` §5–9): their Table 5 isotherms as SCI (multicomponent Sips) with uncentred van 't Hoff
      laws and the new `Arctan` law for the CO₂ exponent, ΔH and k_i of Tables 5–6, total porosity ε_t as fluid
      fraction with ρ_s = ρ_b/(1 − ε_t), molar `GasHeatCapacity`, h_fs from Wakao–Kaguei per run (scaled by
      (1 − ε_b)/(1 − ε_t)), ambient = feed temperature, `InletPressureInletVelocity` for all runs. The model curves of
      Marx et al. are digitized (`digitize_model.py`) and drawn dotted. All 24 runs take about 27 s.
      - t50 RUPTURA/experiment at 5–20 bar: 0.93–0.99 (15/80/5), 0.88–0.96 (40/30/30), i.e. 1–12 % early;
        Marx model/experiment 0.91–1.02; RUPTURA/Marx model 0.88–1.02;
      - 1 bar (T1, T7, T13, T19): 1.17, 1.25, 1.00, 1.06 (Marx model 1.01, 1.02, 1.00, 0.88);
      - largest temperature at 1.1 m: adsorbent 0.2–5.9 K and gas 0.5–7.1 K below the thermocouple; Marx model
        −0.4 to 4.7 K below; T_s − T_g at the peak ≤ 1.8 K;
      - the H2 front of the Marx model arrives 3–21 s after RUPTURA's at 10–20 bar (fitted velocities, downstream
        piping).
- [ ] Optional: feed-velocity fit per run (as Marx et al.), to separate the velocity from the model differences.
