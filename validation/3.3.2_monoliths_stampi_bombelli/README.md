# 3.3.2 Monoliths: Stampi-Bombelli et al. (2024)

Reference: V. Stampi-Bombelli et al., On Comparing Packed Beds and Monoliths for CO₂ Capture from Air
Through Experiments, Theory, and Modeling, Ind. Eng. Chem. Res. 63 (2024) 11637–11653,
doi:10.1021/acs.iecr.4c01392.
Source: `stampi-bombelli/`, with `stampi_bombelli_validation.ipynb`, `dl_from_literature/` (D_L from
correlations), `dl_from_stampibombelli/` (D_L fitted by the authors), and `figures/`.
Status: in the paper (panels 4b, 4d, 5a, 5b).

## Points

- [x] Rewrite the explanation in the Monoliths paragraph. Every run uses `FixedVelocity`, so no
      pressure-drop model is active. RUPTURA's monolith closure is itself Hagen–Poiseuille when the
      Forchheimer coefficient is zero.
- [x] Compare `dl_from_literature/` with `dl_from_stampibombelli/` to test the claim that the steeper
      fronts come from the velocity-independent dispersion coefficient.
      2026-11-03: with the Danckwerts inlet (ROS34PW2), the fitted D_L gives 61–133% of the measured
      10–90% widths and t50 within −15% to +24%. The literature correlations give 11–76% of the widths, and t50
      within 4% of the fitted runs. The paper now uses the fitted D_L. With the earlier Dirichlet inlet, the
      fitted D_L let dispersion admit several times the feed, so it looked unusable.
- [ ] Add the Table 6 comparison ($t_{50}$, $t_{90}$, $q_{90}$) from the notebook.
- [x] Decide on the unused panels (4a, 4c, 5c–5f); move them to the SI or drop.
      2026-11-02: 4c and 5c in the SI (`apx:stampi_bombelli`), 5d–5f dropped; 4a unused (no simulation of the
      400 ppm packed bed).
- [x] Fill the empty `\cite{}` for Stampi-Bombelli in `main.tex` (two places). (None left, 2026-11-02.)
