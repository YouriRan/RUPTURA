# 3.3.4 Humid direct air capture: Liu et al. (2024)

Reference: H. Liu, H. Lin, S. Dai, D.-e. Jiang, Minimal Kinetic Model of Direct Air Capture of CO₂ by
Supported Amine Sorbents in Dry and Humid Conditions, Ind. Eng. Chem. Res. 63 (2024) 5871–5879,
doi:10.1021/acs.iecr.3c04535.
Source: `liu/` (57 simulations, `generate.py`, `run.py`, `liu_paper_figures.ipynb`, `conversions.md`,
`reference_digitized/`).
Status: runs done; figures exist for 1a, 1b, 6a, 6c, and 7.

## Points

- [ ] Paper figures: 1b (dry benchmark), 4 (humid CO₂/H₂O), 6 (dry DAC at 308–328 K), and 7 (dry vs humid).
      Figure 5 overlaps with 3.2.3.
- [ ] Figures 2–3 (parameter sweeps) belong to 3.4.3 (`../3.4.3_sensitivity`).
- [ ] The humid Toth law is frozen at the GAB loading of the feed humidity, because humidity-dependent CO₂
      isotherms are not implemented (`../../TODO.md`). Disclose this, or rerun once they are implemented.
- [x] Figure 4 is less converged (loose tolerances were needed for tractable run times; see
      `liu/README.md`). Report the worst backward step.
      2026-11-03: moot. ROS34PW2 at its default tolerances runs Figure 4 in under a second, and the worst
      backward step is about 1e-16.
- [ ] The bed densities for Figs. 1, 6, and 7 are inferred from the printed time axes. Disclose this.
- [x] State that the absolute CVODE concentration tolerance is lowered to 1e-6 mol m⁻³ for these dilute
      feeds. The default gives visible non-monotonic steps.
      2026-11-03: moot. The runs use ROS34PW2 at its default tolerances.
