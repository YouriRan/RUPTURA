# 3.3.7 Multistacked columns: Park et al. (1998), layered beds

Reference: J.-H. Park, J.-N. Kim, S.-H. Cho, J.-D. Kim, R. T. Yang, Adsorber dynamics and optimal design
of layered beds for multicomponent gas adsorption, Chem. Eng. Sci. 53 (1998) 3951–3963,
doi:10.1016/S0009-2509(98)00196-1. Beds layered with activated carbon followed by zeolite 5A separate H₂
from an H₂/CO₂/CH₄/CO mixture; experiments and simulations.
Backup: J. Yang, C.-H. Lee, Adsorption dynamics of a layered bed PSA for H₂ recovery from coke oven gas,
AIChE J. 44 (1998) 1325–1334, doi:10.1002/aic.690440610. Layered activated carbon/zeolite 5A with
five-component coke-oven gas; breakthrough experiments and two-bed PSA.
Status: blank case directory; PDFs still to obtain.

## Points

- [ ] Obtain the paper. Confirm the adsorbents, the gas system, and which measured data are available
      (breakthrough curves, temperatures).
- [ ] Extract the isotherms, mass-transfer coefficients, bed properties, and layer lengths.
- [ ] Digitize the experiments into `reference_digitized/`.
- [ ] Set up the layered column (`Adsorbents` + `ColumnSections`), with a sharp and a graded interface.
- [ ] `generate.py`, `run.py`, and `conversions.md` as in the other cases.
- [ ] Demonstrations from `../../examples/Multibed/`: layer order, graded interface, and uniform mix vs
      layers.
- [ ] Verification counterpart: a multibed column with one adsorbent equals the single bed (3.1,
      `../../benchmark/multibed/`).
