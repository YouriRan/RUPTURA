# 3.2.3 Water isotherm (GAB)

Source: `../3.3.3_chemisorption_elfving/elfving/equilibrium/` and `../3.3.4_humid_dac_liu/liu/`.

## Points

- [ ] Elfving Fig. 2: GAB water isotherm (`equilibrium_reproduction.svg`; redo in the paper style).
- [ ] Liu Fig. 5a: water isotherm at 308 K. MixturePrediction cannot evaluate a GAB site, so the points
      come from a short breakthrough run (`../3.3.4_humid_dac_liu/liu/README.md`).
- [ ] State in the paper that the GAB saturation pressure is a fixed Antoine equation for water, valid
      from 274 to 373 K (`src/isotherm.h:222`).
