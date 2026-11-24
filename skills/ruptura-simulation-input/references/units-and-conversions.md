# Units and conversions

Ruptura uses SI units throughout:

| Quantity | Unit |
|---|---|
| temperature | K |
| pressure | Pa |
| length | m |
| time | s |
| velocity | m/s |
| viscosity | Pa s |
| density | kg/m3 |
| molar mass | kg/mol |
| loading | mol/kg |
| gas concentration | not an input |
| liquid concentration | mol/m3 |
| energy | J/mol |
| LDF and rate constants | 1/s (unless noted) |
| dispersion and diffusivity | m2/s |
| heat-transfer coefficient | W/(m2 K) |
| thermal conductivity | W/(m K) |
| heat capacity | J/(kg K) |

Pressures are absolute. Ask when a source pressure may be gauge (barg, psig, an overpressure).

## Common factors

```text
bar x 1e5 = Pa;  kPa x 1e3;  MPa x 1e6;  mbar x 100;  atm x 101325;  Torr = mmHg x 133.322368
degC + 273.15 = K;  cm / 100 = m;  mm / 1000 = m;  min x 60 = s;  h x 3600 = s
g/mol / 1000 = kg/mol;  kJ/mol x 1000 = J/mol;  g/cm3 x 1000 = kg/m3;  cP x 1e-3 = Pa s
mmol/g = mol/kg (same number);  mg/g / M[g/mol] = mol/kg
cm3(STP)/g / 22.414 = mol/kg (STP = 273.15 K, 101325 Pa; ask if the source's STP differs)
mmol/L = mol/m3 (same number);  mol/L x 1000 = mol/m3;  mg/L / M[g/mol] = mol/m3
ppm (mole) x 1e-6 = mole fraction;  1/min / 60 = 1/s
```

## Isotherm parameters

Convert each parameter by substituting the source variable into Ruptura's equation
([isotherms.md](isotherms.md)). Do not just divide every constant by 1e5. For a gas,
x[Pa] = 1e5 x[bar]:

- a parameter that multiplies x (b, K_H, K_R, a, k1, k2, ...): value[1/Pa] = value[1/bar] / 1e5
- b in Langmuir-Freundlich multiplies x^nu: b[Pa^-nu] = b[bar^-nu] x (1e5)^(-nu)
- K_F in Freundlich multiplies x^(1/n): K_F[Pa] = K_F[bar] x (1e5)^(-1/n)
- a_R in Redlich-Peterson multiplies x^g: a_R[Pa] = a_R[bar] x (1e5)^(-g)
- c in Quadratic multiplies x^2: c[1/Pa^2] = c[1/bar^2] / 1e10
- loadings convert as loadings; saturation capacities in mmol/g keep their number

For a liquid, x is concentration in mol/m3. Examples: b[m3/mol] = b[L/mol] / 1000, and
b[m3/mol] = b[L/mg] x M[g/mol].

**Exponent conventions.** Ruptura's Sips and Freundlich parameters are the INVERSE of the exponent
([isotherms.md](isotherms.md)). If a source writes Sips as q = q_sat (b p)^n / (1 + (b p)^n), enter
1/n. Langmuir-Freundlich and Toth take the exponent itself. Always rewrite the source equation in
Ruptura's form before mapping its parameters. Never map parameters by model name alone.

**Heats of adsorption.** Enter a POSITIVE magnitude. If the source gives Delta H = -25 kJ/mol, enter
25000.0.

## Flow rates

`ColumnEntranceVelocity` is the INTERSTITIAL inlet velocity v, not the superficial velocity u.

- Packed bed: v = Q / (epsilon_b A), where A = pi d_in^2 / 4 is the empty-tube area, epsilon_b is
  the bed (interparticle) void fraction, and Q is the volumetric flow at the inlet pressure and
  temperature. A superficial velocity converts as v = u / epsilon_b.
- Monolith: v = Q / (N_ch A_ch), the mean channel velocity (A_ch from [geometry.md](geometry.md)).
- Flow at standard conditions: Q = Q_std (T / T_std) (P_std / P). A molar flow converts as
  Q = n_dot R T / P, and a mass flow as Q = m_dot / rho with rho = P M_mix / (R T), where
  R = 8.314462618 J/(mol K).

Ask for the column diameter, the void fraction, the flow reference state, and whether a reported
velocity is superficial or interstitial whenever they are not stated.

## Compositions

`GasPhaseMolFraction` is the FEED mole fraction, and the carrier is included in the sum of 1.
Partial pressures convert as y_j = p_j / P_total. If the source fractions do not sum to 1, normalize
them yourself and report it.

## Densities and porosities

- `ParticleDensity` is the particle (pellet) density, i.e. mass per particle volume including the
  pores. A bulk (packing) density converts as rho_p = rho_bulk / (1 - epsilon_b). Crystal and
  skeletal densities are different quantities, so ask if it is unclear.
- `ColumnVoidFraction` is the interparticle bed void fraction epsilon_b, not the total porosity.

## Transport

- `MassTransferCoefficient` is the linear-driving-force constant k_LDF in 1/s. If only an effective
  diffusivity and particle radius are given, Glueckauf's k_LDF = 15 D_e / r_p^2 is a common
  estimate. Treat it as an estimate (workflow step 5 of SKILL.md: only when the user asks).
- `AxialDispersionCoefficient` D_ax is in m2/s. A column Peclet number converts as D_ax = v L / Pe.
