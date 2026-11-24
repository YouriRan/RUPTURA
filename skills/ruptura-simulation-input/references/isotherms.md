# Isotherms, temperature laws, and mixture methods

Contents: site form | isotherm types | humidity-dependent isotherms | temperature dependence |
mixture methods | humidity inputs | pH inputs

## Site form

`"PhysisorptionSites"` (and the equilibrium isotherm of chemisorption sites) is an array of sites:
`{"Type": "<Type>", "Parameters": [values in the order below]}`. Sites are summed. Each value is a
number, or, where marked [T] below, a temperature law (see
[Temperature dependence](#temperature-dependence)).

Physisorption sites may also contain:

- `"RateEquation": "FirstOrder"` (the only physisorption option), and
  `"RateEquationParameters": {"MassTransferCoefficient": k_LDF}`,
- `"HeatOfAdsorption": value`.

Both overrides must be identical for every site of a component, so set them on the component
instead.

## Isotherm types

Notation: q in mol/kg; x = partial pressure [Pa] for a gas, concentration [mol/m3] for a liquid. The
[T] mark means the type has temperature-scalable parameters.

```text
Henry         [T]   [K_H]                q = K_H x
Freundlich    [T]   [K_F, n]             q = K_F x^(1/n)
Langmuir      [T]   [q_sat, b]           q = q_sat b x / (1 + b x)
Langmuir-Freundlich [T] [q_sat, b, nu]   q = q_sat b x^nu / (1 + b x^nu)
Sips          [T]   [q_sat, b, n]        q = q_sat (b x)^(1/n) / (1 + (b x)^(1/n))
Toth          [T]   [q_sat, b, t]        q = q_sat b x / (1 + (b x)^t)^(1/t)
Anti-Langmuir       [a, b]               q = a x / (1 - b x), valid for b x < 1
Redlich-Peterson [T] [K_R, a_R, g]       q = K_R x / (1 + a_R x^g)
BET                 [q_m, b, c]          q = q_m b x / ((1 - c x)(1 - c + b x))
Quadratic           [q_sat, b, c]        q = q_sat (b x + 2 c x^2) / (1 + b x + c x^2)
Temkin              [q_sat, b, c]        q = q_sat [th + c th^2 (th - 1)], th = b x/(1 + b x)
OBrien&Myers        [q_sat, b, s]        q = q_sat [b x/(1 + b x) + s^2 b x (1 - b x)/(1 + b x)^3]
Bingel&Walton       [q_sat, a, b]        q = q_sat (1 - E) / (1 + (b/a) E), E = exp(-(a + b) x)
Unilan              [q_sat, b, s]        DO NOT USE: the current implementation returns zero.
pH-Langmuir   [T]   [q_sat, b, pH_0]     b_pH = b / (1 + 10^(pH_0 - pH)),
                                         q = q_sat b_pH x / (1 + b_pH x). Liquid only, SPI only.
GAB           [T]   [q_m, C, K]          Water: x_r = p / p_sat(T), with p_sat from
                                         WaterSaturationPressureModel (humidity inputs below) at
                                         the column temperature;
                                         q = q_m C K x_r / ((1 - K x_r)(1 + (C - 1) K x_r)),
                                         valid while K x_r < 1. Gas only, SPI only. A
                                         temperature-dependent C = C_0 exp(H_C/(R T)) is the law
                                         {"Value": C_0, "Law": "VantHoff", "Energy": H_C}
                                         without ReferenceTemperature; likewise K.
```

The BET and O'Brien-Myers forms above are the implemented ones and differ from some textbook
versions. Map a source only after rewriting it in this form, and ask if the forms do not match.

## Humidity-dependent isotherms

Gas only, SPI only, not on the Water component, not as chemisorption sites. RH and q_w [mol/kg] are
those of [Humidity inputs](#humidity-inputs-gas-only); each type needs a Water component or a
`RelativeHumidity`, and the GAB-coadsorption and GAB-modified types need a Water component.

```text
Humidity-enhanced-Langmuir [T] [q_sat, b, alpha, k_w]
                                b_eff = b [1 + alpha k_w RH / (1 + k_w RH)], Langmuir
Moisture-swing-Langmuir    [T] [q_sat, b, n_w]
                                b_eff = b RH^(-n_w) (b at RH = 1), Langmuir; q -> q_sat at RH = 0
GAB-coadsorption-Langmuir  [T] [q_sat, b, alpha, omega]
                                b_eff = b (1 + alpha q_w), q_eff = q_sat max(0, 1 - omega q_w);
                                alpha, omega in kg/mol, 0 disables either modifier
GAB-coadsorption-Sips      [T] [q_sat, b, n, alpha, omega]    Sips with the same modifiers
GAB-coadsorption-Toth      [T] [q_sat, b, t, alpha, omega]    Toth with the same modifiers
GAB-modified-Toth          [T] [q_sat, b, t, gamma, beta]
                                q_eff = q_sat / (1 - gamma q_w), b_eff = b (1 + beta q_w), Toth
                                (Stampi-Bombelli et al. 2020; Liu et al. 2024); gamma q_w < 1
Exponential-humidity-Toth  [T] [q_sat, b, t, a]
                                q_eff = q_sat exp(a RH), Toth (Liu et al. 2024, eq. 8)
```

## Temperature dependence

A plain number is constant. A temperature law replaces the number:
`{"Value": V, "Law": "<law>", "<coefficient key>": c, "ReferenceTemperature": T_ref}` (T_ref in K,
optional).

| Law | p(T) | Coefficient key | Unit |
|---|---|---|---|
| `VantHoff` | V exp[Q/R (1/T - 1/T_ref)] | `"Energy"` Q | J/mol (Q > 0: grows as T drops) |
| `Arrhenius` | V exp[-E/R (1/T - 1/T_ref)] | `"ActivationEnergy"` E | J/mol (E > 0: grows with T) |
| `Exponential` | V exp[chi (1 - T/T_ref)] | `"Chi"` chi | - |
| `Linear` | V + a (T - T_ref) | `"Slope"` a | unit of V per K |
| `Inverse` | V + alpha (1 - T_ref/T) | `"Alpha"` alpha | unit of V |

V is the value at T_ref. A law without `ReferenceTemperature` takes the component's
`ReferenceTemperature`. Without any T_ref, VantHoff and Arrhenius are V exp[Q/(R T)] and
V exp[-E/(R T)] (V is then a pre-exponential factor); the other laws need a T_ref.

Where each law is allowed (any other position must be a number):

- q_sat or q_m (Langmuir, pH-Langmuir, Sips, Langmuir-Freundlich, Toth, GAB, and the [T] humidity
  types): `Exponential`, `Linear`
- affinity b or k, K_H, K_F (Henry, Freundlich, Langmuir, pH-Langmuir, Sips, Langmuir-Freundlich,
  Toth, and the [T] humidity types): `VantHoff`, which may omit `"Energy"` to take
  `HeatOfAdsorption`
- exponent n, nu, t (Freundlich, Sips, Langmuir-Freundlich, Toth, and the Sips and Toth humidity
  types) and g (Redlich-Peterson): `Linear`, `Inverse`
- K_R, a_R (Redlich-Peterson) and C, K (GAB): `VantHoff` with `"Energy"` given

Heat: a given `HeatOfAdsorption` is the energy-balance heat; affinity laws with their own `Energy`
keep it (it may differ). Without `HeatOfAdsorption`, the heat is the `Energy` of the affinity laws,
which must then agree when `energyBalance` is on. A chemisorption site's own `HeatOfAdsorption`,
else the `Energy` of its affinity law, else the component's `HeatOfAdsorption` is that site's heat.

If the source gives a temperature-dependent isotherm, map each parameter to the law of the same
form, and say so; if no law fits, ask.

## Mixture methods

`"MixturePredictionMethod"` (default `IAST`):

- `IAST`: ideal adsorbed solution theory. The general default.
- `SIAST`: IAST applied per site index. Use it when sites represent distinct adsorption
  environments.
- `SPI`: sites evaluated independently, without competition. Required for GAB, pH-Langmuir, and the
  humidity-dependent isotherms.
- `SCI`: segregated competitive isotherm. Each site index must use ONE model across all components,
  chosen from Langmuir, Anti-Langmuir, Sips, Langmuir-Freundlich, Redlich-Peterson, Toth.
- `EI`: explicit size-effect Langmuir. Every physisorption site must be Langmuir; no chemisorption.
- `SEI`: segregated EI. Langmuir only, for physisorption and chemisorption sites alike.
- `MPD`: macrostate particle distribution ([optional-blocks.md](optional-blocks.md)). Gas only.

`"IASTMethod"` for IAST and SIAST: `"FastIAST"` (default) or `"NestedLoopBisection"` (more robust,
slower).

If the source uses a mixture model, use that one. Otherwise use IAST and say so.

## Humidity inputs (gas only)

For the humidity-dependent isotherms above. The relative humidity is RH = p_H2O / p_sat(T), with
p_H2O the partial pressure of the component with `"Water": true` and T the gas temperature, at every
grid point. `"RelativeHumidity"`: a number >= 0 fixes RH instead (the Water component then only
supplies q_w, at p_H2O = RH p_sat(T)). The water loading q_w is the loading of the Water component's
isotherm at that p_H2O.

`"WaterSaturationPressureModel"`: `"Antoine"` (default, 274-373 K), `"Wagner"` (IAPWS equation of
Wagner and Pruss, 273.16-647.096 K), or `"Constant"` with `"WaterSaturationPressure"` [Pa]. It also
sets the p_sat of GAB.

## pH inputs (liquid only)

- `"pHMode": "Fixed"` (default) uses `"pHValue"` (default 7).
- `"pHMode": "HPlus"` or `"OHMinus"` computes the pH from the transported component selected by
  `"pHComponent"`, a zero-based integer component identifier in `Components` (for example, 0 selects
  the first component). Component names are not accepted. `"pKw"` (default 14) is used with
  `OHMinus`.
- The H+ or OH- component needs a disabled isotherm (`{"Type": "Henry", "Parameters": [0.0]}`).
