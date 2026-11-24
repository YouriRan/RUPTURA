# Chemisorption and reactions

Add these blocks only when the user's model contains chemisorption or reactions.

## Chemisorption sites

`"ChemisorptionSites"` is an array inside a component. Each site is

```text
{
  "Type": "<isotherm type>", "Parameters": [...],   equilibrium isotherm q*_chem (isotherms.md)
  "MaximumLoading": q_max [mol/kg, > 0],             required
  "HeatOfAdsorption": Q_chem [J/mol],                optional; default: component value
  "RateEquation": "<model>",                         optional; default "FirstOrder"
  "RateEquationParameters": { ... }                  see below; keys are CASE-SENSITIVE
}
```

A chemisorbing component still needs a `PhysisorptionSites` entry (SKILL.md, Components).

## Rate equations

The rate coefficients (`MassTransferCoefficient`, `RateCoefficient`, `AdsorptionRateCoefficient`,
`DesorptionRateCoefficient`, `Alpha`) are numbers or Arrhenius laws ([isotherms.md](isotherms.md)),
k(T) = k0 exp[-E/R (1/T - 1/T_ref)].

**FirstOrder**: `{"MassTransferCoefficient": k [1/s]}`. Optional; default: the component's k_LDF.
Optional desorption branch (all three together): `"DesorptionMassTransferCoefficient"` k1 [1/s] > 0,
`"DesorptionCoverageEnergy"` beta [J/mol], `"CoverageReferenceLoading"` q_ref [mol/kg] > 0; while
q > q*: dq/dt = k1 exp(beta q / (q_ref R T)) (q* - q).

**PseudoNth**: `{"RateCoefficient": k_n [(mol/kg)^(1-n) / s], "Order": n >= 0}`; optional
`"DesorptionRateCoefficient"` and `"DesorptionOrder"` for the release (q > q*), each a number or a
temperature law (e.g. Arrhenius with its own `ActivationEnergy`).
dq/dt = k_n (q* - q) |q* - q|^(n-1).

**Avrami**: `{"RateCoefficient": k_A [1/s], "Order": n >= 0}`.
Uptake: dq/dt = n k_A s^(n-1) (q* - q), s = max([-ln(1 - q/q*)]^(1/n), 1e-3), so a clean sorbent at
constant q* follows q/q* = 1 - exp(-(k_A t)^n) with t counted from its own start (no global clock);
release (q > q*): dq/dt = k_A (q* - q).

**General**: ALL of the following are required: `"AdsorptionRateCoefficient"`,
`"DesorptionRateCoefficient"`, `"PoreConcentrationOrder"`, `"CapacityOrder"`, `"DesorptionOrder"`
(non-negative integers), `"FilmMassTransferCoefficient"` [m/s], `"PoreDiffusivity"` [m2/s],
`"UsePoreSurfaceTransport"` (bool).
dq/dt = k+(T) c_p^m (q_max - q)^n+ - k-(T) q^n-; give an activation energy as an Arrhenius law of
k+ or k-.

**Elovich**: ALL of `"Alpha"` [m3/(kg s)], `"Beta"` [kg/mol], `"FilmMassTransferCoefficient"`,
`"PoreDiffusivity"`, `"UsePoreSurfaceTransport"`. dq/dt = Alpha c_p exp(-Beta q).

The units of the General rate coefficients depend on the orders. Never guess them; ask.

## Particle transport (physisorption)

The root `"ParticleTransport": {"Model": "HSDM", "RadialCells": 12}` replaces the linear driving
force of physisorption by surface diffusion in a spherical particle, with
D_s / R^2 = MassTransferCoefficient / 15; not with `Reactions` or CVODE KLU. Leave it out unless
asked.

## Reactions

Root `"Reactions"` array; each reaction:

- `Phase`: `"Physisorbed"` | `"Chemisorbed"` | `"PoreConcentration"` | `"BulkConcentration"`
  (BulkConcentration acts on the bulk-fluid concentrations, per fluid volume, with the catalyst
  amount folded into k+; it requires `energyBalance` false).
- `Style`: `"GeneralPowerLaw"` | `"LangmuirHinshelwood"` | `"LangmuirHinshelwoodHougenWatson"` (LH
  and LHHW require Phase `"PoreConcentration"` or `"BulkConcentration"`).
- `Reactants`, `Products`: non-empty arrays of component names (preferred) or zero-based indices; no
  component may appear twice.
- `Stoichiometry`: optional positive numbers, either one array (reactants then products) or
  `{"Reactants": [...], "Products": [...]}`; default all 1.
- `ForwardOrders`, `BackwardOrders`: optional positive arrays; default = stoichiometry.
- `Site`: optional zero-based site index (default 0). A Chemisorbed reaction needs that
  chemisorption site on every participant in every adsorbent.
- `RateLimitTime`: optional [s] (default 1e-4). This is the shortest time over which a reactant may
  be depleted or a product site filled.
- `Kinetics` (required):
  - `"ForwardRateCoefficient"` k+ (>= 0; units follow the orders; a number or an Arrhenius law),
  - `"EquilibriumConstant"` K (> 0; a number or a VantHoff law),
  - optional `"HeatOfReaction"` Q_r [J/mol] (= -dH_r, positive when exothermic) and
    `"ReferenceTemperature"` [K] for the laws.

  Ruptura uses k- = k+ / K. A VantHoff K without `"Energy"` takes `HeatOfReaction`; without
  `HeatOfReaction` the heat is the `Energy` of K (0 if K is constant). To describe an effectively
  irreversible reaction, make K very large and report that.
- LH and LHHW also need `"AdsorptionConstants": {"<component name>": K_j, ...}` in `Kinetics`, one
  K_j [m3/mol] per adsorbing component (a number or a VantHoff law with `"Energy"`; every reactant
  needs K_j > 0; components left out do not adsorb). They set the coverages
  theta_j = K_j c_j / (1 + sum_i K_i c_i). LH is irreversible, r = k+ prod theta_r^order; LHHW is
  r = k+ prod theta_r^order - k- prod theta_p^order.
