# Shared chemisorption for gas and liquid

Revision 5, 12 November 2026. `Components[].ChemisorptionSites` is the public
chemisorption interface for gas/liquid breakthrough and cycles, using either
FDM or FVM. The loading-rate kernel is independent of fluid phase. The column
supplies equilibrium loadings using gas partial pressures (Pa) or liquid
concentrations (mol/m³), and accounts for the resulting material and heat terms.
General reaction networks remain excluded from cycles.

## Rate equations

`RateEquation` accepts `FirstOrder` (`PFO` alias), `PseudoNth` (`PSO` means order 2),
`Avrami`, `General`, `Elovich`, and `DualNth` (`dual_nth` alias).
Existing General/Elovich surface/pore transport controls are retained. They are
separate from the nonreactive physisorption HSDM model; reactive HSDM is not added.

FirstOrder, PseudoNth, Avrami and DualNth accept rate coefficients at
`ReferenceTemperature` (default 298.15 K), separate adsorption/desorption
activation energies, and independent desorption parameters. Their Arrhenius
factor is `exp(E/R * (1/Tref - 1/T))`. General retains its original pre-exponential
coefficient convention `k(T) = k0 exp(-E/(RT))`; it is not silently reinterpreted.

Pseudo-nth rates preserve the sign of the loading difference:

`dq/dt = k_ads (q* - q)^n_ads` for adsorption,
`dq/dt = -k_des (q - q*)^n_des` for desorption.

This corrects the previous raw-power desorption behavior for even or noninteger
orders. Regeneration results using those older settings can therefore change.
Sub-unity orders use the configured `DrivingRegularizer` near zero driving force.
For `RateBasis: Loading`, coefficient units depend on order:
`(mol/kg)^(1-n) / s`. For `RateBasis: Fractional`, rates use the explicitly stated
`ReferenceLoading`, or MaximumLoading if omitted, and coefficients have units 1/s.
Avrami retains a loading basis.

## DualNth

A DualNth entry is ONE physical equilibrium site with TWO kinetic loading
states. Branch targets are `f_b q*`, fractions sum to one for each component,
and total loading is `q_1 + q_2`. Equilibrium capacity is not counted twice.
Each branch has its own adsorption/desorption coefficient, activation energy
and positive order. Optional branch heats participate in both stored energy and
uptake/desorption energy rates.

```json
{
  "Type": "Langmuir",
  "Parameters": [1.0, 0.5],
  "MaximumLoading": 1.0,
  "HeatOfAdsorption": 25000.0,
  "RateEquation": "DualNth",
  "RateEquationParameters": {
    "ReferenceTemperature": 298.15,
    "RateBasis": "Loading",
    "Branches": [
      {"Fraction": 0.35, "AdsorptionRateCoefficient": 0.4,
       "DesorptionRateCoefficient": 0.2, "AdsorptionOrder": 1.0,
       "DesorptionOrder": 1.0},
      {"Fraction": 0.65, "AdsorptionRateCoefficient": 0.1,
       "DesorptionRateCoefficient": 0.06, "AdsorptionOrder": 1.5,
       "DesorptionOrder": 1.5}
    ]
  }
}
```

Here the affinity 0.5 has units m³/mol for a liquid. For the matched gas
example at 298.15 K, use `0.5/(R*298.15)` in Pa⁻¹. Arbitrary gas and liquid
parameter sets are not automatically interchangeable.

## Avrami histories

Existing inputs default to `HistoryPolicy: ElapsedTime`, retaining
`dq/dt = n k^n t^(n-1)(q* - q)` and the existing finite initialization convention.
The shared column integrates this age continuously across phases.

`HistoryPolicy: ContinuousErasure` retains the previous donor model's separate,
dimensionless adsorption/desorption ages. During adsorption, adsorption age
increases at k_ads while desorption age decays; during desorption the roles
reverse. Optional `MemoryErasureRateAds` and `MemoryErasureRateDes` override the
opposite-direction coefficient used for decay. `AgeRegularizer` is dimensionless
in this convention. Stage changes do not reset either history. Their states enter
checkpoints and CSS checks.

`HistoryDrivingRegularizer` is an optional loading scale (mol/kg), default zero.
For a positive value epsilon, age growth/erasure is multiplied by
`drive^2 / (drive^2 + epsilon^2)`. This removes the discontinuity of history
switching at an almost empty/equilibrated node, without changing the instantaneous
loading-rate equation. It changes subsequent age evolution near zero driving
force and must be reported. The moving-front synthetic examples explicitly use
1e-7 mol/kg; legacy imports retain zero unless edited. AgeRegularizer and this
loading-scale regularizer have different meanings and units.

## Equilibrium and coexisting physical uptake

A chemical site's Type/Parameters specify its equilibrium independently of its
kinetic equation. Optional site `ReferenceTemperature` and `EquilibriumHeat`
set its affinity temperature dependence. `AdsorptionReferenceTemperature` defines
the reference for the stated adsorption heat independently; legacy conversion
retains the previous energy reference. Optional `MixturePredictionMethod`
selects a common chemical mixture method; nonempty overrides must agree across
chemical sites/components. Otherwise the column's method is inherited.

For migration of a shared physical/chemical capacity, `EquilibriumSource:
Physisorption` obtains the physical mixture target before its uptake-capacity
factor, and `EquilibriumFraction` partitions that target. The associated
`PhysisorptionFraction` scales physical loading. Independent physical and
chemical sites remain additive. Imported partitions are explicit in canonical
JSON, preventing duplicated capacity.

## Legacy import, GUI and restart

The retired top-level `ChemicalKinetics` syntax is accepted only as an import
format. It is converted to ChemisorptionSites before normal simulation parsing;
the former gas-only guard is removed. The compatibility rate adapter delegates
to the same Chemisorption kernel and has no independent kinetic equations.

```bash
build/ruptura old.json --canonicalize > simulation.json
build/ruptura simulation.json --validate
```

CLI, Python ProcessSimulation, Lab import, JSON editing and cluster export use
the canonical configuration. Lab offers DualNth branch forms and rate/history
controls. The existing original compatibility path still accepts its older scalar
laws; stateful/new chemistry selections use the shared column.

Chemisorption loading/history layouts have changed. Start a fresh calculation
when migrating old chemistry checkpoints; incompatible checkpoints must not be
silently reused. Newly written checkpoints preserve branch and history states.
This revision does not add experimental validation or completed CSS campaigns.

## Reproduce verification

```bash
python -m pytest -q tests
ctest --test-dir build -R '^process_' --output-on-failure
NODE_PATH=/path/to/jsdom/node_modules node tests/chemisorption_editor_dom_test.cjs
```

The verification record is `verification/shared_chemisorption/summary.json`.
All 400 Python checks and 17 shared native groups pass. The unfiltered native
suite still has the same 28 missing-fixture failures recorded in the recovered
checkpoint. All 300 shared example inputs validate, and all 32 new synthetic
chemistry examples complete with plotting disabled during verification. GUI
checks exercise real editor code in a DOM; full screen layout remains unverified.
