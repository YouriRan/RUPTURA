# Haghpanah comparison: model correspondence and validation status

Two references must be kept separate:

* Haghpanah et al. (2013), numerical adsorption/VSA benchmark, DOI
  [10.1021/ie302658y](https://doi.org/10.1021/ie302658y).
* Krishnamurthy et al., including Haghpanah (2014), dry-flue-gas pilot VSA,
  DOI [10.1002/aic.14435](https://doi.org/10.1002/aic.14435).

The supplied cycle package includes source-derived 2013 candidates and a
conditional reconstruction of the 2014 breakthrough, Run 1, Run 6 and Run 7.
Their settings have been converted to main-style JSON. Conversion preserves
reported parameters and explicit donor assumptions. It does not make a conditional
reconstruction an exact reproduction of the experiment.

## Earlier evidence, not new results from this build

The 30 September numerical review of the user's separate `numerics.zip` reported
WENO5/121-node purity 79.172% and recovery 36.446%, versus the 2013 reference
88.740% and 35.830%. Mesh/tolerance refinement did not remove the purity gap.
The review associated the principal compositional discrepancy with excess
nitrogen in the evacuated product and identified inlet-density and energy
equation correspondence as unresolved checks. Those figures describe that older
run, not the newly integrated source.

The earlier 2014 validation document reported engineering agreement for a
reconstruction, with explicit uncertainties in run temperature, feed basis,
scheduling and connection conductance. The native reference markers supplied in
the donor are approximate/digitized; retain that provenance when plotting them.

## What changed here

The shared model now exposes `InletDensityBasis` explicitly:

| Choice | Prescribed-velocity gas inlet |
|---|---|
| `Reservoir` (default) | Density from the stated source pressure and temperature |
| `LocalPressure` | Face pressure recovered from local Darcy/Ergun resistance; density follows that pressure |

For a Darcy inlet with inward superficial velocity U, the recovered face pressure
is the adjacent state pressure plus the boundary hydraulic resistance times U.
For Ergun flow, inlet density also enters the quadratic resistance; the scalar
relation is solved algebraically. Incoming species and enthalpy fluxes use the
same recovered concentration. Imposed molar-flow inputs retain their explicit
common-feed basis and reject the local-pressure velocity option.

`examples/new_examples/integrated/validation/haghpanah_2013_equations/dcb_local_pressure/`
provides an explicit 41-node boundary comparison. The source-derived donor input
is retained beside it. A unit check verifies that the recovered-pressure option
changes inlet molar flux in the expected direction. This alone does not quantify
the correction to the old purity discrepancy.

The energy model is stated in the shared manual. Gas storage uses Cv, transported
enthalpy uses Cp, and adsorbed energy uses a component heat/reference temperature
and adsorbed heat capacity. The earlier paper-to-code energy mapping has **not
been demonstrated algebraically equivalent**. No undocumented heat refit is used
to force agreement.

## Acceptance before calling this literature-validated

1. Compare initial/boundary conditions and molar feed, including pressure recovery.
2. Establish the gas/adsorbed energy and heat-reference correspondence.
3. Compare breakthrough outlet composition, velocity, temperature and integrated
   inventory/flux balances with the same output-time resolution.
4. Compare end-of-step profiles and CSS product balances for the full four-step
   cycle; track evacuation nitrogen as well as carbon dioxide.
5. Separate mesh sensitivity, time tolerance, CSS periodicity, conservation and
   literature agreement. None is a substitute for another.

New integration runs, where performed, are recorded in
`verification/integration/reference_runs.json`. A completed short test or a
passed material-law regression is reported as such, without relabeling it as
experimental validation. The 162-case MFI migration is independent of this
reference-model question.

## New breakthrough calculations in this integration

The paired 41-node source-derived 2013 cases were run for 5,000 s with the same
recipe and numerics, changing the inlet-density basis. The 2014 adsorption
reconstruction was run for 7,000 s. These are completed transients, not CSS runs.

| Quantity | Reservoir-density inlet | Local-pressure inlet |
|---|---:|---:|
| 2013 first-front half-plateau time [s] | 529.35 | 519.16 |
| 2013 cooling-front temperature midpoint [s] | 3786.50 | 3709.64 |
| Maximum outlet temperature [K] | 359.847 | 359.852 |

The boundary choice advances the thermal transition by about 77 s in this paired
calculation. It does not establish that the remaining energy formulation matches
the paper or that the earlier four-step purity discrepancy is resolved.

For the conditional 2014 adsorption reconstruction, CO2 outlet mole-fraction RMSE
against the donor's approximate manual markers is 0.02208. The overlay shows a
sharper first front and a higher intermediate plateau than those markers. This
is a remaining discrepancy, not a validated experimental fit. The digitization
itself has stated resolution of roughly 50 s and 0.002 in mole fraction.

New normalized integrated mass/energy residuals are approximately 1e-9 to 1.2e-8.
The actual metrics, outlet CSVs and `haghpanah_comparison.png` are retained in
`verification/integration`. They document conservation and model comparisons;
full literature cycle validation remains open.
