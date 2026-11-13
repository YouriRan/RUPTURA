# Shared column and cycle GUI implementation

Requested 11 November 2026. Baseline archive SHA256:
7261597945c9b21e59932f39f35469971d072e9bdf77e6082c75fdbfa337dda1.
Cycle reference archive SHA256:
51aa469c4f60b65c311f1ccc99b4ba00ea585ef097bebc4233ed6f3a767c6cbc.

## Agreed scope

- One shared physical column model for explicit FDM and FVM breakthrough/cycles.
- Main JSON names and units remain the public configuration basis.
- Shared numerical choices: Upwind/TVD/WENO3/WENO5; native BDF1, existing adaptive BDF,
  CVODE, SSP-RK3 and a verified semi-implicit third-order method.
- HSDM sphere/radial cylinder, with/without film, nonreactive fixed beds, for both
  axial methods. Same radial shell model and volume-weighted loading source.
- Breakthrough capability parity, with general reaction networks excluded from cycles.
  MPD now uses the original distribution kernel on both shared gas paths.
- Functional cycle flowsheet and schedule GUI with equipment forms, stopping/CSS,
  JSON import/export and cluster launch bundle, using the same simulation configuration.
- Preserve compatibility of original main inputs, update examples, docs and bindings.

## Work sequence

1. Inspect and preserve baseline; build a portable executable.
2. Shared axial operators, WENO3, common input routing and capability discovery.
3. Shared integrators and verification of time order, events and passive ledgers.
4. Common material/kinetic/reaction/pH adapters, with explicit cycle exclusions.
5. Cycle GUI, structured equipment/schedule editing and cluster export.
6. Native, Python, GUI and matched FDM/FVM checks; update documentation and package.

## Completion criteria

Numerical options must execute, not merely parse. Film/no-film HSDM must use the
same radial closure in both axial modes. Integration must preserve material/
energy accounting and localize cycle events. Exported GUI JSON must validate and
run with the native CLI. All reported limitations and tests must distinguish
short verification runs from completed CSS campaigns or literature validation.

## Progress

- Baseline and donor archives inspected and extracted to an isolated working tree.
- The uploaded checkpoint contained shared spatial/integration work, original-main
  material adapters and a cycle editor, but its status files were stale.
- Resumed on 11 November: corrected carrier-site adapter access, network
  surface/pore totals, reaction scope through cycle column overrides, GUI liquid
  buffers/valve identity, optional KLU build handling and GUI installation.
- Added matched FDM/FVM examples, native physical equivalence/conservation tests,
  executed CLI/Python comparisons and cycle-editor DOM checks.
- Added MPD to both shared axial methods, including local-temperature kernel
  comparisons, all shared integrators, HSDM, Lab JSON and cluster file export.
- See `docs/UNIFIED_COLUMN_STATUS.md`, `verification/unified/summary.json` and
  `verification/mpd_shared/summary.json` for completion evidence and limitations.
  Full browser regression, missing historical fixtures and literature/CSS
  campaigns remain open.

## 12 November: common gas/liquid chemisorption

- Consolidated loading-rate equations in Chemisorption and delegated the retired
  donor adapter to that kernel. ChemisorptionSites is the canonical public input.
- Added DualNth branch states, directional/reference-temperature coefficients,
  signed pseudo-nth regeneration, fractional bases and Avrami history policies.
- Carried branch heats and histories through balances, stage changes, restart,
  CSS, public JSON, Python, browser forms and cluster export.
- Added synthetic examples and dedicated analytic, trajectory, conservation,
  migration, state-history and GUI component checks. Experimental validation
  and complete CSS campaigns remain outside this change.
