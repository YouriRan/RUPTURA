# Integration status — 8 November 2026

This is the historical integration record. The resumed 11 November checkpoint
is described in [UNIFIED_COLUMN_STATUS.md](UNIFIED_COLUMN_STATUS.md). Its explicit
FDM/FVM shared path, material adapters and cycle editor supersede the feature
and remaining-work descriptions below where they differ.

This is a first integrated source build based on `RUPTURA-main_7Nov2026`, with
the supplied standalone breakthrough and cycle C++ packages as feature donors.
The source package deliberately retains the original main calculation path while
the FVM process path is made reviewable. It is not a claim of full numerical
equivalence between all three original simulators.

## Flow/thermal revision

See [the agreed flow and thermal update](FLOW_THERMAL_UPDATE.md). The removed
boundary case is no longer available. The original-main default is FixedVelocity.
FVM adds constant-velocity gas/liquid controls, a liquid prescribed-average-gradient
control, Ergun sphericity, explicit thermal selectors, common fluid/solid temperature
and initial-condition presets. Nine short examples accompany this revision.
General reaction networks are deliberately outside cycle scope. MPD was deferred
at this stage and is now implemented for shared FDM/FVM gas columns; see
[MPD_SHARED_MODEL.md](MPD_SHARED_MODEL.md).
Current verification is under `verification/flow_update/`; earlier integration
records below retain their historical context.

## Implemented

| Area | Shared FVM path |
|---|---|
| Breakthrough / cycles | Same native column, four gas/liquid × fixed-bed/monolith combinations |
| Inputs | Main JSON names, component order, adsorbent overrides, section grids, phases; new optional blocks |
| Numerics | Conservative axial FVM, Upwind/TVD/WENO5, native BDF/NDF, CVODE Dense/preconditioned SPGMR |
| Particles | Nonreactive fixed-bed LDF/HSDM, sphere/radial cylinder, finite film, Arrhenius diffusion |
| Heat | Separate fluid/solid/wall, gas LTE option, distributed transient coolant, thermal phase controls |
| Flow | Gas Ergun/Darcy and channel momentum; incompressible liquid column/network hydraulics |
| Materials | SPI/SCI/IAST/SIAST/SEI, imported humidity models, capacity-temperature law, fixed-pH Langmuir |
| Chemistry | Imported gas PFO/PSO/Avrami/dual-nth model through `ChemicalKinetics`; continuous histories |
| Equipment | Sources, valves, direct bed connections, compressors, gas buffers, variable-volume liquid buffers |
| Controls | Timed/conditional phases, independent unconnected-bed readiness, waits, shared connection-phase stops |
| Runtime | Checkpoints, exact RHS cache, native structured FD controls, compact output, integrated ledgers |
| Interfaces | CLI, Python dispatch, Lab JSON controls, recorded-array plots, portable HAL arrays and build script |
| Input migration | 27 MFI breakthrough + 162 MFI cycle inputs; selected equipment, chemistry, moisture and Haghpanah cases |

The Marx T15 inputs now use FixedVelocity because their former boundary mode
was removed. They are migrated cases rather than unchanged upstream copies.
The upstream blob recorded for the pre-revision input was
`6f12c317b193aa1f2afde8be0fc6f785817a3d83`. Earlier reference curves and timing
records refer to that previous build, not a new validation of these migrated inputs.

## Remaining work and limits

1. In this historical build, main's MPD, full pure-isotherm catalogue, legacy `ChemisorptionSites`, general
   reaction networks and transported pH are not connected to the shared FVM cycle
   path. They remain on the original main path. Unsupported FVM selections must
   be resolved explicitly; chemical particle transport and HSDM monoliths are not
   implemented. General reaction networks in cycles are excluded by the agreed
   scope rather than scheduled for this revision.
   The later [shared-column checkpoint](UNIFIED_COLUMN_STATUS.md) supersedes
   these material-adapter and MPD omissions; particle transport limits remain.
2. The breakthrough donor's thermal formulation is replaced by the shared cycle
   energy model. Converted material, flow and hardware values do not establish
   identical thermal curves. Scientific reuse requires matched comparisons.
3. In the earlier integration build, the original main Marx T15 run completed its full 4,000 s. Its shared-FVM
   counterpart passed a 1 s / 101-node startup check (exact Sips, CVODE order 2),
   but reached only 60 s in the 360 s full-run limit. Its complete breakthrough
   comparison and runtime improvement remain open.
4. Haghpanah literature validation remains open. A separate note records the
   earlier discrepancy, the inlet-density option, references, and what the new
   integration checks do and do not establish.
5. Migrated campaigns retain their original screening schedules and their stated
   source/backflow assumptions. They are not optimized recipes, newly converged
   CSS results, or a new plant-sizing/TEA study.
6. The Lab exposes the new options through structured JSON cards and saved-array
   plots. More dedicated form controls, visual connection editing, and end-to-end
   browser testing remain useful follow-up work. No custom Python RHS or stop
   callable is executed by the native JSON engine.

## Verification record

`verification/integration/` records the actual native tests, input checks,
complete small-example runs, Python/Lab checks and reference-case runs performed
for this source build. A passed input check means the input and initial equations
are accepted; it does not establish a completed simulation or literature agreement.
The two MFI campaigns were not fully rerun to CSS in this integration session.

See [the shared model manual](SHARED_FVM_MANUAL.md) for the usable JSON API and
[Haghpanah status](HAGHPANAH_STATUS.md) for the scientific comparison boundary.
