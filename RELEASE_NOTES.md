# Automatic figures for every new example — 12 November 2026

- Enabled Output.Plots in all 296 new inputs and maintained it in generators.
- Added native CLI automatic plots and --plot-only support for retained legacy breakthrough text outputs.
- Legacy equations, physical parameters and WENO3 choices are retained.
- HAL uses Matplotlib by default; the collector reports PNG/PDF counts and plot status.
- Recompile after extracting this revision and activate the project virtual environment before submission.

# Example organization and WENO3 update — 12 November 2026

- Moved integrated, unified and shared chemisorption examples under `examples/new_examples`.
- Moved 14 inherited material collections under `examples/ruptura1_examples`.
- Removed four identical packed-bed input copies; their canonical inputs remain in unified/fvm.
- All 296 new inputs select WENO3; individual spatial methods, integrators and physical settings are retained.
- Added grouped/uncapped HAL arrays, isolated run outputs, validation, frozen manifests and a result collector.
- Updated generators, execution references, CI paths and documentation.
- Historical numerical/validation records refer to their recorded build and input settings.

# Shared gas/liquid chemisorption — 12 November 2026

ChemisorptionSites is now the canonical interface for both fluids and both axial
methods. DualNth retains two branch loadings under one physical site, including
individual rates, orders, activation energies and optional heats. All kinetic
rate calculations delegate to one fluid-independent Chemisorption kernel.
Legacy ChemicalKinetics is an import format, converted before execution.

PseudoNth now preserves the sign of the driving force during regeneration;
this corrects previous even/noninteger-order behavior. Existing elapsed-time
Avrami and the donor continuous-erasure history policy are explicit choices.
Histories survive stage changes and restart. Chemistry checkpoints from older
builds require a fresh run because the state layout changed. See
`docs/SHARED_CHEMISORPTION.md` and the revision-5 verification records.

# Shared column implementation checkpoint — 11 November 2026

## MPD integration update

Gas MPD is now available for explicit FDM and FVM breakthrough and cycles.
Both share the original distribution kernel, including non-carrier ordering,
particle-to-mol/kg conversion and local-temperature reweighting. The distribution
is loaded once per material, and its content is included in restart compatibility.
LDF, fixed-bed HSDM with/without film and monolith LDF use the same MPD equilibrium.

Lab exposes MPD for these gas columns with the original `MPDSettings` JSON.
Cluster bundles include distribution files automatically and change only their
paths. Two runnable synthetic breakthrough examples and forty Python integration
tests were added, with native comparisons to the original MPD kernel. See
`docs/MPD_SHARED_MODEL.md` and `verification/mpd_shared/` for evidence and scope.

## Shared-column recovery

Explicit FDM and FVM selections use the shared physical column for breakthrough
and adsorption cycles, with separate axial and time-integration selections.
The cycle designer, JSON editor, native CLI, Python process interface and cluster
export use the public main-style simulation configuration.

The supplied checkpoint already contained an Avrami rate correction. It is
retained here: for the integrated law
`q(t) = q_eq [1 - exp(-(k t)^n)]`, the differential rate for positive time is
`dq/dt = n k^n t^(n-1) (q_eq - q)`. The corresponding native derivative test
passes. This correction can change calculations made with older main builds;
the code retains its finite initialization convention at zero age.

This recovery also fixes original-isotherm adapter access to padded carrier
sites, missing surface/pore storage in network inventory totals, cycle reaction
scope through column overrides, GUI-created liquid buffer configuration and
valve renaming across phases. The portable build no longer requires the unused
KLU header and includes a browser-GUI launcher by default.

See `docs/UNIFIED_COLUMN_STATUS.md` and `verification/unified/summary.json` for
actual tests and open work. Full browser
rendering and scientific validation campaigns have not been completed. This
archive is an implementation checkpoint, not a fully validated release.
