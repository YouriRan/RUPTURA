# Revision 5 — 12 November 2026

Shared gas/liquid chemisorption supersedes the earlier chemistry-interface
account below. ChemisorptionSites now includes DualNth, directional rates and
explicit Avrami history policies; ChemicalKinetics is import-only. See
[the common chemistry manual](SHARED_CHEMISORPTION.md) and
`verification/shared_chemisorption/summary.json` for current scope and checks.

# Shared column checkpoint — 11 November 2026

This checkpoint resumes `RUPTURA_unified_implementation_checkpoint(1).zip`.
It contains the implemented shared-column path and the cycle editor, plus fixes
and reproducible checks added after inspecting the interrupted work. It is not
a claim of a fully validated scientific release.

## Requested scope

| Requirement | Implementation and verification boundary |
|---|---|
| One physical column model | Explicit FDM and FVM both use `cycle::Column`; breakthrough and cycles both use `CycleRunner` and `ColumnNetwork`. Spatial operators and time integrators are separate. |
| Same numerical choices | Both operations expose Upwind/TVD/WENO3/WENO5 and BDF1/BDF/CVODE/SSPRK3/SemiImplicitSSPRK3. Forty matched selections execute as breakthrough and as a one-phase cycle with equal saved states. |
| No reaction networks in cycles | Rejected by the public adapter, original compatibility parser, cycle constructor and network scope. Column overrides cannot acquire the breakthrough-only capability. Chemisorption rate laws remain allowed. |
| Original main variables | Public JSON names, units, component order, isotherm parameter order, adsorbent sections and overrides are retained. The shared state uses bulk concentration, site loading, particle histories and fluid/solid/wall temperatures. Files/checkpoints retain component-major ordering; internal solver ordering may differ. |
| HSDM with/without film for both axial methods | Shared radial shell closure, spherical and radial-cylinder particles, gas/liquid fixed beds. Sixteen derivative/conservation combinations and matched finite-time trajectories test the volume-weighted mean loading. |
| Equal breakthrough capabilities for FDM and FVM | Shared material, kinetic, reaction, thermal, flow, output and stopping implementations. Tests exercise gas/liquid, packed beds/monoliths, original isotherm adapters, transported pH, five main chemisorption rate laws and conservative reaction production. Gas MPD now calls the original kernel on both shared paths. |
| Cycle GUI and one configuration | The Cycle designer edits the same simulation object as Input JSON and cluster export. Equipment, valves, schedules, CSS, buffer controls and compressor fields have component tests. Full browser rendering remains unverified in this environment. |

## Corrections made during recovery

- The original material parser pads shorter isotherm lists, including carriers.
  The shared adapter now skips those synthetic entries when reading raw site
  settings, preventing a crash for original-main isotherm models.
- Column inventories already included surface/pore storage in the uploaded
  checkpoint. Network totals did not. Network totals now include the same
  contribution, so material-balance and CSS inventory calculations account for it.
  Closed-system mass and adiabatic energy derivative checks pass.
- Cycle-only reaction exclusion now also applies to compatibility inputs,
  programmatically constructed compatibility cycles and per-column overrides.
- Liquid buffers created in the GUI now use `InitialConcentration` and volume
  controls. Valve renaming preserves the same valve identity across phases.
  Compound `All`/`Any` stop conditions are preserved. Per-column temperature
  coupling controls display the inherited value.
- The native CVODE header no longer includes an unused optional KLU header.
  The installer includes the browser launcher, uses static dependency builds
  for the portable compiler, and disables module scanning and host CPU tuning.
- The Avrami correction already present in the supplied checkpoint is retained
  and documented in `RELEASE_NOTES.md`.

## Reproduce checks

```bash
./install.sh --python --tests
RUPTURA_CLI_EXECUTABLE="$PWD/build/ruptura" python -m pytest -q tests
ctest --test-dir build -R '^process_' --output-on-failure
python tools/generate_unified_examples.py
```

The inherited native suite references 28 tests whose input fixtures are absent
from the supplied archive. Consequently the unfiltered native suite (and
`install.sh --tests`, which runs it) reports these failures. The required
missing files are listed in `verification/unified/summary.json`; no substitute
literature data were fabricated. Shared process tests are separately runnable
with the CTest filter above.

Optional cycle-editor component checks:

```bash
npm install --prefix .build-tools/dom-tests jsdom
NODE_PATH="$PWD/.build-tools/dom-tests/node_modules" node tests/cycle_editor_dom_test.cjs
```

These checks execute the actual editor JavaScript in a DOM test environment.
They verify configuration changes and highlighting classes, not screen layout,
pointer dragging or accessibility in a full browser. Full browser verification
could not be completed because the available browser could not access the local
test server, and its security policy rejected local-file navigation. Existing
browser smoke code remains available for a suitable local environment.

The recovery results and source hashes are retained in `verification/unified/`.
The subsequent MPD integration results are in `verification/mpd_shared/`.
Input validation now covers the existing 236 integrated inputs and 32 matched
FDM/FVM inputs. It does not establish that the full campaigns have run to CSS.

## MPD addition

MPD is available for explicit FDM/FVM gas breakthrough and cycles, including LDF,
fixed-bed HSDM with/without film and monolith LDF. The original MPD loading,
carrier mapping, framework-mass conversion and local-temperature reweighting
are reused. Lab JSON and cluster export carry the same settings; exports include
distribution bytes and change only their file paths. Changed distribution
content invalidates a saved checkpoint.

Forty Python integration tests compare both methods against a synthetic,
exactly equivalent competitive Langmuir model, exercise all five integrators,
HSDM geometries/film options, cycles, CLI/Python and actual cluster-bundle
execution, and check invalid inputs and restart compatibility. Native checks
also compare against the original MPD kernel at different local temperatures
with a conditional-energy column and a carrier between the adsorbing species.
See [MPD integration notes](MPD_SHARED_MODEL.md) for the verification boundary.

## Still outside the completed verification

- HSDM is for nonreactive physisorption in fixed beds; HSDM monoliths and coupled
  particle reaction transport are rejected for both axial methods.
- Full browser regression, the 28 missing-fixture native tests, mesh-convergence
  studies, completed MFI CSS campaigns, and Haghpanah literature validation are
  not completed by this checkpoint.
