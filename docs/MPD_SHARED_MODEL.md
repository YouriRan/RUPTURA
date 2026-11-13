# MPD in the shared FDM/FVM column

Explicit `SpatialDiscretization: "FDM"` and `"FVM"` both support
`MixturePredictionMethod: "MPD"` for gas breakthrough and adsorption cycles.
They use the same physical column and original `MixturePrediction` MPD kernel.
Existing inputs without the explicit spatial selection retain the compatibility
path. MPD remains gas-only.

## Public input

Keep the original `MPDSettings` alongside the ordinary main-style components,
geometry, integration and schedule settings:

```json
"MixturePredictionMethod": "MPD",
"MPDSettings": {
  "FileName": "distribution.data",
  "ReferenceTemperature": 300.0,
  "ReferenceFugacity": 100000.0,
  "ReferenceFrameworkMass": 1.0e-24,
  "ComponentBounds": [
    {"Component": "A", "NMin": 0, "NMax": 10, "DeltaN": 1},
    {"Component": "B", "NMin": 0, "NMax": 20, "DeltaN": 1}
  ]
}
```

The bounds describe non-carrier components in their `Components` order. A
carrier may occur anywhere in that list and has zero MPD loading. Bounds are
inclusive; each extent is `(NMax-NMin)/DeltaN+1`, and the interval must be
divisible by `DeltaN`. The number of distribution rows equals the product of
those extents. Rows use C order, with the last component varying fastest.
Each row supplies a finite probability in [0,1], optionally followed by its
conditional mean Hamiltonian in joules. Blank lines and full-line `#` comments
are allowed. The original loader performs validation.

`ReferenceTemperature` is in K, `ReferenceFugacity` in Pa and
`ReferenceFrameworkMass` in kg for the represented framework. `FrameworkMass`
remains an accepted alias. The kernel returns the mean particle count divided
by Avogadro's constant and that mass, giving mol/kg. An analytic
`PhysisorptionSites` model is not required.

Each cell uses its local gas temperature and component partial pressures. The
original fugacity/temperature reweighting is retained, including the
reference-to-target temperature factor. The optional energy column adds the
first-order Hamiltonian temperature correction; without it that contribution
is omitted. MPD does not infer an adsorption heat for the column energy
balance: retain explicit `HeatOfAdsorption` and thermal settings as appropriate.

## Breakthrough, cycles and particles

MPD changes equilibrium, independently of spatial discretization, convection
and time integration. BDF1, BDF, CVODE, SSPRK3 and SemiImplicitSSPRK3 remain
available on both axial paths. Short explicit runs still require suitable step
sizes for the selected spatial and particle stiffness.

The same equilibrium supplies LDF, fixed-bed spherical/radial-cylinder HSDM
with or without film resistance, and monolith LDF. Existing exclusions remain:
HSDM monoliths and coupled reactive particle transport are unavailable; general
reaction networks are excluded from cycles. For a cycle, use the same MPD block
with `SimulationType: "SwingAdsorption"` and `SwingAdsorptionPhases`.

The distribution is loaded during material construction and cached across
cells and temperature changes; the RHS does not reopen the file. To apply an
edited file, construct a new simulation. Its file content digest participates
in checkpoint compatibility, so an altered distribution cannot silently resume
an old state. Ordinary compatible cycle extensions remain supported.

## Lab, JSON and cluster execution

Lab offers MPD for shared gas breakthrough and cycles. Its Setup panel exposes
the original `MPDSettings` JSON; the main JSON editor and cycle designer retain
the same simulation object. Resolve `FileName` relative to the simulation's
directory, or use an absolute local path.

Cluster export includes each referenced MPD distribution under `mpd/` and
rewrites only the corresponding `FileName` paths. Duplicate references to one
source file are packaged once. The source configuration and distribution are
not modified. Exported `simulation.json` and `run.sh` can therefore execute the
same case on a cluster with the same build. Other external material/grid files
still require inline data before cluster export.

## Verification

`tests/process/test_mpd.cpp` compares shared equilibrium directly with the
original kernel using two adsorbing species, a carrier between them, correlated
occupancies, conditional energies and local temperatures. It also checks
mass conversion, LDF component conservation, copy behavior, invalid dimension
names and cached-file independence.

Forty tests in `tests/mpd_shared_test.py` compare both axial methods against a
synthetic, exactly equivalent competitive Langmuir model: all five integrators,
complete 60 s breakthrough trajectories, repeated cycle phases, both HSDM
geometries and film choices, monolith LDF, CLI/Python agreement, Lab JSON and
executed cluster bundles. They also exercise malformed data and compatible
versus changed-distribution restarts. Results are in `verification/mpd_shared/`.

Runnable synthetic inputs and distributions are in
`examples/new_examples/unified/{fdm,fvm}/gas_mpd_breakthrough/`. These checks verify integration
and conservation; they do not establish experimental agreement, mesh
convergence, cycle optimization or completed CSS campaigns.
