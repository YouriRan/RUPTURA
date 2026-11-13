# Matched shared-column examples

`fdm/` and `fvm/` contain the same physical inputs and schedules. Only
`SpatialDiscretization` differs within each matched pair. These examples cover
gas/liquid breakthrough and cycles, spherical/radial-cylinder HSDM with and
without film resistance, and one-, two- and four-column PSA recipes. The
two-column example includes pressure equalization. The 32 inputs also include
matched gas-phase MPD breakthrough cases with their distribution files.

Regenerate them from the bundled public configurations:

```bash
python tools/generate_unified_examples.py
build/ruptura examples/new_examples/unified/fdm/gas_breakthrough/simulation.json --validate
build/ruptura examples/new_examples/unified/fvm/gas_breakthrough/simulation.json --validate
```

All configurations appear in `manifest.txt`, with paths relative to the project
root. The same JSON can be imported into Lab, executed with the CLI, or loaded
with `ruptura.ProcessSimulation`. Lab's cluster bundle exports that configuration.

The recipes are editable demonstrations, not optimized processes. Input
validation does not establish complete breakthrough, CSS or literature agreement.
Verification uses separate short cases in `tests/shared_column_test.py`.

## Synthetic MPD example

```bash
build/ruptura examples/new_examples/unified/fdm/gas_mpd_breakthrough/simulation.json
build/ruptura examples/new_examples/unified/fvm/gas_mpd_breakthrough/simulation.json
```

Both use 300 K, 2.5 bar, feed fractions A/He/B = 0.2/0.2/0.6 and a clean
carrier-filled bed. `distribution.data` is synthetic single-site verification
data with probabilities 1/6, 2/6, 3/6 and zero for (A,B) = (0,0), (0,1),
(1,0), (1,1). With the illustrative framework mass 1/Avogadro kg, it is exactly
equivalent at 300 K to competitive Langmuir capacities 1 mol/kg and affinities
3e-5 Pa^-1 for A and 2e-5 Pa^-1 for B. It is not an experimental adsorbent.

The 60 s trajectories are checked against that analytic equivalent in
`tests/mpd_shared_test.py`. See `docs/MPD_SHARED_MODEL.md` for MPD input details,
cycle use, HSDM and portable cluster export.

All inputs in this collection select WENO3. See [HAL execution and collection](../README.md).
