# Scope of native verification

The current targeted native run is `current_native_checks.log`: all 13 process
test groups and the original-main Ergun-sphericity test pass. This covers the
implemented flow, equilibrium/transport, integration, heat/coolant, conservation,
input, runtime and output checks. `python_tests.log` records the complete Python
and Lab suite. `complete_runs.json` records complete executions of the nine new
examples; `input_validation.json` records the accepted input inventory.

A broader portable native run, retained as `native_tests_portable.log`, executed
143 tests: 114 passed and 29 failed. It is **not** reported as a fully passing
legacy suite. Twenty-eight failures depend on example or MPD fixture files absent
from the supplied source archive; several consequently parse an empty stream.

The remaining failure is
`ChemisorptionRate.AvramiMatchesDerivativeOfIntegratedModel`: the original rate
function returns 0.0303002110 while the test expects 0.0442383081. This revision
does not change the original chemisorption implementation. Its source SHA256
remains `f8869984561f6d8e40a32a378747b4c5946fd3a10aaff02bce8034f918316a66`.

The old default-boundary CVODE scalar regression has been replaced with a
comparison against independent RK3 integration under the new FixedVelocity
default; `default_velocity_check.log` records that additional passing check.
The portable scalar mass-derivative reference uses the same 1e-12 tolerance as
its independent formula checks, rather than a tighter 1e-14 absolute threshold.

The standalone GCC process suite also passed all 13 groups. Integrated CLI and
Python verification use the portable Zig CPU baseline; automatic host CPU
selection emitted an unsupported AVX-512 instruction during an initial local
build. No compiled binaries are included in the source ZIP. The installer already
selects the portable Zig baseline.

To reproduce the targeted native checks after building with tests enabled:

```bash
ctest --test-dir build -R '^process_|^Geometry.Ergun|^GasBreakthroughRegression.DefaultFixedVelocity' --output-on-failure
RUPTURA_CLI_EXECUTABLE="$PWD/build/ruptura" python -m pytest tests -q
python tools/verify_flow_controls.py --executable build/ruptura --validate-all
```

The new flow checks do not establish literature agreement. The Marx inputs have
been migrated away from the deleted boundary mode, and full T15/Haghpanah
scientific comparisons remain open.
