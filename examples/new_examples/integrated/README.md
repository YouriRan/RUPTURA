# Integrated main-style examples

The `breakthrough/` and `cycles/` folders contain gas/liquid monolith examples.
The four basic packed-bed cases are in [the matched collection](../unified/README.md):
`fvm/gas_breakthrough`, `fvm/liquid_breakthrough`, `fvm/gas_cycle` and
`fvm/liquid_cycle`. All 232 integrated inputs select WENO3. `particle/` covers sphere/radial-cylinder HSDM with finite film;
`thermal/` covers transient coolant. Parameters in these tutorials are illustrative.

The main JSON conventions and units are documented in
[`docs/SHARED_FVM_MANUAL.md`](../../../docs/SHARED_FVM_MANUAL.md). From a case folder:

```bash
./run --validate
./run
```

`campaigns/co2_ch4_breakthrough` contains the 27 supplied MFI breakthrough cases;
`campaigns/co2_ch4_cycles` contains the 162 supplied cycle recipes. Their JSON
inputs are self-contained. Material fit precision, separate equilibrium/energy
heats, molar feed, phase ramps, transient coolant and stream tags are retained.
The breakthrough auto-stop additionally requires loading equilibrium. The shared
energy formulation differs from the standalone breakthrough donor, so do not
substitute the converted cases for archived scientific results without comparing.

The campaign source hashes are stored in each input's `Provenance`. Original
screening notes are retained in `Notes`. These inputs have been migrated, not all
rerun to CSS. `tools/submit_integrated_hal.sh cycles` submits all 162 jobs at once
without a concurrency cap; use `breakthrough` for the 27-case array.

`validation/marx_T15_original` retains the legacy formulation with FixedVelocity
and explicit WENO3; the supplied source input remains under `marx2/ternary/T15`.
`marx_T15_fvm` is the explicit shared-model comparison. Haghpanah inputs retain
conditional/source-derived status; see
[`docs/HAGHPANAH_STATUS.md`](../../../docs/HAGHPANAH_STATUS.md).

`equipment/`, `chemistry/` and `moisture/` contain selected converted tutorials.
Use [the integration status](../../../docs/INTEGRATION_STATUS.md) to check the
remaining differences from the full original-main feature catalogue.

For the complete 296-case HAL array and collector, see [the new example index](../README.md).
