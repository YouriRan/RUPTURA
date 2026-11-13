# Verification for the integration review build

| Check | Observed outcome |
|---|---|
| Build | Main C++23 executable + Python extension; standalone C++20 process library |
| Native tests | 12 groups passed, including unchanged donor gas/equilibrium reference values, chemistry, integrators, boundaries, liquid hydraulics and HSDM |
| Python complete suite | 126 tests passed, including the 115 Lab/plot checks |
| Main-style input checks | 227 passed; includes 27 MFI breakthrough + 162 MFI cycle recipes |
| Complete small simulations | 16 passed: four breakthrough, four cycles, four HSDM, two thermal, automatic saturation and native BDF |
| Liquid tank | Expected variable-volume balance, solute balance and positive pump work passed |
| Python restart | Two-cycle restart and uninterrupted run final states agree exactly in the tested liquid case |
| Original Marx T15 | Full 4,000 s completed |
| Shared FVM Marx T15 | 1 s startup passed at 101 nodes; full run reached 60 s in the 360 s limit. Full comparison remains incomplete |
| Haghpanah breakthrough | Paired 2013 cases completed 5,000 s; conditional 2014 case completed 7,000 s; plots retain discrepancies |
| HPC scripts | Bash syntax and case counts checked; no HAL submission performed |

`smoke_runs.json` records complete simulations, not input-only checks. The native
HSDM mean/shell invariant is also checked with nonzero radial profiles; recorded
trajectory discrepancies are at most 1.3e-13 mol/kg in the four examples.

Full raw rerun directories are excluded from the source archive. Selected outlet
CSVs, plots, actual settings/statuses, and reproducible verification scripts are
included. Run `tools/verify_shared_model.py --executable build/ruptura --smoke`
and optionally `--references`; reference runs may take substantially longer.
Use `tools/analyse_integrated_checks.py` after generating the raw rerun outputs.

These checks do not establish MPD/reaction/transported-pH parity in the FVM cycle
path, scientific equivalence to every donor equation, or literature cycle validation.
See `docs/INTEGRATION_STATUS.md` and `docs/HAGHPANAH_STATUS.md`.
