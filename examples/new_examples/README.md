# New examples

This collection has **296 unique simulation inputs: 92 breakthrough and 204 cycle
inputs**. All explicitly select `ConvectionScheme: WENO3` and
`Output.Plots: true`. Every case generates PNG and PDF figures automatically
after successful execution, including the retained legacy Marx comparison. WENO3 is the axial
convective reconstruction; FDM/FVM and time integration remain separate choices.
Near boundaries and material interfaces the implementation uses supported local
stencils and lower-order fallback. HSDM radial shell balances are unchanged.

| Collection | Inputs | Purpose |
|---|---:|---|
| `integrated/` | 232 | 43 tutorials/reference inputs plus the two MFI campaigns |
| `unified/` | 32 | 16 matched FDM/FVM pairs, including HSDM and MPD |
| `shared_chemisorption/` | 32 | 16 matched gas/liquid FDM/FVM chemistry pairs |

The four basic packed-bed inputs have one canonical home in `unified/fvm/`:

| Former integrated location | Canonical location relative to this folder |
|---|---|
| `breakthrough/gas_fixed_bed` | `unified/fvm/gas_breakthrough` |
| `breakthrough/liquid_fixed_bed` | `unified/fvm/liquid_breakthrough` |
| `cycles/gas_fixed_bed` | `unified/fvm/gas_cycle` |
| `cycles/liquid_fixed_bed` | `unified/fvm/liquid_cycle` |

`manifest.txt`, `cases.json` and `cases.csv` enumerate the canonical inputs using
project-relative paths. Refresh them after adding examples:

```bash
python3 tools/hal_examples.py refresh
```

The remaining FDM/FVM pairs are intentional. Only `SpatialDiscretization` differs
within each pair. Numerical agreement and convergence must be assessed from runs.
Tutorial parameters are illustrative; input validation is not literature validation.
Most tutorials run a few fixed cycles. The MFI cycle campaign requests CSS with
up to 200 cycles. A successful fixed-cycle run is reported separately from CSS.
Changing these inputs to WENO3 requires new calculations before comparing paper
results; historical verification files retain the settings used for those records.

## HAL execution

Build the portable executable with the project installer and CVODE support first.
From the project root, inspect and submit:

```bash
./tools/submit_new_examples_hal.sh all --dry-run
./tools/submit_new_examples_hal.sh all
```

`all` submits **one uncapped array, indices 0–295**. HAL controls actual scheduling.
No jobs are submitted by the dry run. Available selections are:

| Selection | Jobs |
|---|---:|
| `all` | 296 |
| `tutorials` | 43 |
| `mfi_breakthrough` | 27 |
| `mfi_cycles` | 162 |
| `unified` | 32 |
| `chemisorption` | 32 |
| `integrated` | 232 |

Defaults match the existing HAL scripts: partition `serial`, one node/task/CPU,
168 hours, excluding `c171`, and single-thread numerical libraries. Overrides:

```bash
./tools/submit_new_examples_hal.sh unified --partition serial --time 24:00:00
./tools/submit_new_examples_hal.sh all --executable /path/to/build/ruptura
./tools/submit_new_examples_hal.sh all --no-plots
./tools/submit_new_examples_hal.sh all --validate-only
```

`--exclude ''` removes the node exclusion. `--max-parallel N` adds a concurrency
cap only when requested. The environment variables `RUPTURA_EXECUTABLE`,
`RUPTURA_HAL_PARTITION`, `RUPTURA_HAL_TIME` and `RUPTURA_HAL_EXCLUDE` set defaults.
Native CLI options can follow `--`; example: `all -- --no-plots`.

Normal submission validates every selected input before calling `sbatch`.
`--skip-validation` is available for an already checked selection. Automatic plots
are enabled in every configuration. The HAL worker uses Matplotlib by default
and reports a plotting error rather than silently substituting an SVG-only result
when the plotting environment is unavailable. Keep the installer virtual
environment activated when submitting. `--no-plots` is an explicit opt-out. A new timestamped directory
under `hal_runs/` holds the frozen manifest, validation report and scheduler logs.
Each worker stages the input and required local files in its own case directory;
outputs never share a directory with a different case or a previous submission.
The source JSON and executable identity are recorded, and changed inputs/builds
are rejected after submission. Keep that build available until all tasks finish.

The compatibility launcher `tools/submit_integrated_hal.sh cycles` selects the
162 MFI cycle cases; `breakthrough` selects the 27 MFI breakthrough cases.

## Collecting results

Use the run directory printed by the submit command:

```bash
python3 tools/collect_new_examples_hal.py hal_runs/RUN_DIRECTORY
```

This writes `results_summary.json` and `results_summary.csv`, including runtime,
host, execution status, requested CSS, computed convergence, completed cycles,
all available final performance/balance metrics, PNG/PDF file counts and plot status. Missing jobs, failed jobs,
interrupted/running jobs, fixed-cycle completion, CSS and completed non-CSS cases
remain distinct. Figures are in each case’s `results/plots/` directory. The
retained legacy Marx comparison still writes its original text outputs and is
identified separately when it has no shared-model summary; a read-only plotter
automatically converts those results into PNG and PDF figures in the same plot
location. Its physical formulation is retained.

For a paper, use completed transients and verified CSS results as appropriate;
check component/energy balances and selected mesh/time refinements separately.
The catalog is a capability collection, not 296 experimentally validated cases.

## Replotting saved results

Both shared-model cases and the retained legacy breakthrough case support:

```bash
build/ruptura hal_runs/RUN_DIRECTORY/cases/CASE_PATH/simulation.json --plot-only
```

This reads recorded results and regenerates figures without running the solver.
Recompile with `bash install.sh --no-gui --toolchain` after updating to this
archive, then activate `.venv-ruptura` before submission.
