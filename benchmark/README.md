# Ruptura integrator benchmarks

This directory contains an integrator benchmark suite derived from:

- `examples/BEA-alkanes-C7/breakthrough/simulation.json`

The suite is the full cross product of 5 grid sizes, 7 time-step sizes,
3 axial-dispersion coefficients, and 3 integrators (315 simulations total).
Each simulation covers a common 0.1-second physical horizon, with the
step count computed as `ceil(0.1 s / TimeStep)`. This keeps the complete scaling
matrix practical while preserving the expected timestep-dependent work. The
selected axial-dispersion coefficient is applied to every component, including
the carrier gas. The `sirk3` directories are generated for future use, even
though that integrator is not currently runnable for these cases.

Regenerate the simulation directories and `manifest.csv` from the current base
inputs with:

```sh
python3 benchmark/generate.py
```

Run all currently runnable simulations from the repository root with up to
eight simultaneous `build/ruptura` subprocesses:

```sh
python3 benchmark/run.py
```

The runner writes each process's combined output to `run.log` in its simulation
directory and writes timing and exit-code information to `benchmark/results.csv`.
Use `--dry-run` to list selected simulations without starting them, `--case` to
select one base case, or `--jobs` to override the default concurrency. SIRK3 is
excluded by default; use `--include-sirk3` to include it once it is runnable.

Open `analysis.ipynb` to inspect live completion status, runtime scaling, paired
integrator speedups, breakthrough-curve agreement and breakthrough-time shifts.
The notebook can be rerun while a benchmark is in progress.

## Requested comparison figures

`figures/` contains the reproducible simulation matrix for the two new notebook
figures: the BEA/Marx CVODE-versus-RK3 2x2 comparison and the Marx-only
IAST/SIAST/SEI/SPI bar plots. It includes generation, repeated execution,
timing/error summarization, and reusable plotting helpers. See
`figures/README.md` for the exact matrix and commands.

## Single-bed versus multibed matrix

`multibed/` contains the Google Benchmark and correctness package for comparing
the regular column against a one-adsorbent multibed column. Its 120 settings
cover 2–6 components, energy balance on/off, physisorption and direct/transported
chemisorption workloads, and IAST/SIAST/SEI/SPI. See
`multibed/README.md` for build, filtering, and result-output instructions.
