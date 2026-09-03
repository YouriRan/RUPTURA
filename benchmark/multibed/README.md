# Single-bed versus one-adsorbent multibed benchmark

This Google Benchmark executable compares 100 RK3 simulation steps for the
regular `Column` implementation and `MultibedColumn` configured with exactly
one adsorbent. A third benchmark runs both implementations and reports
`max_abs_error`, `max_rel_error`, and a `results_equal` counter. The latter is
one only when state, pressure, and interstitial-velocity results agree within
`1e-6` absolute or relative tolerance.

The benchmark runs the full 120-case product of these settings:

- total fluid components (including the carrier): 2, 3, 4, 5, and 6;
- energy balance: off and on;
- workload:
  - physisorption only;
  - direct chemisorption without surface/pore transport;
  - chemisorption with surface/pore transport;
- mixture prediction: IAST, SIAST, SEI, and SPI.

All cases use 50 grid points and two Langmuir physisorption sites per
adsorbate. Chemisorption cases use one General/Langmuir site per adsorbate.
The generated adsorbates have distinct affinities and capacities, so increasing
the component count exercises real mixture-prediction work rather than adding
inactive species.

The numeric `setting` argument is the index into the shared `settings` array in
`benchmark_cases.h`; every result also carries a descriptive label containing
all four dimensions. Use Google Benchmark's `--benchmark_filter` to select a
subset.

Configure and run it with:

```sh
cmake -S . -B build -DBUILD_BENCHMARKS=ON
cmake --build build --target multibed_single_adsorbent_benchmark -j
./build/benchmark/multibed/multibed_single_adsorbent_benchmark
```

Google Benchmark must be installed with a discoverable CMake package config.
Use `--benchmark_out=results.json --benchmark_out_format=json` to save results.
The repository's executed `../analysis.ipynb` reads
`performance_results.json`; reproduce that file with:

```sh
./build/benchmark/multibed/multibed_single_adsorbent_benchmark \
  --benchmark_min_time=0.05s \
  --benchmark_out=benchmark/multibed/performance_results.json \
  --benchmark_out_format=json
```

To run all 120 correctness simulations and write every final-state comparison
to CSV:

```sh
cmake --build build --target multibed_single_adsorbent_results -j
./build/benchmark/multibed/multibed_single_adsorbent_results
```

To refresh the correctness data used by `../analysis.ipynb`, pass its expected
output directory explicitly:

```sh
./build/benchmark/multibed/multibed_single_adsorbent_results \
  benchmark/multibed/simulations
```

This creates the root-level `multibed/` directory. Each of its 120 setting
directories contains `singlebed/` and `multibed/` subdirectories with normal
breakthrough `.data` output for the initial and final states, reproducible
`simulation.json` inputs, and a detailed `comparison.csv`. The top-level
`summary.csv` includes the component count, energy mode, workload,
chemisorption/transport flags, mixture method, maximum errors, and an equality
flag. The executable returns a nonzero status if any setting falls outside the
documented tolerance.
