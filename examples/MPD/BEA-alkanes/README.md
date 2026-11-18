# BEA C7 alkanes: multinomial MPD versus Langmuir

This example validates macrostate particle-distribution (MPD) reweighting
against the competitive multisite Langmuir equation for the seven C7 alkane
isomers in `examples/BEA-alkanes-C7`.

The reference MPD is constructed as a product of two independent multinomial
site-family distributions. The primary family has three discrete sites and the
secondary family has one. Their framework-mass conversion gives shared
capacities of 1.2 and 0.4 mol/kg. The BEA affinities at 552 K are retained; only
the saturation capacities are made common within each family, as required by a
competitive multinomial site model. The four more highly branched isomers have
zero activity on the secondary family.

For each site family, the generator evaluates

```text
N_s! / [(N_s - sum_i n_i,s)! product_i(n_i,s!)]
    * product_i[(b_i,s f_ref)^(n_i,s)]
```

and normalizes it by `(1 + sum_i b_i,s f_ref)^N_s`. It then convolves the two
site-resolved distributions over every allocation with the same total particle
vector. The resulting dense, C-ordered file has `5^7 = 78,125` entries and is
consumed directly by Ruptura's MPD implementation.

The `langmuir` run uses `MixturePredictionMethod: "SCI"`, Ruptura's direct
segregated competitive-isotherm evaluation. It evaluates exactly

```text
q_i = sum_s q_sat,s b_i,s y_i f / (1 + sum_j b_j,s y_j f).
```

Run the complete example from the repository root:

```sh
python3 examples/MPD/BEA-alkanes/run_example.py
```

Use `--binary /path/to/ruptura` if the executable is not in `build/ruptura`.
The workflow regenerates `mpd/particle_distribution.data`, runs both inputs,
writes `comparison.csv`, optionally writes `comparison.png` when matplotlib is
available, and fails if the maximum scaled discrepancy exceeds `1e-10`.

This is an isothermal validation (`T = T_ref = 552 K`), so the conditional
energy term is identically inactive. A simulation-derived MPD can add the
conditional mean Hamiltonian as the second data column without changing the
input structure.

## Breakthrough and performance benchmark

Open `breakthrough_benchmark.ipynb` to run matched breakthrough simulations for
MPD and competitive Langmuir using CVODE and RungeKutta3, 20/50/100 grid cells,
and inlet velocities of 0.019/0.19 m/s. The notebook checkpoints all 24 runs,
plots outlet curves and grid convergence, compares wall and internal timings,
and tabulates interpolated curve discrepancies. Benchmark output is written to
the individual directories beneath `breakthrough`. Install the notebook
dependencies with `pip install -e '.[lab]'` when the optional plotting
environment is not already available.

All 24 `simulation.json` inputs are materialized so they can also be run without
the notebook. Launch up to eight Ruptura processes simultaneously with:

```sh
python3 examples/MPD/BEA-alkanes/breakthrough/run.py --jobs 8
```

Successful runs are reused by default; pass `--force` to recompute them. Use
`--dry-run` to list the selected directories, or provide exact directory names
after `run.py` to run only a subset. Each run writes `ruptura.log` and
`benchmark.json`; the launcher writes the aggregate `benchmark_results.csv`
consumed by the notebook. Regenerate the inputs after changing the benchmark
matrix with `breakthrough/generate_simulations.py`.
