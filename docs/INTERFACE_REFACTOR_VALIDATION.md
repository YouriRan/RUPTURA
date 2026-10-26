# Interface refactor validation

Validated on 2026-10-06 against the implementation in `Refactor.md`.

Interfaces are separate linear sections. Their lengths add to the total column length, and their grid counts add to the total interval count. The ten new `MultibedInterfaces` tests cover lengths, grids, fractions, integrated inventories, long interfaces, sharp boundaries, minimum grid counts, explicit nodes, single adsorbents, and invalid inputs. They also check that copied columns retain interface grid counts and that layered columns derive their length from their sections.

## Builds and tests

- CLI, Qt application, and C++ tests: `cmake --build build -j 8` succeeded.
- Python bindings: `cmake --build build-python --target ruptura_python -j 8` succeeded. A Python smoke check confirmed the new `adsorbent_interface_grid_points` fields, a 0.15 m column with 76 nodes, and initialization. Additional checks verified upstream ownership on 300 sharp-interface grids and bounded fractions at interface endpoints.
- Full suite: `ctest --test-dir build --output-on-failure -j 8` finished with 142 of 160 tests passing. All ten new interface tests and all 56 parameterized single-adsorbent equivalence tests passed. The same 18 tests failed before and after; no new failures were introduced. Existing failures concern chemisorption inputs and kinetics, swing and monolith inputs, unavailable BEA data, duplicate MPD benchmark settings, and stored gas/liquid regression values.
- Paper: `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` succeeded. The rendered geometry definitions, inventory equation, JSON example, figure caption, and spatial-discretization text were inspected.
- No validation JSON inputs use `ColumnSections`; only the two multibed examples needed migration.

## Breakthrough comparison

All four multibed examples were run before and after with their original physical and solver settings. Output was sampled every 0.1 s. The times below are the first interpolated crossings of normalized outlet partial pressure at 0.5. The gradient column's regenerated explicit grid samples both interfaces, so its changes also include discretization effects.

| Example | Component | Before t50 / s | After t50 / s | Change / % |
|---|---|---:|---:|---:|
| C6-alkanes-layered | nC6 | 433.440 | 462.103 | +6.61 |
| C6-alkanes-layered | 22DMB | 228.619 | 244.765 | +7.06 |
| CO2-C3H8-gradient | CO2 | 92.332 | 112.586 | +21.94 |
| CO2-C3H8-gradient | C3H8 | 125.990 | 146.891 | +16.59 |
| CO2-N2-AB-mixture | CO2 | 199.817 | 199.817 | +0.00 |
| CO2-N2-AB-mixture | N2 | 19.894 | 19.894 | +0.00 |
| CO2-N2-polishing | CO2 | 413.427 | 413.427 | +0.00 |
| CO2-N2-polishing | N2 | 9.494 | 9.494 | +0.00 |

The finite-interface columns break through later, consistent with their added length and material. The sharp-interface polishing and uniform-mixture outlet curves are identical before and after. All output values are finite.

Inputs, logs, raw output, and reproducible `generate.py`, `run.py` (eight workers), and `compare.py` scripts are in `/private/tmp/ruptura-interface-refactor/`. The final four-page comparison PDF is `output/pdf/multibed-breakthrough-comparison.pdf` inside that directory. Test and build logs are in `/private/tmp/ruptura-interface-*.log`. Original example output files were preserved.
