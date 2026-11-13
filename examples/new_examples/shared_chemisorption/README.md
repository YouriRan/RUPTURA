# Shared gas/liquid chemisorption examples

These 32 synthetic examples cover FirstOrder, PseudoNth, Avrami, General,
Elovich and DualNth in gas/liquid columns with FDM/FVM. Avrami/DualNth also
include two-cycle adsorption/countercurrent-regeneration schedules. Only the
spatial selector differs between matched FDM/FVM pairs. Rates and material
properties are illustrative; these are not fitted materials or validated recipes.

Run `build/ruptura PATH/simulation.json --validate`, then the same command without
`--validate`. The same JSON imports into Lab or ProcessSimulation. General and
Elovich examples disable additional surface/pore transport to isolate their
kinetic law. See docs/SHARED_CHEMISORPTION.md for units and scope. Regenerate with
`python tools/generate_chemisorption_examples.py`.

All inputs in this collection select WENO3. See [HAL execution and collection](../README.md).
