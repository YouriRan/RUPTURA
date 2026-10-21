# Paper simulations

These self-contained inputs are bundled with Ruptura Lab, so the installed app
does not depend on a paper checkout or local study directories.

| Preset | Source input | Conditions |
| --- | --- | --- |
| T15 | `marx2/ternary/T15/simulation.json` | Marx ternary experiment T15, CO₂/N₂/H₂ 40/30/30, 298 K, 10 bar, energy balance, 4,000 s |
| BEA alkanes | `examples/BEA-alkanes-C7/breakthrough/simulation.json` | Seven C₇ isomers at 1% each in helium, 552 K, 1.2 bar, IAST, 1,200 s |
| Monolith, 400 ppm | `stampi-bombelli/dl_from_literature/monolith_400ppm_2mmol/simulation.json` | 2 mmol/s, 295.15 K, 144 square channels, 24,000 s |

The sources correspond to the Marx energy-balance validation, BEA benchmark,
and Stampi-Bombelli monolith figure in `ruptura-2.0-paper/main.tex`. The monolith
validation notebook selects `dl_from_literature` as its base model; its fitted
dispersion alternative is not used here.

Physical properties, geometry, adsorption, and transport parameters are retained.
T15's feed percentages are expressed as equivalent mole fractions, and display
names identify the experiments. CVODE uses SPGMR and returns its adaptive solution
every 2 s, saving every return for T15/BEA and every sixth return (12 s) for the
monoliths. Print intervals are 100 s for T15/BEA and 600 s for the monoliths.
The original T15 and monolith horizons are retained; BEA uses the paper's
1,200 s complete-breakthrough horizon. All inputs remain editable in the Lab.

The initial examples and paper presets have numbered names, while the generic
Quick start choices keep their original names; creating a simulation from any
preset assigns a new `untitled_i` name.
The `legacy/` inputs are retained only to recognize and remove untouched monolith
defaults from earlier workspaces; edited or previously run cases are preserved
as personal simulations.
