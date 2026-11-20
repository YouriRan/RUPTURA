# Reference settings of Hassan and Chung for T1-T24 (frozen)

The inputs of the Marx et al. (2014) ternary runs as used until 2026-11-19: IAST with the single-site Sips
isotherms of the reference settings, fluid fraction 0.403, particle density 456.6 kg/m3, `heatTransferGasSolid`
20 W/(m2 K), `heatCapacityGas` per unit mass, and an ambient temperature of 298 K.

They are kept unchanged because the C++ unit tests read them (the user's choice). Since 2026-11-20 every paper study of
T15 (verification, cycles, self-driving laboratory, benchmark, convection, boundary-condition, and sensitivity studies)
uses the inputs of the model of Marx et al. in `../ternary/` (TODO V-PHY-12).

The Marx validation itself uses `../ternary/T*/simulation.json`, written by `../generate.py` with the inputs of
the model of Marx et al. Do not regenerate the files here.
