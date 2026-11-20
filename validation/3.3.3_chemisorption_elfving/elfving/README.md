# Elfving and Sainio (2021): humid-air CO2 adsorption

This directory maps the equilibrium and fixed-bed models from J. Elfving and T. Sainio, *Kinetic approach to modelling CO2 adsorption from humid air using amine-functionalized resin: Equilibrium isotherms and column dynamics*, Chemical Engineering Science 246 (2021) 116885, https://doi.org/10.1016/j.ces.2021.116885, to runnable Ruptura examples.

## What is reproduced

- `generate.py` evaluates the paper equations directly and writes CSV data plus `equilibrium/equilibrium_reproduction.svg` for the H2O GAB model (Fig. 2), the 5-parameter CO2 model (Fig. 3), the bounded 7-parameter model (Fig. S3), and the supplementary Toth model (Fig. S5). The relaxed 7-parameter result is written separately as `equilibrium/figure4_relaxed_reproduction.svg`.
- Each case directory is named `<model>_<temperature>c_<nominal H2O>pct_<flow>mlmin`.
- `five_25c_0p2pct_500mlmin`-`five_35c_2pct_500mlmin` (runs 1-5) reproduce the five 8.33e-6 m3/s fixed-bed fits made with the 5-parameter model in Table 4 and Figs. 5-6.
- `seven_bounded_25c_0p2pct_500mlmin`-`seven_bounded_35c_2pct_500mlmin` (runs 6-9) reproduce the four reported 8.33e-6 m3/s fits made with the bounded 7-parameter model in Table 4 and Fig. 7.
- `five_25c_2pct_100mlmin` and `seven_bounded_25c_2pct_100mlmin` (runs 10-11) reproduce the 1.67e-6 m3/s, 298 K, 0.020 H2O mole-fraction comparison in Fig. 8.
- `toth_ldf_25c_2pct_500mlmin` and `toth_ldf_35c_2pct_500mlmin` (runs 12-13) use the paper's supplementary Toth+LDF parameters directly for the 298 K and 308 K, 0.020 H2O mole-fraction comparison in Fig. S6.

## Paper-to-Ruptura mapping

The paper's CO2 model has two reaction pathways coupled through one shared free-amine balance, and the humid pathway depends explicitly on H2O partial pressure. Ruptura's current native chemisorption interface does not expose that cross-component kinetic state. The dynamic approximation uses a reversible `General` chemisorption rate proportional to `c_CO2 * (q_m/2 - q)^2 - k_d q`. The `q_m/2` ceiling represents two amine sites consumed per adsorbed CO2, while the case-specific forward and reverse coefficients reproduce the paper model's combined clean-bed feed rate and feed loading. This retains the requested squared free-surface dependence, although the single component-local loading cannot reproduce the paper's separate dry and humid products or their exact shared `q_m - 2 q_1 - q_2` balance. The Toth projections are retained for the Figure 4 equilibrium comparison and projection-error diagnostics.

The two supplementary Toth+LDF cases retain the Table S1 feed loading and Table 4 clean-bed rate when mapped to the same second-order `General` chemisorption form. They remain supplementary checks rather than distinct main-paper breakthrough panels.

Other mappings are:

- Bed length 0.01746 m, inner radius 0.0045 m, void fraction 0.375, particle diameter 0.6 mm, and particle density 720 kg/m3 come from Table 1.
- Interstitial velocities are calculated from standard volumetric flows of 8.33e-6 or 1.67e-6 m3/s at the stated bed cross-section and void fraction.
- Axial dispersion, H2O LDF, effective axial thermal conductivity, and bed-wall heat-transfer coefficients come from Table 4. H2O LDF values above 0.25 1/s are capped at 0.25 1/s in the runnable projection because the water front is already effectively instantaneous on the plotted timescale and the larger fitted values make the native energy-coupled system unnecessarily stiff. Both values are recorded in `paper_mapping.json`; remove `WATER_LDF_CAP` in `generate.py` for a literal high-stiffness run.
- Every case uses CVODE with the `SPGMR` Krylov linear solver (`CVODELinearSolver`).
- H2O uses Ruptura's native, temperature-dependent GAB isotherm with the five Table 2 parameters (`q_m,mono`, `C0`, `K0`, `Delta H_C`, `Delta H_K`). The GAB contribution and the second Langmuir empty-column correction are both represented as `FirstOrder` chemisorption sites with rates proportional to `(q_eq - q)`. They use the reported H2O LDF coefficient, subject to the existing 0.25 1/s runtime cap.

  In a Ruptura input, use `{"Type": "GAB", "Parameters": [q_mono, C0, K0, delta_H_C, delta_H_K]}` together with `"nonIsothermal": true` and `"HeatOfAdsorption": 1.0`, without a reference temperature. Loading is in mol/kg, the two prefactors are dimensionless, and both enthalpies are in J/mol. The component scale is then `exp(1 / (R T))`; GAB evaluates `C = C0 * scale^delta_H_C` and `K = K0 * scale^delta_H_K`.

- A large gas-solid heat-transfer coefficient approximates the paper's gas-solid thermal-equilibrium assumption. The wall is held near the initial temperature to represent the paper's prescribed wall temperature.
- The paper used 200 finite volumes. These examples use 30 grid cells for the 8.33e-6 m3/s cases and 15 for the long 1.67e-6 m3/s cases so the complete suite remains practical; increase `NumberOfGridPoints` in generated inputs for a grid-convergence study.

The paper publishes plotted experimental histories but not their underlying numerical samples. The main validator therefore uses model equations and reported landmarks. The notebook additionally uses explicitly labeled, pixel-digitized experimental controls; these do not claim the precision of the undisclosed source samples.

## Run everything

Build Ruptura in the repository root, then run:

```sh
python3 elfving/generate.py
python3 elfving/run.py
```

`generate.py` writes all 13 case directories' JSON inputs and recreates the equilibrium curves. `run.py` runs all 13 Ruptura cases and writes `validation_summary.csv` and `breakthrough_summary.svg`.

Open `elfving_numbered_examples.ipynb` for sections matching paper Figures 4-8. The notebook contains its complete equilibrium-model, fitting, parsing, tabulation, and Matplotlib plotting code. Each Figure 5-8 section contains separate 8 x 6 cells for normalized CO$_2$/H$_2$O breakthrough, midpoint solid temperature in kelvin, and grid-averaged CO$_2$ loading. The breakthrough plots compare Ruptura with experimental CO$_2$ controls recovered from the paper's red pixels; Figure 4 compares the equilibrium equations on a logarithmic CO$_2$-pressure scale, all three temperatures overlaid in one panel using light/normal/dark shades of the C0/C1/C2 humidity colors. Run `python3 digitize_breakthrough.py` only when the CSV controls need to be regenerated.

Run a single already-generated case directly, for example:

```sh
cd validation/3.3.3_chemisorption_elfving/elfving/five_25c_2pct_500mlmin && ../../../../build/ruptura
```

`generate.py` is the source of truth for repeated input values, including each case's directory name. Re-run it after changing a shared modelling assumption.

## Expected comparisons

- At 298 K, raising the H2O mole fraction from 0.002 to 0.020 increases the CO2 equilibrium loading and delays the midpoint breakthrough.
- H2O breaks through before CO2 and dominates the initial temperature peak.
- At 35 C the CO2 capacity is lower than at 25 C.
- At 1.67e-6 m3/s the first CO2 breakthrough occurs on the roughly 1e4 s scale reported for Fig. 8.
- The exact equilibrium evaluator shows that the bounded 7-parameter model fits more sharply than the 5-parameter model, while the relaxed 7-parameter model produces the paper's very high apparent reaction orders.
