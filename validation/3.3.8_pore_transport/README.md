# 3.3.8 Pore transport: CO₂ from air on PEI/fumed silica, four particle sizes

Validates the pore transport of RUPTURA 2.0 (film and pore reservoirs, `UsePoreSurfaceTransport`) with `General`
chemisorption (TODO V-PHY-2). Data: Goeppert et al., ChemSusChem 7 (2014) 1386, as digitized by Auton et al.,
Appl. Math. Model. 130 (2024) 827 (arXiv:2311.10161), whose dual-porosity model and fitted parameters the runs use
unchanged. All numbers and unit conversions are in `conversions.md`.

## Why this case and not Shi et al. (2024)

Shi et al., Chem. Eng. J. 496 (2024) 154090 (full text read 2026-11-09, without its SI) does not map onto RUPTURA:

- its amine kinetics are Toth-type, ∂q/∂t = k_T {c_p RT [1 − (q/q_sat)^τ]^{1/τ} − q/(b q_sat)} with τ(T) from
  0.25 to 0.29 between 298 and 318 K (Table 5), which `General` (integer powers of q_max − q) cannot write;
- the runs are non-isothermal with an Ergun pressure drop, and the correlations (Table S5) and the kinetic fits
  (Note S3) are only in the SI;
- its film and macropore resistances are lumped into one Biot-number coefficient (Eqs. 8–10), which would map, but
  the curves of its Fig. 5 are a raster image.

Auton's model is isothermal at constant velocity, its kinetics are `General` with orders 1, 2, 2, and its code and the
digitized data are public (https://github.com/aguareles/Intra-particle-diffusion).

## Case

CO₂ at 405 ppm in air, 298.15 K, 101325 Pa; column 8 mm × 0.1 m, 2.72 g sorbent, ρ_b = 541.1 kg/m³, ϕ = 0.48,
v = 0.2314 m/s (interstitial), `FixedVelocity`. Four sieve fractions, R = 1.25e-4, 1.875e-4, 5.5e-4, 8.5e-4 m.

- Kinetics in the pores: ∂q/∂t = k₊ c_p (q_max − q)² − k₋ q², `General` with orders 1, 2, 2 (Table 8 k₊, k₋, Table 5
  q_max). The site is a Langmuir isotherm with q_sat = q_max and b = 10³ Pa⁻¹ (as in `elfving/`), because the `General`
  law takes its free capacity from the site isotherm at the bulk gas, min(q*, q_max) − q. Since 2026-11-20 the paper
  writes the law this way too (TODO I-PHY-14: the code is kept).
- Exchange: Auton's single k_p (Table 8) as RUPTURA's film and pore resistances in series, 1/k_p = 1/k_f + R/(5 D_p),
  half each by default. Particle porosity: Auton's ε_p = 0.416 only multiplies the pore holdup, which is negligible
  (α = 4e-6), so RUPTURA's ε_p = 1 (C-PHY-1) does not matter here; the split study confirms it.
- The smallest size has β = 5.9e5 (k_p = 2.33e3 m/s), no intra-particle limit; it takes up directly from the bulk.

## Runs

`generate.py` writes 44 runs, `run.py` runs them with 8 workers (`python run.py --binary PATH`). CVODE, rtol 1e-7,
atol 1e-9 mol/m³ and 1e-8 mol/kg, `CVODEMaximumTimeStep` 0. Up to 200 cells the dense solver; above, SPGMR with the
`SorptionJacobian` preconditioner (TODO P-ITG-9, P-ITG-10). WENO3, Danckwerts inlet. The runs of 2026-11-20
(`run_2026-11-20_cvode.log`) used a binary built from the source of commit ba573b9 (2026-11-19) and frozen in a scratch
directory: 820 s of run time summed over the 44 runs (the longest 108 s; 2.5 min of wall time with 8 workers), and no
outlet curve decreases anywhere. Every output is byte-identical to that of the first runs, on 2026-11-09
(`run_2026-11-09_cvode.log`: 1146 s summed, with the binary of 2026-11-09 17:08).

| Study | Column | D_ax | Cells | What |
|---|---|---|---|---|
| grid | L_fit | 0 | 50–800 | convergence |
| split | L_fit | 0 | 100 | 99.9% of 1/k_p in the film, in the pores, or half each; direct uptake |
| beta | L_fit | 0 | 100 | smallest size: β = 100 and 1000 against direct uptake |
| affinity | L_fit | 0 | 100 | site affinity 10⁶ against 10³ Pa⁻¹ (0.5–1.7 mm) |
| wave | 70 𝓛 | 0 | 200 | 14 column profiles, convergence to Eq. (62) |
| experiment | 0.1 m and L_fit | Table 5 | 100 | against the data |

L_fit (0.0855, 0.0970, 0.0932, 0.1163 m) is the column whose stoichiometric time is that of Auton's fitted curve
(conversions.md §6); 𝓛 is their length scale (2.3e-3 to 1.35e-2 m).

## Results

`analysis.ipynb` writes the panels to `figures/`: `verification_*`, `travelling_wave_*`, `experiment_*` (one per sieve
fraction), and `grid_convergence.pdf`.

**Mapping checks** (RMS of the outlet c/c_in):

- The split of 1/k_p does not matter: film-only and pore-only against half and half differ by 1.1e-6 or less.
- β → ∞: β = 1000 is 2.2e-5 from direct uptake, β = 100 is 1.5e-4.
- Site affinity 10⁶ against 10³ Pa⁻¹: 1.3e-5.
- Grid: second order. 100 cells are 4.2e-4 or less from 800 cells, with t_1/2 within 13 s. The dense (≤ 200 cells) and
  SPGMR (400, 800) runs converge together.
- The pore transport is exercised: direct uptake from the bulk differs from it by 2.2e-3, 3.3e-2 (up to 0.12), and
  7.4e-3 for the three larger sizes.

**Verification against Auton's travelling wave, Eq. (62)** (no axial dispersion). Eq. (62) is the l = L/𝓛 → ∞ limit of
the same model. RUPTURA converges to it as the front develops: the time at which a profile of the wave column reaches C
deviates from Eq. (62) by

| Sieve fraction | β | l | C = 0.05 | 0.25 | 0.75 | 0.95 |
|---|---|---|---|---|---|---|
| < 0.25 mm | 5.9e5 | 62 | −0.009 ks | −0.004 ks | −0.014 ks | −0.29 ks |
| 0.25–0.5 mm | 9.7 | 57 | −0.023 | 0.001 | −0.038 | −0.81 |
| 0.5–1.7 mm | 0.56 | 52 | 0.000 | 0.005 | −0.075 | −1.35 |
| > 1.7 mm | 5.3 | 57 | −0.049 | −0.024 | −0.157 | −2.79 |

(the C → 1 tail converges last; a 1 m test column of the > 1.7 mm case at 1000 cells reached −1.7 ks at l = 65).
At the experimental length (l = 38, 22, 19, 8.6) the front of the larger particles has not reached that shape: at the
outlet of L_fit, RUPTURA is +0.8, +1.5, +9.1 ks off at C = 0.05 and −3.0, −4.2, −27 ks at C = 0.95 for the three larger
sizes, while the smallest is within 0.11 ks up to C = 0.75. So the β and k_p that Auton et al. fitted with Eq. (62) are
those of a developed wave; for the > 1.7 mm particles the full model with these values gives a much steeper upper half
than the data in a 0.1 m column.

**Experiment** (Table 5 axial dispersion, no fitting). RMSE in c/c_in on the data points, and the half times:

| Sieve fraction | t_1/2 data | RUPTURA 0.1 m | RUPTURA L_fit | RMSE 0.1 m | RMSE L_fit | RMSE Eq. (62) |
|---|---|---|---|---|---|---|
| < 0.25 mm | 38.6 ks | 45.6 ks | 38.8 ks | 0.309 | 0.042 | 0.041 |
| 0.25–0.5 mm | 43.8 | 45.9 | 44.5 | 0.111 | 0.064 | 0.047 |
| 0.5–1.7 mm | 43.5 | 47.6 | 44.3 | 0.064 | 0.025 | 0.030 |
| > 1.7 mm | 50.4 | 46.9 | 55.4 | 0.056 | 0.068 | 0.021 |

- With the published column the half times are 18%, 5%, 10% late and 7% early. Table 5's q_e in 0.1 m does not
  hold the CO₂ the data take up: the outlet curves give 1.418, 1.645, 1.667, and ≥ 2.028 mol/kg against q_e = 1.614,
  1.675, 1.700, 1.860 mol/kg (−12%, −2%, −2%, ≥ +9%). Auton's Eq. (62) does not see this, because it takes t_1/2 from
  the data. These q_e are the totals that Goeppert et al. measured (their Fig. 2: 71.0, 73.7, 74.8, 81.7 mg/g), each
  integrated with the sample mass and inlet of its own run. Goeppert et al. report neither of these, nor the bed length
  ("typically 3 g" in an 8 mm U-tube, 400–420 ppm; conversions.md §1 and §7). So the mismatch lies in the column, which
  is not fully known, and L_fit stands for it: its equivalent sample masses are 2.33, 2.64, 2.54, and 3.16 g.
- With L_fit, RUPTURA matches the data as well as Auton's own curve for the smallest (no transport limit) and the
  0.5–1.7 mm particles, and worse for > 1.7 mm, where the front is not developed (above).

Other findings: the printed Eq. (62) has −β(τ − τ_1/2) on the left; the authors' code, which reproduces their Fig. 9,
has +β. The 0.5–1.7 mm figure of Auton et al. leaves out 19 data points above C = 1.

## Paper

Since 2026-11-20 in `ruptura-2.0-paper/main.tex`: Results `sec:results_chemisorption` (Fig. `fig:goeppert`, the four
`experiment_*` panels), a row of the verification table, and SI `apx:goeppert` (Fig. `fig:goeppert_verification`:
`verification_0p5-1p7mm`, `grid_convergence`, `travelling_wave_0p5-1p7mm`, `travelling_wave_gt1p7mm`). The panels are
copied to `ruptura-2.0-paper/figures/goeppert/`; copy them again after rerunning the notebook. In the `experiment_*`
panels the travelling wave of Auton et al. is dotted, like the other authors' models in the paper.

## Open

- A refit of k_p with the full model for the two larger sizes would show how much of their k_p (not monotonic in
  the size) comes from the travelling-wave assumption; not done here.
