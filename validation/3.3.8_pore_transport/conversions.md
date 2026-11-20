# Unit conversions: Goeppert et al. (2014) with the model of Auton et al. (2024)

Sources:

- Experiments: A. Goeppert, H. Zhang, M. Czaun, R. B. May, G. K. S. Prakash, G. A. Olah, S. R. Narayanan,
  "Easily regenerable solid adsorbents based on polyamines for carbon dioxide capture from the air",
  ChemSusChem 7 (2014) 1386–1397, doi:10.1002/cssc.201301114 (read 2026-11-20; FS-PEI-50, 50 wt% PEI of
  M_w 25000 on Aerosil 380 fumed silica, Figs. 1 and 2 and the Experimental Section).
- Model and parameters: L. C. Auton, M. Aguareles, A. Valverde, T. G. Myers, M. Calvo-Schwarzwalder, "An analytical
  investigation into solute transport and sorption via intra-particle diffusion in the dual-porosity limit",
  Appl. Math. Model. 130 (2024) 827–851, doi:10.1016/j.apm.2024.03.023; arXiv:2311.10161v1, Tables 5 and 8.
- Their MATLAB code and the digitized data: https://github.com/aguareles/Intra-particle-diffusion
  (`Intra-Particle-codes.zip`, `Figure_9_driver_Table_8.m`, `Figure_8_9_data.m`). Where Table 5 rounds, the numbers
  below are those of the code.

`generate.py` evaluates every expression below; the numbers here are its output.

## 1. Feed and column

The code of Auton et al. sets 405 ppm CO₂ at 25 °C and 1 atm, 335 mL/min, 2.72 g of sorbent in a column of 8 mm
diameter and 0.1 m length, and an inter-particle porosity ϕ = 0.48 (estimated by them from the column-to-particle
diameter ratio). The temperature and the humidity of the experiments are not given in their paper.

What Goeppert et al. report (Experimental Section and the caption of Fig. 1): 25 °C in an oil bath; 335 mL/min of air
with 400–420 ppm CO₂, mixed by mass flow controllers from 1900 ppm CO₂ in air and CO₂-free air, without a reference
state for the flow; "typically 3 g" of sorbent in a U-tube of 8 mm inner diameter, weighed after 3 h at 85 °C under
vacuum (weight loss 2–10%), and that weight is the basis of the capacities. They give neither the bed length nor the
mass, flow, or inlet concentration of each run, and no humidity (the gases come from cylinders). So the 0.1 m and
2.72 g of Auton et al. are their own choices (2.72 g is 3 g less 9%), and so is 405 ppm. The bypass readings at the
start of Fig. 1 and the plateaus of the curves (about 405, 390, 400, and ≥ 380 ppm) show that the inlet differed from
run to run. Particles below 0.25 mm built up 0.10 bar upstream of the bed; `FixedVelocity` neglects that.

$$c_\mathrm{in} = \frac{y\,p}{R T} = \frac{405\times10^{-6}\cdot 101325}{8.314463\cdot 298.15}
  = 0.016554\ \mathrm{mol\,m^{-3}}\quad(\text{Table 5: } 0.0166)$$

$$A = \frac{\pi d^2}{4} = 5.0265\times10^{-5}\ \mathrm{m^2},\qquad
  v = \frac{Q}{\phi A} = \frac{335\times10^{-6}/60}{0.48\cdot 5.0265\times10^{-5}} = 0.23141\ \mathrm{m\,s^{-1}}$$

$v$ is interstitial, as `ColumnEntranceVelocity`.

$$\rho_b = \frac{m}{A L} = \frac{2.72\times10^{-3}}{5.0265\times10^{-5}\cdot 0.1} = 541.13\ \mathrm{kg\,m^{-3}},\qquad
  \rho_p = \frac{\rho_b}{1-\phi} = 1040.6\ \mathrm{kg\,m^{-3}}$$

RUPTURA's `ParticleDensity` is the mass per particle volume, so that the bulk sink factor
$(1-\varepsilon)\rho_p/\varepsilon = \rho_b/\phi$ is that of Auton's Eq. (17). Auton's $\rho_b/(1-\phi)$ in the pore
balance is the same $\rho_p$.

```python
R = 8.31446261815324
c_in = 405e-6 * 101325.0 / (R * 298.15)            # 0.016554 mol/m3
area = math.pi * 8.0e-3**2 / 4.0                  # 5.0265e-5 m2
v = (335e-6 / 60.0) / (0.48 * area)               # 0.23141 m/s, interstitial
rho_b = 2.72e-3 / (area * 0.1)                    # 541.13 kg/m3
rho_p = rho_b / (1.0 - 0.48)                      # 1040.6 kg/m3 -> "ParticleDensity"
```

`FixedVelocity` keeps $v$ constant (trace CO₂, as in Auton's model) and computes no pressure drop, so the viscosity
(1.8e-5 Pa s) does not enter.

## 2. Axial dispersion

Table 5 gives $D$ with the unit "m/s²"; it is a diffusivity, m²/s, as in Auton's Eq. (2a) and their code:
0.6e-4, 0.9e-4, 2.6e-4, 4.1e-4 m²/s for the four sizes (`AxialDispersionCoefficient`). They take them from a chart in
Levenspiel (1999). Their travelling-wave solution, Eq. (62), neglects dispersion ($\mathrm{Pe}^{-1} = 0$ although
Table 8 gives 0.09–0.23), so the runs that are compared with Eq. (62) use $D = 0$; the runs against the data use
Table 5.

## 3. Particle sizes

Goeppert et al. sieved four fractions: < 0.25, 0.25–0.5, 0.5–1.7, and > 1.7 mm. Auton et al. use the mean diameter of
the two inner fractions and the bounds for the outer two: R = 1.25e-4, 1.875e-4, 5.5e-4, 8.5e-4 m.
`ParticleDiameter` = 2R (both at the top level and in `Geometry`). RUPTURA's pore-diffusion length is d_p/2 = R
(`src/geometry.cpp`).

## 4. Kinetics: Sips (a = 1, b = 2) as `General`

Auton's Eq. (18) with a = 1, b = 2,

$$\frac{\partial \bar m_p}{\partial t} = k_+\,\bar c_p\,(\bar m_\max - \bar m_p)^2 - k_-\,\bar m_p^2,$$

is the `General` rate with `PoreConcentrationOrder` 1, `CapacityOrder` 2, `DesorptionOrder` 2. Units: $\bar m$ in mol/kg,
$c_p$ in mol/m³, $k_+$ in m³ kg mol⁻² s⁻¹, $k_-$ in kg mol⁻¹ s⁻¹; all SI already.

| R / m | q_max / mol kg⁻¹ | k₊ / m³ kg mol⁻² s⁻¹ | k₋ / kg mol⁻¹ s⁻¹ |
|---|---|---|---|
| 1.25e-4 | 1.662 | 0.0348 | 5.13e-7 |
| 1.875e-4 | 1.725 | 0.0166 | 2.45e-7 |
| 5.5e-4 | 1.751 | 0.0145 | 2.14e-7 |
| 8.5e-4 | 1.915 | 0.00440 | 6.49e-8 |

The code multiplies $k_+$ by $(q^* - q)^2$, with $q^* = \min(q^*_\mathrm{site}, q_\max)$ and $q^*_\mathrm{site}$ the
site isotherm at the bulk gas, not by $(q_\max - q)^2$ as the paper's Eq. `general_reversible_chemisorption` writes
(`src/chemisorption.cpp`, `Chemisorption::rate`). To get Auton's law, the site is a Langmuir isotherm with
$q_\mathrm{sat} = q_\max$ and $b = 10^3$ Pa⁻¹, as in `elfving/`: at the feed, $b p = 4.1\times10^4$ and
$q^*_\mathrm{site} = q_\max(1 - 2.4\times10^{-5})$. The `affinity` study shows that $b = 10^6$ Pa⁻¹ gives the same
curve.

The equilibrium of the law at the feed, $k_+ c_\mathrm{in}(q_\max - q_e)^2 = k_- q_e^2$, gives Auton's Eq. (19):

$$\mu = \frac{q_\max}{q_e} = 1 + \sqrt{\frac{k_-}{k_+ c_\mathrm{in}}} = 1.0298\quad\text{for all sizes}$$

(Auton et al. fix μ = 1.02985 from the batch values 0.0759 and 0.0737 g/g and derive $k_-$ from it.)
$q_e = q_\max/\mu$ = 1.614, 1.675, 1.700, 1.860 mol/kg.

```python
mu = 1.0 + math.sqrt(k_minus / (k_plus * c_in))
q_e = q_max / mu
```

## 5. Bulk–particle exchange: k_p as RUPTURA's film and pore resistances

Auton's bulk and pore balances (Eqs. 15 and 17), with spheres ($|\partial\omega|/|\omega| = 3/R$):

$$\frac{\partial c}{\partial t} + v\frac{\partial c}{\partial z} - D\frac{\partial^2 c}{\partial z^2}
  = -\frac{1-\phi}{\phi}\,\frac{3}{R}\,k_p\,(c - c_p),\qquad
  \phi_p\frac{\partial c_p}{\partial t} = \frac{3}{R}\,k_p\,(c - c_p) - \rho_p\frac{\partial \bar m}{\partial t}.$$

RUPTURA with `UsePoreSurfaceTransport` has two reservoirs per particle, the surface $c_s$ and the pore $c_p$, both per
particle volume (`src/integrators/transport.cpp`):

$$\text{bulk sink} = \frac{1-\varepsilon}{\varepsilon}\,\frac{6}{d_p}\,k_f\,(c - c_s),\qquad
  \frac{\partial c_s}{\partial t} = \frac{6}{d_p}k_f(c - c_s) - \frac{15 D_p}{R^2}(c_s - c_p),\qquad
  \frac{\partial c_p}{\partial t} = \frac{15 D_p}{R^2}(c_s - c_p) - \rho_p\frac{\partial q}{\partial t}.$$

The bulk sink has Auton's form with $6/d_p = 3/R$. The holdups of the reservoirs are negligible: Auton's
$\alpha = \phi_p\,\mathrm{Da}/\varphi = 4\times10^{-6}$ (Table 8), and RUPTURA's reservoirs (no particle porosity, the
case $\varepsilon_p = 1$ of TODO C-PHY-1) have $\alpha/\phi_p = 10^{-5}$. With quasi-steady reservoirs the flux per
particle area is the same through the film and the pores, $k_f(c - c_s) = (5D_p/R)(c_s - c_p)$, so the two act as one
coefficient

$$\frac{1}{k_p} = \frac{1}{k_f} + \frac{R}{5 D_p}.$$

Any split of $1/k_p$ gives the same column; `generate.py` uses half and half by default,

$$k_f = 2k_p,\qquad D_p = \frac{2 R k_p}{5},$$

and the `split` study puts 99.9% of the resistance in the film or in the pores (`SPLITS`, share $s$ of the film:
$k_f = k_p/s$, $D_p = R k_p/(5(1-s))$).

| R / m | k_p / m s⁻¹ (Table 8) | β | k_f / m s⁻¹ | D_p / m² s⁻¹ |
|---|---|---|---|---|
| 1.25e-4 | 2.33e3 | 5.93e5 | direct uptake | direct uptake |
| 1.875e-4 | 0.0294 | 9.71 | 0.0588 | 2.205e-6 |
| 5.5e-4 | 0.0045 | 0.563 | 0.0090 | 9.90e-7 |
| 8.5e-4 | 0.0237 | 5.28 | 0.0474 | 8.058e-6 |

For the smallest size the published $k_p$ makes the film and pore exchange about $10^8$ s⁻¹, which CVODE does not get
through in 2 min of CPU. Auton et al. read its β = O(10⁵) as no intra-particle limit, so the run takes up directly from
the bulk (`UsePoreSurfaceTransport` false, the limit β → ∞); the `beta` study runs β = 100 and 1000 with the
equal split ($k_p$ scaled by β/5.93e5) against it.

```python
share = 0.5                                      # share of 1/k_p in the film
k_f = k_p / share                                # "FilmMassTransferCoefficient"
D_p = R * k_p / (5.0 * (1.0 - share))            # "PoreDiffusivity"
k_p_beta = k_p * beta_target / beta_published    # beta study, smallest size
```

## 6. Dimensionless groups and the column length of Auton's fit

Auton's Eqs. (21) and (24):

$$T = \frac{1}{q_e k_+ c_\mathrm{in}},\qquad \mathcal{L} = \frac{v T \phi c_\mathrm{in}}{\rho_b q_e},\qquad
  \mathrm{Da} = \frac{\mathcal{L}}{vT},\qquad \beta = k_p\frac{3}{R}\,T\,\mathrm{Da}\,\frac{1-\phi}{\phi}.$$

`generate.py` reproduces Table 8 to 0.5% (T = 1076, 2173, 2450, 7383 s; β = 5.93e5, 9.71, 0.563, 5.28;
$\mathcal{L}$ = 2.265e-3, 4.408e-3, 4.897e-3, 1.349e-2 m).

Eq. (62), the travelling wave for a = 1, b = 2, with α = Da = Pe⁻¹ = 0, as in `Analytical_calculator.m`:

$$t - t_{1/2} = \frac{T}{\beta}\bigl[G(C) - G(\tfrac12)\bigr],\qquad
  G(C) = \Bigl(1 + \frac{\beta}{\mu^2}\Bigr)\ln C - \frac{\beta + (\mu-1)^2}{\mu^2 - 1}\ln(1-C)
       + \frac{\beta + \mu^2(\mu-1)^2}{\mu^2(\mu^2-1)}\ln(\mu^2 - C).$$

The printed Eq. (62) has $-\beta(\tau - \tau_{1/2})$ on the left; with that sign the curve runs backwards in time
(RMS 2.7–15 h from their Fig. 9), while the code's sign reproduces the Fig. 9 curves extracted from the PDF to 0.02 h.

In Eq. (62) the half time $t_{1/2}$ is taken from the data, not from a mass balance. A front that moves at
$w = \mathcal{L}/T = v\phi c_\mathrm{in}/(\rho_b q_e)$ has its stoichiometric time $t_\mathrm{st}$ at the mean of
$t(C)$ over $0 < C < 1$,

$$t_\mathrm{st} = t_{1/2} + \delta,\qquad
  \delta = \frac{T}{\beta}\Bigl[-A + B + E\bigl(\mu^2\ln\mu^2 - (\mu^2-1)\ln(\mu^2-1) - 1\bigr) - G(\tfrac12)\Bigr],$$

with $A$, $B$, $E$ the coefficients of $G$ ($\int_0^1 \ln C\,\mathrm{d}C = \int_0^1\ln(1-C)\,\mathrm{d}C = -1$). The column
length that has this $t_\mathrm{st}$ is $L_\mathrm{fit} = w\,t_\mathrm{st}$ (the gas holdup changes $w$ by 1e-5):

| R / m | t₁/₂ / h (Table 5) | δ / h | L_fit / m | t_st at L = 0.1 m / h |
|---|---|---|---|---|
| 1.25e-4 | 10.730 | 0.548 | 0.0855 | 13.19 |
| 1.875e-4 | 12.188 | 1.088 | 0.0970 | 13.69 |
| 5.5e-4 | 12.080 | 0.879 | 0.0932 | 13.90 |
| 8.5e-4 | 14.028 | 3.644 | 0.1163 | 15.20 |

So the curves of Auton's Fig. 9 do not hold the q_e of Table 5 in a 0.1 m column: their capacities are −14%, −3%,
−7%, and +16% of $q_e$ (= L_fit/0.1 − 1). The runs compared with Eq. (62) use L_fit; the runs against the data use both.

```python
T = 1.0 / (q_e * k_plus * c_in)
L_ref = v * T * 0.48 * c_in / (rho_b * q_e)
beta = k_p * (3.0 / R) * T * (L_ref / (v * T)) * (1.0 - 0.48) / 0.48
L_fit = (t_half + delta) * L_ref / T
```

## 7. Reference data

`reference_digitized/extract_reference.py` writes the four arrays of `Figure_8_9_data.m` (time in h, outlet CO₂ in ppm,
digitized from Goeppert et al. by A. Valverde) to CSV. The notebook normalizes them as Auton et al. do,
$C = \mathrm{ppm}/405$, and converts h to s. The arrays match the markers of Auton's Fig. 9 (read from the vector PDF)
to 0.001 h and 1e-4 in C, except that the 0.5–1.7 mm figure leaves out the 19 points above C = 1 (up to 1.045, between
38 and 44.5 h).

The capacity held by a data set follows from the outlet alone,

$$q_\mathrm{data} = \frac{Q\,c_\mathrm{in}}{m}\int_0^{t_\mathrm{end}}(1 - C)\,\mathrm{d}t,$$

a lower bound when $C(t_\mathrm{end}) < 1$ (the notebook prints it next to $q_e$).

The $q_e$ of Table 5 are the total uptakes that Goeppert et al. measured: their Fig. 2 (a raster figure, read by eye;
the 0.25–0.5 mm value is the 73.7 mg/g of their Table 1) gives about 71.0, 73.7, 74.8, and 81.7 mg/g, i.e. 1.613,
1.675, 1.700, and 1.856 mol/kg with M = 44.0095 g/mol, within 0.2% of $q_e$ = 1.614, 1.675, 1.700, 1.860 mol/kg.
Goeppert et al. integrated each run with its own sample mass and inlet, which they do not report, while
$q_\mathrm{data}$ above uses 2.72 g, 405 ppm, and 335 mL/min at 25 °C for all four. That is why the data curves hold
−12%, −2%, −2%, and ≥ +9% of $q_e$ at 0.1 m. The equivalent sample masses of the $L_\mathrm{fit}$ columns,
$2.72\,\mathrm{g}\cdot L_\mathrm{fit}/0.1\,\mathrm{m}$, are 2.33, 2.64, 2.54, and 3.16 g (Goeppert: "typically 3 g").
Their Fig. 2 also gives the uptake before the outlet reaches 10 ppm: about 57.6, 51.8, 41.3, and 24.0 mg/g.
