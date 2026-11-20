# Conversions for the Marx et al. (2014) ternary runs

`generate.py` writes `ternary/T1-T24/simulation.json` with the inputs of the model of Marx et al. (2014). The
feed flow rate, pressures, and temperatures are converted as in
`../../3.4.1_boundary_conditions_marx_10bar/conversions.md` (sections 1 and 2):
`ColumnEntranceVelocity` is $v = Q/(A\,\varepsilon_t)$ with $\varepsilon_t = 0.742096$. Sections 1-4 document the
gas–wall heat-transfer coefficient `heatTransferGasWall` (`heat_transfer.py`, imported by `generate.py`); sections
5-9 the gas heat capacity, the total-porosity mapping, the isotherms, the gas–solid coefficient, and the ambient.

Until 2026-11-19 the inputs followed the reference settings of Hassan & Chung (IAST, their Sips parameters,
$\varepsilon_b$ as fluid fraction, $h_{fs} = 20$ W m⁻² K⁻¹, ambient 298 K). Those inputs are frozen in
`reference_settings_hassan_chung/` because the tests and every other T15 study read them; the parts of sections 5
and 6 marked "reference settings" describe them.

## 1. Correlation of Marx et al.

Marx et al. (2014), Eq. 3, with the parameters of Casas et al. (2012) (Marx et al., Table 6):

$$\mathrm{Nu} = \frac{h_L d_i}{K_L} = \eta_1\,\mathrm{Re}_0^{\eta_2}\exp\!\left(-\frac{6 d_p}{d_i}\right),
\qquad \eta_1 = 41.13,\quad \eta_2 = 0.32,$$

$$\mathrm{Re}_0 = \frac{\rho_0\,u\,d_p}{\mu}.$$

- $\rho_0$ is the initial fluid density in the column. The column is filled with helium before the feed is
  switched on, so $\rho_0 = P_\mathrm{in} M_\mathrm{He} / (R\,T_\mathrm{in})$ and $\mathrm{Re}_0$ is constant
  for each experiment (Casas, Diss. ETH 20681, 2012, Sect. 6.4.2; Marx, Diss. ETH 22996, 2015, Eq. 3.3:
  "initial Reynolds number").
- $u$ is the superficial feed velocity, $u = Q/A = v\,\varepsilon_t$ (Marx et al., Table 3, write the mass
  balance with $\partial(uc_i)/\partial z$ and the Ergun equation with $\varepsilon_b^3$).
- $\mu = 1.9\times10^{-5}$ Pa s, $d_p = 3\times10^{-3}$ m, $d_i = 0.025$ m (Marx et al., Table 4).
- $h_L$ is the bed–wall coefficient of the single-temperature model of Marx et al.; it is used as
  RUPTURA's gas–wall coefficient $h_{fw}$ (`heatTransferGasWall`). The wall–ambient coefficient stays
  $h_w = 5$ W m⁻² K⁻¹ (Marx et al., Table 6).

## 2. Axial thermal conductivity $K_L$

Marx et al. use one representative value of $K_L$ for all experiments (Marx, Diss. 2015, Sect. 3.3.2) but
do not report it. Casas (2012, Sect. 3.3.4) defines $K_L = D_L C_g$. The value is recovered from the
correlation curve in Fig. 6 of Marx et al., which plots $h_L$ against $\mathrm{Re}_0$ (0–25) from Eq. 3.
The curve was digitized from a 600 dpi rendering of page 504 (page 12 of `s10450-013-9593-5-1.pdf`),
calibrated on the dotted grid lines (20, 40, 60, 80 W m⁻² K⁻¹ and $\mathrm{Re}_0$ = 5, 10, 15, 20):

| Fit to 1194 curve points ($1 < \mathrm{Re}_0 < 24.5$) | Value |
|---|---|
| $h_L = C\,\mathrm{Re}_0^{0.32}$ | $C = 32.0$ W m⁻² K⁻¹, rms residual 0.06 W m⁻² K⁻¹ |
| free exponent | 0.3205 (Eq. 3: 0.32) |

$$K_L = \frac{C\,d_i}{\eta_1 \exp(-6 d_p/d_i)} = \frac{32.0 \times 0.025}{41.13 \times 0.4868} = 0.0400\ \mathrm{W\,m^{-1}\,K^{-1}}.$$

The fitted points in Fig. 6 (B31, B8, B30, B22, B5) lie within ±35% of $\mathrm{Re}_0$ computed this way;
Marx et al. fitted the feed velocity of every experiment, which explains the spread.

## 3. Resulting values

```python
import math
R, M_he, mu, d_p, d_i, K_L = 8.31446261815324, 0.004, 1.9e-5, 3.0e-3, 0.025, 0.0400
def h_L(P_in, T_in, u):
    reynolds = P_in * M_he / (R * T_in) * u * d_p / mu
    return 41.13 * reynolds**0.32 * math.exp(-6.0 * d_p / d_i) * K_L / d_i
h_L(1.0e6, 298.0, 0.02037)   # T15: Re_0 = 5.19, h_L = 54.3 W/(m2 K)
```

| Runs | $P_\mathrm{in}$ / Pa | $T_\mathrm{in}$ / K | $Q$ / m³ s⁻¹ | $\mathrm{Re}_0$ | $h_L$ / W m⁻² K⁻¹ |
|---|---|---|---|---|---|
| T1, T13 | 1e5 | 298 | 1e-5 | 0.519 | 25.97 |
| T2, T14 | 5e5 | 298 | 1e-5 | 2.596 | 43.47 |
| T3, T15 | 1e6 | 298 | 1e-5 | 5.193 | 54.26 |
| T4, T16 | 1e6 | 298 | 2e-5 | 10.386 | 67.74 |
| T5, T17 | 1e6 | 298 | 3e-5 | 15.579 | 77.13 |
| T6, T18 | 2e6 | 298 | 1e-5 | 10.386 | 67.74 |
| T7, T19 | 1e5 | 318 | 1e-5 | 0.487 | 25.44 |
| T8, T20 | 5e5 | 318 | 1e-5 | 2.433 | 42.58 |
| T9, T21 | 1e6 | 318 | 1e-5 | 4.866 | 53.15 |
| T10, T22 | 1e6 | 318 | 2e-5 | 9.733 | 66.35 |
| T11, T23 | 1e6 | 318 | 3e-5 | 14.599 | 75.54 |
| T12, T24 | 2e6 | 318 | 1e-5 | 9.733 | 66.35 |

Before 2026-10-20, every run used 25.6 W m⁻² K⁻¹, the Hassan & Chung value for T21.

## 4. Mapping onto RUPTURA's energy balance

RUPTURA resolves the gas, adsorbent, and wall temperatures. The wall balance receives
$h_{fw}\,\pi d_i\,(T_g - T_w)$ per unit length, as in Marx et al. The gas balance, written per unit fluid
volume, loses $h_{fw}\,4/(\varepsilon_b d_i)\,(T_g - T_w)$ per unit fluid volume, i.e. the same heat per unit
length (`src/geometry.cpp`, `fluidWallPerFluidVolume`). With $h_{fw} = h_L$, both sides match Marx et al.

Until 2026-10-26 the gas balance used $4/d_i$ per unit fluid volume and so removed only $\varepsilon_b$ times the
heat the wall gained; the runs before that date were made with this version.

## 5. Gas heat capacity

**Marx inputs (`ternary/`).** Every component carries the molar heat capacity of Marx et al. (2014), Table 4,
`GasHeatCapacity` = 29.51 J mol⁻¹ K⁻¹ (helium included), and `heatCapacityGas` is removed. RUPTURA then forms
$\rho_g C_{p,g} = \sum_j c_j C_{p,j} = c_{\mathrm{tot}}\,C_{p,g}^{\mathrm{mol}}$ for every composition, which is the
$C_g$ (per unit volume) of Marx et al.

**Reference settings (`reference_settings_hassan_chung/`).** RUPTURA's gas energy balance uses the volumetric heat
capacity $\varepsilon\,\rho_g C_{p,g}$ with the mass density $\rho_g = \sum_j M_j c_j$, so `heatCapacityGas` is a
specific heat capacity in J kg⁻¹ K⁻¹. The molar value was used directly until 2026-10-20, which gave the gas
$1/\overline{M}$ (about 38) times too little heat capacity. This hardly changed the temperatures, because the adsorbent
and the wall hold almost all of the heat capacity, but it made the gas energy balance stiff: explicit RK3 failed for
T15 at $\Delta t = 10^{-3}$ s (it ran at $7\times10^{-4}$ s), and runs at $3\times10^{-2}$ s with the value below.
`heat_capacity.py` writes the specific heat capacity of the feed,

$$C_{p,g} = \frac{C_{p,g}^{\mathrm{mol}}}{\overline{M}},\qquad \overline{M} = \sum_j y_{j,\mathrm{in}} M_j,$$

```python
M = {"CO2": 0.044, "N2": 0.028, "H2": 0.002016}           # kg/mol, as in the inputs
def specific_heat_capacity(y):                           # feed mole fractions
    return 29.51 / sum(y[j] * M[j] for j in y)
specific_heat_capacity({"CO2": 0.15, "N2": 0.80, "H2": 0.05})  # T1-T12:  M = 0.029101 kg/mol, 1014.1 J/(kg K)
specific_heat_capacity({"CO2": 0.40, "N2": 0.30, "H2": 0.30})  # T13-T24: M = 0.026605 kg/mol, 1109.2 J/(kg K)
```

In the helium that fills the column at the start, this gives about 1/7 of the feed value per mole; that gas is
displaced within the first column volume.

## 6. Total porosity: `ColumnVoidFraction`, `ParticleDensity`, `ColumnEntranceVelocity`

Marx et al. (2014), Table 3, write the mass and energy balances with the total porosity $\varepsilon_t$, the
superficial velocity $u = Q/A$, and the bulk density $\rho_b = 502$ kg m⁻³ (Table 4):

$$\varepsilon_t\frac{\partial c}{\partial t} + \frac{\partial (u c)}{\partial z} = -\rho_b\frac{\partial q}{\partial t},
\qquad
\left(\varepsilon_t C_g + \rho_b C_s + \rho_b C_{\mathrm{ads}}\right)\frac{\partial T}{\partial t} + \dots$$

**Marx inputs (`ternary/`).** RUPTURA writes the balances per unit fluid volume with the fluid fraction
$\varepsilon$, the interstitial velocity $v$, and the adsorbent $(1-\varepsilon)\rho_s$ per unit column volume.
`generate.py` sets

$$\varepsilon = \varepsilon_t = 0.742096,\qquad \rho_s = \frac{\rho_b}{1-\varepsilon_t} = 1946.46\ \mathrm{kg\,m^{-3}},
\qquad v = \frac{Q}{A\,\varepsilon_t},$$

so that the adsorbent per unit column volume is $\rho_b$, the gas and solid heat capacities per unit column volume are
$\varepsilon_t C_g$ and $\rho_b C_s$, and the gas–wall heat per unit length is that of section 4 (it does not depend on
$\varepsilon$). Not in RUPTURA: the adsorbed-phase heat capacity $\rho_b C_{\mathrm{ads}}$.

```python
bulk_density, total_porosity = 502.0, 0.403 + (1.0 - 0.403) * 0.568   # 0.742096
particle_density = bulk_density / (1.0 - total_porosity)               # 1946.46 kg/m3
```

Two terms also see $\varepsilon$:

- the Ergun pressure drop, now evaluated with $\varepsilon_t$ instead of $\varepsilon_b = 0.403$. With the Ergun
  equation of Marx et al. (Table 3, $\varepsilon_b$, $\mu = 1.9\times10^{-5}$ Pa s) it is below
  $4.6\times10^{-4}$ of the inlet pressure in every run (T1, 1 bar), so the change does not matter here;
- the external particle surface per unit column volume, $(1-\varepsilon)\,6/d_p$; see section 8.

**Reference settings (`reference_settings_hassan_chung/`).** These keep $\varepsilon = \varepsilon_b$ and
$v = Q/(A\varepsilon_t)$ and match the adsorbent per unit fluid volume, $\rho_p(1-\varepsilon_b)/\varepsilon_b =
\rho_b/\varepsilon_t$:

$$\rho_p = \frac{\rho_b\,\varepsilon_b}{(1-\varepsilon_b)\,\varepsilon_t}
= \frac{502 \times 0.403}{0.597 \times 0.742096} = 456.6\ \mathrm{kg\,m^{-3}}$$

(`particle_density.py`; 502 was used until 2026-10-26, 10% more adsorbent than the reference model).

## 7. Isotherms

Marx et al. (2014), Table 5 (Schell et al. 2011), Sips isotherm per component,

$$q_i = q_{\mathrm{sat},i}\frac{(K_i p_i)^{s_i}}{1 + (K_i p_i)^{s_i}},\quad
q_{\mathrm{sat},i} = a_i e^{-b_i/RT},\quad K_i = A_i e^{-B_i/RT},\quad
s_i = \alpha_i \arctan[\beta_i (T - T_{\mathrm{ref},i})] + s_{\mathrm{ref},i},$$

and the multicomponent Sips equation for the mixture (Sect. 4.2.1), i.e. RUPTURA's `SCI` with Sips sites,
$q_i = q_{\mathrm{sat},i}(K_ip_i)^{s_i}/[1 + \sum_j (K_jp_j)^{s_j}]$. RUPTURA's Sips site is
$q_{\mathrm{sat}}(kp)^{1/\nu}/[1+(kp)^{1/\nu}]$, so (`sips_site()` in `generate.py`):

| Marx | RUPTURA `Parameters` | Law |
|---|---|---|
| $a$ [mmol/g = mol/kg], $b$ [J/mol] | 1: `{"Value": a, "Law": "VantHoff", "Energy": -b}` (a number if $b = 0$) | uncentred van 't Hoff, no `ReferenceTemperature` |
| $A$ [MPa⁻¹], $B$ [J/mol] | 2: `{"Value": A * 1e-6, "Law": "VantHoff", "Energy": -B}` [Pa⁻¹] | uncentred van 't Hoff |
| $s_{\mathrm{ref}}$, $\alpha$, $\beta$ [K⁻¹], $T_{\mathrm{ref}}$ [K] | 3: `{"Value": 1/s_ref, "Law": "Arctan", "Alpha": α, "Beta": β, "ReferenceTemperature": T_ref}` (CO₂); `1/s_ref` (N₂, H₂, $\alpha = 0$) | $1/\nu = 1/\nu_{\mathrm{ref}} + \alpha\arctan[\beta(T-T_{\mathrm{ref}})]$ |

The component `ReferenceTemperature` is removed, so the van 't Hoff laws stay uncentred. Because every law has its own
energy, `HeatOfAdsorption` (= $-\Delta H$ of Table 5: 26000, 15600, 9800 J/mol) enters only the energy balance.
Isotherms are evaluated at the adsorbent temperature.

## 8. Gas–solid heat-transfer coefficient `heatTransferGasSolid`

Marx et al. assume gas and adsorbent at one temperature; RUPTURA resolves both. The coefficient follows Wakao,
Kaguei & Funazkri, Chem. Eng. Sci. 34 (1979) 325,

$$\mathrm{Nu} = \frac{h_{fs} d_p}{\lambda_g} = 2 + 1.1\,\mathrm{Pr}^{1/3}\mathrm{Re}_p^{0.6},\qquad
\mathrm{Re}_p = \frac{\rho_g u d_p}{\mu},\qquad \mathrm{Pr} = \frac{C_{p,g}^{\mathrm{mol}}}{\overline{M}}\frac{\mu}{\lambda_g},$$

with the feed at $P_{\mathrm{in}}$, $T_{\mathrm{in}}$ ($\rho_g = P\overline{M}/RT$), $u = Q/A$, $\mu = 1.92\times10^{-5}$ Pa s
(the input), $\lambda_g = 0.028$ W m⁻¹ K⁻¹, $d_p = 3$ mm (`gas_solid_heat_transfer()` in `generate.py`). RUPTURA's
transfer area per unit column volume is $(1-\varepsilon)\,6/d_p$; with $\varepsilon = \varepsilon_t$ the input is scaled
to keep the packed-bed area $(1-\varepsilon_b)\,6/d_p$:

$$h_{fs}^{\mathrm{input}} = h_{fs}\,\frac{1-\varepsilon_b}{1-\varepsilon_t} = 2.3148\,h_{fs}.$$

| Run | $h_{fs}$ / W m⁻² K⁻¹ | input |
|---|---|---|
| T1, T13 | 38.7, 38.3 | 89.66, 88.56 |
| T15 | 96.7 | 223.75 |
| T5, T17 (highest) | 173.1, 169.4 | 400.68, 392.24 |
| T19 (lowest) | 37.5 | 86.83 |

`python generate.py --dry-run` prints all 24 input values. The reference settings use 20 W m⁻² K⁻¹ for every run.

## 9. Ambient temperature and boundary condition

The column of Marx et al. is heated to the experimental temperature and held there by an electrical heater controlled
by a thermocouple on the outer wall, so `Temperature` (ambient and initial) equals `InfluxTemperature` (298 or 318 K).
The reference settings use 298 K for all runs, which cooled the 318 K columns below the measured temperatures once
the heat of adsorption had left. All runs use `InletPressureInletVelocity` (T8 and T21 used
`InletVelocityOutletPressure` in the reference settings), and `SimulationTime` is 4000, 2000, and 1400 s for 1, 2,
and 3 × 10⁻⁵ m³/s.
