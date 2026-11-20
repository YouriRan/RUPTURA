# Unit conversions: Liu et al. (2024) to RUPTURA

Every number in `generate.py` is derived here. The article mixes atm, kPa, mL min$^{-1}$
and relative humidity; RUPTURA takes SI throughout. Symbols follow the paper:
$n_s$ saturation loading, $b$ affinity, $t_T$ Toth exponent, $k$ linear-driving-force
rate constant, $\varphi$ relative humidity.

Reference: H. Liu, H. Lin, S. Dai, D.-e. Jiang, *Ind. Eng. Chem. Res.* **63** (2024)
5871-5879, DOI 10.1021/acs.iecr.3c04535.

## Constants

| Symbol | Value | Meaning |
| --- | --- | --- |
| $R$ | 8.31446261815324 J mol$^{-1}$ K$^{-1}$ | gas constant |
| $p^\circ$ | 100000 Pa | standard pressure, and the total pressure of Figures 2-7 |
| $T^\circ$ | 273.15 K | standard temperature for the volumetric flows |
| 1 atm | 101325 Pa | used by the Table 1 and Table 2 affinities |

## 1. Toth affinity: atm$^{-1}$ and kPa$^{-1}$ to Pa$^{-1}$

RUPTURA evaluates its Toth site as $q = n_s\,bp / \left(1 + (bp)^{t_T}\right)^{1/t_T}$
with $p$ in Pa, so $b$ must be in Pa$^{-1}$.

$$b\,[\mathrm{Pa^{-1}}] = \frac{b\,[\mathrm{atm^{-1}}]}{101325}, \qquad
  b\,[\mathrm{Pa^{-1}}] = \frac{b\,[\mathrm{kPa^{-1}}]}{1000}$$

```python
FLUE_GAS_TOTH["affinity0"] = 8.3e5 / 101_325.0   # 8.19146 Pa-1   (Table 1, Figure 1)
BASELINE["affinity"]       = 1.0e5 / 101_325.0   # 0.98692 Pa-1   (Table 2, Figures 2-4)
FIBER_TOTH["affinity0"]    = 6.2e3 / 1000.0      # 6.2     Pa-1   (Table 3, Figures 5-7)
```

The Figure 3b affinity sweep is tabulated in atm$^{-1}$ and converted the same way:

| $b$ / atm$^{-1}$ | $b$ / Pa$^{-1}$ |
| --- | --- |
| 1e4 | 0.0986923 |
| 5e4 | 0.4934616 |
| 1e5 | 0.9869233 |
| 2e5 | 1.9738465 |
| 3e5 | 2.9607698 |

## 2. Temperature dependence of the Toth parameters (equations 4-6)

$$n_s(T) = n_{s0}\exp\!\left[\chi\left(1 - \frac{T}{T_0}\right)\right], \quad
  b(T) = b_0\exp\!\left[\frac{\Delta H}{R T_0}\left(\frac{T_0}{T} - 1\right)\right], \quad
  t_T(T) = t_{T0} + \alpha\left(1 - \frac{T_0}{T}\right)$$

These are evaluated once per run, at that run's temperature, because every
breakthrough here is isothermal (`energyBalance: false`). The results are written
straight into the site's `Parameters`.

Flue-gas sorbent, $T_0 = 298.15$ K, $n_{s0} = 0.97$, $b_0 = 8.19146$ Pa$^{-1}$,
$\Delta H = 84000$ J mol$^{-1}$, $t_{T0} = 0.21$, $\chi = \alpha = 0$:

| $T$ / K | $n_s$ / mol kg$^{-1}$ | $b$ / Pa$^{-1}$ | $t_T$ / [-] |
| --- | --- | --- | --- |
| 298.15 | 0.970000 | 8.191463 | 0.21 |
| 318.15 | 0.970000 | 0.973314 | 0.21 |

PEI/silica fiber sorbent, $T_0 = 308$ K, $n_{s0} = 0.81$, $b_0 = 6.2$ Pa$^{-1}$,
$\Delta H = 210000$ J mol$^{-1}$, $t_{T0} = 0.40$, $\chi = 6.6$, $\alpha = 10.8$:

| $T$ / K | $n_s$ / mol kg$^{-1}$ | $b$ / Pa$^{-1}$ | $t_T$ / [-] |
| --- | --- | --- | --- |
| 308 | 0.810000 | 6.200000 | 0.400000 |
| 318 | 0.653765 | 0.470392 | 0.739623 |
| 328 | 0.527666 | 0.041765 | 1.058537 |

## 3. Feed concentration to mole fraction

The sweeps of Figures 2-4 specify a CO$_2$ concentration; RUPTURA takes a mole
fraction at the inlet pressure.

$$y_i = \frac{c_i R T}{p_\mathrm{total}}$$

$$y_{\mathrm{CO_2}} = \frac{1.6\times10^{-2} \times 8.31446 \times 298.15}{100000}
  = 3.96633\times10^{-4} \quad (p_{\mathrm{CO_2}} = 39.6633\ \mathrm{Pa})$$

Figures 6 and 7 give 380 ppm directly, so $y_{\mathrm{CO_2}} = 380\times10^{-6}$.
Figure 1 gives 0.1 atm of CO$_2$ in a 1 atm stream, so $y_{\mathrm{CO_2}} = 0.1$
with `InletPressure` 101325 Pa.

## 4. Bulk density to particle density

The paper writes its mass balance with the bulk density $\rho_b$ (kilograms of
sorbent per cubic metre of bed). RUPTURA takes the particle density $\rho_p$ and
the void fraction $\varepsilon$ separately and forms $\rho_b = (1-\varepsilon)\rho_p$
internally, so

$$\rho_p = \frac{\rho_b}{1 - \varepsilon}$$

| Figures | $\rho_b$ / kg m$^{-3}$ | $\varepsilon$ / [-] | $\rho_p$ / kg m$^{-3}$ | Source |
| --- | --- | --- | --- | --- |
| 1-4 | 500 | 0.5 | 1000.00 | Table 2 |
| 5, 6 | 95.7 | 0.43 | 167.89 | calibrated, see README.md |
| 7 | 160.9 | 0.43 | 282.28 | calibrated, see README.md |

The two fiber densities are not published values. They are inferred from the printed time
axes by scaling a trial density with the ratio of the digitized to the simulated 50%
breakthrough time, and so are specific to this solver.

## 5. Superficial to interstitial velocity

`ColumnEntranceVelocity` is the interstitial velocity, which is what `column.data`
reports. The paper tabulates superficial velocities.

$$u_\mathrm{int} = \frac{u_\mathrm{sup}}{\varepsilon}$$

The Figure 2b sweep of 0.002 to 0.010 m s$^{-1}$ superficial at $\varepsilon = 0.5$
therefore enters RUPTURA as 0.004 to 0.020 m s$^{-1}$; the Table 2 baseline
0.008 m s$^{-1}$ enters as 0.016 m s$^{-1}$.

## 6. Standard volumetric flow to velocity

Figures 6 and 7 give standard volumetric flows referenced to $T^\circ = 273.15$ K
and $p^\circ = 100000$ Pa, through a fiber module of radius 0.002 m.

$$A = \pi r^2 = \pi (0.002)^2 = 1.256637\times10^{-5}\ \mathrm{m^2}$$
$$u_\mathrm{sup} = \frac{Q^\circ}{A}\cdot\frac{T}{T^\circ}\cdot\frac{p^\circ}{p}, \qquad
  u_\mathrm{int} = \frac{u_\mathrm{sup}}{\varepsilon}, \quad \varepsilon = 0.43$$

The runs are at $p = p^\circ$, so the pressure ratio is 1.

| $T$ / K | $Q^\circ$ / m$^3$ s$^{-1}$ | $u_\mathrm{sup}$ / m s$^{-1}$ | $u_\mathrm{int}$ / m s$^{-1}$ |
| --- | --- | --- | --- |
| 308 | 6.6667e-7 | 0.0598203 | 0.139117 |
| 308 | 1.0000e-6 | 0.0897304 | 0.208675 |
| 308 | 1.5000e-6 | 0.1345960 | 0.313013 |
| 318 | 6.6667e-7 | 0.0617625 | 0.143634 |
| 318 | 1.0000e-6 | 0.0926437 | 0.215451 |
| 318 | 1.5000e-6 | 0.1389660 | 0.323176 |
| 328 | 6.6667e-7 | 0.0637047 | 0.148150 |
| 328 | 1.0000e-6 | 0.0955571 | 0.222226 |
| 328 | 1.5000e-6 | 0.1433360 | 0.333339 |
| 308 | 3.3333e-6 (Figure 7) | 0.2991010 | 0.695585 |

## 7. Relative humidity to water mole fraction

The paper prints the Buck equation (equation 9) for the water saturation pressure.
RUPTURA's GAB site, and the relative humidity that the humidity-dependent CO$_2$ sites
see, use the default `WaterSaturationPressureModel`, a Kelvin Antoine correlation
converted from mmHg to Pa:

$$p_\mathrm{sat} = 10^{\,8.07131 - 1730.63/(T - 39.724)} \times 133.32236842105263$$

Because RUPTURA's GAB site uses this correlation and not Buck, `generate.py` and the
notebook must use it too, or the humidity axis would not agree with the loadings.
The two differ by well under a percent over 274-373 K, but mixing them would be
inconsistent rather than merely imprecise.

$$y_{\mathrm{H_2O}} = \frac{\varphi\,p_\mathrm{sat}(T)}{p_\mathrm{total}}$$

| $T$ / K | $p_\mathrm{sat}$ / Pa | $\varphi$ | $y_{\mathrm{H_2O}}$ |
| --- | --- | --- | --- |
| 298.15 | 3158.09 | 0.20 | 0.006316 |
| 298.15 | 3158.09 | 0.30 | 0.009474 |
| 298.15 | 3158.09 | 0.40 | 0.012632 |
| 298.15 | 3158.09 | 0.50 | 0.015790 |
| 308.00 | 5562.65 | 0.85 | 0.047283 |

Going the other way, the notebook recovers the outlet humidity of Figure 4b as
$\varphi = p_{\mathrm{H_2O}}/p_\mathrm{sat}$ from RUPTURA's reported partial pressure,
with $p_\mathrm{sat}$ from the same Antoine expression (stored per case in
`simulations/index.json` as `water_saturation_pressure_pa`).

## 8. GAB parameters

The paper's GAB form (equation 7) is

$$q_{\mathrm{H_2O}} = \frac{c_m c_g k \varphi}{(1 - k\varphi)\left(1 + (c_g - 1)k\varphi\right)}$$

RUPTURA's GAB site takes `[q_mono, C, K]` with $\varphi = p/p_\mathrm{sat}(T)$, and gets $T$
from the column (`setTemperature`); $C$ and $K$ may carry VantHoff laws, which
these constant parameters do not need. The parameters map one to one:

```python
SWEEP_GAB = (35.0, 0.15, 0.51)   # Figure 4  -> Parameters [35.0, 0.15, 0.51]
FIBER_GAB = (26.0, 0.22, 0.84)   # Figure 5a -> Parameters [26.0, 0.22, 0.84]
```

The H$_2$O component carries `"Water": true`, so that its partial pressure sets the
relative humidity of the humidity-dependent CO$_2$ sites (sections 9 and 10) and its
GAB isotherm their water loading.

## 9. Humid enhancement of the CO$_2$ isotherm (equation 8)

For the fiber sorbent the paper multiplies the whole Toth loading by an empirical
humid factor, which is RUPTURA's `Exponential-humidity-Toth` site
`[n_s, b, t_T, 1.11]`, evaluated at the local relative humidity:

$$q_{\mathrm{CO_2}}(\varphi) = q_{\mathrm{CO_2}}(0)\,e^{1.11\varphi}
  \quad\Longleftrightarrow\quad n_s(\varphi) = n_s(0)\,e^{1.11\varphi}$$

At the Figure 7 feed humidity the factor is

$$e^{1.11 \times 0.85} = 2.568957 \quad\Longrightarrow\quad
  n_s(308\ \mathrm{K},\ \varphi = 0.85) = 0.81 \times 2.568957 = 2.080855\ \mathrm{mol\ kg^{-1}}$$

Figure 5b's humid isotherm is a MixturePrediction with `"RelativeHumidity": 0.85`;
Figure 7 takes $\varphi$ from the water partial pressure in each grid cell, so the
enhancement rises as the water front passes.

## 10. The Figure 4 modified Toth law

For the sweep sorbent the paper makes the Toth parameters functions of the
adsorbed water loading $q_{\mathrm{H_2O}}$:

$$n_s = \frac{n_{s0}}{1 - 0.006\,q_{\mathrm{H_2O}}}, \qquad b = b_0\left(1 + 25\,q_{\mathrm{H_2O}}\right)$$

This is RUPTURA's `GAB-modified-Toth` site `[n_s0, b_0, t_T, 0.006, 25]`, with
$q_{\mathrm{H_2O}}$ the GAB loading at the local water partial pressure. Until
2026-11-04 the runs froze the law at the GAB loading of the feed humidity; the coupled
law changes the CO$_2$ breakthrough by at most $1.2\times10^{-4}\,c_0$, because water
saturates the bed in roughly 2000 to 9000 s while the CO$_2$ front needs 90000 s.
The dry case keeps the plain Toth site. The parameters at the feed humidity:

| $\varphi$ | $q_{\mathrm{H_2O}}$ / mol kg$^{-1}$ | $n_s$ / mol kg$^{-1}$ | $b$ / Pa$^{-1}$ |
| --- | --- | --- | --- |
| 0.00 | 0.0000 | 1.0000 | 0.9869 |
| 0.20 | 0.6529 | 1.0039 | 17.0968 |
| 0.30 | 1.0901 | 1.0066 | 27.8835 |
| 0.40 | 1.6277 | 1.0099 | 41.1479 |
| 0.50 | 2.2943 | 1.0140 | 57.5934 |

## 11. Rate constants

The paper gives the fiber uptake rate constants per minute and the sweep constants
per second:

$$k\,[\mathrm{s^{-1}}] = \frac{k\,[\mathrm{min^{-1}}]}{60}$$

| Case | $k$ as printed | `MassTransferCoefficient` / s$^{-1}$ |
| --- | --- | --- |
| Figures 2-4 CO$_2$ baseline | 0.02 s$^{-1}$ | 0.02 |
| Figure 4 H$_2$O | 0.03 s$^{-1}$ | 0.03 |
| Figure 6, 308 K | 0.040 min$^{-1}$ | 6.6667e-4 |
| Figure 6, 318 K | 0.048 min$^{-1}$ | 8.0000e-4 |
| Figure 6, 328 K | 0.100 min$^{-1}$ | 1.6667e-3 |
| Figure 7 humid CO$_2$ | 0.020 min$^{-1}$ | 3.3333e-4 |
| Figure 7 H$_2$O | 0.360 min$^{-1}$ | 6.0000e-3 |

## 12. Solver tolerances and the dilute feed

CVODE's absolute tolerance on concentration is an absolute number in mol m$^{-3}$, so its
meaning depends entirely on how dilute the feed is. These DAC streams are far more dilute
than a typical separation:

| Figures | Feed CO$_2$ | Default 1e-5 as a fraction of feed |
| --- | --- | --- |
| 1 | 4.088 mol m$^{-3}$ | 0.00024% |
| 2-4 | 0.016 mol m$^{-3}$ | 0.063% |
| 6, 7 | 0.0148 mol m$^{-3}$ | 0.068% |

For Figure 1 the default is ample. For everything else it allows errors of order a tenth of
a percent of the feed, which is enough to break the monotonicity of a breakthrough curve
visibly. `generate.py` therefore sets it to 1e-6 everywhere, except Figure 4, where 1e-6 is
not tractable. See README.md for the measurements behind both choices.

## 13. Output columns

`component_<i>_<name>.data` is read by the notebook as follows. Column 13 is the
outlet breakthrough curve $c/c_0$; it normalises by the *inlet* partial pressure,
which is what the paper plots.

| Column | Quantity | Used for |
| --- | --- | --- |
| 2 | time / min | multiplied by 60 for the time axis |
| 3 | position / m | last value is the outlet |
| 4 | concentration / mol m$^{-3}$ | Figure 4a |
| 7 | loading / mol kg$^{-1}$ | averaged over the bed for Figures 4c, 4d |
| 11 | partial pressure / Pa | Figure 4b humidity, Figure 5a isotherm |
| 12 | equilibrium loading / mol kg$^{-1}$ | Figure 5a isotherm |
| 13 | $p_i / (p_t y_{i,0})$ / [-] | Figures 1b, 2, 3, 6, 7 |
