# Unit conversions: Stadler et al. (2022) to RUPTURA

Reference: T. J. Stadler, L. J. Bender, P. Pfeifer, Front. Chem. Eng. 4 (2022) 1000064
(`fceng-04-1000064.pdf`). Flow rates of the breakthrough experiments are from Stadler et al.,
Energies 14 (2021) 355, Table 2.

RUPTURA works in SI units: Pa, K, m, s, mol, kg. All numbers below are computed in `generate.py`
(simulation inputs) or `digitize_reference.py` (reference data).

```python
R = 8.31446261815324   # J/(mol K)
T_STP = 273.15         # K
P_STP = 101325.0       # Pa
```

## 1. Temperature

$T/\mathrm{K} = \vartheta/^\circ\mathrm{C} + 273.15$: 225, 250, 275 °C → 498.15, 523.15, 548.15 K.

## 2. Flow at STP to molar flow and velocity

The paper gives volumetric flows at STP (mass flow controllers). STP is taken as 273.15 K and
101325 Pa (DIN 1343); the paper does not say which reference state it uses.

$$\dot n = \frac{p_{\mathrm{STP}}\,\dot V_{\mathrm{STP}}}{R\,T_{\mathrm{STP}}},\qquad
u = \frac{\dot n\,R\,T}{p\,A},\qquad v = \frac{u}{\varepsilon_b}$$

$u$ is the superficial velocity (the $u$ of the paper's Eq. 15). $v$ is the interstitial velocity,
which is RUPTURA's `ColumnEntranceVelocity`. The cross-section is two slits of 2 mm × 50 mm:
$A = 2 \times 0.002 \times 0.050 = 2.0\times10^{-4}$ m².

```python
A = 2 * 0.002 * 0.050                              # m2
n = P_STP * 2000e-6 / 60 / (R * T_STP)             # 2000 mL/min -> 1.4872e-3 mol/s
u = n * R * 523.15 / (8e5 * A)                     # 8 bar, 523.15 K -> 0.04043 m/s
v = u / 0.4                                        # -> 0.10107 m/s
```

| step | $\dot V_{\mathrm{STP}}$ / mL min⁻¹ | $p$ / bar | $\dot n$ / mmol s⁻¹ | $u$ / m s⁻¹ | $v$ / m s⁻¹ |
|---|---|---|---|---|---|
| adsorption (Figs. 6–8), WGS activity (Fig. 4) | 2000 | 8 | 1.4872 | 0.04043 | 0.10107 |
| regeneration (Figs. 6–8) | 1000 | 1 | 0.7436 | 0.16172 | 0.40430 |

The gas residence time in the bed is $L/v \approx 1$ s.

## 3. Bed and sorbent

The slits are 100 mm long in the flow direction, the bed void fraction is 0.4, and the particle
radius is 100 µm (Table 2). The sorbent masses come from Section 2.5. RUPTURA stores the sorbent
through `ParticleDensity`, so that $\rho_p (1-\varepsilon_b) A L = m_{\mathrm{sorb}}$:

```python
rho_p = m_sorbent / ((1 - 0.4) * A * 0.100)        # 9.70 g -> 808.3 kg/m3; 9.35 g -> 779.2 kg/m3
```

The resulting bulk density is 485 kg m⁻³. Table 2 lists 1096 kg m⁻³, but that table describes the
pilot plant. The laboratory runs use the measured sorbent mass.

`Geometry` must be given for a breakthrough run, but RUPTURA has no slit geometry. A `PackedBed` with
the diameter of a tube of equal area is used, $d = \sqrt{4A/\pi}$ = 15.96 mm. Only the wall
heat-transfer model uses this diameter, and the runs are isothermal.

## 4. Freundlich isotherm (Eq. 11, Table 3)

Paper: $q = k_{\mathrm{Fr}}\,p^{1/n_{\mathrm{Fr}}}$, with $q$ in mmol g⁻¹ (= mol kg⁻¹), $p$ in bar,
$k_{\mathrm{Fr}}$ = 0.177 mmol g⁻¹ bar⁻¹ᐟⁿ, and $1/n_{\mathrm{Fr}}$ = 0.220.
RUPTURA: `Freundlich` with `Parameters` $[a, \nu]$ means $q = a\,p^{1/\nu}$ with $p$ in Pa.

$$a = k_{\mathrm{Fr}}\,(10^5)^{-1/n_{\mathrm{Fr}}},\qquad \nu = n_{\mathrm{Fr}}$$

```python
a  = 0.177 * 1e5 ** -0.220     # 0.0140596 mol kg-1 Pa^(-1/nu)
nu = 1 / 0.220                 # 4.54545
```

Check: $q_{\mathrm{eq}}$(0.4 bar) = 0.1447 mol kg⁻¹, which is the CO₂ partial pressure of the Fig. 6
feed at 8 bar.

## 5. Sorption kinetics (Eqs. 6–10, 13)

- Adsorption: $\mathrm{d}q/\mathrm{d}t = k_{\mathrm{ads}}(q_{\mathrm{eq}} - q)$, with $k_{\mathrm{ads}}$ = 0.1125 s⁻¹. This
  maps directly onto RUPTURA's LDF `MassTransferCoefficient`, in s⁻¹.
- The bulk-to-pore LDF of Eq. 13 is not modelled separately. With $D_{\mathrm{eff}} \sim 10^{-6}$ m² s⁻¹,
  $k_{\mathrm{LDF}} = 15 D_{\mathrm{eff}}/r_p^2 \sim 10^3$ s⁻¹, four orders of magnitude above
  $k_{\mathrm{ads}}$. The pore gas therefore follows the bulk gas.
- Desorption (q > q*): $k_{\mathrm{des}} = k^1_{\mathrm{des}} \exp[\beta_{\mathrm{des}}(q/q_{\max})/RT]$, with
  $k^1_{\mathrm{des}}$ = 5.38×10⁻⁴ s⁻¹ and $\beta_{\mathrm{des}}$ = 74,000 J mol⁻¹.
  Here $q_{\max} = q_{\mathrm{eq}}$(20 bar) = 0.177·20^0.22 = 0.3421 mol kg⁻¹, and
  $\beta_{\mathrm{des}}/RT$ = 17.01 at 523.15 K.
- In the cycles, CO₂ sits on a chemisorption site: Freundlich with `MaximumLoading` = q_max, which never
  caps q* at 8 bar or less, and the `FirstOrder` law with `MassTransferCoefficient` = k_ads,
  `DesorptionMassTransferCoefficient` = k¹_des, `DesorptionCoverageEnergy` = β_des, and
  `CoverageReferenceLoading` = q_max.

## 6. Axial dispersion (Eq. 12)

$$D_{\mathrm{ax},i} = D_{\mathrm{mix},i}\sqrt{\varepsilon_b} + u\,r_p$$

Here $u$ is the superficial velocity. RUPTURA's `AxialDispersionCoefficient` is $D_{\mathrm{ax}}$
itself: dividing Eq. 15 by $\varepsilon_b$ gives
$\partial c/\partial t = -v\,\partial c/\partial z + D_{\mathrm{ax}}\,\partial^2 c/\partial z^2 + \dots$.
For $D_{\mathrm{mix},i}$, the binary diffusivity of $i$ in N₂ is used, from the Fuller equation
(Poling et al., *The Properties of Gases and Liquids*, 5th ed., Eq. 11-4.4). N₂ is 65–75 % of every
feed.

$$D_{AB} = \frac{0.00143\,T^{1.75}}{p\,M_{AB}^{1/2}\left[(\Sigma_v)_A^{1/3} + (\Sigma_v)_B^{1/3}\right]^2}
\ \mathrm{cm^2\,s^{-1}},\qquad M_{AB} = \frac{2}{1/M_A + 1/M_B}$$

Here $T$ is in K, $p$ in bar, and $M$ in g mol⁻¹. The diffusion volumes are N₂ 18.5, CO 18.0,
CO₂ 26.9, H₂ 6.12, and H₂O 13.1.

```python
M_AB = 2 / (1 / M_i + 1 / M_N2)                                   # g/mol
D = 1e-4 * 0.00143 * T**1.75 / (p_bar * M_AB**0.5 * (V_i**(1/3) + V_N2**(1/3))**2)   # m2/s
D_ax = D * 0.4**0.5 + u * 100e-6
```

| component | $D_{\mathrm{mix}}$ (8 bar) / m² s⁻¹ | $D_{\mathrm{ax}}$ (8 bar) / m² s⁻¹ | $D_{\mathrm{ax}}$ (1 bar, regeneration) / m² s⁻¹ |
|---|---|---|---|
| CO  | 6.97×10⁻⁶ | 8.45×10⁻⁶ | 5.14×10⁻⁵ |
| CO₂ | 5.49×10⁻⁶ | 7.52×10⁻⁶ | 4.40×10⁻⁵ |
| H₂  | 2.64×10⁻⁵ | 2.07×10⁻⁵ | 1.50×10⁻⁴ |
| H₂O | 8.73×10⁻⁶ | 9.56×10⁻⁶ | 6.03×10⁻⁵ |

The column Péclet number in the adsorption step is $vL/D_{\mathrm{ax}} \approx 1300$, which is close to plug flow.
RUPTURA holds one $D_{\mathrm{ax}}$ per component, so the cycles use the 8 bar values in both steps. This has
no effect, because the upwind convection adds a numerical dispersion $v\Delta z/2$ that is larger still
(see the Fig. 4 grid study in the notebook).

## 7. Dry-basis outlet fractions

The analyser measures CO, CO₂ and H₂ after the steam is condensed. RUPTURA reports wet mole fractions:

$$y_{i,\mathrm{dry}} = \frac{y_i}{1 - y_{\mathrm{H_2O}}}$$

The extra-column mixing is applied to these dry fractions (see the notebook), because the CO tracer
that measures it is itself a dry-basis signal.

## 8. Digitised reference data

The data are in Figure 5 axes: mmHg → Pa (× 133.322387415) and mmol g⁻¹ → mol kg⁻¹ (× 1). All other
figures are already in s and [-].

## 9. WGS rate (Eq. 4, Table 3) to a RUPTURA `BulkConcentration` reaction (Figs. 4, 7, 8B)

Paper: $r = k_\infty e^{-E_a/RT}\,(p_{\mathrm{CO}}p_{\mathrm{H_2O}} - p_{\mathrm{CO_2}}p_{\mathrm{H_2}}/K_{\mathrm{eq}})$.
Here $r$ is in mol g_cat⁻¹ h⁻¹ and $p$ in bar, with $k_\infty$ = 2.01×10⁵ mol bar⁻² g⁻¹ h⁻¹ and
$E_a$ = 51,845 J mol⁻¹. This maps onto `GeneralPowerLaw`, with orders 1/1 on CO and H₂O (forward)
and on CO₂ and H₂ (backward); RUPTURA sets $k_b = k_f/K$.

The reaction acts on the bulk gas: the authors' Fig. 4 model is an ideal plug-flow reactor without
pores (reproduced to 1e-4). RUPTURA's `BulkConcentration` rate is per unit **fluid** volume, with
concentrations in mol m⁻³ ($p_i/\mathrm{bar} = c_i R T/10^5$). With $\rho_{\mathrm{cat}}$ the catalyst
mass per **bed** volume and $\varepsilon_b$ the bed void fraction:

$$R_f = \underbrace{k_\infty\left(\frac{RT}{10^5}\right)^2 \frac{10^3}{3600}\,\frac{\rho_{\mathrm{cat}}}{\varepsilon_b}}_{\texttt{forwardRateCoefficient}}
e^{-E_a/RT}\left(c_{\mathrm{CO}}c_{\mathrm{H_2O}} - \frac{c_{\mathrm{CO_2}}c_{\mathrm{H_2}}}{K}\right)$$

The factor $10^3$ converts g to kg, and $1/3600$ converts h to s. RUPTURA's Arrhenius term has no
reference temperature ($k = k_0 e^{-E_a/RT}$), so the $(RT)^2$ factor makes `forwardRateCoefficient`
specific to each isothermal run.

```python
rho_cat = m_cat / (A * L_catalyst)                                           # kg per m3 of bed
k0 = 2.01e5 * (R * T / 1e5) ** 2 * 1e3 / 3600 * rho_cat / eps_b              # forwardRateCoefficient
Ea = 51845.0                                                                 # ActivationEnergy of the Arrhenius law, J/mol
```

| case | catalyst | $\rho_{\mathrm{cat}}$ / kg m⁻³ of bed | $k_0$ / m³ mol⁻¹ s⁻¹ |
|---|---|---|---|
| Fig. 4, 225 / 250 / 275 °C | 0.50 g in the first 25 mm (with SiC) | 100.0 | 2.395×10⁴ / 2.641×10⁴ / 2.899×10⁴ |
| Figs. 7, 8B, 250 °C | 0.51 g in the first 25 mm (with sorbent) | 102.0 | 2.694×10⁴ |

In the SE-WGS runs, the reaction carries `"Adsorbents": ["catalyst-zone"]`, so it acts only in the first
quarter of the two-section multibed column.

Moe's $K_{\mathrm{eq}} = \exp(4577.8/T - 4.33)$ (Eq. 5) is dimensionless, because the reaction conserves
moles. It maps exactly onto a VantHoff law of RUPTURA's equilibrium constant, $K(T) = K_0\exp[Q_r/(RT)]$ (no
reference temperature), whose Energy is the heat of reaction $Q_r = -\Delta H_r$:

```python
"EquilibriumConstant": {"Value": math.exp(-4.33), "Law": "VantHoff"}   # K_0 = 0.0131675
"HeatOfReaction": 4577.8 * R                                          # Q_r = 38061.9 J/mol
```

$K_{\mathrm{eq}}$ = 128.97, 83.13 and 55.78 at 498.15, 523.15 and 548.15 K. The digitised equilibrium
lines of Fig. 4 match $X^2 = K(1-X)(S/G - X)$ to four decimal places, for a feed of 10 % CO and
$S/G = \dot n_{\mathrm{H_2O}}/\dot n_{\mathrm{CO}}$.
