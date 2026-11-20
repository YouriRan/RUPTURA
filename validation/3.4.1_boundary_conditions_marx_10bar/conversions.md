# Conversions for the boundary-condition check

All inputs are taken from `../../3.3.1_energy_balance_marx_ternary/marx2/ternary/<run>/simulation.json`. This file documents how the
paper's operating conditions map onto those inputs, and how `generate.py` derives the outlet pressure
and pressure gradient for the pressure-controlled boundary modes. All quantities are SI.

## 1. Volumetric flow rate to `ColumnEntranceVelocity` (inherited from marx2)

Marx et al. report a feed flow rate of 10 cm³/s, i.e. $Q = 1\times10^{-5}$ m³/s. marx2 converts it with the
total porosity of the reference settings,

$$\varepsilon_t = \varepsilon_b + (1-\varepsilon_b)\,\varepsilon_p = 0.403 + 0.597 \times 0.568 = 0.742096,$$

$$v = \frac{Q}{A\,\varepsilon_t}, \qquad A = \frac{\pi d_\mathrm{in}^2}{4} = 4.908739\times10^{-4}\ \mathrm{m^2},$$

which gives $v = 0.0274517$ m/s, the value in every base input at $1\times10^{-5}$ m³/s.

```python
import math

area = math.pi * 0.025**2 / 4.0          # m2
total_porosity = 0.403 + (1.0 - 0.403) * 0.568
velocity = 1.0e-5 / (area * total_porosity)  # 0.0274517 m/s
```

RUPTURA reads `ColumnEntranceVelocity` as the interstitial velocity of a bed with void fraction
$\varepsilon_b = 0.403$, so the superficial velocity in the simulation is
$\varepsilon_b v = 0.01106$ m/s, whereas $Q/A = 0.02037$ m/s. Together with `ParticleDensity` = 502 kg/m³,
which Marx et al. report as the bulk density (particle density 850 kg/m³), this mapping gives a solid
capacity per unit fluid volume that is about 10% above that of the reference model. This is listed in
`../../TODO.md`; it is the same for every boundary mode and does not affect the comparison between modes.

`generate.py` inverts this relation only to label the directories ($Q = v A \varepsilon_t$).

## 2. Pressure and temperature

| Paper | SI | Input key |
|---|---|---|
| 10 bar | $1\times10^{6}$ Pa | `InletPressure` |
| 25 °C | 298 K | `InfluxTemperature` (T3, T15) |
| 45 °C | 318 K | `InfluxTemperature` (T9, T21) |

`Temperature` is 298 K in every base input and sets the ambient temperature (see the SI note on the
318 K experiments).

## 3. Pressure drop of a feed-filled column

RUPTURA evaluates the Ergun law in terms of the interstitial velocity with
$r = (1-\varepsilon_b)/\varepsilon_b$ (`src/geometry.cpp`, `src/geometry.h`):

$$G = -\frac{\partial P}{\partial z} = \frac{150\,r^2}{d_p^2}\,\mu\,v + \frac{1.75\,r}{d_p}\,\rho\,v\,|v|.$$

For the column filled with the feed at the feed temperature,

$$\bar M = \sum_j y_j M_j, \qquad \rho_\mathrm{feed} = \frac{P_\mathrm{in}\,\bar M}{R\,T_\mathrm{feed}},$$

and `generate.py` sets

$$P_\mathrm{out} = P_\mathrm{in} - G\,L, \qquad \texttt{PressureGradient} = -G.$$

With these values every mode has $P(0) = P_\mathrm{in}$, $P(L) = P_\mathrm{out}$, and $v(0) = v$ once the
column is filled with the feed.

```python
R = 8.31446261815324                      # J mol-1 K-1, src/utils.h
mu, void, d_p, length = 1.92e-5, 0.403, 3.0e-3, 1.2
r = (1.0 - void) / void
molar_mass = 0.40 * 0.044 + 0.30 * 0.028 + 0.30 * 0.002016   # T15 feed, kg/mol
density = 1.0e6 * molar_mass / (R * 298.0)                   # kg/m3
gradient = 150.0 * r**2 / d_p**2 * mu * velocity + 1.75 * r / d_p * density * velocity**2
outlet_pressure = 1.0e6 - gradient * length
```

## 4. Velocity for a prescribed gradient

`InletPressureOutletPressure` (shooting) and `FixedPressureInletVelocity` (local Ergun balance) obtain the
velocity from the gradient. With $a = 150\,r^2\mu/d_p^2$ and $b = 1.75\,r\,\rho/d_p$ the stable root is
(`src/geometry.h`)

$$v = \frac{2G}{a + \sqrt{a^2 + 4\,b\,G}}.$$

In the initially helium-filled column ($\bar M = 0.004$ kg/mol) the same gradient drives a higher velocity
than in the feed. The 2 s smoke test of T15 reproduced $v_\mathrm{He}$ below exactly, in both modes.

```python
a = 150.0 * r**2 * mu / d_p**2
b = 1.75 * r * (1.0e6 * 0.004 / (R * 298.0)) / d_p
helium_velocity = 2.0 * gradient / (a + math.sqrt(a * a + 4.0 * b * gradient))
```

## 5. Resulting values

| Run | Feed CO₂/N₂/H₂ | $T_\mathrm{feed}$ / K | $\bar M$ / g mol⁻¹ | $\rho_\mathrm{feed}$ / kg m⁻³ | $G$ / Pa m⁻¹ | $G L$ / Pa | $P_\mathrm{out}$ / Pa | $v_\mathrm{He}$ / m s⁻¹ |
|---|---|---|---|---|---|---|---|---|
| T3  | 15/80/5  | 298 | 29.101 | 11.745 | 26.926 | 32.312 | 999967.688 | 0.03580 |
| T9  | 15/80/5  | 318 | 29.101 | 11.006 | 26.445 | 31.734 | 999968.266 | 0.03533 |
| T15 | 40/30/30 | 298 | 26.605 | 10.738 | 26.270 | 31.524 | 999968.476 | 0.03498 |
| T21 | 40/30/30 | 318 | 26.605 | 10.062 | 25.831 | 30.997 | 999969.003 | 0.03456 |

The pressure drop is about $3\times10^{-5}$ of the operating pressure, so any difference between the modes
comes from how they treat the velocity, not from the pressure level.
