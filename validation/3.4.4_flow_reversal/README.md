# 3.4.4 Flow reversal: co-current and counter-current regeneration against equilibrium theory

Verifies the flow reversal of RUPTURA 2.0 (`FlowDirection`, TODO F-BC-1) against an exact solution, on the question
that motivates counter-current purge steps in pressure-swing cycles (Skarstrom 1960): is a partially loaded column
regenerated faster with or against the feed direction?

## Case

A trace adsorbate A, mole fraction $y_\mathrm{F} = 10^{-3}$ in helium, at 298.15 K and $10^5$ Pa, with a Langmuir
isotherm $q^* = q_s b p/(1 + b p)$, $q_s = 1$ mol/kg and $b = 0.04$ Pa⁻¹, so that $b p_\mathrm{F} = 4$. Column:
$L = 0.1$ m, $\varepsilon = 0.4$, $\rho_p = 1000$ kg/m³, $v = 0.05$ m/s (interstitial). $k_\mathrm{LDF} = 1$ s⁻¹ (the
fronts move on time scales of $10^4$ s, so the solid is at equilibrium), no axial dispersion, Danckwerts inlet.

Each run is a two-phase `SwingAdsorption` run with ROS34PW2:

1. Loading, 30 000 s: feed into the helium-filled column (forward). The shock front reaches $z_f = 0.504\,L$.
2. Purge, 320 000 s, with helium: co-current (`FlowDirection` `Forward`, purge at $z = 0$) or counter-current
   (`Reverse`, purge at $z = L$).

`generate.py` writes both directions with WENO3 and Upwind on 25 to 400 cells (20 runs); `run.py` runs them with
eight workers. The runs use `FixedVelocity`: in the trace limit the velocity is constant, which is what the theory
assumes. (With `InletPressureInletVelocity` and ROS34PW2 this very stiff trace case shows the drift of
$\sum_j c_j$ logged in TODO.md; the continuity closures in reverse flow are covered by the unit tests
`FlowReversal.*`.)

## Equilibrium theory

With local equilibrium, constant $v$, and no dispersion (Rhee, Aris, and Amundson 1986; Ruthven 1984),

$$\left(1 + F\frac{\mathrm{d}q^*}{\mathrm{d}c}\right)\frac{\partial c}{\partial t} + v\frac{\partial c}{\partial z} = 0,
\qquad F = \frac{(1-\varepsilon)\rho_p}{\varepsilon},$$

so a concentration $c$ travels at $w(c) = v/(1 + F\,\mathrm{d}q^*/\mathrm{d}c)$, and a jump from $c_l$ to $c_r$ at
$w_s = v/(1 + F\,[q^*(c_l) - q^*(c_r)]/(c_l - c_r))$. For the Langmuir isotherm, with $K = bRT$,

$$w(c) = \frac{v}{1 + F q_s K/(1 + Kc)^2},
\qquad c(w) = \frac{1}{K}\left(\sqrt{\frac{F q_s K}{v/w - 1}} - 1\right).$$

- **Loading:** a shock at $w_s = v/(1 + F q^*(c_\mathrm{F})/c_\mathrm{F})$; stoichiometric time
  $t_\mathrm{st} = L/w_s = 59\,497$ s; front $z_f = w_s t_\mathrm{load}$.
- **Counter-current purge** (time $t$ from the start of the purge): a desorption wave centred at the front runs back to
  $z = 0$, $c = c(w)$ with $w = (z_f - z)/t$ for $w(0)t < z_f - z < w(c_\mathrm{F})t$, $c_\mathrm{F}$ closer to $z = 0$,
  and 0 beyond the front. The outlet at $z = 0$ holds $c_\mathrm{F}$ until $z_f/w(c_\mathrm{F}) = 0.101\,t_\mathrm{st}$
  and is clean after $z_f/w(0) = 2.52\,t_\mathrm{st}$.
- **Co-current purge:** the wave is centred at $z = 0$, $c = c(z/t)$, with the plateau $c_\mathrm{F}$ up to the shock at
  $z_f + w_s t$. Its fastest characteristic overtakes the shock at $t^* = z_f/(w(c_\mathrm{F}) - w_s) = 0.126\,t_\mathrm{st}$;
  afterwards the weakened shock follows $\mathrm{d}z_s/\mathrm{d}t = v/(1 + F q^*(c_l)/c_l)$ with $c_l = c(z_s/t)$,
  integrated with RK4 (the notebook). It leaves the column at $0.666\,t_\mathrm{st}$, and the column is clean after
  $L/w(0) = 5.00\,t_\mathrm{st}$.

Error measure: $E_1 = \int |c_\mathrm{out} - c_\mathrm{out}^\mathrm{ref}|\,\mathrm{d}t / \int c_\mathrm{out}^\mathrm{ref}\,\mathrm{d}t$
over the purge, with $c_\mathrm{out}$ the first cell (counter-current, the gas leaving at $z = 0$) or the last cell
(co-current). Regeneration: amount left $n(t) = \sum_i h_i(c_i + F q_i)$ relative to the end of loading.

## Results (2026-11-19, binary with the B-ITG-17 fix)

| | Counter-current | Co-current |
| --- | --- | --- |
| $E_1$, WENO3, $N_z$ = 25 / 100 / 400 | 5.0 / 1.24 / 0.34 % | 5.0 / 1.28 / 0.33 % |
| $E_1$, Upwind, $N_z$ = 25 / 100 / 400 | 12.2 / 4.1 / 1.4 % | 12.3 / 4.2 / 1.3 % |
| 99 % regenerated, theory | 2.090 $t_\mathrm{st}$ | 4.384 $t_\mathrm{st}$ |
| 99 % regenerated, WENO3 100 / 400 cells | 2.181 / 2.112 $t_\mathrm{st}$ | 4.465 / 4.405 $t_\mathrm{st}$ |

The error decreases in proportion to the grid spacing (first order at the shocks and corners) with both schemes. The
counter-current purge regenerates the column in half the time and keeps its clean part clean. After loading, the
column holds the fed amount to within $2\times10^{-5}$ with WENO3 (25 cells; $3\times10^{-6}$ on 100) and $2\times10^{-8}$
with Upwind. The runs of 2026-11-06 held $3\times10^{-4}$ more with WENO3: ROS34PW2 clamped the WENO3 undershoots ahead
of the front at every step (TODO.md B-ITG-17).

Figures: `analysis.ipynb` writes `ruptura-2.0-paper/figures/flow_reversal/` (outlet histories, profiles of both purges,
grid convergence, amount left).

## References

- H.-K. Rhee, R. Aris, N. R. Amundson, First-Order Partial Differential Equations, Vol. I, Prentice-Hall (1986).
- D. M. Ruthven, Principles of Adsorption and Adsorption Processes, Wiley (1984).
- C. W. Skarstrom, US Patent 2,944,627 (1960).
