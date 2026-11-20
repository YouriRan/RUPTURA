# Unit conversions: Sharma et al. 2023, CO₂/propane in MOR

Paper: S. Sharma et al., *Modelling of adsorbate-size dependent explicit isotherms using a segregated approach to
account for surface heterogeneities*, Mol. Phys. 121 (2023) e2183721.

## Simulation box and framework mass

The TMMC box is 2 × 2 × 4 unit cells of all-silica MOR (Si₄₈O₉₆ per unit cell, RASPA3 pseudo-atom masses
Si 28.0855 and O 15.9994 g/mol).

$$
M_\mathrm{uc} = 48 \cdot 28.0855 + 96 \cdot 15.9994 = 2884.04\ \mathrm{g\,mol^{-1}}, \qquad
m_\mathrm{fw} = \frac{16\,M_\mathrm{uc}}{1000\,N_\mathrm{A}} = 7.6625 \times 10^{-23}\ \mathrm{kg}
$$

`ReferenceFrameworkMass` in `MPDSettings` is $m_\mathrm{fw}$ (written by `mpd/build_distribution.py`). A loading
$q$ in mol/kg corresponds to $q\, m_\mathrm{fw} N_\mathrm{A} = 46.145\, q$ molecules in the box.

```python
AVOGADRO = 6.02214076e23
FRAMEWORK_MASS = 16 * (48 * 28.0855 + 96 * 15.9994) * 1.0e-3 / AVOGADRO   # kg
molecules = loading * FRAMEWORK_MASS * AVOGADRO                            # loading in mol/kg
```

The crystal density of the CIF cell ($V_\mathrm{uc} = 2798.9$ Å³), $M_\mathrm{uc} / (N_\mathrm{A} V_\mathrm{uc}) =
1711$ kg/m³, is the particle density of Table 1 (1711.06 kg/m³).

## Macrostate ranges

From SIAST with the Table 3 fits (equimolar, 300 K): $10^6$ Pa gives 1.67 mol/kg CO₂ and 1.05 mol/kg propane
(77 and 48 molecules), and $10^8$ Pa gives 3.81 and 0.40 mol/kg (176 and 18 molecules). The grid is
$n_\mathrm{CO_2} + n_\mathrm{C_3H_8} \le 200$ with $n_\mathrm{C_3H_8} \le 100$. Pure CO₂ reaches 200 molecules
(4.33 mol/kg) at about $6 \times 10^6$ Pa (Table 3 fit), so the pure-CO₂ MPD is valid up to a few MPa only.

## Isotherm parameters (Table 3)

Table 3 gives $q_\mathrm{max}$ in mol/kg and $k$ in 1/Pa. RUPTURA's Langmuir parameters are $[q_\mathrm{sat}, b]$
in the same units, so they are used unchanged:

| | site 1 (small pockets) | site 2 (large pockets) |
| --- | --- | --- |
| CO₂ | 1.60 mol/kg, 2.617 × 10⁻⁴ Pa⁻¹ | 3.62 mol/kg, 4.859 × 10⁻⁷ Pa⁻¹ |
| C₃H₈ | 1.73 mol/kg, 6.506 × 10⁻⁹ Pa⁻¹ | 1.09 mol/kg, 2.376 × 10⁻⁴ Pa⁻¹ |

`examples/MOR-CO2-C3H8` carries different values (e.g. CO₂ 4.4 and 10.0 with 2.91 × 10⁻⁴ and 6.96 × 10⁻⁷); they are
not the Table 3 fits and are not used here.

## Breakthrough (Sec. 3, Tables 1 and 2)

| Paper | Value | RUPTURA key |
| --- | --- | --- |
| $\rho_p$ | 1711.06 kg/m³ | `ParticleDensity` |
| $\varepsilon_b$ | 0.128 | `ColumnVoidFraction` |
| $v_\mathrm{in}$ (interstitial) | 0.1 m/s | `ColumnEntranceVelocity` (interstitial) |
| $p_\mathrm{tot}$ | 1.0 bar = 1 × 10⁵ Pa | `InletPressure` |
| $L$ | 0.1 m | `ColumnLength` |
| $k_M$ | 0.06 s⁻¹ | `MassTransferCoefficient` |
| feed | 10 % CO₂, 10 % C₃H₈, 80 % He | `GasPhaseMolFraction` |

The paper's grid ($\Delta z = 0.005$ m, 20 cells, first-order upwind, SSP-RK3 at $\Delta t = 10^{-3}$ s) is not
reproduced; all methods here use the same RUPTURA numerics.

## Fugacity

Ideal gas, fugacity coefficient 1 (paper Sec. 3; RASPA3 `FugacityCoefficient` 1). The TMMC runs report ln Π at
the fugacity `REFERENCE_PRESSURE` = 10⁵ Pa of their sampled component. `build_distribution.py` shifts the
distribution to a common reference fugacity $f_\mathrm{ref}$,

$$
\ln \Pi(n, m; f_\mathrm{ref}) = \ln \Pi(n, m; 10^5\ \mathrm{Pa}) + (n + m) \ln (f_\mathrm{ref} / 10^5\ \mathrm{Pa}),
$$

chosen such that the populated states span less than 700 in ln Π (doubles underflow below $e^{-708}$).

The N-d TMMC distribution `lnpi_co2_propane_mor_300K_nd_tmmc_walks.txt` is given at
$f_\mathrm{CO_2} = f_\mathrm{C_3H_8} = 10^5$ Pa; `mpd/build_nd_tmmc_distribution.py` keeps that reference fugacity
(ln Π spans 433 there) and only renormalizes.

## Energies

RASPA3 writes $U/k_\mathrm{B}$ in K. RUPTURA's mean Hamiltonian $\langle H \rangle$ is in J:
$\langle H \rangle = k_\mathrm{B} \langle U - U_\mathrm{intra} \rangle$ with
$k_\mathrm{B} = 1.380649 \times 10^{-23}$ J/K.

For the N-d TMMC distribution, $\langle H \rangle = k_\mathrm{B} (\langle U \rangle - \tfrac{3}{2} T_\mathrm{ref}\,
n_\mathrm{C_3H_8})$, as its header prescribes: the ideal-gas mean intramolecular energy of propane (two bonds and one
bend, harmonic) is subtracted instead of the sampled one. The two differ by up to 3.5 % of
$\langle U - U_\mathrm{intra} \rangle$ at the highest propane loadings, by 0.02 % at the breakthrough feed.

```python
hamiltonian = BOLTZMANN * (U_mean - 1.5 * 300.0 * n_C3H8)   # J; U_mean in K
```
