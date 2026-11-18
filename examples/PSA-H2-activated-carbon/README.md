# H2 purification by a Skarstrom PSA cycle on activated carbon

A four-step pressure swing adsorption cycle for the 40/30/30 CO2/N2/H2 feed of experiment T15 of Marx et al., with the
sorbent, column, and heat-transfer data of that experiment (`marx2/ternary/T15`). The cycle itself is an illustration,
not an experiment. The column starts filled with helium at 1 bar; about 10 of the 20 cycles reach the cyclic steady
state.

| Step | `ClosedEnd` | Pressure | Duration | Gas that enters |
| --- | --- | --- | --- | --- |
| Pressurization | `ProductEnd` | 1 to 10 bar, linear | 40 s | feed at z = 0 |
| Adsorption | | 10 bar | 60 s | feed at z = 0, 0.03 m/s |
| Blowdown | `ProductEnd` | 10 to 1 bar, linear | 40 s | none; gas leaves at z = 0 |
| Purge | | 1 bar | 60 s | H2 at z = L, 0.1 m/s |

## Purity and recovery

`phases.data` (or `SwingAdsorption.phase_records` in Python) holds the amount of each component in the gas and on the
solid at the start and end of every phase, per unit column cross-section. The change in a phase is what entered minus
what left, so the product of the adsorption step follows from its known feed, without integrating the outlet:

```python
import numpy as np
import ruptura

simulation = ruptura.load_simulation("simulation.json")
simulation.run()
pressurization, adsorption, blowdown, purge = simulation.phase_records[-4:]

change = lambda phase: (np.array(phase.gas_at_end) + phase.adsorbed_at_end
                        - phase.gas_at_start - phase.adsorbed_at_start)
R, T, void_fraction = 8.31446261815324, 298.0, 0.403
feed_fractions = np.array([0.0, 0.4, 0.3, 0.3])  # He, CO2, N2, H2
fed = void_fraction * 0.03 * 1.0e6 / (R * T) * 60.0 * feed_fractions
purge_gas = void_fraction * 0.1 * 1.0e5 / (R * T) * 60.0  # H2

product = fed - change(adsorption)
fed_h2 = fed[3] + change(pressurization)[3]
print("H2 purity", product[3] / product.sum())       # 0.9991
print("H2 recovery", (product[3] - purge_gas) / fed_h2)  # 0.162
```

N2 limits the purity: it adsorbs weakly on activated carbon and breaks through soon after H2, while CO2 stays in the
column. A longer adsorption step or less purge gas raises the recovery at the cost of the purity.
`examples/Self-driving-lab-H2-PSA` searches this trade-off.

## Model limits

- The whole column follows the pressure of the open end, and the time derivative of the pressure drop is left out. The
  pressure drop must stay small against the lowest pressure: the 40 s ramps keep `Max |R T sum(c) / P - 1|` in
  `run.log`, and the largest one at a phase end in column 5 of `phases.data`, below 0.008. Faster ramps raise it: with
  exponential ramps of 5 s it reaches 0.02, and the faster pressure change also lowers the purity, to 0.96, as the
  feed then enters faster than mass transfer takes it up.
- CVODE; ROS34PW2 gives the same purity (0.99913) and recovery (0.1618) in 13.5 s of CPU time instead of 21.5 s,
  with `"BreakthroughIntegrator": "ROS34PW2"` (since 2026-11-11, when its stages were projected on the closure).
- `heatCapacityGas` is a constant per unit mass, here that of the feed, while the gas swings between H2-rich and
  CO2-rich. The solid holds most of the heat, so that the result depends little on it: 14300 J/(kg K), the value of H2,
  in place of 1109 J/(kg K) gives a purity of 0.9991 either way, a recovery of 0.161 instead of 0.162, and gas
  temperatures within 2 K.
