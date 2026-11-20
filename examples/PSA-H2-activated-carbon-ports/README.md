# H2 purification by a PSA cycle with ports (MomentumBalance)

The cycle of `examples/PSA-H2-activated-carbon` with `BoundaryCondition` `MomentumBalance`: the pressure of every cell
is that of its gas, and the gas flows between cells and through the ends by the pressure-drop law of the bed. Each
phase sets the two ends with ports.

| Step | `FeedEnd` (z = 0) | `ProductEnd` (z = L) | Duration |
| --- | --- | --- | --- |
| Pressurization | pressure, 1 to 10 bar, linear ramp | closed | 40 s |
| Adsorption | velocity, 0.03 m/s of feed at 10 bar | pressure, 10 bar | 60 s |
| Blowdown | pressure, 10 to 1 bar, linear ramp | closed | 40 s |
| Purge | pressure, 1 bar | velocity, 0.1 m/s at 1 bar | 60 s |

The purge gas is the product: `"ProductEndGasPhaseMolFractions": "Adsorption"` takes the mean gas that left the column
in the adsorption step. `CyclicSteadyState` ends the run once a cycle ends within 1e-3 of the one before (20 cycles
here; `NumberOfCycles` 40 is the most), and `cycles.data` lists the changes per cycle. Every component has a
`GasHeatCapacity`, which the energy balance of this boundary condition needs. With a purge of pure H2
(`"ProductEndGasPhaseMolFractions": {"H2": 1.0}`) the cycle gives the purity and recovery of
`examples/PSA-H2-activated-carbon`, 0.9991 and 0.162.

## Purity and recovery

`phases.data` (or `SwingAdsorption.phase_records[i].stream_amounts` in Python) holds what entered and left at each end
in every phase, per unit column cross-section:

```python
import numpy as np

n = 4  # He, CO2, N2, H2
table = np.loadtxt("phases.data")
streams = table[:, 5 + 4 * n:].reshape(len(table), 4, n)  # entered z = 0, left z = 0, entered z = L, left z = L
pressurization, adsorption, blowdown, purge = streams[-4:]

product = adsorption[3]
fed = pressurization[0] + adsorption[0]
print("H2 purity", product[3] / product.sum())                  # 0.9988
print("H2 recovery", (product[3] - purge[2][3]) / fed[3])       # 0.163
```
