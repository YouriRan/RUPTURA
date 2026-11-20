# H2 purification by a two-bed PSA cycle with pressure equalization

The beds and the feed of `examples/PSA-H2-activated-carbon-ports`, as two columns A and B that run
half a cycle apart (`Columns`). Valves between the product ends and a product tank (`Connections`,
`Buffers`) take the place of the set product and purge streams of the one-column cycle.

| Step | A | B | Connections | Duration |
| --- | --- | --- | --- | --- |
| Pressurization A, blowdown B | feed end: pressure, to 10 bar, linear ramp | feed end: pressure, to 1 bar, linear ramp | | 40 s |
| Adsorption A, purge B | feed end: velocity, 0.03 m/s of feed at 10 bar | feed end: pressure, 1 bar | A to tank (check valve), A to B | 40 s |
| Equalization A to B | closed | closed | A to B | 20 s |

The three steps then run with A and B exchanged. All connections are at the product ends and
are linear valves, n = k ΔP. The tank holds 10 bar (`MaximumPressure`): what it takes in above
that pressure leaves it as the product. The equalization passes the gas of the bed that has
adsorbed to the bed that is pressurized next, which the one-column cycle loses in its blowdown.
`CyclicSteadyState` ends the run after 23 cycles; the two halves of the last cycle then agree
within 1e-4 of the largest stream.

## Purity and recovery

`buffers.data` holds, per phase, what left the tank, and `A/phases.data` and `B/phases.data`
what entered each column at its feed end per unit cross-section:

```python
import numpy as np

n = 4  # He, CO2, N2, H2
area = np.pi * 0.025**2 / 4
product = np.loadtxt("buffers.data")[-6:, 2 + n:2 + 2 * n].sum(axis=0)
fed = sum(np.loadtxt(f"{bed}/phases.data")[-6:, 5 + 4 * n:5 + 5 * n].sum(axis=0) for bed in "AB") * area
print("H2 purity", product[3] / product.sum())   # 0.9991
print("H2 recovery", product[3] / fed[3])        # 0.363
```

The one-column cycle gives 0.99993 and 0.170 with an adsorption step of 60 s: the equalization
returns hydrogen, and with it some of the nitrogen ahead of the front, to the other bed.
