# Figure benchmark matrices

This package generates the data for the notebook's CVODE/RK3 comparison and
Marx-only mixture-method bar plots. The integrator matrix covers a common
0.1-second physical horizon for both base cases, seven stable time steps, the requested
grid sweep, and the requested adsorbing-component-count sweep. Carrier gas is
not counted in `ncomp`; when adsorbates are removed, their inlet fraction is
returned to the carrier.

BEA uses `1e-5`–`1e-2` s. Marx uses `1e-6`–`1e-3` s because its explicit RK3
solution becomes unstable above `1e-3` s for this matrix, while CVODE remains
stable.

The method matrix contains 32 Marx simulations: IAST/SIAST/SEI/SPI crossed with
RK3/CVODE, energy off/on, and transported chemisorption off/on. Transported
chemisorption adds a General kinetic site to the Marx CO2 component while
retaining its original physisorption site.

Because SEI requires single-site Langmuir isotherms, the method matrix projects
each Marx Sips site to its Langmuir limit by retaining the original saturation
capacity and affinity and removing only the Sips exponent. This projection is
applied consistently to all four methods and only to the method matrix.

Run the complete workflow from the repository root:

```sh
python3 benchmark/figures/generate.py
python3 benchmark/figures/run.py --suite cvode-rk3 --repeats 3
python3 benchmark/figures/run.py --suite methods --repeats 3
python3 benchmark/figures/summarize.py
```
