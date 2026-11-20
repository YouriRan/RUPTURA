# Liu et al. (2024) minimal DAC kinetic model

Reproduction of Figures 1-7 of:

> H. Liu, H. Lin, S. Dai, and D.-e. Jiang, "Minimal Kinetic Model of Direct Air Capture of
> CO$_2$ by Supported Amine Sorbents in Dry and Humid Conditions," *Ind. Eng. Chem. Res.*
> **63** (2024) 5871-5879. DOI: 10.1021/acs.iecr.3c04535.

Every model curve in `liu_paper_figures.ipynb` is read from a RUPTURA simulation under
`simulations/`. The notebook computes no adsorption physics of its own and loads no stored
result table; it opens `component_*.data` files and plots them. Every reference curve is a
digitization of a published panel, kept in `reference_digitized/`.

## Running it

```sh
python generate.py          # write simulations/<case>/simulation.json and simulations/index.json
python run.py               # run them all, 8 workers, skipping any that already finished
jupyter lab liu_paper_figures.ipynb
```

`generate.py` is idempotent and does not delete existing RUPTURA output; pass `--clean` to
start from an empty `simulations/`. `run.py` skips finished directories unless given
`--rerun`, so an interrupted queue can simply be restarted.

Expect the full queue to take on the order of 100 core-hours. The cost is inherent: the
Toth isotherm with $t_T < 1$ is only once differentiable at $p \to 0$, so CVODE drops to
low order at the leading edge and needs tens of internal steps per simulated second, while
the figures run out to 132000 s. Raising the axial dispersion would remove the stiffness
but moves the 50% breakthrough time by tens of percent, so it is not used.

## Case map

| Figure | Cases | Simulation type | Simulated time |
| --- | --- | --- | --- |
| 1a | 2 | MixturePrediction | - |
| 1b | 1 | Breakthrough | 1800 s |
| 2a, 2b, 2c | 5 + 6 + 5 | Breakthrough | 25200 / 114000 / 114000 s |
| 3a, 3b, 3c | 6 + 5 + 5 | Breakthrough | 132000 / 43200 / 36000 s |
| 4 | 5 | Breakthrough, CO$_2$ + H$_2$O | 90000 s |
| 5a | 1 | Breakthrough, H$_2$O | 600 s |
| 5b | 4 | MixturePrediction | - |
| 6a, 6b, 6c | 3 + 3 + 3 | Breakthrough | 18000 s |
| 7 | 2 | Breakthrough, one with H$_2$O | 18000 s |

`simulations/index.json` maps each directory to its figure, panel, legend label and the
parameters needed to match it against the digitized reference. The notebook looks cases up
through that file rather than hard-coding directory names.

## Numerical settings

All breakthrough runs use the ROS34PW2 integrator at its default tolerances (relative 1e-6, absolute
1e-7 mol m$^{-3}$) with the Danckwerts inlet, an axial dispersion coefficient of 1e-5 m$^2$ s$^{-1}$, and
`energyBalance: false`, matching the paper's isothermal treatment. Thirty grid points agree with eighty to
within 0.6% on the breakthrough times of Figure 1b; Figure 4 keeps the 15 grid points it had with CVODE.
All 50 runs take about 3 s in total, and every curve is monotonic to rounding (worst backward step about
1e-16 in the Figure 4 cells).

Until 2026-11-03 the runs used CVODE. These DAC feeds needed an uncapped internal step and an absolute
concentration tolerance of 1e-6 instead of the 1e-5 default, which put non-monotonic steps of a few percent
into the curves. Figure 4 (affinity up to 57.6 Pa$^{-1}$ against a 0.016 mol m$^{-3}$ feed, a partition
coefficient near 1e8) was tractable only on 15 grid points with an absolute tolerance of 0.6% of the feed,
about 8.5 h per case. 36 of the 50 stored CVODE runs stopped with `mxstep` (TODO B-EX-10).

## Fidelity limits

- **The humid CO$_2$ laws follow the local humidity (since 2026-11-04).** Figure 4 uses
  RUPTURA's `GAB-modified-Toth` site, whose $n_s$ and $b$ follow the GAB water loading at
  the local water partial pressure, the equilibrium loading rather than the transported
  one. Figures 5b and 7 use the `Exponential-humidity-Toth` site for equation 8. The H$_2$O
  component is marked `"Water": true`. Before, each run froze the law at its feed humidity;
  the coupled Figure 4 curves differ from the frozen ones by at most $1.2\times10^{-4}\,c_0$
  with the same 50% breakthrough times, since water saturates the bed in 2000-9000 s while
  the CO$_2$ front needs 90000 s (`conversions.md`, sections 9 and 10).
- **Figure 5a's water isotherm comes from a short breakthrough, not MixturePrediction.**
  A 600 s water breakthrough into a dry module reports, at every grid point, the local
  water partial pressure beside RUPTURA's own GAB loading at that pressure, and the
  developing front sweeps the whole 0 to 0.85 humidity range. The points plotted are
  therefore still RUPTURA's own isotherm evaluations. (Until 2026-11-04 a MixturePrediction
  run could not write a GAB site; it can now.)
- **Reference traces are graphical reconstructions.** The article publishes no
  machine-readable series. The curves and markers in `reference_digitized/` were recovered
  from the printed panels and are visual comparison controls, not the authors' data. No
  RUPTURA input is fitted against them.
- **Some apparatus dimensions are calibrated against RUPTURA.** Liu et al. do not repeat
  the Figure 1 benchmark dimensions from their reference 21, nor the loaded module masses
  for Figures 6-7, so the corresponding bed densities are inferred from the printed time
  axes. Because they are inferred, they depend on the solver they were inferred with, and
  the values inherited from the earlier Python implementation do not carry over. They were
  recalibrated here by scaling a trial density with the ratio of the digitized to the
  simulated 50% breakthrough time:

  | Case | Trial | Ratio | Calibrated |
  | --- | --- | --- | --- |
  | Figure 6, geometric mean of all nine curves | 80 kg m$^{-3}$ | 1.197 (spread 1.043-1.295) | 95.7 kg m$^{-3}$ |
  | Figure 7, dry curve | 110 kg m$^{-3}$ | 1.568, refined once | 160.9 kg m$^{-3}$ |

  Breakthrough time is not exactly proportional to the density, because the
  linear-driving-force spreading does not scale with it, so Figure 7 needed a second pass.
  Residuals after calibration are 0.993 for Figure 6 as a geometric mean over its nine
  curves, with individual curves within 14%, and 1.010 for the Figure 7 dry curve. The
  Figure 1b length and density needed no adjustment: RUPTURA puts its 50% breakthrough at
  294 s against roughly 300 s digitized.

- **Figure 7's dry and humid curves do not share one density.** Calibrating on the dry
  curve gives 160.9 kg m$^{-3}$; calibrating on the humid curve alone would give roughly
  105 kg m$^{-3}$. RUPTURA's humid-to-dry ratio of 50% breakthrough times is 2.60 against
  1.84 digitized, because equation 8's $e^{1.11\varphi}$ raises the CO$_2$ capacity by a
  factor 2.569 wherever the gas is at the feed humidity $\varphi = 0.85$. Letting the
  enhancement rise with the local humidity as the water front passes, instead of freezing
  it from $t = 0$, moved the humid 50% breakthrough by only $-1.5\%$ (ratio 2.64 before),
  so the gap is in equation 8, not in its coupling. The dry calibration is used because it
  involves no humid-model assumption, so the humid curve is left visibly late rather than
  split the difference and bias both. A further RUPTURA limitation pushes the same way: the
  bed cannot be started pre-loaded with water the way the experiment's pre-humidification
  protocol does.
- **The discretization is not the paper's.** Liu et al. used MATLAB `ode15s` with
  second-order centred differences. RUPTURA solves the full mass balance including gas
  accumulation, which the paper drops in its quasi-steady plug-flow limit. Breakthrough
  times agree closely; the approach to $c/c_0 = 1$ is slightly sharper here.

Unit conversions between the article and the RUPTURA inputs are derived in `conversions.md`.

## Superseded files

The following predate this pipeline and are no longer referenced by the notebook. They
implemented the paper's equations directly in Python and wrote their results to
`runNN_figureN/results.csv`, which is what made the plotted data untraceable. They are kept
only until the RUPTURA results have been checked against them, and can be deleted after
that.

`model.py`, `reproduce.py`, `generate_runs.py`, `validate_all.py`, `run_all.sh`,
`requirements.txt`, `run01_figure1/` through `run07_figure7/`.
