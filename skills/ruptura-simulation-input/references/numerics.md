# Numerical settings

You may choose these yourself; report every choice under "Assumptions and conversions".

## Integrator and time step

`BreakthroughIntegrator`:

- `"CVODE"`: adaptive BDF; robust and the recommended default.
- `"ROS34PW2"`: adaptive, linearly implicit Rosenbrock-W method with the analytical sorption
  Jacobian; often much faster than CVODE, all column features. The integrator for
  `MomentumBalance` (CVODE works there only with `CVODELinearSolver` `"KLU"`).
- `"RungeKutta3"`: explicit; the parser default.

`TimeStep` [s] (`RungeKutta3` only):

- The integration step. It must satisfy roughly TimeStep <= 0.5 dz / v, TimeStep <= dz^2 / (2 D_ax),
  and TimeStep <= 2 / k for the largest `MassTransferCoefficient` k, where
  dz = ColumnLength / NumberOfGridPoints. With the energy balance also
  TimeStep <= dz^2 rho c_p / (2 lambda), with the gas density rho, `heatCapacityGas` c_p, and
  `gasThermalConductivity` lambda.
- The step is not checked: a larger `TimeStep` gives wrong curves without an error.
- CVODE and ROS34PW2 choose their own internal steps; an input with either and `TimeStep` is
  rejected.

## Grid and discretization

- `NumberOfGridPoints`: 100 is a good default; use at least 50.
- `ConvectionScheme`: `"WENO3"` (default; third-order face values, little numerical dispersion),
  `"WENO5"` (fifth order on a uniform grid, WENO3 next to the column ends; slightly sharper fronts,
  more steps with ROS34PW2), `"TVD"` (second order, minmod limiter), or `"Upwind"` (first order;
  adds a numerical dispersion of about v dz / 2). Leave it out.
- `InletCondition`: `"Danckwerts"` (default; the total inflow equals the feed flux) or
  `"Dirichlet"` (the inlet face holds the feed state, so axial dispersion adds to the feed flux; it
  overfeeds when ColumnEntranceVelocity x ColumnLength / AxialDispersionCoefficient is below about
  100). Leave it out.

## Run length and output

- `SimulationTime` [s]:
  - `"auto"` runs until every outlet mole fraction is within 1% of the feed (with reactions: until
    the column states change by less than 1e-4 of their value per second), then 10% longer. Use it
    for adsorption breakthrough.
  - Use a time for desorption, purge, a fixed duration, or any run whose outlet never reaches the
    feed.
  - SwingAdsorption takes no root `SimulationTime`: the phase durations set it.
- `WriteInterval` [s]: one column profile and one progress line every `WriteInterval` (default 5);
  keep it well below the breakthrough time, aiming for a few hundred output rows.
- `WriteOutletOnly` (true/false, default false): write only the outlet row of every output (z = L,
  or the first cell in a reverse-flow phase) instead of the whole column; use it for long runs whose
  column profiles are not needed. Leave it out otherwise.
- `NumberOfInitTimeSteps`: RungeKutta3 ramp-up steps (default 0). Leave it out.
- `NumberOfTimeSteps`, `WriteEvery`, and `PrintEvery` (RUPTURA 2) no longer exist; the reader
  rejects them and prints the equivalent times.

## Solver tuning (add only on request)

CVODE (defaults in parentheses):

- `CVODERelativeTolerance` (1e-5); `CVODEAbsoluteToleranceConcentration` [mol/m3] (1e-6);
  `CVODEAbsoluteToleranceLoading` [mol/kg] (1e-7); `CVODEAbsoluteToleranceTemperature` [K] (1e-4).
- `CVODEMaximumTimeStep` [s] (-1 = 4 dz/v automatic, 0 = unlimited); `CVODEMaximumOrder` (5, BDF
  order 1-5; 1 if a steep front stalls).
- `CVODELinearSolver` `"SPGMR"` (default), `"Dense"` or `"KLU"`; `CVODEKrylovDimension` (30, SPGMR
  only).
- `CVODEPreconditioner` `"None"` (default) or `"SorptionJacobian"` (banded I - gamma J; SPGMR only,
  not with MPD).

ROS34PW2:

- `ROS34PW2RelativeTolerance` (1e-6); `ROS34PW2AbsoluteToleranceConcentration` [mol/m3] (1e-7);
  `ROS34PW2AbsoluteToleranceLoading` [mol/kg] (1e-8); `ROS34PW2AbsoluteToleranceTemperature` [K]
  (1e-5); `ROS34PW2MaximumTimeStep` [s] (0 = unlimited).

For dilute feeds (ppm level), lower `CVODEAbsoluteToleranceConcentration` well below the feed
concentration, and report it.
