#include "rk3.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <vector>

void updateStateRK(Column& column, Column& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.adsorption.size(); i++)
  {
    newColumn.adsorption[i] =
        alpha * column.adsorption[i] + beta * (newColumn.adsorption[i] + timeStep * newColumn.adsorptionDot[i]);
    newColumn.moleFraction[i] =
        std::max(0.0, alpha * column.moleFraction[i] +
                          beta * (newColumn.moleFraction[i] + timeStep * newColumn.moleFractionDot[i]));
  }
  if (column.energyBalance)
  {
    for (size_t grid = 0; grid < column.numberOfGridPoints + 1; grid++)
    {
      newColumn.gasTemperature[grid] =
          alpha * column.gasTemperature[grid] +
          beta * (newColumn.gasTemperature[grid] + timeStep * newColumn.gasTemperatureDot[grid]);
      newColumn.solidTemperature[grid] =
          alpha * column.solidTemperature[grid] +
          beta * (newColumn.solidTemperature[grid] + timeStep * newColumn.solidTemperatureDot[grid]);
      newColumn.wallTemperature[grid] =
          alpha * column.wallTemperature[grid] +
          beta * (newColumn.wallTemperature[grid] + timeStep * newColumn.wallTemperatureDot[grid]);
    }
  }
}

bool RungeKutta3::propagate(Column& column, size_t step, Timing& timings)
{
  auto totalTimer = timings.scoped(timings.total);

  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  if (autoNumberOfSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.moleFraction[numberOfGridPoints * numberOfComponents + j] /
                                                column.components[j].initialGasMoleFraction) -
                                               1.0));
    }

    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    if (tolerance < 0.01)
    {
      std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoNumberOfSteps = false;
    }
  }

  // column.writeJSON(std::format("{}_0.json", step));

  // SSP-RK Step 1
  // ======================================================================

  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P
  timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
  Column newColumn(column);

  updateStateRK(column, newColumn, 0.0, 1.0, timeStep);
  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });

  // Dqdt and Dp/dt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  timings.measure(timings.computeDerivatives, [&] { computeDerivatives(newColumn); });
  updateStateRK(column, newColumn, 0.75, 0.25, timeStep);
  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });

  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  // newColumn.writeJSON(std::format("{}_1.json", step));
  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // newColumn.writeJSON(std::format("{}_2.json", step));
  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  timings.measure(timings.computeDerivatives, [&] { computeDerivatives(newColumn); });
  updateStateRK(column, newColumn, (1.0 / 3.0), (2.0 / 3.0), timeStep);
  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });

  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // newColumn.writeJSON(std::format("{}_2.json", step));
  // update to the new time step
  column = newColumn;
  // enforceBoundaryCondition(column);
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}
