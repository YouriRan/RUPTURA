#include "rk3.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <vector>

void updateStateRK(Column& column, Column& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.physisorption.size(); i++)
  {
    newColumn.physisorption[i] =
        alpha * column.physisorption[i] +
        beta * (newColumn.physisorption[i] + timeStep * newColumn.physisorptionDot[i]);
    newColumn.concentration[i] =
        std::max(0.0, alpha * column.concentration[i] +
                          beta * (newColumn.concentration[i] + timeStep * newColumn.concentrationDot[i]));
  }

  for (size_t i = 0; i < column.chemisorption.size(); ++i)
  {
    newColumn.chemisorption[i] =
        alpha * column.chemisorption[i] +
        beta * (newColumn.chemisorption[i] + timeStep * newColumn.chemisorptionDot[i]);
  }

  for (size_t i = 0; i < column.surfaceConcentration.size(); ++i)
  {
    newColumn.surfaceConcentration[i] =
        std::max(0.0, alpha * column.surfaceConcentration[i] +
                          beta * (newColumn.surfaceConcentration[i] +
                                  timeStep * newColumn.surfaceConcentrationDot[i]));
    newColumn.poreConcentration[i] =
        std::max(0.0, alpha * column.poreConcentration[i] +
                          beta * (newColumn.poreConcentration[i] + timeStep * newColumn.poreConcentrationDot[i]));
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
      if (column.components[j].initialGasMoleFraction <= 0.0) continue;

      const size_t outlet = numberOfGridPoints * numberOfComponents + j;
      const double feed = column.components[j].initialGasMoleFraction;
      tolerance = std::max(tolerance, std::abs((column.moleFraction[outlet] / feed) - 1.0));
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
  timings.measure(timings.computeDerivatives,
                  [&]
                  {
                    computeSorptionDerivatives(column);
                    computeDerivatives(column);
                  });
  Column newColumn(column);

  updateStateRK(column, newColumn, 0.0, 1.0, timeStep);
  timings.measure(timings.updateVelocityAndPressure, [&] { updateVelocityAndPressure(newColumn); });

  // Dqdt and Dp/dt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  timings.measure(timings.computeDerivatives,
                  [&]
                  {
                    computeSorptionDerivatives(newColumn);
                    computeDerivatives(newColumn);
                  });
  updateStateRK(column, newColumn, 0.75, 0.25, timeStep);
  timings.measure(timings.updateVelocityAndPressure, [&] { updateVelocityAndPressure(newColumn); });

  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  // newColumn.writeJSON(std::format("{}_2.json", step));
  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  timings.measure(timings.computeDerivatives,
                  [&]
                  {
                    computeSorptionDerivatives(newColumn);
                    computeDerivatives(newColumn);
                  });
  updateStateRK(column, newColumn, (1.0 / 3.0), (2.0 / 3.0), timeStep);
  timings.measure(timings.updateVelocityAndPressure, [&] { updateVelocityAndPressure(newColumn); });

  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });

  // newColumn.writeJSON(std::format("{}_2.json", step));
  // update to the new time step
  column = newColumn;
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}
