#include "rk3.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <vector>

void updateStateRK(Column& column, Column& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.adsorption.size(); i++)
  {
    newColumn.adsorption[i] =
        alpha * column.adsorption[i] + beta * (newColumn.adsorption[i] + timeStep * newColumn.adsorptionDot[i]);
    newColumn.partialPressure[i] = alpha * column.partialPressure[i] +
                                   beta * (newColumn.partialPressure[i] + timeStep * newColumn.partialPressureDot[i]);
  }
  if (column.energyBalance)
  {
    for (size_t grid = 0; grid < column.Ngrid + 1; grid++)
    {
      // newColumn.totalPressure[grid] =
      //     alpha * column.totalPressure[grid] +
      //     beta * (newColumn.totalPressure[grid] + timeStep * newColumn.totalPressureDot[grid]);
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

bool RungeKutta3::propagate(Column& column, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.partialPressure[Ngrid * Ncomp + j] /
                                                (column.exitPressure * column.components[j].Yi0)) -
                                               1.0));
    }

    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  // SSP-RK Step 1
  // ======================================================================

  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P
  computeFirstDerivatives(column);
  Column newColumn(column);

  updateStateRK(column, newColumn, 0.0, 1.0, timeStep);

  // Dqdt and Dpdt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme
  computeEquilibriumLoadings(newColumn);

  computeVelocity(newColumn);

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(newColumn);
  updateStateRK(column, newColumn, 0.75, 0.25, timeStep);

  computeEquilibriumLoadings(newColumn);

  computeVelocity(newColumn);

  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(newColumn);
  updateStateRK(column, newColumn, (1.0 / 3.0), (2.0 / 3.0), timeStep);

  computeEquilibriumLoadings(newColumn);

  computeVelocity(newColumn);

  // update to the new time step
  column = newColumn;

  // pulse boundary condition
  for (size_t j = 0; j < Ncomp; ++j)
  {
    if (column.pulse && t > column.pulseTime)
    {
      if (j == column.carrierGasComponent)
      {
        column.partialPressure[0 * Ncomp + j] = column.externalPressure;
      }
      else
      {
        column.partialPressure[0 * Ncomp + j] = 0.0;
      }
    }
    else
    {
      column.partialPressure[0 * Ncomp + j] = column.externalPressure * column.components[j].Yi0;
    }
  }

  return (!autoSteps && step >= numberOfSteps - 1);
}