#include "rk3.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <vector>

bool RungeKutta3::propagate(BreakthroughState& state, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  size_t Ngrid = state.Ngrid;
  size_t Ncomp = state.Ncomp;
  BreakthroughState newState(state);

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(
          tolerance,
          std::abs((state.partialPressure[Ngrid * Ncomp + j] / (state.exitPressure * state.components[j].Yi0)) - 1.0));
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
  computeFirstDerivatives(state);

  // Dqdt and Dpdt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme
  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    newState.adsorption[i] = state.adsorption[i] + timeStep * state.adsorptionDot[i];
    newState.partialPressure[i] = state.partialPressure[i] + timeStep * state.pressureDot[i];
  }

  computeEquilibriumLoadings(newState);

  computeVelocity(newState);

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(newState);

  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    newState.adsorption[i] =
        0.75 * state.adsorption[i] + 0.25 * (newState.adsorption[i] + timeStep * newState.adsorptionDot[i]);
    newState.partialPressure[i] =
        0.75 * state.partialPressure[i] + 0.25 * (newState.partialPressure[i] + timeStep * newState.pressureDot[i]);
  }

  computeEquilibriumLoadings(newState);

  computeVelocity(newState);

  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep
  computeFirstDerivatives(newState);

  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    newState.adsorption[i] = (1.0 / 3.0) * state.adsorption[i] +
                             (2.0 / 3.0) * (newState.adsorption[i] + timeStep * newState.adsorptionDot[i]);
    newState.partialPressure[i] = (1.0 / 3.0) * state.partialPressure[i] +
                                  (2.0 / 3.0) * (newState.partialPressure[i] + timeStep * newState.pressureDot[i]);
  }

  computeEquilibriumLoadings(newState);

  computeVelocity(newState);

  // update to the new time step
  state = newState;

  // pulse boundary condition
  if (state.pulse)
  {
    if (t > state.pulseTime)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        if (j == state.carrierGasComponent)
        {
          state.partialPressure[0 * Ncomp + j] = state.externalPressure;
        }
        else
        {
          state.partialPressure[0 * Ncomp + j] = 0.0;
        }
      }
    }
  }

  return (!autoSteps && step >= numberOfSteps - 1);
}