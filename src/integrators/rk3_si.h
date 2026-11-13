#pragma once

#include "column.h"
#include "rk3.h"
#include "timing.h"

/**
 * \brief Legacy split semi-implicit column integrator, selected by SIRK3.
 *
 * Advances the physisorbed loadings (LDF, implicit) and the gas concentrations (upwind convection and axial
 * dispersion, implicit) of a single adsorbent; the input reader rejects everything else (chemisorption, reactions,
 * surface-pore transport, the energy balance, WENO3). Not third order even for constant-equilibrium scalar LDF
 * relaxation.
 */
struct SemiImplicitRungeKutta3 : FixedStepIntegrator
{
  SemiImplicitRungeKutta3(const InputReader& inputReader)
      : FixedStepIntegrator{inputReader.timeStep, inputReader.numberOfInitTimeSteps}
  {
  }

  SemiImplicitRungeKutta3(double timeStep, size_t numberOfInitTimeSteps = 0)
      : FixedStepIntegrator{timeStep, numberOfInitTimeSteps}
  {
  }

  /**
   * \brief Advances the column by one semi-implicit RK3 step of size dt in s, starting at startTime in s.
   */
  void step(Column& column, double startTime, double dt, Timing& timings);
};
