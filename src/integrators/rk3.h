#pragma once

#include "column.h"
#include "compute.h"
#include "timing.h"

/**
 * \brief Third-order strong-stability-preserving Runge-Kutta integrator.
 *
 * Stores integration settings and advances a Column by one explicit RK3 step.
 */
struct RungeKutta3
{
  /**
   * \brief Constructs the integrator from parsed input settings.
   */
  RungeKutta3(const InputReader& inputReader)
      : timeStep(inputReader.timeStep),
        autoNumberOfSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps) {};

  /**
   * \brief Constructs the integrator from explicit time-step settings.
   */
  RungeKutta3(double timeStep, bool autoNumberOfSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoNumberOfSteps(autoNumberOfSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;         ///< Integration time step in s.
  bool autoNumberOfSteps;  ///< Continue until breakthrough criterion when true.
  size_t numberOfSteps;    ///< Requested number of integration steps.

  /**
   * \brief Advances the column by one RK3 step.
   */
  bool propagate(Column& column, size_t step, Timing& timings);
};
