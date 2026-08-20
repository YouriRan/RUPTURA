#pragma once

#include "column.h"
#include "timing.h"

struct MultibedColumn;

void computeDerivatives(Column& column);
void computeDerivatives(MultibedColumn& column);
void precompute(Column& column, Timing& timings);
void precompute(MultibedColumn& column, Timing& timings);
bool reactionAutoStopReached(const Column& column, double timeStep) noexcept;
bool reactionAutoStopReached(const MultibedColumn& column, double timeStep) noexcept;

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
  template <typename ColumnType>
  bool propagate(ColumnType& column, size_t step, Timing& timings);
};

extern template bool RungeKutta3::propagate(Column& column, size_t step, Timing& timings);
extern template bool RungeKutta3::propagate(MultibedColumn& column, size_t step, Timing& timings);
