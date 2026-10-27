#pragma once

#include "column.h"
#include "timing.h"

struct MultibedColumn;

void computeDerivatives(Column& column, double elapsedTime = 0.0);
void computeDerivatives(MultibedColumn& column, double elapsedTime = 0.0);
void precompute(Column& column, Timing& timings);
void precompute(MultibedColumn& column, Timing& timings);
bool reactionAutoStopReached(const Column& column) noexcept;
bool reactionAutoStopReached(const MultibedColumn& column) noexcept;

/**
 * \brief Third-order strong-stability-preserving Runge-Kutta integrator.
 *
 * Stores integration settings and advances a Column over a time interval in explicit RK3 steps.
 */
struct RungeKutta3
{
  /**
   * \brief Constructs the integrator from parsed input settings.
   */
  RungeKutta3(const InputReader& inputReader)
      : timeStep(inputReader.timeStep), numberOfInitTimeSteps(inputReader.numberOfInitTimeSteps) {};

  /**
   * \brief Constructs the integrator from explicit time-step settings.
   */
  RungeKutta3(double timeStep, size_t numberOfInitTimeSteps = 0)
      : timeStep(timeStep), numberOfInitTimeSteps(numberOfInitTimeSteps) {};

  double timeStep;               ///< Integration time step in s.
  size_t numberOfInitTimeSteps;  ///< Number of initial steps over which the step ramps up from 1e-4 timeStep.
  size_t numberOfStepsTaken{0};  ///< Steps taken so far; drives the initial ramp.

  /**
   * \brief Advances the column from startTime to endTime in s.
   *
   * Takes steps of timeStep (ramped during the first numberOfInitTimeSteps steps); the last step is shortened to
   * land on endTime unless the interval is a whole number of steps.
   */
  template <typename ColumnType>
  void advance(ColumnType& column, double startTime, double endTime, Timing& timings);

  /**
   * \brief Advances the column by one RK3 step of size dt in s, starting at startTime in s.
   */
  template <typename ColumnType>
  void step(ColumnType& column, double startTime, double dt, Timing& timings);
};

extern template void RungeKutta3::advance(Column& column, double startTime, double endTime, Timing& timings);
extern template void RungeKutta3::advance(MultibedColumn& column, double startTime, double endTime, Timing& timings);
extern template void RungeKutta3::step(Column& column, double startTime, double dt, Timing& timings);
extern template void RungeKutta3::step(MultibedColumn& column, double startTime, double dt, Timing& timings);
