#pragma once

#include "column.h"
#include "timing.h"

struct MultibedColumn;

namespace RK3Helpers
{
void updateVelocityAndPressure(Column& column);
void computeEquilibriumLoadings(Column& column);
void computeDerivatives(Column& column);
void computeSorptionDerivatives(Column& column);
void computePhysisorption(Column& column);
void computeChemisorption(Column& column);
void computeChemisorptionTransportDerivatives(Column& column);
void computeBulkSpeciesSink(Column& column);
void computeReactionDerivatives(Column& column);
void computeMassDerivatives(Column& column);
void computeEnergyDerivatives(Column& column);
bool reactionAutoStopReached(const Column& column, double timeStep) noexcept;
}  // namespace RK3Helpers

namespace RK3MultibedHelpers
{
void updateVelocityAndPressure(MultibedColumn& column);
void computeEquilibriumLoadings(MultibedColumn& column);
void computeDerivatives(MultibedColumn& column);
void computeSorptionDerivatives(MultibedColumn& column);
void computePhysisorption(MultibedColumn& column);
void computeChemisorption(MultibedColumn& column);
void computeChemisorptionTransportDerivatives(MultibedColumn& column);
void computeBulkSpeciesSink(MultibedColumn& column);
void computeReactionDerivatives(MultibedColumn& column);
void computeMassDerivatives(MultibedColumn& column);
void computeEnergyDerivatives(MultibedColumn& column);
bool reactionAutoStopReached(const MultibedColumn& column, double timeStep) noexcept;
}  // namespace RK3MultibedHelpers

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

  /**
   * \brief Advances a multibed column by one RK3 step.
   */
  bool propagate(MultibedColumn& column, size_t step, Timing& timings);
};
