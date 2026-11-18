#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "column.h"
#include "timing.h"

/**
 * \brief Evaluates the state derivatives (sorption rates, pore transport, reactions, bulk sink, mass and energy
 *        balances) of the state in views, with the fields that precompute set in column. The time enters through the
 *        pressure fields that Column::advancePressureRamp sets.
 *
 * views holds the state and receives its derivatives: the column's own, or a trial state of CVODE. With sorptionOnly
 * the evaluation stops after the bulk sink, without the velocity rebuild and the mass and energy balances.
 */
void computeDerivatives(ColumnData& column, const ColumnStateViews& views, bool sorptionOnly = false);

/**
 * \brief computeDerivatives of the column's own state.
 */
void computeDerivatives(Column& column);

/**
 * \brief Updates the fields that depend on the state in views only: the bulk sink (from the derivatives in views),
 *        velocity, pressure, mole fractions and equilibrium loadings.
 */
void precompute(ColumnData& column, const ColumnStateViews& views, Timing& timings);

/**
 * \brief precompute of the column's own state.
 */
void precompute(Column& column, Timing& timings);

/**
 * \brief True when a column with reactions has reached a steady state (every state rate is small).
 */
bool reactionAutoStopReached(const Column& column) noexcept;

/**
 * \brief Third-order strong-stability-preserving Runge-Kutta integrator (explicit), with a fixed step.
 *
 * During the first numberOfInitTimeSteps steps the step grows smoothly from 1e-4 timeStep to timeStep. Afterwards
 * whole steps of timeStep are taken, and the last step is shortened to land on endTime unless the interval is a whole
 * number of steps up to 1e-6. Step start times are computed, not accumulated, so they do not drift.
 *
 * The step is not checked for stability. An explicit step is stable below the time in which the fluid crosses the
 * smallest grid cell, dz / v, below about 2.5 / k for the largest mass-transfer coefficient k, and below about
 * 0.6 dz^2 / D for the largest axial dispersion coefficient D; a larger step gives wrong results without an error.
 */
struct RungeKutta3
{
  RungeKutta3(const InputReader& inputReader)
      : timeStep(inputReader.timeStep), numberOfInitTimeSteps(inputReader.numberOfInitTimeSteps)
  {
  }

  RungeKutta3(double timeStep, size_t numberOfInitTimeSteps = 0)
      : timeStep(timeStep), numberOfInitTimeSteps(numberOfInitTimeSteps)
  {
  }

  double timeStep;               ///< Integration time step in s.
  size_t numberOfInitTimeSteps;  ///< Number of initial steps over which the step ramps up from 1e-4 timeStep.
  size_t numberOfStepsTaken{0};  ///< Steps taken so far; drives the initial ramp.

  /**
   * \brief Advances the column from startTime to endTime in s.
   */
  void advance(Column& column, double startTime, double endTime, Timing& timings)
  {
    // Timed per call, not per step, to keep the clock calls out of the step loop.
    auto totalTimer = timings.scoped(timings.total);
    constexpr double initialTimeStepFraction = 1.0e-4;

    double time = startTime;
    while (numberOfInitTimeSteps > 0 && numberOfStepsTaken <= numberOfInitTimeSteps)
    {
      // The first two steps both use the smallest step.
      const double progress = static_cast<double>(std::max<size_t>(numberOfStepsTaken, 1) - 1) /
                              static_cast<double>(numberOfInitTimeSteps);
      const double rampTimeStep =
          timeStep * (initialTimeStepFraction + (1.0 - initialTimeStepFraction) *
                                                    (3.0 * progress * progress - 2.0 * progress * progress * progress));
      const double remaining = endTime - time;
      if (remaining <= 0.0) return;
      step(column, std::min(rampTimeStep, remaining), timings);
      ++numberOfStepsTaken;
      if (rampTimeStep >= remaining) return;
      time += rampTimeStep;
    }

    const double numberOfSteps = (endTime - time) / timeStep;
    if (numberOfSteps <= 0.0) return;
    const bool wholeNumberOfSteps = std::abs(numberOfSteps - std::round(numberOfSteps)) <= 1.0e-6;
    const size_t fullSteps =
        static_cast<size_t>(wholeNumberOfSteps ? std::round(numberOfSteps) : std::floor(numberOfSteps));
    for (size_t i = 0; i < fullSteps; ++i)
    {
      step(column, timeStep, timings);
    }
    numberOfStepsTaken += fullSteps;
    if (!wholeNumberOfSteps)
    {
      const double lastStartTime = time + static_cast<double>(fullSteps) * timeStep;
      step(column, endTime - lastStartTime, timings);
      ++numberOfStepsTaken;
    }
  }

  /**
   * \brief Advances the column by one RK3 step of size dt in s.
   */
  void step(Column& column, double dt, Timing& timings);
};
