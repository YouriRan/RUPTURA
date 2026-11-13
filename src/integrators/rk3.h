#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "column.h"
#include "timing.h"

/**
 * \brief Evaluates the state derivatives (sorption rates, pore transport, reactions, bulk sink, mass and energy
 *        balances) from the current state and the fields set by precompute; elapsedTime in s.
 */
void computeDerivatives(Column& column, double elapsedTime = 0.0);

/**
 * \brief Updates the fields that depend on the state only: the bulk sink, velocity, pressure, mole fractions and
 *        equilibrium loadings.
 */
void precompute(Column& column, Timing& timings);

/**
 * \brief True when a column with reactions has reached a steady state (every state rate is small).
 */
bool reactionAutoStopReached(const Column& column) noexcept;

/**
 * \brief Fixed-step integration over a time interval, shared by RungeKutta3 and SemiImplicitRungeKutta3.
 *
 * During the first numberOfInitTimeSteps steps the step grows smoothly from 1e-4 timeStep to timeStep. Afterwards
 * whole steps of timeStep are taken, and the last step is shortened to land on endTime unless the interval is a whole
 * number of steps up to 1e-6. Step start times are computed, not accumulated, so they do not drift.
 */
struct FixedStepIntegrator
{
  double timeStep;               ///< Integration time step in s.
  size_t numberOfInitTimeSteps;  ///< Number of initial steps over which the step ramps up from 1e-4 timeStep.
  size_t numberOfStepsTaken{0};  ///< Steps taken so far; drives the initial ramp.

  /**
   * \brief Advances the column from startTime to endTime in s with self.step(column, start, dt, timings).
   */
  void advance(this auto& self, Column& column, double startTime, double endTime, Timing& timings)
  {
    constexpr double initialTimeStepFraction = 1.0e-4;

    double time = startTime;
    while (self.numberOfInitTimeSteps > 0 && self.numberOfStepsTaken <= self.numberOfInitTimeSteps)
    {
      // The first two steps both use the smallest step.
      const double progress = static_cast<double>(std::max<size_t>(self.numberOfStepsTaken, 1) - 1) /
                              static_cast<double>(self.numberOfInitTimeSteps);
      const double rampTimeStep =
          self.timeStep * (initialTimeStepFraction + (1.0 - initialTimeStepFraction) *
                                                         (3.0 * progress * progress - 2.0 * progress * progress * progress));
      const double remaining = endTime - time;
      if (remaining <= 0.0) return;
      self.step(column, time, std::min(rampTimeStep, remaining), timings);
      ++self.numberOfStepsTaken;
      if (rampTimeStep >= remaining) return;
      time += rampTimeStep;
    }

    const double numberOfSteps = (endTime - time) / self.timeStep;
    if (numberOfSteps <= 0.0) return;
    const bool wholeNumberOfSteps = std::abs(numberOfSteps - std::round(numberOfSteps)) <= 1.0e-6;
    const size_t fullSteps =
        static_cast<size_t>(wholeNumberOfSteps ? std::round(numberOfSteps) : std::floor(numberOfSteps));
    for (size_t i = 0; i < fullSteps; ++i)
    {
      self.step(column, time + static_cast<double>(i) * self.timeStep, self.timeStep, timings);
    }
    self.numberOfStepsTaken += fullSteps;
    if (!wholeNumberOfSteps)
    {
      const double lastStartTime = time + static_cast<double>(fullSteps) * self.timeStep;
      self.step(column, lastStartTime, endTime - lastStartTime, timings);
      ++self.numberOfStepsTaken;
    }
  }
};

/**
 * \brief Third-order strong-stability-preserving Runge-Kutta integrator (explicit).
 */
struct RungeKutta3 : FixedStepIntegrator
{
  RungeKutta3(const InputReader& inputReader)
      : FixedStepIntegrator{inputReader.timeStep, inputReader.numberOfInitTimeSteps}
  {
  }

  RungeKutta3(double timeStep, size_t numberOfInitTimeSteps = 0)
      : FixedStepIntegrator{timeStep, numberOfInitTimeSteps}
  {
  }

  /**
   * \brief Advances the column by one RK3 step of size dt in s, starting at startTime in s.
   */
  void step(Column& column, double startTime, double dt, Timing& timings);
};
