#pragma once

#include "column.h"
#include "compute.h"
#include "timing.h"

/**
 * \brief Legacy split semi-implicit column integrator, selected by SIRK3.
 *
 * This column path predates the verified ROS3 kernel in ros3.h and is not third order
 * even for constant-equilibrium scalar LDF relaxation. Migration requires the shared
 * column RHS and its analytical Jacobian; see docs/SIRK3_ANALYTIC_JACOBIAN_PLAN.md.
 */
struct SemiImplicitRungeKutta3
{
  /**
   * \brief Constructs the integrator from parsed input settings.
   */
  SemiImplicitRungeKutta3(const InputReader& inputReader)
      : timeStep(inputReader.timeStep), numberOfInitTimeSteps(inputReader.numberOfInitTimeSteps) {};

  /**
   * \brief Constructs the integrator from explicit time-step settings.
   */
  SemiImplicitRungeKutta3(double timeStep, size_t numberOfInitTimeSteps = 0)
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
  void advance(Column& column, double startTime, double endTime, Timing& timings);

  /**
   * \brief Advances the column by one semi-implicit RK3 step of size dt in s, starting at startTime in s.
   */
  void step(Column& column, double startTime, double dt, Timing& timings);
};

/**
 * \brief Builds and solves the semi-implicit concentration update matrix.
 */
void computeConcentrationUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved);

/**
 * \brief Builds and solves the final-stage semi-implicit concentration update matrix.
 */
void computeConcentrationUpdateMatrixFinal(Column& column, double timeStep, std::vector<double>& solved);

/**
 * \brief Builds and solves the semi-implicit concentration update matrix with energy-balance coupling.
 */
void computeConcentrationUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved);

/**
 * \brief Builds and solves the final-stage concentration update matrix with energy-balance coupling.
 */
void computeConcentrationUpdateMatrixEnergyBalanceFinal(Column& column, double timeStep, std::vector<double>& solved);
