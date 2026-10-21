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
      : timeStep(inputReader.timeStep),
        autoNumberOfSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps) {};

  /**
   * \brief Constructs the integrator from explicit time-step settings.
   */
  SemiImplicitRungeKutta3(double timeStep, bool autoNumberOfSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoNumberOfSteps(autoNumberOfSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;         ///< Integration time step in s.
  bool autoNumberOfSteps;  ///< Continue until breakthrough criterion when true.
  size_t numberOfSteps;    ///< Requested number of integration steps.

  /**
   * \brief Advances the column by one semi-implicit RK3 step.
   */
  bool propagate(Column& column, size_t step, Timing& timings);
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
