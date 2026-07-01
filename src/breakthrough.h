#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "column.h"
#include "component.h"
#include "cvode.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "rk3_si.h"
#include "timing.h"

/**
 * \brief Simulates a breakthrough process in an adsorption column.
 *
 * The Breakthrough struct encapsulates the parameters and methods required to simulate
 * the breakthrough of gases in an adsorption column. It handles the initialization of
 * simulation parameters, computation of time steps, and generation of output scripts
 * for plotting and visualization.
 */
struct Breakthrough
{
  /**
   * \brief Time-integration scheme used for breakthrough propagation.
   */
  enum class IntegrationScheme
  {
    SSP_RK = 0,     ///< Strong Stability Preserving Runge-Kutta method.
    CVODE = 1,      ///< CVODE integration.
    Iterative = 2,  ///< Iterative integration scheme.
    SIRK3 = 3,      ///< Semi-implicit third-order Runge-Kutta method.
  };

  /**
   * \brief Constructs a Breakthrough simulation from explicit simulation objects and settings.
   *
   * This constructor is the argument-based alternative to constructing from an InputReader.
   */
  Breakthrough(std::string displayName, size_t carrierGasComponent, size_t numberOfComponents,
               size_t numberOfGridPoints, size_t printEvery, size_t writeEvery, double timeStep,
               size_t numberOfInitTimeSteps, size_t numberOfTimeSteps, bool autoNumberOfTimeSteps,
               size_t maxIsothermTerms, Column column, RungeKutta3 rk3, SemiImplicitRungeKutta3 sirk3, CVODE cvode,
               IntegrationScheme integrationScheme, std::optional<std::string> readColumnFile = std::nullopt);

  /**
   * \brief Constructs a Breakthrough simulation using an InputReader.
   *
   * Initializes a Breakthrough object with parameters specified in the InputReader.
   *
   * \param inputreader Reference to an InputReader containing simulation parameters.
   */
  Breakthrough(const InputReader& inputreader);

  /**
   * \brief Prints the representation of the Breakthrough object to the console.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the Breakthrough object.
   *
   * \return A string representing the Breakthrough object.
   */
  std::string repr() const;

  /**
   * \brief Runs the Breakthrough simulation.
   *
   * Executes the simulation over the specified number of time steps.
   */
  void run();

  /**
   * \brief Computes a single simulation step.
   *
   * Advances the simulation by one time step.
   *
   * \param step The current time step index.
   */
  void computeStep(size_t step);

  const std::string displayName;  ///< Name of the simulation for display purposes.
  size_t carrierGasComponent{0};  ///< Index of the carrier gas component.
  size_t numberOfComponents;      ///< Number of components.
  size_t numberOfGridPoints;      ///< Number of grid intervals; node count is numberOfGridPoints + 1.

  size_t printEvery;  ///< Frequency of printing time steps to the screen.
  size_t writeEvery;  ///< Frequency of writing data to files.

  double timeStep;               ///< Time step for integration in s.
  size_t numberOfInitTimeSteps;  ///< Number of initial ramp-up steps.
  size_t numberOfSteps;          ///< Total number of steps.
  bool autoNumberOfSteps;        ///< Flag to use automatic number of steps.
  size_t maxIsothermTerms;       ///< Maximum number of isotherm terms.

  Column column;                        ///< Column model advanced by the simulation.
  RungeKutta3 rk3;                      ///< Explicit RK3 integrator instance.
  SemiImplicitRungeKutta3 sirk3;        ///< Semi-implicit RK3 integrator instance.
  CVODE cvode;                          ///< CVODE integrator instance.
  IntegrationScheme integrationScheme;  ///< Selected integration scheme.
  Timing timings;                       ///< Accumulated timing diagnostics.
};
