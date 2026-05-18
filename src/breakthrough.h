#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

#include "column.h"
#include "component.h"
#include "cvode.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "rk3_si.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

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
  enum class IntegrationScheme
  {
    SSP_RK = 0,     ///< Strong Stability Preserving Runge-Kutta method.
    CVODE = 1,      ///< CVODE integration
    Iterative = 2,  ///< Iterative integration scheme.
    SIRK3 = 3,
  };

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
  size_t Ncomp;                   ///< Number of components.
  size_t Ngrid;                   ///< Number of grid points.

  size_t printEvery;  ///< Frequency of printing time steps to the screen.
  size_t writeEvery;  ///< Frequency of writing data to files.

  double dt;                     ///< Time step for integration.
  size_t numberOfInitTimeSteps;  ///< Ramping up the timestep.
  size_t Nsteps;                 ///< Total number of steps.
  bool autoSteps;                ///< Flag to use automatic number of steps.
  size_t maxIsothermTerms;       ///< Maximum number of isotherm terms.

  Column column;
  RungeKutta3 rk3;
  SemiImplicitRungeKutta3 sirk3;
  CVODE cvode;
  IntegrationScheme integrationScheme;
};
