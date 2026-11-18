#pragma once

#include <cstddef>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "column.h"
#include "cvode.h"
#include "inputreader.h"
#include "rk3.h"
#include "ros34pw2.h"
#include "timing.h"

/**
 * \brief Time-integration scheme used for breakthrough propagation.
 */
enum class BreakthroughIntegrationScheme
{
  SSP_RK = 0,    ///< Strong Stability Preserving Runge-Kutta method.
  CVODE = 1,     ///< CVODE integration.
  ROS34PW2 = 2,  ///< Adaptive four-stage Rosenbrock-W method with the analytical sorption Jacobian.
};

/**
 * \brief Simulates a breakthrough process: feeds the column and writes its state at every output time.
 */
struct Breakthrough
{
  explicit Breakthrough(const InputReader& inputReader);

  void print() const;
  [[nodiscard]] std::string repr() const;
  void run();

  /**
   * \brief Integrates the column from the current time to endTime in s with the selected integrator.
   */
  void advanceTo(double endTime, Timing& accumulatedTimings);

  /**
   * \brief Writes column.json and one snapshot to column.data and the component files, which the first call starts,
   *        and prints the progress lines. Throws when a file cannot be written.
   */
  void writeOutput();

  /**
   * \brief Integrates to the next output time and applies the automatic stop criteria.
   *
   * \return True once the simulation time is reached.
   */
  bool advanceToNextOutput();

  const std::string displayName;  ///< Name of the simulation for display purposes.
  size_t carrierGasComponent{0};  ///< Index of the carrier gas component.
  size_t numberOfComponents;      ///< Number of components.
  size_t numberOfGridPoints;      ///< Number of grid intervals; node count is numberOfGridPoints + 1.

  double writeInterval;      ///< Interval between outputs (and progress lines) in s.
  double simulationTime;     ///< Simulated time in s; set once the stop criterion is met when autoSimulationTime.
  bool autoSimulationTime;   ///< Run until the breakthrough or reaction steady-state criterion is met when true.
  double time{0.0};          ///< Current simulated time in s.
  size_t numberOfOutputs{0};  ///< Output intervals completed; the next output is at (numberOfOutputs + 1) writeInterval.
  size_t maxIsothermTerms;   ///< Maximum number of isotherm terms.
  /// run() continues the output files and leaves out the snapshot at its start, which they hold already: after
  /// ReadColumnFile without ResetClock, or once the files are started.
  bool appendOutput;
  bool warnedIdealGasDeviation{false};          ///< The warning of writeOutput is printed once.
  std::vector<std::ofstream> componentStreams;  ///< Output file of each component, opened by the first writeOutput.
  std::ofstream columnStream;                   ///< column.data, opened by the first writeOutput.

  Column column;                                    ///< Column model advanced by the simulation.
  RungeKutta3 rk3;                                  ///< Explicit RK3 integrator instance.
  ROS34PW2 ros34pw2;                         ///< ROS34PW2 integrator instance.
  CVODE cvode;                                      ///< CVODE integrator instance.
  BreakthroughIntegrationScheme integrationScheme;  ///< Selected integration scheme.
  Timing timings;                                   ///< Accumulated timing diagnostics.
  std::function<void()> afterOutput;  ///< Called by run() after every output; Python checks for Ctrl-C here.
};
