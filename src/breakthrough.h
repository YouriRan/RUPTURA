#pragma once

#include <cstddef>
#include <string>
#include <type_traits>

#include "column.h"
#include "column_multibed.h"
#include "cvode.h"
#include "inputreader.h"
#include "rk3.h"
#include "rk3_si.h"
#include "timing.h"

/**
 * \brief Time-integration scheme used for breakthrough propagation.
 */
enum class BreakthroughIntegrationScheme
{
  SSP_RK = 0,  ///< Strong Stability Preserving Runge-Kutta method.
  CVODE = 1,   ///< CVODE integration.
  SIRK3 = 2,   ///< Semi-implicit third-order Runge-Kutta method.
};

/**
 * \brief Simulates a breakthrough process for a concrete column model.
 *
 * Column and MultibedColumn share this implementation at compile time. The
 * input boundary selects the concrete instantiation.
 */
template <typename ColumnType>
struct Breakthrough
{
  static_assert(std::is_same_v<ColumnType, Column> || std::is_same_v<ColumnType, MultibedColumn>,
                "Breakthrough supports Column and MultibedColumn");

  explicit Breakthrough(const InputReader& inputReader);

  void print() const;
  [[nodiscard]] std::string repr() const;
  void run();

  /**
   * \brief Integrates the column from the current time to endTime in s with the selected integrator.
   */
  void advanceTo(double endTime, Timing& accumulatedTimings);

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

  ColumnType column;                                ///< Concrete column model advanced by the simulation.
  RungeKutta3 rk3;                                  ///< Explicit RK3 integrator instance.
  SemiImplicitRungeKutta3 sirk3;                    ///< Semi-implicit RK3 integrator instance.
  CVODE cvode;                                      ///< CVODE integrator instance.
  BreakthroughIntegrationScheme integrationScheme;  ///< Selected integration scheme.
  Timing timings;                                   ///< Accumulated timing diagnostics.
};

extern template struct Breakthrough<Column>;
extern template struct Breakthrough<MultibedColumn>;
