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

  ColumnType column;                                ///< Concrete column model advanced by the simulation.
  RungeKutta3 rk3;                                  ///< Explicit RK3 integrator instance.
  SemiImplicitRungeKutta3 sirk3;                    ///< Semi-implicit RK3 integrator instance.
  CVODE cvode;                                      ///< CVODE integrator instance.
  BreakthroughIntegrationScheme integrationScheme;  ///< Selected integration scheme.
  Timing timings;                                   ///< Accumulated timing diagnostics.
};

extern template struct Breakthrough<Column>;
extern template struct Breakthrough<MultibedColumn>;
