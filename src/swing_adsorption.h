#pragma once

#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "breakthrough.h"
#include "column.h"
#include "column_multibed.h"
#include "inputreader.h"
#include "timing.h"

/**
 * \brief Operating conditions for one swing-adsorption sub-stage.
 */
struct SwingAdsorptionSubStage
{
  std::string name;         ///< Sub-stage label.
  double temperature{0.0};  ///< Sub-stage temperature in K.
  double pressure{0.0};     ///< Sub-stage inlet pressure in Pa.
  size_t numberOfSteps{0};  ///< Number of time steps in this sub-stage.
};

/**
 * \brief Coordinates swing-adsorption stages for a concrete column model.
 */
template <typename ColumnType>
struct SwingAdsorption
{
  static_assert(std::is_same_v<ColumnType, Column> || std::is_same_v<ColumnType, MultibedColumn>,
                "SwingAdsorption supports Column and MultibedColumn");

  explicit SwingAdsorption(const InputReader& inputReader);

  void run();
  void print() const;
  [[nodiscard]] std::string repr() const;

  Breakthrough<ColumnType> breakthrough;           ///< Breakthrough simulation reused across sub-stages.
  std::vector<SwingAdsorptionSubStage> subStages;  ///< Ordered swing-adsorption stages.
  Timing timings;                                  ///< Accumulated timing diagnostics.
};

extern template struct SwingAdsorption<Column>;
extern template struct SwingAdsorption<MultibedColumn>;
