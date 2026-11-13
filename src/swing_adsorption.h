#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "breakthrough.h"
#include "column.h"
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
  double duration{0.0};     ///< Sub-stage duration in s.
  std::vector<double> moleFractions{};  ///< Sub-stage feed mole fractions per component; empty = Components feed.
  double velocity{0.0};                 ///< Sub-stage inlet (interstitial) velocity in m/s.
};

/**
 * \brief Runs a breakthrough column through a sequence of swing-adsorption stages.
 */
struct SwingAdsorption
{
  explicit SwingAdsorption(const InputReader& inputReader);

  void run();
  void print() const;
  [[nodiscard]] std::string repr() const;

  Breakthrough breakthrough;                       ///< Breakthrough simulation reused across sub-stages.
  std::vector<SwingAdsorptionSubStage> subStages;  ///< Ordered swing-adsorption stages.
  Timing timings;                                  ///< Accumulated timing diagnostics.
};
