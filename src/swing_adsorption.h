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
  Column::FlowDirection flowDirection{Column::FlowDirection::Forward};  ///< End at which the feed enters.
  /// Mole fractions of the gas that flows back in at z = L in a forward sub-stage; empty = the gas of the last cell.
  std::vector<double> productEndMoleFractions{};
  bool closedEnd{false};  ///< The end opposite the one where the feed enters is closed.
  Column::PressureRamp::Shape rampShape{Column::PressureRamp::Shape::None};  ///< Ramp from the previous pressure.
  double rampTimeConstant{0.0};                                              ///< Time constant of the ramp in s.
};

/**
 * \brief Times, pressures, and amounts of one finished swing-adsorption sub-stage.
 *
 * The amounts are per component and unit column cross-section in mol/m^2, as Column::gasInventory and
 * Column::adsorbedInventory, at the start of the sub-stage (after a pressure step has rescaled the gas) and at its end.
 * Their change is what entered minus what left during the sub-stage, in the boundary flux eps v c.
 */
struct SwingAdsorptionPhaseRecord
{
  std::string name;                ///< Sub-stage label.
  double startTime{0.0};           ///< In s.
  double endTime{0.0};             ///< In s.
  double pressureAtZero{0.0};      ///< Pressure at z = 0 at the end, in Pa.
  double pressureAtLength{0.0};    ///< Pressure at z = L at the end, in Pa.
  double gasDeviation{0.0};        ///< Column::idealGasDeviation at the end.
  std::vector<double> gasAtStart{};
  std::vector<double> adsorbedAtStart{};
  std::vector<double> gasAtEnd{};
  std::vector<double> adsorbedAtEnd{};
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
  std::vector<SwingAdsorptionPhaseRecord> phaseRecords;  ///< One per finished sub-stage; also written to phases.data.
  Timing timings;                                  ///< Accumulated timing diagnostics.
};
