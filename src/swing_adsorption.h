#pragma once

#include <array>
#include <cstddef>
#include <deque>
#include <optional>
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
  /// The feed, or the product-end gas, is the mean gas that left the column in this phase of the cycle the last time
  /// it ran; until then the mole fractions above apply.
  std::optional<size_t> feedFromPhase{};
  std::optional<size_t> productEndFromPhase{};
  /// The sub-stage ends before its duration once these hold (all, or any with stopOnAny) after minimumDuration.
  std::vector<InputReader::StopCriterion> stopCriteria{};
  bool stopOnAny{false};
  double minimumDuration{0.0};               ///< In s.
  bool errorOnTimeout{false};                ///< Reaching the duration with the stop criteria unmet is an error.
  std::array<InputReader::Port, 2> ports{};  ///< MomentumBalance: the ports at z = 0 and z = L.
  double jacketTemperature{0.0};  ///< Temperature beyond the outer wall in K: the sub-stage one without ThermalControl.
  double jacketHeatTransfer{0.0};                      ///< Heat-transfer coefficient of the outer wall in W/(m^2 K).
  double coolantMassFlow{0.0};                         ///< Mass flow of the Coolant in kg/s.
  std::vector<InputReader::Connection> connections{};  ///< Of the phase; on the sub-stage of the first column.
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
  /// Amount of each component that entered and that left at z = 0, then at z = L, in the sub-stage, in mol/m^2: four
  /// blocks of one value per component, from the net amounts of Column::streamAmount, so that a component that
  /// entered and left at one end shows their difference. Empty when the column does not integrate them.
  std::vector<double> streamAmounts{};
  /// Heat that entered through the outer wall (negative: cooling), and work that brings the gas that left through
  /// pressure ports below the exhaust pressure to it, in the sub-stage, in J/m^2; with the stream amounts.
  double wallHeat{0.0};
  double exhaustWork{0.0};
};

/**
 * \brief Runs a breakthrough column through a sequence of swing-adsorption stages.
 *
 * A stage ends at its duration, or earlier at the time at which its stop criteria turn true. With a cyclic-steady-state
 * tolerance the run ends once, for steadyStateCycles cycles in a row, the end of a cycle differs from the end of the
 * cycle before by less than the tolerance in three relative measures, which cycles.data lists per cycle: every state
 * value, the amount of every component in the column, and every stream amount of every stage. Each change is taken
 * relative to the size of the value plus the largest of its kind (its state field, the amounts, the stream amounts of
 * the cycle), so that values near zero, such as a carrier gas that has left the column, are measured against that
 * largest one.
 *
 * With columnNames the run holds that many identical columns, and subStages and phaseRecords one entry per phase and
 * column, in the order of the names. The connections of a phase carry gas between the ends of the columns and the
 * buffers. Each column keeps its own integrator, and the connections are advanced in coupling steps: the flows follow
 * from the pressures at the start of a step, the columns then take the step one after the other with their connected
 * ends at these flows, a column that gives gas before those that take it, so that what enters a column is what the
 * stream amounts of the others say left them. A buffer is a mixed vessel at a constant temperature. The stop criteria
 * are then looked at after each coupling step.
 */
struct SwingAdsorption
{
  explicit SwingAdsorption(const InputReader& inputReader);

  void run();
  void print() const;
  [[nodiscard]] std::string repr() const;

  Breakthrough breakthrough;              ///< Breakthrough simulation reused across sub-stages; the first column.
  std::deque<Breakthrough> otherColumns;  ///< Those of the columns after the first.
  std::vector<std::string>
      columnNames;  ///< Names of the columns; empty: one column, which writes to the run directory.
  std::vector<InputReader::Buffer> buffers;        ///< Gas vessels of the connections.
  std::vector<std::vector<double>> bufferAmounts;  ///< Amount of each component in each buffer in mol.
  std::vector<SwingAdsorptionSubStage> subStages;  ///< Ordered swing-adsorption stages.
  size_t phasesPerCycle;                           ///< Phases of one cycle: subStages repeats them.
  double steadyStateTolerance;                     ///< Cyclic-steady-state tolerance; 0 runs every sub-stage.
  size_t steadyStateCycles;                        ///< Cycles in a row within the tolerance that end the run.
  std::vector<SwingAdsorptionPhaseRecord> phaseRecords;  ///< One per finished sub-stage; also written to phases.data.
  Timing timings;                                  ///< Accumulated timing diagnostics.
};
