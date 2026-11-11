#include "swing_adsorption.h"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <format>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "rk3.h"
#include "utils.h"

SwingAdsorption::SwingAdsorption(const InputReader& inputReader) : breakthrough(inputReader)
{
  if (inputReader.swingAdsorptionPhases.empty())
  {
    throw std::invalid_argument("SwingAdsorption: at least one swing adsorption phase must be provided.");
  }

  subStages.reserve(inputReader.swingAdsorptionPhases.size());
  for (const InputReader::SwingAdsorptionPhase& phase : inputReader.swingAdsorptionPhases)
  {
    const double temperature = phase.temperature.value_or(inputReader.temperature);
    const double pressure = phase.inletPressure.value_or(inputReader.inletPressure);
    if (!(phase.duration > 0.0))
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has a non-positive duration.");
    }
    if (temperature < 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid temperature.");
    }
    if (pressure <= 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid inlet pressure.");
    }
    const double velocity = phase.columnEntranceVelocity.value_or(inputReader.columnEntranceVelocity);
    subStages.emplace_back(SwingAdsorptionSubStage{
        phase.name, temperature, pressure, phase.duration, phase.gasPhaseMolFractions, velocity,
        Column::FlowDirection(phase.flowDirection.value_or(inputReader.flowDirection)),
        phase.productEndGasPhaseMolFractions, phase.closedEnd, Column::PressureRamp::Shape(phase.pressureRampShape),
        phase.pressureRampTimeConstant});
  }
}

void SwingAdsorption::run()
{
  Column& column = breakthrough.column;
  const size_t n = column.numberOfComponents;
  const size_t lastNode = column.numberOfGridPoints;

  std::vector<std::ofstream> streams;
  for (size_t comp = 0; comp < breakthrough.numberOfComponents; ++comp)
  {
    const std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    streams.emplace_back(std::ofstream{fileName, std::ios_base::app});
  }
  std::ofstream columnStream("column.data", std::ios_base::app);
  std::ofstream phaseStream("phases.data");
  std::print(phaseStream,
             "# column 1: Phase start time [s]\n# column 2: Phase end time [s]\n"
             "# column 3: Pressure at z = 0 at the end of the phase [Pa]\n"
             "# column 4: Pressure at z = L at the end of the phase [Pa]\n"
             "# column 5: Max |R T sum(c) / P - 1| at the end of the phase [-]\n");
  for (size_t comp = 0; comp < n; ++comp)
  {
    std::print(phaseStream,
               "# columns {} to {}: {}: gas and adsorbed at the phase start, gas and adsorbed at its end [mol/m^2]\n",
               6 + 4 * comp, 9 + 4 * comp, column.components[comp].name);
  }
  phaseRecords.clear();

  auto applyStageConditions = [&](const SwingAdsorptionSubStage& stage)
  {
    const bool reverse = stage.flowDirection == Column::FlowDirection::Reverse;
    using Boundary = Column::BoundaryCondition;
    // The gas takes the pressure profile of the new phase at its composition, where a continuity closure would keep the
    // old amount: after a pressure step without a ramp, a temperature step without the energy balance, and when the end
    // where the pressure is set changes, which moves the profile by the pressure drop. Left alone, that shift would add up
    // over the cycles. FixedVelocity and a fixed pressure profile let the feed displace the gas, and
    // InletVelocityOutletPressure does not follow the phase pressure.
    const bool ramped = stage.rampShape != Column::PressureRamp::Shape::None;
    const bool rescaleGas =
        column.fluidPhase == Column::FluidPhase::Gas && column.boundaryCondition != Boundary::FixedVelocity &&
        column.boundaryCondition != Boundary::FixedPressureInletVelocity &&
        column.boundaryCondition != Boundary::InletVelocityOutletPressure &&
        ((!ramped && stage.pressure != column.inletPressure) || stage.flowDirection != column.flowDirection ||
         (!column.energyBalance && stage.temperature != column.gasTemperature[1]));
    column.pressureRamp = {stage.rampShape, stage.rampTimeConstant, breakthrough.time,
                           column.totalPressure[reverse ? lastNode : 0], stage.pressure, breakthrough.time};
    column.pressureRate = 0.0;
    column.relativePressureRate = 0.0;
    column.flowDirection = stage.flowDirection;
    column.closedEnd = stage.closedEnd;
    column.setTemperature(stage.temperature);
    column.inletPressure = stage.pressure;
    column.advancePressureRamp(breakthrough.time);
    // Boundary conditions that derive the velocity from the pressures recompute it in precompute.
    column.columnEntranceVelocity = stage.velocity;

    // The feed enters from node 0 in a forward phase and from the product-end gas in a reverse phase. There node 0,
    // which gas that flows back in at z = 0 does not use, takes the state of the first cell, so that the output at
    // z = 0 starts from the gas that leaves there.
    std::vector<double> feed(column.concentration.begin(), column.concentration.begin() + static_cast<std::ptrdiff_t>(n));
    if (column.fluidPhase == decltype(column.fluidPhase)::Liquid)
    {
      for (size_t comp = 0; comp < n; ++comp) feed[comp] = column.components[comp].inletLiquidConcentration;
    }
    else if (column.inletPressure > 0.0)
    {
      const double inletConcentration = column.inletPressure / (R * std::max(1e-10, stage.temperature));
      for (size_t comp = 0; comp < n; ++comp)
      {
        const double feedFraction = stage.moleFractions.empty() ? column.components[comp].initialGasMoleFraction
                                                                 : stage.moleFractions[comp];
        feed[comp] = feedFraction * inletConcentration;
      }
    }
    if (reverse)
    {
      column.productEndConcentration = feed;
      column.productEndTemperature = stage.temperature;
      std::copy_n(column.concentration.begin() + static_cast<std::ptrdiff_t>(n), n, column.concentration.begin());
      if (!column.gasTemperature.empty()) column.gasTemperature[0] = column.gasTemperature[1];
    }
    else
    {
      column.productEndConcentration.clear();
      column.productEndTemperature = 0.0;
      std::copy(feed.begin(), feed.end(), column.concentration.begin());
      if (!column.gasTemperature.empty()) column.gasTemperature[0] = stage.temperature;
    }
    // Derived fields (velocity, pressure, equilibrium), then the rates and the sink they feed back into.
    if (rescaleGas) std::fill(column.stateDot.begin(), column.stateDot.end(), 0.0);  // the rates of the old state
    precompute(column, timings);
    // The second pass takes the pressure drop of the rescaled gas.
    for (size_t pass = 0; rescaleGas && pass < 2; ++pass)
    {
      for (size_t grid = reverse ? 0 : 1; grid <= lastNode; ++grid)
      {
        const double factor = column.totalPressure[grid] /
                              std::max(1e-30, R * column.gasTemperature[grid] * column.totalConcentration[grid]);
        for (size_t comp = 0; comp < n; ++comp) column.concentration[grid * n + comp] *= factor;
      }
      precompute(column, timings);
    }
    if (!reverse && !stage.productEndMoleFractions.empty())
    {
      // The gas that flows back in at z = L, at the pressure there and the phase temperature.
      const double productEndTotal =
          column.totalPressure[column.numberOfGridPoints] / (R * std::max(1e-10, stage.temperature));
      column.productEndConcentration.resize(n);
      for (size_t comp = 0; comp < n; ++comp)
      {
        column.productEndConcentration[comp] = stage.productEndMoleFractions[comp] * productEndTotal;
      }
      column.productEndTemperature = stage.temperature;
      precompute(column, timings);
    }
    computeDerivatives(column);
    precompute(column, timings);
  };

  const double writeInterval = breakthrough.writeInterval;
  bool warnedIdealGasDeviation = false;
  auto writeOutput = [&]
  {
    column.writeJSON("column.json");
    column.writeOutput(streams, columnStream, breakthrough.time);
    const double averageIastSteps = column.iastPerformance.second == 0
                                        ? 0.0
                                        : static_cast<double>(column.iastPerformance.first) /
                                              static_cast<double>(column.iastPerformance.second);
    std::print("Output {}, time: {:6.5f} [s]\n", breakthrough.numberOfOutputs, breakthrough.time);
    std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageIastSteps);
    // As in Breakthrough::run: the continuity closures keep sum_j c_j at P/RT.
    if (const double deviation = column.idealGasDeviation(); deviation > 0.0)
    {
      std::print("    Max |R T sum(c) / P - 1|: {:.3e}\n", deviation);
      if (deviation > 1.0e-2 && !warnedIdealGasDeviation)
      {
        std::print("    Warning: R T sum(c) departs from the closure pressure by more than 1 %; the closures do "
                   "not follow changes of the pressure drop\n");
        warnedIdealGasDeviation = true;
      }
    }
    std::fflush(stdout);
  };

  try
  {
    // An output that falls on the start of a phase is written once that phase's conditions are applied, as is the
    // initial state; the state at the end of the last phase is always written.
    bool outputPending = true;
    double phaseEndTime = 0.0;
    for (size_t stageIndex = 0; stageIndex < subStages.size(); ++stageIndex)
    {
      const SwingAdsorptionSubStage& stage = subStages[stageIndex];
      try
      {
        applyStageConditions(stage);
      }
      catch (const std::runtime_error& error)
      {
        throw std::runtime_error(std::format("{}at the start of swing phase '{}'; a large pressure step can need a "
                                             "PressureRamp or a ClosedEnd\n",
                                             error.what(), stage.name));
      }
      if (breakthrough.integrationScheme == BreakthroughIntegrationScheme::CVODE)
      {
        // Restart the multistep history and re-derive the 4 dz/v step cap from this phase's velocity.
        breakthrough.cvode.reinitialize(column);
      }
      // The switched boundary conditions start a new transient, which ROS34PW2 starts with a small step.
      if (breakthrough.integrationScheme == BreakthroughIntegrationScheme::ROS34PW2) breakthrough.ros34pw2.reinitialize();
      const bool reverse = stage.flowDirection == Column::FlowDirection::Reverse;
      std::print("\nSwing phase {}: {} (T={} K, inlet pressure={} Pa, inlet velocity={} m/s, duration={} s{}{}{})\n",
                 stageIndex + 1, stage.name, stage.temperature, stage.pressure, column.columnEntranceVelocity,
                 stage.duration, reverse && !stage.closedEnd ? ", reverse flow" : "",
                 !stage.closedEnd ? "" : reverse ? ", closed at z = 0" : ", closed at z = L",
                 stage.rampShape == Column::PressureRamp::Shape::None ? "" : ", pressure ramp");
      SwingAdsorptionPhaseRecord record{.name = stage.name,
                                        .startTime = breakthrough.time,
                                        .gasAtStart = column.gasInventory(),
                                        .adsorbedAtStart = column.adsorbedInventory()};
      if (outputPending)
      {
        writeOutput();
        outputPending = false;
      }

      // Integrate to each output time inside the phase and then to the phase end, where the conditions switch.
      phaseEndTime += stage.duration;
      const double tolerance = 1.0e-9 * writeInterval;
      while (breakthrough.time < phaseEndTime - tolerance)
      {
        // Output times are computed from their index, so that they do not drift over many phases.
        const double nextOutputTime = static_cast<double>(breakthrough.numberOfOutputs + 1) * writeInterval;
        const bool outputBeforePhaseEnd = nextOutputTime < phaseEndTime - tolerance;
        const bool outputAtPhaseEnd = !outputBeforePhaseEnd && nextOutputTime <= phaseEndTime + tolerance;
        breakthrough.advanceTo(outputBeforePhaseEnd ? nextOutputTime : phaseEndTime, timings);
        if (outputBeforePhaseEnd || outputAtPhaseEnd) ++breakthrough.numberOfOutputs;
        if (outputBeforePhaseEnd) writeOutput();
        outputPending = outputAtPhaseEnd;
      }

      record.endTime = breakthrough.time;
      record.pressureAtZero = column.totalPressure[0];
      record.pressureAtLength = column.totalPressure[lastNode];
      record.gasDeviation = column.idealGasDeviation();
      record.gasAtEnd = column.gasInventory();
      record.adsorbedAtEnd = column.adsorbedInventory();
      std::print(phaseStream, "{} {} {} {} {}", record.startTime, record.endTime, record.pressureAtZero,
                 record.pressureAtLength, record.gasDeviation);
      for (size_t comp = 0; comp < n; ++comp)
      {
        std::print(phaseStream, " {} {} {} {}", record.gasAtStart[comp], record.adsorbedAtStart[comp],
                   record.gasAtEnd[comp], record.adsorbedAtEnd[comp]);
      }
      std::print(phaseStream, "  # {}\n", stage.name);
      phaseStream.flush();
      phaseRecords.push_back(std::move(record));
    }
    writeOutput();
  }
  catch (const std::exception&)
  {
    column.writeJSON("failed_state.json");
    throw;
  }

  std::print("Final time: {:6.5f} [s]\n\n", breakthrough.time);
  timings.print();
}

void SwingAdsorption::print() const
{
  std::print("{}", repr());
}

std::string SwingAdsorption::repr() const
{
  std::string result = breakthrough.repr();
  result += std::format(
      "Swing adsorption phases\n"
      "=======================================================\n");
  for (size_t i = 0; i < subStages.size(); ++i)
  {
    const SwingAdsorptionSubStage& stage = subStages[i];
    const bool reverse = stage.flowDirection == Column::FlowDirection::Reverse;
    result += std::format("{}: {}  T={} [K], inlet pressure={} [Pa], inlet velocity={} [m/s], duration={} [s]{}{}{}\n",
                          i + 1, stage.name, stage.temperature, stage.pressure, stage.velocity, stage.duration,
                          reverse ? ", reverse flow (feed at z = L)" : "",
                          !stage.closedEnd ? "" : reverse ? ", closed at z = 0" : ", closed at z = L",
                          stage.rampShape == Column::PressureRamp::Shape::None
                              ? ""
                              : std::format(", pressure ramp over {} [s]", stage.rampTimeConstant));
    for (const auto& [label, fractions] : {std::pair{"feed", &stage.moleFractions},
                                            std::pair{"product-end gas", &stage.productEndMoleFractions}})
    {
      if (fractions->empty()) continue;
      result += std::format("    {} mole fractions:", label);
      for (size_t comp = 0; comp < fractions->size(); ++comp)
      {
        result += std::format(" {}={}", breakthrough.column.components[comp].name, (*fractions)[comp]);
      }
      result += "\n";
    }
  }
  result += "\n\n";
  return result;
}
