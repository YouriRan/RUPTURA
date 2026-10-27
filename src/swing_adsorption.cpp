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
    subStages.emplace_back(SwingAdsorptionSubStage{phase.name, temperature, pressure, phase.duration,
                                                   phase.gasPhaseMolFractions, velocity});
  }
}

void SwingAdsorption::run()
{
  Column& column = breakthrough.column;

  std::vector<std::ofstream> streams;
  for (size_t comp = 0; comp < breakthrough.numberOfComponents; ++comp)
  {
    const std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    streams.emplace_back(std::ofstream{fileName, std::ios_base::app});
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  auto applyStageConditions = [&](const SwingAdsorptionSubStage& stage)
  {
    column.setTemperature(stage.temperature);
    column.inletPressure = stage.pressure;
    // Boundary conditions that derive the velocity from the pressures recompute it in precompute.
    column.columnEntranceVelocity = stage.velocity;
    if (!column.gasTemperature.empty())
    {
      column.gasTemperature[0] = stage.temperature;
    }
    if (column.fluidPhase == decltype(column.fluidPhase)::Liquid)
    {
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        column.concentration[comp] = column.components[comp].inletLiquidConcentration;
      }
    }
    else if (column.inletPressure > 0.0)
    {
      const double inletTemperature =
          std::max(1e-10, column.gasTemperature.empty() ? stage.temperature : column.gasTemperature[0]);
      const double inletConcentration = column.inletPressure / (R * inletTemperature);
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        const double feedFraction = stage.moleFractions.empty() ? column.components[comp].initialGasMoleFraction
                                                                 : stage.moleFractions[comp];
        column.concentration[comp] = feedFraction * inletConcentration;
      }
    }
    // Derived fields (velocity, pressure, equilibrium), then the rates and the sink they feed back into.
    precompute(column, timings);
    computeDerivatives(column);
    precompute(column, timings);
  };

  const double writeInterval = breakthrough.writeInterval;
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
      applyStageConditions(stage);
      if (breakthrough.integrationScheme == BreakthroughIntegrationScheme::CVODE)
      {
        // Restart the multistep history and re-derive the dz/v step cap from this phase's velocity.
        breakthrough.cvode.reinitialize(column);
      }
      std::print("\nSwing phase {}: {} (T={} K, inlet pressure={} Pa, inlet velocity={} m/s, duration={} s)\n",
                 stageIndex + 1, stage.name, stage.temperature, stage.pressure, column.columnEntranceVelocity,
                 stage.duration);
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
    }
    writeOutput();
  }
  catch (const std::runtime_error&)
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
    result += std::format("{}: {}  T={} [K], inlet pressure={} [Pa], inlet velocity={} [m/s], duration={} [s]\n",
                          i + 1, stage.name, stage.temperature, stage.pressure, stage.velocity, stage.duration);
    if (!stage.moleFractions.empty())
    {
      result += "    feed mole fractions:";
      for (size_t comp = 0; comp < stage.moleFractions.size(); ++comp)
      {
        result += std::format(" {}={}", breakthrough.column.components[comp].name, stage.moleFractions[comp]);
      }
      result += "\n";
    }
  }
  result += "\n\n";
  return result;
}
