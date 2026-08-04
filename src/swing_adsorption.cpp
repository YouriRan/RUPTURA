#include "swing_adsorption.h"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <format>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "breakthrough.h"
#include "column.h"
#include "rk3.h"
#include "sorption.h"
#include "utils.h"

SwingAdsorption::SwingAdsorption(const InputReader& inputReader) : breakthrough(inputReader)
{
  if (inputReader.swingAdsorptionPhases.empty())
  {
    throw std::invalid_argument("SwingAdsorption: at least one swing adsorption phase must be provided.");
  }

  subStages.reserve(inputReader.swingAdsorptionPhases.size());
  for (size_t i = 0; i < inputReader.swingAdsorptionPhases.size(); ++i)
  {
    const InputReader::SwingAdsorptionPhase& phase = inputReader.swingAdsorptionPhases[i];
    const double temperature = phase.temperature.value_or(inputReader.temperature);
    const double pressure = phase.inletPressure.value_or(inputReader.inletPressure);
    if (phase.numberOfSteps == 0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has zero time steps.");
    }
    if (temperature < 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid temperature.");
    }
    if (pressure <= 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid inlet pressure.");
    }
    subStages.emplace_back(SwingAdsorption::SubStage{phase.name, temperature, pressure, phase.numberOfSteps});
  }
}

void SwingAdsorption::run()
{
  Column& column = breakthrough.column;

  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < breakthrough.numberOfComponents; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + column.components[i].name + ".data";
    streams.emplace_back(std::ofstream{fileName, std::ios_base::app});
  }
  std::ofstream movieStream("column.data", std::ios_base::app);

  size_t step = 0;
  double realTime = 0.0;
  double xi = 1e-4;

  const size_t totalNumberOfSteps = std::accumulate(
      subStages.begin(), subStages.end(), size_t{0},
      [](size_t total, const SubStage& stage) { return total + stage.numberOfSteps; });

  breakthrough.rk3.autoNumberOfSteps = false;
  breakthrough.sirk3.autoNumberOfSteps = false;
  breakthrough.cvode.autoNumberOfSteps = false;
  breakthrough.rk3.numberOfSteps = totalNumberOfSteps;
  breakthrough.sirk3.numberOfSteps = totalNumberOfSteps;
  breakthrough.cvode.numberOfSteps = totalNumberOfSteps;

  if (breakthrough.numberOfInitTimeSteps > 0)
  {
    breakthrough.rk3.timeStep = breakthrough.timeStep * xi;
    breakthrough.sirk3.timeStep = breakthrough.timeStep * xi;
    breakthrough.cvode.timeStep = breakthrough.timeStep * xi;
  }

  auto applyStageConditions = [&](const SubStage& stage)
  {
    column.setTemperature(stage.temperature);
    column.inletPressure = stage.pressure;
    if (!column.gasTemperature.empty())
    {
      column.gasTemperature[0] = stage.temperature;
    }
    if (column.inletPressure > 0.0)
    {
      const double inletTemperature = std::max(1e-10, column.gasTemperature.empty() ? stage.temperature
                                                                                    : column.gasTemperature[0]);
      const double inletConcentration = column.inletPressure / (R * inletTemperature);
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        column.concentration[comp] = column.components[comp].initialGasMoleFraction * inletConcentration;
      }
    }
    RK3Helpers::updateVelocityAndPressure(column);
    RK3Helpers::computeEquilibriumLoadings(column);
    RK3Helpers::computeSorptionDerivatives(column);
    RK3Helpers::updateVelocityAndPressure(column);
    RK3Helpers::computeEquilibriumLoadings(column);
  };

  try
  {
    for (size_t stageIndex = 0; stageIndex < subStages.size(); ++stageIndex)
    {
      const SubStage& stage = subStages[stageIndex];
      applyStageConditions(stage);
      std::print("\nSwing phase {}: {} (T={} K, inlet pressure={} Pa, steps={})\n",
                 stageIndex + 1, stage.name, stage.temperature, stage.pressure, stage.numberOfSteps);

      for (size_t subStep = 0; subStep < stage.numberOfSteps; subStep++)
      {
        const double stepTime = realTime;

        if (step % breakthrough.writeEvery == 0)
        {
          const std::string outputFile = std::format("column.json", step);
          column.writeJSON(outputFile);
        }

        switch (breakthrough.integrationScheme)
        {
          case Breakthrough::IntegrationScheme::SSP_RK:
          {
            static_cast<void>(breakthrough.rk3.propagate(column, step, timings));
            realTime += breakthrough.rk3.timeStep;
            break;
          }
          case Breakthrough::IntegrationScheme::CVODE:
          {
            static_cast<void>(breakthrough.cvode.propagate(column, step, timings));
            realTime += breakthrough.cvode.timeStep;
            break;
          }
          case Breakthrough::IntegrationScheme::SIRK3:
          {
            static_cast<void>(breakthrough.sirk3.propagate(column, step, timings));
            realTime += breakthrough.sirk3.timeStep;
            break;
          }
          default:
            break;
        }

        if (step < breakthrough.numberOfInitTimeSteps)
        {
          double i = static_cast<double>(step) / static_cast<double>(breakthrough.numberOfInitTimeSteps);
          double nextTime = breakthrough.timeStep * (xi + (1 - xi) * (3 * i * i - 2 * i * i * i));
          breakthrough.rk3.timeStep = nextTime;
          breakthrough.sirk3.timeStep = nextTime;
          breakthrough.cvode.timeStep = nextTime;
        }
        else if (step == breakthrough.numberOfInitTimeSteps)
        {
          breakthrough.rk3.timeStep = breakthrough.timeStep;
          breakthrough.sirk3.timeStep = breakthrough.timeStep;
          breakthrough.cvode.timeStep = breakthrough.timeStep;
        }

        if (step % breakthrough.writeEvery == 0)
        {
          column.writeOutput(streams, movieStream, stepTime);
        }
        if (step % breakthrough.printEvery == 0)
        {
          const double averageIastSteps =
              column.iastPerformance.second == 0
                  ? 0.0
                  : static_cast<double>(column.iastPerformance.first) /
                        static_cast<double>(column.iastPerformance.second);
          std::print("Timestep {}, time: {:6.5f} [s]\n", step, stepTime);
          std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageIastSteps);
          std::fflush(stdout);
        }
        step++;
      }
    }
  }
  catch (const std::runtime_error&)
  {
    column.writeJSON("failed_state.json");
    throw;
  }

  std::print("Final timestep {}, time: {:6.5f} [s]\n\n", step, realTime);
  timings.print();
}

void SwingAdsorption::print() const { std::print("{}", repr()); }
std::string SwingAdsorption::repr() const
{
  std::string s = breakthrough.repr();
  s += std::format(
      "Swing adsorption phases\n"
      "=======================================================\n");
  for (size_t i = 0; i < subStages.size(); ++i)
  {
    const SubStage& stage = subStages[i];
    s += std::format("{}: {}  T={} [K], inlet pressure={} [Pa], steps={}\n",
                     i + 1, stage.name, stage.temperature, stage.pressure, stage.numberOfSteps);
  }
  s += "\n\n";
  return s;
}
