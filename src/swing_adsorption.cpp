#include "swing_adsorption.h"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "breakthrough.h"
#include "column.h"

SwingAdsorption::SwingAdsorption(const InputReader& inputReader) : breakthrough(inputReader)
{
  const size_t numTemperatures = inputReader.swingTemperatures.size();
  const size_t numPressures = inputReader.swingPressures.size();
  const size_t numSteps = inputReader.swingSteps.size();

  if (numTemperatures == 0 && numPressures == 0 && numSteps == 0)
  {
    throw std::invalid_argument(
        "SwingAdsorption: at least one of swingTemperatures, swingPressures, or swingSteps must be provided.");
  }

  size_t numSubStages = 0;

  if (numTemperatures != 0)
  {
    numSubStages = numTemperatures;
  }

  if (numPressures != 0)
  {
    if (numSubStages != 0 && numPressures != numSubStages)
    {
      throw std::invalid_argument(
          "SwingAdsorption: swingPressures must have the same size as the other swing input vectors.");
    }
    numSubStages = numPressures;
  }

  if (numSteps != 0)
  {
    if (numSubStages != 0 && numSteps != numSubStages)
    {
      throw std::invalid_argument(
          "SwingAdsorption: swingSteps must have the same size as the other swing input vectors.");
    }
    numSubStages = numSteps;
  }

  subStages.reserve(numSubStages);

  for (size_t i = 0; i < numSubStages; ++i)
  {
    const double temperature = (numTemperatures != 0) ? inputReader.swingTemperatures[i] : inputReader.temperature;
    const double pressure = (numPressures != 0) ? inputReader.swingPressures[i] : inputReader.totalPressure;
    const size_t numberOfSteps = (numSteps != 0) ? inputReader.swingSteps[i] : inputReader.numberOfTimeSteps;
    subStages.emplace_back(SwingAdsorption::SubStage{temperature, pressure, numberOfSteps});
  }
}

void SwingAdsorption::run()

{
  Column& column = breakthrough.column;

  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < breakthrough.Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + column.components[i].name + ".data";
    streams.emplace_back(std::ofstream{fileName, std::ios_base::app});
  }
  std::ofstream movieStream("column.data", std::ios_base::app);

  bool finished = false;
  size_t step = 0;
  double realTime = 0.0;

  for (auto& stage : subStages)
  {
    column.externalTemperature = stage.temperature;
    column.externalPressure = stage.pressure;

    for (size_t subStep = 0; subStep < stage.numberOfSteps; subStep++)
    {
      if (step % breakthrough.writeEvery == 0)
      {
        const std::string outputFile = std::format("column.json", step);
        column.writeJSON(outputFile);
      }

      switch (breakthrough.integrationScheme)
      {
        case Breakthrough::IntegrationScheme::SSP_RK:
        {
          finished = breakthrough.rk3.propagate(column, step);
          realTime += breakthrough.rk3.timeStep;
          break;
        }
        case Breakthrough::IntegrationScheme::CVODE:
        {
          finished = breakthrough.cvode.propagate(column, step);
          realTime += breakthrough.cvode.timeStep;
          break;
        }
        case Breakthrough::IntegrationScheme::SIRK3:
        {
          finished = breakthrough.sirk3.propagate(column, step);
          realTime += breakthrough.sirk3.timeStep;
          break;
        }
        default:
          break;
      }

      if (step % breakthrough.writeEvery == 0)
      {
        column.writeOutput(streams, movieStream, realTime);
      }
      if (step % breakthrough.printEvery == 0)
      {
        std::print("Timestep {}, time: {:6.5f} [s]\n", step, realTime);
        std::print(
            "    Average number of mixture-prediction steps: {:6.5f}\n",
            static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second));
        std::cout << std::flush;
      }
      step++;
    }
  }
}

void SwingAdsorption::print() const { std::print("{}", repr()); }
std::string SwingAdsorption::repr() const { return breakthrough.repr(); }
