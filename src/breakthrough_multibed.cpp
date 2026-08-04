#include "breakthrough_multibed.h"

#include <algorithm>
#include <cstdio>
#include <format>
#include <fstream>
#include <print>
#include <stdexcept>
#include <vector>

BreakthroughMultibed::BreakthroughMultibed(const InputReader& inputReader)
    : displayName(inputReader.displayName),
      numberOfComponents(inputReader.components.size()),
      numberOfGridPoints(inputReader.numberOfGridPoints),
      printEvery(inputReader.printEvery),
      writeEvery(inputReader.writeEvery),
      timeStep(inputReader.timeStep),
      numberOfInitTimeSteps(inputReader.numberOfInitTimeSteps),
      numberOfSteps(inputReader.numberOfTimeSteps),
      autoNumberOfSteps(inputReader.autoNumberOfTimeSteps),
      maxIsothermTerms(inputReader.maxIsothermTerms),
      column(inputReader),
      rk3(inputReader)
{
  if (inputReader.adsorbentComponents.size() < 2)
  {
    throw std::runtime_error("Error: multibed breakthrough requires at least two adsorbents");
  }
  if (inputReader.breakthroughIntegrator != 0)
  {
    throw std::runtime_error("Error: multibed breakthrough currently supports the RungeKutta3 integrator");
  }
  if (!inputReader.reactions.empty())
  {
    throw std::runtime_error("Error: reactions are not yet supported in multibed breakthrough simulations");
  }
  for (const std::vector<Component>& components : inputReader.adsorbentComponents)
  {
    if (std::ranges::any_of(components,
                            [](const Component& component) { return component.chemisorption.numberOfSites > 0; }))
    {
      throw std::runtime_error("Error: chemisorption is not yet supported in multibed breakthrough simulations");
    }
  }

  column.initialize();
  if (inputReader.readColumnFile.has_value())
  {
    column.readJSON(*inputReader.readColumnFile);
  }

  const auto openMode = inputReader.readColumnFile.has_value() ? std::ios_base::app : std::ios_base::trunc;
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    componentStreams[comp] = std::ofstream{fileName, openMode};
  }
  std::ofstream columnStream("column.data", openMode);
  if (!inputReader.readColumnFile.has_value())
  {
    column.writeOutputHeader(componentStreams, columnStream);
  }
}

void BreakthroughMultibed::run()
{
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    componentStreams[comp] = std::ofstream{fileName, std::ios_base::app};
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  bool finished = false;
  size_t step = 0;
  double realTime = 0.0;
  constexpr double initialTimeStepFraction = 1.0e-4;

  if (numberOfInitTimeSteps > 0)
  {
    rk3.autoNumberOfSteps = false;
    rk3.timeStep = timeStep * initialTimeStepFraction;
  }

  try
  {
    while (!finished)
    {
      const double stepTime = realTime;
      if (step % writeEvery == 0)
      {
        column.writeJSON("column.json");
      }

      finished = rk3.propagate(column, step, timings);
      realTime += rk3.timeStep;

      if (step < numberOfInitTimeSteps)
      {
        const double progress = static_cast<double>(step) / static_cast<double>(numberOfInitTimeSteps);
        rk3.timeStep = timeStep * (initialTimeStepFraction +
                                   (1.0 - initialTimeStepFraction) *
                                       (3.0 * progress * progress - 2.0 * progress * progress * progress));
      }
      else if (step == numberOfInitTimeSteps)
      {
        rk3.autoNumberOfSteps = autoNumberOfSteps;
        rk3.timeStep = timeStep;
      }

      if (step % writeEvery == 0)
      {
        column.writeOutput(componentStreams, columnStream, stepTime);
      }
      if (step % printEvery == 0)
      {
        const double averageMixtureSteps = column.iastPerformance.second == 0
                                               ? 0.0
                                               : static_cast<double>(column.iastPerformance.first) /
                                                     static_cast<double>(column.iastPerformance.second);
        std::print("Timestep {}, time: {:6.5f} [s]\n", step, stepTime);
        std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageMixtureSteps);
        std::fflush(stdout);
      }

      ++step;
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

void BreakthroughMultibed::print() const { std::print("{}", repr()); }

std::string BreakthroughMultibed::repr() const
{
  std::string result = std::format(
      "Column properties\n"
      "=======================================================\n"
      "Display-name:                          {}\n",
      displayName);
  result += column.repr();
  result += std::format(
      "Breakthrough settings\n"
      "=======================================================\n"
      "Number of adsorbents:          {}\n"
      "Number of time steps:          {}\n"
      "Print every step:              {}\n"
      "Write data every step:         {}\n\n\n"
      "Integration details\n"
      "=======================================================\n"
      "Time step:                     {} [s]\n"
      "Number of column grid points:  {}\n\n\n"
      "Adsorbent component data\n"
      "=======================================================\n"
      "maximum isotherm terms:        {}\n",
      column.numberOfAdsorbents, numberOfSteps, printEvery, writeEvery, timeStep, numberOfGridPoints, maxIsothermTerms);

  for (size_t ads = 0; ads < column.numberOfAdsorbents; ++ads)
  {
    result += std::format("Adsorbent {}\n", ads);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      result += std::format("{}\n", column.physisorptionMixtures[ads].components[comp].repr());
    }
  }
  return result;
}
