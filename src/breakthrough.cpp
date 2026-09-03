#include "breakthrough.h"

#include <cstdio>
#include <format>
#include <fstream>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace
{
void validateMultibedInput(const InputReader& inputReader)
{
  if (inputReader.adsorbentComponents.empty())
  {
    throw std::runtime_error("Error: MultibedColumn requires at least one adsorbent");
  }
  if (inputReader.adsorbentComponents.size() == 1 && !inputReader.debugForceMultibed)
  {
    throw std::runtime_error(
        "Error: MultibedColumn requires at least two adsorbents unless DebugForceMultibed is enabled");
  }
  if (inputReader.breakthroughIntegrator == 2)
  {
    throw std::runtime_error("Error: MultibedColumn does not support the SIRK3 integrator");
  }
}

template <typename ColumnType>
std::string componentFileName(const ColumnType& column, size_t component)
{
  return "component_" + std::to_string(component) + "_" + column.components[component].name + ".data";
}

template <typename ColumnType>
double averageMixturePredictionSteps(const ColumnType& column)
{
  return column.iastPerformance.second == 0
             ? 0.0
             : static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second);
}
}  // namespace

template <typename ColumnType>
Breakthrough<ColumnType>::Breakthrough(const InputReader& inputReader)
    : displayName(inputReader.displayName),
      carrierGasComponent(inputReader.carrierGasComponent),
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
      rk3(inputReader),
      sirk3(inputReader),
      cvode(inputReader),
      integrationScheme(BreakthroughIntegrationScheme(inputReader.breakthroughIntegrator))
{
  if constexpr (std::is_same_v<ColumnType, MultibedColumn>)
  {
    validateMultibedInput(inputReader);
  }
  else if (inputReader.adsorbentComponents.size() > 1)
  {
    throw std::runtime_error("Error: multiple adsorbents require Breakthrough<MultibedColumn>");
  }

  column.initialize();
  if (inputReader.readColumnFile.has_value())
  {
    column.readJSON(*inputReader.readColumnFile);
  }

  if (integrationScheme == BreakthroughIntegrationScheme::CVODE) cvode.initialize(column);

  const auto openMode = inputReader.readColumnFile.has_value() ? std::ios_base::app : std::ios_base::trunc;
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    componentStreams[comp] = std::ofstream{componentFileName(column, comp), openMode};
  }
  std::ofstream columnStream("column.data", openMode);
  if (!inputReader.readColumnFile.has_value())
  {
    column.writeOutputHeader(componentStreams, columnStream);
  }
}

template <typename ColumnType>
void Breakthrough<ColumnType>::run()
{
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    componentStreams[comp] = std::ofstream{componentFileName(column, comp), std::ios_base::app};
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  bool finished = false;
  size_t step = 0;
  double realTime = 0.0;
  constexpr double initialTimeStepFraction = 1.0e-4;

  if (numberOfInitTimeSteps > 0)
  {
    rk3.autoNumberOfSteps = false;
    sirk3.autoNumberOfSteps = false;
    cvode.autoNumberOfSteps = false;
    rk3.timeStep = timeStep * initialTimeStepFraction;
    sirk3.timeStep = timeStep * initialTimeStepFraction;
    cvode.timeStep = timeStep * initialTimeStepFraction;
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

      switch (integrationScheme)
      {
        case BreakthroughIntegrationScheme::SSP_RK:
        {
          finished = rk3.propagate(column, step, timings);
          realTime += rk3.timeStep;
          break;
        }
        case BreakthroughIntegrationScheme::CVODE:
        {
          finished = cvode.propagate(column, step, timings);
          realTime += cvode.timeStep;
          break;
        }
        case BreakthroughIntegrationScheme::SIRK3:
        {
          if constexpr (std::is_same_v<ColumnType, Column>)
          {
            finished = sirk3.propagate(column, step, timings);
            realTime += sirk3.timeStep;
          }
          else
          {
            throw std::runtime_error("Error: SIRK3 is not supported by MultibedColumn");
          }
          break;
        }
      }

      if (step < numberOfInitTimeSteps)
      {
        const double progress = static_cast<double>(step) / static_cast<double>(numberOfInitTimeSteps);
        const double nextTime = timeStep * (initialTimeStepFraction +
                                            (1.0 - initialTimeStepFraction) *
                                                (3.0 * progress * progress - 2.0 * progress * progress * progress));
        rk3.timeStep = nextTime;
        sirk3.timeStep = nextTime;
        cvode.timeStep = nextTime;
      }
      else if (step == numberOfInitTimeSteps)
      {
        rk3.autoNumberOfSteps = autoNumberOfSteps;
        sirk3.autoNumberOfSteps = autoNumberOfSteps;
        cvode.autoNumberOfSteps = autoNumberOfSteps;
        rk3.timeStep = timeStep;
        sirk3.timeStep = timeStep;
        cvode.timeStep = timeStep;
      }

      if (step % writeEvery == 0)
      {
        column.writeOutput(componentStreams, columnStream, stepTime);
      }
      if (step % printEvery == 0)
      {
        std::print("Timestep {}, time: {:6.5f} [s]\n", step, stepTime);
        std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageMixturePredictionSteps(column));
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
  if (integrationScheme == BreakthroughIntegrationScheme::CVODE) cvode.printStatistics();
}

template <typename ColumnType>
void Breakthrough<ColumnType>::computeStep(size_t step)
{
  switch (integrationScheme)
  {
    case BreakthroughIntegrationScheme::SSP_RK:
      static_cast<void>(rk3.propagate(column, step, timings));
      break;
    case BreakthroughIntegrationScheme::CVODE:
      static_cast<void>(cvode.propagate(column, step, timings));
      break;
    case BreakthroughIntegrationScheme::SIRK3:
      if constexpr (std::is_same_v<ColumnType, Column>)
      {
        static_cast<void>(sirk3.propagate(column, step, timings));
      }
      else
      {
        throw std::runtime_error("Error: SIRK3 is not supported by MultibedColumn");
      }
      break;
  }
}

template <typename ColumnType>
void Breakthrough<ColumnType>::print() const
{
  std::print("{}", repr());
}

template <typename ColumnType>
std::string Breakthrough<ColumnType>::repr() const
{
  std::string result = std::format(
      "Column properties\n"
      "=======================================================\n"
      "Display-name:                          {}\n",
      displayName);
  result += column.repr();

  result += std::format(
      "Breakthrough settings\n"
      "=======================================================\n");
  if constexpr (std::is_same_v<ColumnType, MultibedColumn>)
  {
    result += std::format("Number of adsorbents:          {}\n", column.numberOfAdsorbents);
  }
  result += std::format(
      "Number of time steps:          {}\n"
      "Print every step:              {}\n"
      "Write data every step:         {}\n\n\n"
      "Integration details\n"
      "=======================================================\n"
      "Time step:                     {} [s]\n"
      "Number of column grid points:  {}\n"
      "Column spacing:                {} [m]\n\n\n",
      numberOfSteps, printEvery, writeEvery, timeStep, numberOfGridPoints, column.resolution);

  if constexpr (std::is_same_v<ColumnType, MultibedColumn>)
  {
    result += std::format(
        "Adsorbent component data\n"
        "=======================================================\n"
        "maximum isotherm terms:        {}\n",
        maxIsothermTerms);
    for (size_t ads = 0; ads < column.numberOfAdsorbents; ++ads)
    {
      result += std::format("Adsorbent {}\n", ads);
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        result += std::format("{}\n", column.physisorptionMixtures[ads].components[comp].repr());
      }
    }
  }
  else
  {
    result += std::format(
        "Component data\n"
        "=======================================================\n"
        "maximum isotherm terms:        {}\n",
        maxIsothermTerms);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      result += std::format("{}\n", column.components[comp].repr());
    }
  }

  return result;
}

template struct Breakthrough<Column>;
template struct Breakthrough<MultibedColumn>;
