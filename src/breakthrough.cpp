#include "breakthrough.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <fstream>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
std::string componentFileName(const Column& column, size_t component)
{
  return "component_" + std::to_string(component) + "_" + column.components[component].name + ".data";
}

double averageMixturePredictionSteps(const Column& column)
{
  return column.iastPerformance.second == 0
             ? 0.0
             : static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second);
}

/**
 * \brief True when every fed component leaves the column within 1% of its feed value, which is visibly
 * indistinguishable from a converged breakthrough.
 */
bool breakthroughConverged(const Column& column)
{
  double tolerance = 0.0;
  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    const bool gasPhase = static_cast<size_t>(column.fluidPhase) == 0;
    const double feed = gasPhase ? column.components[comp].initialGasMoleFraction
                                 : column.components[comp].inletLiquidConcentration;
    if (feed <= 0.0) continue;

    const size_t outlet = column.numberOfGridPoints * column.numberOfComponents + comp;
    const double outletValue = gasPhase ? column.moleFraction[outlet] : column.concentration[outlet];
    tolerance = std::max(tolerance, std::abs((outletValue / feed) - 1.0));
  }
  return tolerance < 0.01;
}
}  // namespace

Breakthrough::Breakthrough(const InputReader& inputReader)
    : displayName(inputReader.displayName),
      carrierGasComponent(inputReader.carrierGasComponent),
      numberOfComponents(inputReader.components.size()),
      numberOfGridPoints(inputReader.numberOfGridPoints),
      writeInterval(inputReader.writeInterval),
      simulationTime(inputReader.simulationTime),
      autoSimulationTime(inputReader.autoSimulationTime),
      maxIsothermTerms(inputReader.maxIsothermTerms),
      column(inputReader),
      rk3(inputReader),
      sirk3(inputReader),
      cvode(inputReader),
      integrationScheme(BreakthroughIntegrationScheme(inputReader.breakthroughIntegrator))
{
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

void Breakthrough::run()
{
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    componentStreams[comp] = std::ofstream{componentFileName(column, comp), std::ios_base::app};
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  try
  {
    bool finished = false;
    while (true)
    {
      column.writeJSON("column.json");
      column.writeOutput(componentStreams, columnStream, time);
      std::print("Output {}, time: {:6.5f} [s]\n", numberOfOutputs, time);
      std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageMixturePredictionSteps(column));
      std::fflush(stdout);

      if (finished) break;
      finished = advanceToNextOutput();
    }
  }
  catch (const std::runtime_error&)
  {
    column.writeJSON("failed_state.json");
    throw;
  }

  std::print("Final time: {:6.5f} [s]\n\n", time);
  timings.print();
  if (integrationScheme == BreakthroughIntegrationScheme::CVODE) cvode.printStatistics();
}

void Breakthrough::advanceTo(double endTime, Timing& accumulatedTimings)
{
  switch (integrationScheme)
  {
    case BreakthroughIntegrationScheme::SSP_RK:
      rk3.advance(column, time, endTime, accumulatedTimings);
      break;
    case BreakthroughIntegrationScheme::CVODE:
      cvode.advance(column, endTime, accumulatedTimings);
      break;
    case BreakthroughIntegrationScheme::SIRK3:
      sirk3.advance(column, time, endTime, accumulatedTimings);
      break;
  }
  time = endTime;
}

bool Breakthrough::advanceToNextOutput()
{
  // Output times are computed from their index, so that they do not drift over a long run.
  const double nextOutputTime = static_cast<double>(numberOfOutputs + 1) * writeInterval;
  const bool reachesEnd = !autoSimulationTime && nextOutputTime >= simulationTime - 1.0e-9 * writeInterval;
  advanceTo(reachesEnd ? simulationTime : nextOutputTime, timings);
  ++numberOfOutputs;
  if (reachesEnd) return true;

  if (autoSimulationTime)
  {
    bool converged = false;
    if (column.reactions.empty())
    {
      converged = breakthroughConverged(column);
      if (converged) std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
    }
    else
    {
      // RK3 leaves the rates of its last stage behind; CVODE has already evaluated them at the output state.
      if (integrationScheme != BreakthroughIntegrationScheme::CVODE) computeDerivatives(column, time);
      converged = reactionAutoStopReached(column);
      if (converged) std::print("\nReaction convergence criteria reached, running 10% longer\n\n\n");
    }
    if (converged)
    {
      // 10% longer for display purposes, rounded up to a whole output interval.
      simulationTime = static_cast<double>(numberOfOutputs + (numberOfOutputs + 9) / 10) * writeInterval;
      autoSimulationTime = false;
    }
  }
  return false;
}

void Breakthrough::print() const
{
  std::print("{}", repr());
}

std::string Breakthrough::repr() const
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
  result += std::format("Number of adsorbents:          {}\n", column.numberOfAdsorbents);
  result += std::format(
      "Simulation time:               {}\n"
      "Write interval:                {} [s]\n\n\n"
      "Integration details\n"
      "=======================================================\n",
      autoSimulationTime ? std::string{"auto"} : std::format("{} [s]", simulationTime), writeInterval);
  switch (integrationScheme)
  {
    case BreakthroughIntegrationScheme::SSP_RK:
      result += std::format("Time step:                     {} [s]\n", rk3.timeStep);
      break;
    case BreakthroughIntegrationScheme::CVODE:
      result += "Time step:                     adaptive (CVODE)\n";
      break;
    case BreakthroughIntegrationScheme::SIRK3:
      result += std::format("Time step:                     {} [s]\n", sirk3.timeStep);
      break;
  }
  result += std::format(
      "Number of column grid points:  {}\n"
      "Column spacing:                {} [m]\n\n\n",
      numberOfGridPoints, column.resolution);

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

  return result;
}
