#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#if __cplusplus >= 201703L && __has_include(<filesystem>)
#include <filesystem>
#elif __cplusplus >= 201703L && __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#else
#include <sys/stat.h>
#endif

#include "breakthrough.h"
#include "compute.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "utils.h"

Breakthrough::Breakthrough(const InputReader& inputReader)
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
      integrationScheme(IntegrationScheme(inputReader.breakthroughIntegrator))
{
  column.initialize();
  if (inputReader.readColumnFile.has_value())
  {
    column.readJSON(inputReader.readColumnFile.value());
  }

  if (integrationScheme == IntegrationScheme::CVODE) cvode.initialize(column);

  // open the output files in append mode, unless there is a
  // create the output files

  auto openMode = (inputReader.readColumnFile.has_value()) ? std::ios_base::app : std::ios_base::trunc;
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; comp++)
  {
    std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    componentStreams[comp] = std::ofstream{fileName, openMode};
  }
  std::ofstream columnStream("column.data", openMode);
  if (!inputReader.readColumnFile.has_value()) column.writeOutputHeader(componentStreams, columnStream);
}

void Breakthrough::run()
{
  std::vector<std::ofstream> componentStreams(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; comp++)
  {
    std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    componentStreams[comp] = std::ofstream{fileName, std::ios_base::app};
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  bool finished = false;
  size_t step = 0;
  double realTime = 0.0;
  double xi = 1e-4;

  if (numberOfInitTimeSteps > 0)
  {
    rk3.autoNumberOfSteps = false;
    cvode.autoNumberOfSteps = false;
    sirk3.autoNumberOfSteps = false;
    rk3.timeStep = timeStep * xi;
    cvode.timeStep = timeStep * xi;
    sirk3.timeStep = timeStep * xi;
  }

  try
  {
    while (!finished)
    {
      if (step % writeEvery == 0)
      {
        const std::string outputFile = std::format("column.json", step);
        column.writeJSON(outputFile);
      }

      switch (integrationScheme)
      {
        case IntegrationScheme::SSP_RK:
        {
          finished = rk3.propagate(column, step, timings);
          realTime += rk3.timeStep;
          break;
        }
        case IntegrationScheme::CVODE:
        {
          finished = cvode.propagate(column, step, timings);
          realTime += cvode.timeStep;
          break;
        }
        case IntegrationScheme::SIRK3:
        {
          finished = sirk3.propagate(column, step, timings);
          realTime += sirk3.timeStep;
          break;
        }
        default:
          break;
      }

      if (step < numberOfInitTimeSteps)
      {
        double i = static_cast<double>(step) / numberOfInitTimeSteps;
        double nextTime = timeStep * (xi + (1 - xi) * (3 * i * i - 2 * i * i * i));
        rk3.timeStep = nextTime;
        cvode.timeStep = nextTime;
        sirk3.timeStep = nextTime;
      }
      else if (step == numberOfInitTimeSteps)
      {
        rk3.autoNumberOfSteps = autoNumberOfSteps;
        cvode.autoNumberOfSteps = autoNumberOfSteps;
        sirk3.autoNumberOfSteps = autoNumberOfSteps;
        rk3.timeStep = timeStep;
        cvode.timeStep = timeStep;
        sirk3.timeStep = timeStep;
      }

      if (step % writeEvery == 0)
      {
        column.writeOutput(componentStreams, columnStream, realTime);
      }
      if (step % printEvery == 0)
      {
        std::print("Timestep {}, time: {:6.5f} [s]\n", step, realTime);
        std::print(
            "    Average number of mixture-prediction steps: {:6.5f}\n",
            static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second));
        std::fflush(stdout);
      }

      step++;
    }
  }
  catch (const std::runtime_error&)
  {
    column.writeJSON("failed_state.json");
    throw;
  }

  std::print("Final timestep {}, time: {:6.5f} [s]\n\n", step, timeStep * static_cast<double>(step));

  timings.print();
}

void Breakthrough::print() const { std::print("{}", repr()); }

std::string Breakthrough::repr() const
{
  std::string s;

  // Column properties
  s += std::format(
      "Column properties\n"
      "=======================================================\n"
      "Display-name:                          {}\n",
      displayName);

  s += column.repr();

  // Breakthrough settings
  s += std::format(
      "Breakthrough settings\n"
      "=======================================================\n"
      "Number of time steps:          {}\n"
      "Print every step:              {}\n"
      "Write data every step:         {}\n"
      "\n\n",
      numberOfSteps, printEvery, writeEvery);

  // Integration details
  s += std::format(
      "Integration details\n"
      "=======================================================\n"
      "Time step:                     {} [s]\n"
      "Number of column grid points:  {}\n"
      "Column spacing:                {} [m]\n"
      "\n\n",
      timeStep, numberOfGridPoints, column.resolution);

  // Component data header
  s += std::format(
      "Component data\n"
      "=======================================================\n"
      "maximum isotherm terms:        {}\n",
      maxIsothermTerms);

  // Append each component’s repr()
  for (std::size_t i = 0; i < numberOfComponents; ++i)
  {
    s += std::format("{}\n", column.components[i].repr());
  }

  return s;
}
