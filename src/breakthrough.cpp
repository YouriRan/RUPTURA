#include <algorithm>
#include <cmath>
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
#include "mixture_prediction.h"
#include "rk3.h"
#include "utils.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

Breakthrough::Breakthrough(const InputReader& inputReader)
    : displayName(inputReader.displayName),
      carrierGasComponent(inputReader.carrierGasComponent),
      Ncomp(inputReader.components.size()),
      Ngrid(inputReader.numberOfGridPoints),
      printEvery(inputReader.printEvery),
      writeEvery(inputReader.writeEvery),
      dt(inputReader.timeStep),
      numberOfInitTimeSteps(inputReader.numberOfInitTimeSteps),
      Nsteps(inputReader.numberOfTimeSteps),
      autoSteps(inputReader.autoNumberOfTimeSteps),
      pulse(inputReader.pulseBreakthrough),
      tpulse(inputReader.pulseTime),
      maxIsothermTerms(inputReader.maxIsothermTerms),
      column(inputReader),
      rk3(dt, inputReader.autoNumberOfTimeSteps, inputReader.numberOfTimeSteps),
      sirk3(dt, inputReader.autoNumberOfTimeSteps, inputReader.numberOfTimeSteps),
      cvode(dt, inputReader.autoNumberOfTimeSteps, inputReader.numberOfTimeSteps),
      integrationScheme(IntegrationScheme(inputReader.breakthroughIntegrator))

{
  column.initialize();
  if (inputReader.readColumnFile.has_value())
  {
    column.readJSON(inputReader.readColumnFile.value());
  }

  if (integrationScheme == IntegrationScheme::CVODE) cvode.initialize(column);
}

// Breakthrough::Breakthrough(std::string _displayName, std::vector<Component> _components, size_t _carrierGasComponent,
//                            size_t _numberOfGridPoints, size_t _printEvery, size_t _writeEvery, double _temperature,
//                            double _externalPressure, double _columnVoidFraction, double _pressureGradient,
//                            double _particleDensity, double _columnEntranceVelocity, double _columnLength,
//                            double _timeStep, size_t _numberOfTimeSteps, bool _autoSteps, bool _pulse, double
//                            _pulseTime, double _particleDiameter, double _dynamicViscosity, const MixturePrediction
//                            _mixture, size_t _breakthroughIntegrator, size_t _velocityProfile,
//                            std::optional<std::string> readColumnFile)
//     : displayName(_displayName),
//       carrierGasComponent(_carrierGasComponent),
//       Ncomp(_components.size()),
//       Ngrid(_numberOfGridPoints),
//       printEvery(_printEvery),
//       writeEvery(_writeEvery),
//       dt(_timeStep),
//       Nsteps(_numberOfTimeSteps),
//       autoSteps(_autoSteps),
//       pulse(_pulse),
//       tpulse(_pulseTime),
//       maxIsothermTerms(_mixture.getMaxIsothermTerms()),
//       column(_mixture, _components, Ngrid, Ncomp, maxIsothermTerms, _temperature, _externalPressure,
//       _pressureGradient,
//              _columnVoidFraction, _particleDensity, _columnEntranceVelocity, _columnLength, _pulse, _pulseTime,
//              _particleDiameter, _dynamicViscosity, _carrierGasComponent, _velocityProfile),
//       rk3(dt, _autoSteps, _numberOfTimeSteps),
//       sirk3(dt, _autoSteps, _numberOfTimeSteps),
//       cvode(dt, _autoSteps, _numberOfTimeSteps),
//       integrationScheme(IntegrationScheme(_breakthroughIntegrator))
// {
//   column.initialize();
//   if (readColumnFile.has_value())
//   {
//     column.readJSON(readColumnFile.value());
//   }

//   if (integrationScheme == IntegrationScheme::CVODE) cvode.initialize(column);
// }

void Breakthrough::run()
{
  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < Ncomp; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + column.components[i].name + ".data";
    streams.emplace_back(std::ofstream{fileName});
  }

  std::ofstream movieStream("column.data");

  size_t column_nr = 1;
  std::print(movieStream, "# column {}: z (column position)\n", column_nr++);
  std::print(movieStream, "# column {}: V  (velocity)\n", column_nr++);
  std::print(movieStream, "# column {}: Pt (total pressure)\n", column_nr++);
  for (size_t j = 0; j < Ncomp; ++j)
  {
    std::print(movieStream, "# column {}: component {} Q     (loading)\n", column_nr++, j);
    std::print(movieStream, "# column {}: component {} Qeq   (equlibrium loading)\n", column_nr++, j);
    std::print(movieStream, "# column {}: component {} P     (partial pressure)\n", column_nr++, j);
    std::print(movieStream, "# column {}: component {} Pnorm (normalized partial pressure)\n", column_nr++, j);
    std::print(movieStream, "# column {}: component {} Dpdt  (derivative P with t)\n", column_nr++, j);
    std::print(movieStream, "# column {}: component {} Dqdt  (derivative Q with t)\n", column_nr++, j);
  }

  bool finished = false;
  size_t step = 0;
  double realTime = 0.0;
  double xi = 1e-4;

  if (numberOfInitTimeSteps > 0)
  {
    rk3.autoSteps = false;
    cvode.autoSteps = false;
    sirk3.autoSteps = false;
    rk3.timeStep = dt * xi;
    cvode.timeStep = dt * xi;
    sirk3.timeStep = dt * xi;
  }

  while (!finished)
  {
    // compute new step
    // computeStep(step);

    if (step % writeEvery == 0)
    {
      const std::string outputFile = std::format("column.json", step);
      column.writeJSON(outputFile);
    }

    switch (integrationScheme)
    {
      case IntegrationScheme::SSP_RK:
      {
        finished = rk3.propagate(column, step);
        realTime += rk3.timeStep;
        break;
      }
      case IntegrationScheme::CVODE:
      {
        finished = cvode.propagate(column, step);
        realTime += cvode.timeStep;
        break;
      }
      case IntegrationScheme::SIRK3:
      {
        finished = sirk3.propagate(column, step);
        realTime += sirk3.timeStep;
        break;
      }
      default:
        break;
    }

    if (step < numberOfInitTimeSteps)
    {
      double i = static_cast<double>(step) / numberOfInitTimeSteps;

      rk3.timeStep = dt * (xi + (1 - xi) * (3 * i * i - 2 * i * i * i));
      cvode.timeStep = dt * (xi + (1 - xi) * (3 * i * i - 2 * i * i * i));
      sirk3.timeStep = dt * (xi + (1 - xi) * (3 * i * i - 2 * i * i * i));
    }
    else if (step == numberOfInitTimeSteps)
    {
      rk3.autoSteps = autoSteps;
      cvode.autoSteps = autoSteps;
      sirk3.autoSteps = autoSteps;
      rk3.timeStep = dt;
      cvode.timeStep = dt;
      sirk3.timeStep = dt;
    }

    if (step % writeEvery == 0)
    {
      column.writeOutput(streams, movieStream, realTime);
    }
    if (step % printEvery == 0)
    {
      std::print("Timestep {}, time: {:6.5f} [s]\n", step, realTime);
      std::print(
          "    Average number of mixture-prediction steps: {:6.5f}\n",
          static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second));
      std::cout << std::flush;
    }

    step++;
  }

  std::print("Final timestep {}, time: {:6.5f} [s]\n", step, dt * static_cast<double>(step));
}

#ifdef PYBUILD

py::array_t<double> Breakthrough::compute()
{
  size_t colsize = 6 * Ncomp + 5;
  std::vector<std::vector<std::vector<double>>> brk;

  // loop can quit early if autoSteps
  for (size_t step = 0; (step < Nsteps || autoSteps); ++step)
  {
    // check for error from python side (keyboard interrupt)
    if (PyErr_CheckSignals() != 0)
    {
      throw py::error_already_set();
    }

    computeStep(step);
    double t = static_cast<double>(step) * dt;
    if (step % writeEvery == 0)
    {
      std::vector<std::vector<double>> t_brk(Ngrid + 1, std::vector<double>(colsize));
      for (size_t i = 0; i < Ngrid + 1; ++i)
      {
        t_brk[i][0] = t * v_in / L;
        t_brk[i][1] = t / 60.0;
        t_brk[i][2] = static_cast<double>(i) * dx;
        t_brk[i][3] = V[i];
        t_brk[i][4] = Pt[i];

        for (size_t j = 0; j < Ncomp; ++j)
        {
          t_brk[i][5 + 6 * j] = Q[i * Ncomp + j];
          t_brk[i][6 + 6 * j] = Qeq[i * Ncomp + j];
          t_brk[i][7 + 6 * j] = P[i * Ncomp + j];
          t_brk[i][8 + 6 * j] = P[i * Ncomp + j] / (Pt[i] * column.components[j].Yi0);
          t_brk[i][9 + 6 * j] = Dpdt[i * Ncomp + j];
          t_brk[i][10 + 6 * j] = Dqdt[i * Ncomp + j];
        }
      }
      brk.push_back(t_brk);
    }
    if (step % printEvery == 0)
    {
      std::cout << "Timestep " + std::to_string(step) + ", time: " + std::to_string(t) + " [s]" << std::endl;
      std::cout << "    Average number of mixture-prediction steps: " +
                       std::to_string(static_cast<double>(iastPerformance.first) /
                                      static_cast<double>(iastPerformance.second))
                << std::endl;
    }
  }
  std::cout << "Final timestep " + std::to_string(Nsteps) +
                   ", time: " + std::to_string(dt * static_cast<double>(Nsteps)) + " [s]"
            << std::endl;

  std::vector<double> buffer;
  buffer.reserve(brk.size() * (Ngrid + 1) * colsize);
  for (const auto& vec1 : brk)
  {
    for (const auto& vec2 : vec1)
    {
      buffer.insert(buffer.end(), vec2.begin(), vec2.end());
    }
  }
  std::array<size_t, 3> shape{{brk.size(), Ngrid + 1, colsize}};
  py::array_t<double> py_breakthrough(shape, buffer.data());
  py_breakthrough.resize(shape);

  return py_breakthrough;
}

void Breakthrough::setComponentsParameters(std::vector<double> molfracs, std::vector<double> params)
{
  size_t index = 0;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    column.components[i].Yi0 = molfracs[i];
    size_t n_params = column.components[i].isotherm.numberOfParameters;
    std::vector<double> slicedVec(params.begin() + index, params.begin() + index + n_params);
    index = index + n_params;
    column.components[i].isotherm.setParameters(slicedVec);
  }

  // also set for mixture
  mixture.setComponentsParameters(molfracs, params);
}

std::vector<double> Breakthrough::getComponentsParameters()
{
  std::vector<double> params;
  for (size_t i = 0; i < Ncomp; ++i)
  {
    std::vector<double> compParams = column.components[i].isotherm.getParameters();
    params.insert(params.end(), compParams.begin(), compParams.end());
  }
  return params;
}
#endif  // PYBUILD

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
      Nsteps, printEvery, writeEvery);

  // Integration details
  s += std::format(
      "Integration details\n"
      "=======================================================\n"
      "Time step:                     {} [s]\n"
      "Number of column grid points:  {}\n"
      "Column spacing:                {} [m]\n"
      "\n\n",
      dt, Ngrid, column.resolution);

  // Component data header
  s += std::format(
      "Component data\n"
      "=======================================================\n"
      "maximum isotherm terms:        {}\n",
      maxIsothermTerms);

  // Append each component’s repr()
  for (std::size_t i = 0; i < Ncomp; ++i)
  {
    s += std::format("{}\n", column.components[i].repr());
  }

  return s;
}