#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <print>
#include <sstream>
#include <string>

#include "breakthrough.h"
#include "component.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "utils.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

void BreakthroughState::initialize()
{
  // precomputed factor for mass transfer
  for (size_t j = 0; j < Ncomp; ++j)
  {
    prefactorMassTransfer[j] =
        R * externalTemperature * ((1.0 - voidFraction) / voidFraction) * particleDensity * components[j].Kl;
  }

  // set P and Q to zero
  std::fill(partialPressure.begin(), partialPressure.end(), 0.0);
  std::fill(adsorption.begin(), adsorption.end(), 0.0);

  // initial pressure along the column
  std::vector<double> initialPressure(Ngrid + 1);

  // set the initial total pressure along the column assuming the pressure gradient is constant
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    initialPressure[i] = externalPressure + pressureGradient * static_cast<double>(i) * resolution;
  }

  // initialize the interstitial gas velocity in the column
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    interstitialGasVelocity[i] = columnEntranceVelocity * externalPressure / initialPressure[i];
  }

  // set the partial pressure of the carrier gas to the total initial pressure
  // for the column except for the entrance (i=0)
  for (size_t i = 1; i < Ngrid + 1; ++i)
  {
    partialPressure[i * Ncomp + carrierGasComponent] = initialPressure[i];
  }

  // at the column entrance, the mol-fractions of the components in the gas phase are fixed
  // the partial pressures of the components at the entrance are the mol-fractions times the
  // total pressure
  for (size_t j = 0; j < Ncomp; ++j)
  {
    partialPressure[0 * Ncomp + j] = externalPressure * components[j].Yi0;
  }

  // at the entrance: mol-fractions Yi are the gas-phase mol-fractions
  // for the column: the initial mol-fraction of the carrier-gas is 1, and 0 for the other components
  //
  // the K of the carrier gas is chosen as zero
  // so Qeq is zero for all components in the column after the entrance
  // only the values for Yi at the entrance are effected by adsorption
  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    double sum = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      idealGasMolFractions[j] = std::max(partialPressure[i * Ncomp + j] / initialPressure[i], 0.0);
      sum += idealGasMolFractions[j];
    }
    for (size_t j = 0; j < Ncomp; ++j)
    {
      idealGasMolFractions[j] /= sum;
    }

    iastPerformance += mixture.predictMixture(idealGasMolFractions, initialPressure[i], adsorbedMolFractions,
                                              numberOfMolecules, &cachedPressure[i * Ncomp * maxIsothermTerms],
                                              &cachedGrandPotential[i * maxIsothermTerms]);

    for (size_t j = 0; j < Ncomp; ++j)
    {
      equilibriumAdsorption[i * Ncomp + j] = numberOfMolecules[j];
    }
  }

  for (size_t i = 0; i < Ngrid + 1; ++i)
  {
    totalPressure[i] = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      totalPressure[i] += std::max(0.0, partialPressure[i * Ncomp + j]);
    }
  }
}

void BreakthroughState::writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& movieStream,
                                    double time)
{
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    // write breakthrough output to files
    // column 1: dimensionless time
    // column 2: time [minutes]
    // column 3: normalized partial pressure
    double normalizedPressure = partialPressure[Ngrid * Ncomp + comp] / (exitPressure * components[comp].Yi0);
    std::print(componentStreams[comp], "{} {} {}\n", time * timeNormalizationFactor, time / 60.0, normalizedPressure);
  }

  for (size_t grid = 0; grid < Ngrid; ++grid)
  {
    // write breakthrough output to files
    // column 1: column position
    // column 2: velocity
    // column 3: total pressure
    std::print(movieStream, "{} {} {} ", static_cast<double>(grid) * resolution, interstitialGasVelocity[grid],
               totalPressure[grid]);
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      // write breakthrough output to files
      // column 4 + 6 * comp: loading
      // column 5 + 6 * comp: equlibrium loading
      // column 6 + 6 * comp: partial pressure
      // column 7 + 6 * comp: normalized partial pressure
      // column 8 + 6 * comp: pressure time derivative
      // column 9 + 6 * comp: adsorption time derivative
      std::print(movieStream, "{} {} {} {} {} {} ", adsorption[grid * Ncomp + comp],
                 equilibriumAdsorption[grid * Ncomp + comp], partialPressure[grid * Ncomp + comp],
                 partialPressure[grid * Ncomp + comp] / (totalPressure[grid] * components[comp].Yi0),
                 partialPressureDot[grid * Ncomp + comp], adsorptionDot[grid * Ncomp + comp]);
    }
    std::print(movieStream, "\n");
  }
  std::print(movieStream, "\n\n");
}

std::string BreakthroughState::repr() const
{
  std::string s;

  // Column properties
  return std::format(
      "Temperature:                           {} [K]\n"
      "Column length:                         {} [m]\n"
      "Column void-fraction:                  {} [-]\n"
      "Particle density:                      {} [kg/m^3]\n"
      "Total pressure:                        {} [Pa]\n"
      "Pressure gradient:                     {} [Pa/m]\n"
      "Column entrance interstitial velocity: {} [m/s]\n"
      "\n\n",
      externalTemperature, columnLength, voidFraction, particleDensity, externalPressure, pressureGradient,
      columnEntranceVelocity);
}