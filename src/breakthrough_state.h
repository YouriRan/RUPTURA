#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>

#include "component.h"
#include "inputreader.h"
#include "mixture_prediction.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

struct BreakthroughState
{
  BreakthroughState(MixturePrediction mixture, std::vector<Component> components, size_t Ngrid, size_t Ncomp,
                    size_t maxIsothermTerms, double externalTemperature, double externalPressure,
                    double pressureGradient, double voidFraction, double particleDensity, double columnEntranceVelocity,
                    double columnLength, bool pulse, double pulseTime, size_t carrierGasComponent)
      : mixture(mixture),
        components(components),
        Ngrid(Ngrid),
        Ncomp(Ncomp),
        maxIsothermTerms(maxIsothermTerms),
        externalTemperature(externalTemperature),
        externalPressure(externalPressure),
        pressureGradient(pressureGradient),
        voidFraction(voidFraction),
        particleDensity(particleDensity),
        columnEntranceVelocity(columnEntranceVelocity),
        columnLength(columnLength),
        pulse(pulse),
        pulseTime(pulseTime),
        carrierGasComponent(carrierGasComponent),
        resolution(columnLength / static_cast<double>(Ngrid)),
        timeNormalizationFactor(columnEntranceVelocity / columnLength),
        exitPressure(externalPressure + pressureGradient * columnLength),
        prefactorMassTransfer(Ncomp),
        idealGasMolFractions(Ncomp),
        adsorbedMolFractions(Ncomp),
        numberOfMolecules(Ncomp),
        interstitialGasVelocity(Ngrid + 1),
        totalPressure(Ngrid + 1),
        totalPressureDot(Ngrid + 1),
        temperature(Ngrid + 1),
        temperatureDot(Ngrid + 1),
        wallTemperature(Ngrid + 1),
        wallTemperatureDot(Ngrid + 1),
        partialPressure((Ngrid + 1) * Ncomp),
        pressureDot((Ngrid + 1) * Ncomp),
        adsorption((Ngrid + 1) * Ncomp),
        adsorptionDot((Ngrid + 1) * Ncomp),
        equilibriumAdsorption((Ngrid + 1) * Ncomp),
        moleFraction((Ngrid + 1) * Ncomp),
        moleFractionDot((Ngrid + 1) * Ncomp),
        cachedPressure((Ngrid + 1) * Ncomp * maxIsothermTerms),
        cachedGrandPotential((Ngrid + 1) * maxIsothermTerms) {};

  MixturePrediction mixture;
  std::vector<Component> components;

  size_t Ngrid;
  size_t Ncomp;
  size_t maxIsothermTerms;
  size_t numCalls;

  double externalTemperature;
  double externalPressure;
  double pressureGradient;
  double voidFraction;
  double particleDensity;
  double columnEntranceVelocity;
  double columnLength;
  bool pulse;
  double pulseTime;
  size_t carrierGasComponent;

  // derived quantities
  double resolution;
  double timeNormalizationFactor;
  double exitPressure;

  std::pair<size_t, size_t> iastPerformance{
      0,
      0,
  };

  // vector of size 'Ncomp'
  std::vector<double> prefactorMassTransfer;  ///< Precomputed factors for mass transfer.
  std::vector<double> idealGasMolFractions;   ///< Ideal gas mole fractions for each component.
  std::vector<double> adsorbedMolFractions;   ///< Adsorbed mole fractions for each component.
  std::vector<double> numberOfMolecules;      ///< Number of molecules for each component.

  // vector of size '(Ngrid + 1)'
  std::vector<double> interstitialGasVelocity;  ///< Interstitial gas velocity along the column.
  std::vector<double> totalPressure;            ///< Total pressure along the column.
  std::vector<double> totalPressureDot;         ///< Derivative wall temperature along the column wrt time.
  std::vector<double> temperature;              ///< Total temperature along the column.
  std::vector<double> temperatureDot;           ///< Derivative total temperature along the column wrt time.
  std::vector<double> wallTemperature;          ///< Wall temperature along the column.
  std::vector<double> wallTemperatureDot;       ///< Derivative wall temperature along the column wrt time.

  // vector of size '(Ngrid + 1) * Ncomp', for each grid point, data per component (contiguous)
  std::vector<double> partialPressure;  ///< Partial pressure at every grid point for each component.
  std::vector<double> pressureDot;      ///< Derivative of partial pressure wrt time.
  std::vector<double> adsorption;       ///< Volume-averaged adsorption amount at every grid point for each component.
  std::vector<double> adsorptionDot;    ///< Derivative of adsorption wrt time.
  std::vector<double> equilibriumAdsorption;  ///< Equilibrium adsorption amount at every grid point for each component.
  std::vector<double> moleFraction;           ///< Mole fraction of all components along column
  std::vector<double> moleFractionDot;        ///< Derivative mole fraction of all components along column wrt time.

  // vector of size '(Ngrid + 1) * Ngrid * maxIsothermTerms' for caching pressure
  std::vector<double> cachedPressure;  ///< Cached hypothetical pressure.
  // vector of size '(Ngrid + 1) * maxIsothermTerms' for caching grand potential
  std::vector<double> cachedGrandPotential;  ///< Cached reduced grand potential over the column.

  void initialize();
  void writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& movieStream, double time);
  std::string repr() const;
};