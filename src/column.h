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

struct Column
{
  enum class VelocityProfile
  {
    FixedPressureGradient = 0,  ///< Calculate velocities given a fixed gradient dPdx in the input file.
    Ergun = 1,                  ///< Calculate velocities from actual pressure gradient using Ergun equation.
    FixedVelocity = 2           ///< Ignore pressure gradient.
  };

  Column(const InputReader& inputReader)
      : mixture(MixturePrediction(inputReader)),
        components(inputReader.components),
        velocityProfile(VelocityProfile(inputReader.velocityProfile)),
        Ngrid(inputReader.numberOfGridPoints),
        Ncomp(inputReader.components.size()),
        maxIsothermTerms(inputReader.maxIsothermTerms),
        carrierGasComponent(inputReader.carrierGasComponent),
        externalTemperature(inputReader.temperature),
        externalPressure(inputReader.totalPressure),
        pressureGradient(inputReader.pressureGradient),
        voidFraction(inputReader.columnVoidFraction),
        particleDensity(inputReader.particleDensity),
        columnEntranceVelocity(inputReader.columnEntranceVelocity),
        columnLength(inputReader.columnLength),
        dynamicViscosity(inputReader.dynamicViscosity),
        particleDiameter(inputReader.particleDiameter),
        pulse(inputReader.pulseBreakthrough),
        pulseTime(inputReader.pulseTime),
        influxTemperature(inputReader.influxTemperature),
        internalDiameter(inputReader.internalDiameter),
        outerDiameter(inputReader.outerDiameter),
        wallDensity(inputReader.wallDensity),
        gasThermalConductivity(inputReader.gasThermalConductivity),
        wallThermalConductivity(inputReader.wallThermalConductivity),
        heatTransferGasSolid(inputReader.heatTransferGasSolid),
        heatTransferGasWall(inputReader.heatTransferGasWall),
        heatTransferWallExternal(inputReader.heatTransferWallExternal),
        heatCapacityGas(inputReader.heatCapacityGas),
        heatCapacitySolid(inputReader.heatCapacitySolid),
        heatCapacityWall(inputReader.heatCapacityWall),
        energyBalance(inputReader.energyBalance),
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
        gasTemperature(Ngrid + 1),
        gasTemperatureDot(Ngrid + 1),
        solidTemperature(Ngrid + 1),
        solidTemperatureDot(Ngrid + 1),
        wallTemperature(Ngrid + 1),
        wallTemperatureDot(Ngrid + 1),
        partialPressure((Ngrid + 1) * Ncomp),
        partialPressureDot((Ngrid + 1) * Ncomp),
        adsorption((Ngrid + 1) * Ncomp),
        adsorptionDot((Ngrid + 1) * Ncomp),
        equilibriumAdsorption((Ngrid + 1) * Ncomp),
        moleFraction((Ngrid + 1) * Ncomp),
        moleFractionDot((Ngrid + 1) * Ncomp),
        cachedPressure((Ngrid + 1) * Ncomp * maxIsothermTerms),
        cachedGrandPotential((Ngrid + 1) * maxIsothermTerms),
        gasDensity(Ngrid + 1),
        coeffGasGas(Ngrid + 1),
        coeffGasSolid(Ngrid + 1),
        coeffGasWall(Ngrid + 1),
        coeffDiffusion(Ngrid + 1),
        facePressures(Ngrid),
        energyFlow((Ngrid + 1) * Ncomp)

  {};

  MixturePrediction mixture;
  std::vector<Component> components;
  VelocityProfile velocityProfile;
  bool energyBalance;

  size_t Ngrid;
  size_t Ncomp;
  size_t maxIsothermTerms;
  size_t numCalls;
  size_t carrierGasComponent;

  double externalTemperature;
  double externalPressure;
  double pressureGradient;
  double voidFraction;
  double particleDensity;
  double columnEntranceVelocity;
  double columnLength;

  // ergun coefficients
  double dynamicViscosity;
  double particleDiameter;

  // coefficients for energy balance
  double influxTemperature;
  double internalDiameter;
  double outerDiameter;
  double wallDensity;
  double gasThermalConductivity;
  double wallThermalConductivity;
  double heatTransferGasSolid;
  double heatTransferGasWall;
  double heatTransferWallExternal;
  double heatCapacityGas;
  double heatCapacitySolid;
  double heatCapacityWall;

  // pulse vaars
  bool pulse;
  double pulseTime;

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
  std::vector<double> gasDensity;
  std::vector<double> totalPressure;        ///< Total pressure along the column.
  std::vector<double> totalPressureDot;     ///< Derivative wall temperature along the column wrt time.
  std::vector<double> gasTemperature;       ///< Total temperature along the column.
  std::vector<double> gasTemperatureDot;    ///< Derivative total temperature along the column wrt time.
  std::vector<double> solidTemperature;     ///< Total temperature along the column.
  std::vector<double> solidTemperatureDot;  ///< Derivative total temperature along the column wrt time.
  std::vector<double> wallTemperature;      ///< Wall temperature along the column.
  std::vector<double> wallTemperatureDot;   ///< Derivative wall temperature along the column wrt time.

  // vector of size '(Ngrid + 1) * Ncomp', for each grid point, data per component (contiguous)
  std::vector<double> partialPressure;     ///< Partial pressure at every grid point for each component.
  std::vector<double> partialPressureDot;  ///< Derivative of partial pressure wrt time.
  std::vector<double> adsorption;     ///< Volume-averaged adsorption amount at every grid point for each component.
  std::vector<double> adsorptionDot;  ///< Derivative of adsorption wrt time.
  std::vector<double> equilibriumAdsorption;  ///< Equilibrium adsorption amount at every grid point for each component.
  std::vector<double> moleFraction;           ///< Mole fraction of all components along column
  std::vector<double> moleFractionDot;        ///< Derivative mole fraction of all components along column wrt time.

  // vector of size '(Ngrid + 1) * Ngrid * maxIsothermTerms' for caching pressure
  std::vector<double> cachedPressure;  ///< Cached hypothetical pressure.
  // vector of size '(Ngrid + 1) * maxIsothermTerms' for caching grand potential
  std::vector<double> cachedGrandPotential;  ///< Cached reduced grand potential over the column.

  // scratch space vectors of size (Ngrid + 1) to prevent reallocation
  std::vector<double> coeffGasGas;
  std::vector<double> coeffGasSolid;
  std::vector<double> coeffGasWall;
  std::vector<double> coeffDiffusion;
  std::vector<double> facePressures;
  std::vector<double> energyFlow;

  void initialize();
  void writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& movieStream, double time);
  std::string repr() const;

  void writeJSON(const std::string& fileName) const;
  void readJSON(const std::string& fileName);
};
