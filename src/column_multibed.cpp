#include "column_multibed.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "json.h"
#include "mixture_prediction.h"
#include "utils.h"


namespace
{
std::vector<double> toVector(std::span<const double> s) { return std::vector<double>(s.begin(), s.end()); }

std::string vectorToString(const std::vector<double>& values)
{
  std::ostringstream stream;
  stream << "[";
  for (size_t i = 0; i < values.size(); ++i)
  {
    if (i != 0) stream << ", ";
    stream << values[i];
  }
  stream << "]";
  return stream.str();
}

void validateAdsorbentVectors(const ColumnMultibed& column)
{
  if (column.numberOfAdsorbents == 0)
  {
    throw std::runtime_error("Error: multibed column requires at least one adsorbent");
  }
  if (column.adsorbentLengths.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: AdsorbentLengths size must equal numberOfAdsorbents");
  }
  if (column.adsorbentInterfaceLengths.size() != column.numberOfAdsorbents - 1)
  {
    throw std::runtime_error("Error: AdsorbentInterfaceLengths size must equal numberOfAdsorbents - 1");
  }
  if (column.adsorbentVoidFractions.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: AdsorbentVoidFractions size must equal numberOfAdsorbents");
  }
  if (column.particleDensities.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: ParticleDensities size must equal numberOfAdsorbents");
  }
  if (column.particleDiameters.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: ParticleDiameters size must equal numberOfAdsorbents");
  }
  if (column.columnLength <= 0.0)
  {
    throw std::runtime_error("Error: sum of AdsorbentLengths must be positive");
  }
}
}  // namespace

size_t ColumnMultibed::stateSize() const noexcept { return (2 * numberOfComponents + 3) * (numberOfGridPoints + 1); }

void ColumnMultibed::bindStateViews() noexcept
{
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;

  double* base = state.data();
  double* baseDot = stateDot.data();

  concentration = std::span<double>(base + 0 * big, big);
  physisorption = std::span<double>(base + 1 * big, big);
  gasTemperature = std::span<double>(base + 2 * big + 0 * small, small);
  solidTemperature = std::span<double>(base + 2 * big + 1 * small, small);
  wallTemperature = std::span<double>(base + 2 * big + 2 * small, small);

  concentrationDot = std::span<double>(baseDot + 0 * big, big);
  physisorptionDot = std::span<double>(baseDot + 1 * big, big);
  gasTemperatureDot = std::span<double>(baseDot + 2 * big + 0 * small, small);
  solidTemperatureDot = std::span<double>(baseDot + 2 * big + 1 * small, small);
  wallTemperatureDot = std::span<double>(baseDot + 2 * big + 2 * small, small);
}

ColumnMultibed::ColumnMultibed(const ColumnMultibed& other)
    : mixture(other.mixture),
      components(other.components),
      boundaryCondition(other.boundaryCondition),
      energyBalance(other.energyBalance),
      numberOfGridPoints(other.numberOfGridPoints),
      numberOfComponents(other.numberOfComponents),
      numberOfAdsorbents(other.numberOfAdsorbents),
      maxIsothermTerms(other.maxIsothermTerms),
      numberOfCalls(other.numberOfCalls),
      carrierGasComponent(other.carrierGasComponent),
      externalTemperature(other.externalTemperature),
      inletPressure(other.inletPressure),
      outletPressure(other.outletPressure),
      pressureGradient(other.pressureGradient),
      adsorbentLengths(other.adsorbentLengths),
      adsorbentInterfaceLengths(other.adsorbentInterfaceLengths),
      adsorbentGridPoints(other.adsorbentGridPoints),
      adsorbentVoidFractions(other.adsorbentVoidFractions),
      particleDensities(other.particleDensities),
      columnEntranceVelocity(other.columnEntranceVelocity),
      columnLength(other.columnLength),
      dynamicViscosity(other.dynamicViscosity),
      columnDistances(other.columnDistances),
      particleDiameters(other.particleDiameters),
      influxTemperature(other.influxTemperature),
      internalDiameter(other.internalDiameter),
      outerDiameter(other.outerDiameter),
      wallDensity(other.wallDensity),
      gasThermalConductivity(other.gasThermalConductivity),
      wallThermalConductivity(other.wallThermalConductivity),
      heatTransferGasSolid(other.heatTransferGasSolid),
      heatTransferGasWall(other.heatTransferGasWall),
      heatTransferWallExternal(other.heatTransferWallExternal),
      heatCapacityGas(other.heatCapacityGas),
      heatCapacitySolid(other.heatCapacitySolid),
      heatCapacityWall(other.heatCapacityWall),
      resolution(other.resolution),
      timeNormalizationFactor(other.timeNormalizationFactor),
      iastPerformance(other.iastPerformance),
      prefactorMassTransfer(other.prefactorMassTransfer),
      idealGasMolFractions(other.idealGasMolFractions),
      adsorbedMolFractions(other.adsorbedMolFractions),
      numberOfMolecules(other.numberOfMolecules),
      interstitialGasVelocity(other.interstitialGasVelocity),
      gasDensity(other.gasDensity),
      totalConcentration(other.totalConcentration),
      totalPressure(other.totalPressure),
      totalVoidFraction(other.totalVoidFraction),
      particleDensity(other.particleDensity),
      moleFraction(other.moleFraction),
      partialPressure(other.partialPressure),
      equilibriumAdsorption(other.equilibriumAdsorption),
      fractionOfAdsorbent(other.fractionOfAdsorbent),
      hasAdsorbentOfType(other.hasAdsorbentOfType),
      adsorbentScaledVoidFraction(other.adsorbentScaledVoidFraction),
      cachedPressure(other.cachedPressure),
      cachedGrandPotential(other.cachedGrandPotential),
      coeffDiffusion(other.coeffDiffusion),
      facePressures(other.facePressures),
      massFlux(other.massFlux),
      state(other.state),
      stateDot(other.stateDot)
{
  bindStateViews();
}

ColumnMultibed& ColumnMultibed::operator=(const ColumnMultibed& other)
{
  if (this == &other) return *this;

  mixture = other.mixture;
  components = other.components;
  boundaryCondition = other.boundaryCondition;
  energyBalance = other.energyBalance;
  numberOfGridPoints = other.numberOfGridPoints;
  numberOfComponents = other.numberOfComponents;
  numberOfAdsorbents = other.numberOfAdsorbents;
  maxIsothermTerms = other.maxIsothermTerms;
  numberOfCalls = other.numberOfCalls;
  carrierGasComponent = other.carrierGasComponent;
  externalTemperature = other.externalTemperature;
  inletPressure = other.inletPressure;
  outletPressure = other.outletPressure;
  pressureGradient = other.pressureGradient;
  adsorbentLengths = other.adsorbentLengths;
  adsorbentInterfaceLengths = other.adsorbentInterfaceLengths;
  adsorbentGridPoints = other.adsorbentGridPoints;
  adsorbentVoidFractions = other.adsorbentVoidFractions;
  particleDensities = other.particleDensities;
  columnEntranceVelocity = other.columnEntranceVelocity;
  columnLength = other.columnLength;
  dynamicViscosity = other.dynamicViscosity;
  columnDistances = other.columnDistances;
  particleDiameters = other.particleDiameters;
  influxTemperature = other.influxTemperature;
  internalDiameter = other.internalDiameter;
  outerDiameter = other.outerDiameter;
  wallDensity = other.wallDensity;
  gasThermalConductivity = other.gasThermalConductivity;
  wallThermalConductivity = other.wallThermalConductivity;
  heatTransferGasSolid = other.heatTransferGasSolid;
  heatTransferGasWall = other.heatTransferGasWall;
  heatTransferWallExternal = other.heatTransferWallExternal;
  heatCapacityGas = other.heatCapacityGas;
  heatCapacitySolid = other.heatCapacitySolid;
  heatCapacityWall = other.heatCapacityWall;
  resolution = other.resolution;
  timeNormalizationFactor = other.timeNormalizationFactor;
  iastPerformance = other.iastPerformance;
  prefactorMassTransfer = other.prefactorMassTransfer;
  idealGasMolFractions = other.idealGasMolFractions;
  adsorbedMolFractions = other.adsorbedMolFractions;
  numberOfMolecules = other.numberOfMolecules;
  interstitialGasVelocity = other.interstitialGasVelocity;
  gasDensity = other.gasDensity;
  totalConcentration = other.totalConcentration;
  totalPressure = other.totalPressure;
  totalVoidFraction = other.totalVoidFraction;
  particleDensity = other.particleDensity;
  moleFraction = other.moleFraction;
  partialPressure = other.partialPressure;
  equilibriumAdsorption = other.equilibriumAdsorption;
  fractionOfAdsorbent = other.fractionOfAdsorbent;
  hasAdsorbentOfType = other.hasAdsorbentOfType;
  adsorbentScaledVoidFraction = other.adsorbentScaledVoidFraction;
  cachedPressure = other.cachedPressure;
  cachedGrandPotential = other.cachedGrandPotential;
  coeffDiffusion = other.coeffDiffusion;
  facePressures = other.facePressures;
  massFlux = other.massFlux;
  state = other.state;
  stateDot = other.stateDot;

  bindStateViews();
  return *this;
}

void ColumnMultibed::initialize()
{
  validateAdsorbentVectors(*this);

  std::vector<double> boundaries(numberOfAdsorbents > 1 ? numberOfAdsorbents - 1 : 0, 0.0);
  double cumulativeLength = 0.0;
  for (size_t ads = 0; ads + 1 < numberOfAdsorbents; ++ads)
  {
    cumulativeLength += adsorbentLengths[ads];
    boundaries[ads] = cumulativeLength;
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const double z = columnDistances[grid];
    const size_t base = grid * numberOfAdsorbents;

    size_t region = 0;
    while (region + 1 < numberOfAdsorbents && z > boundaries[region])
    {
      ++region;
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      fractionOfAdsorbent[base + ads] = ads == region ? 1.0 : 0.0;
    }

    for (size_t interface = 0; interface < adsorbentInterfaceLengths.size(); ++interface)
    {
      const double interfaceLength = adsorbentInterfaceLengths[interface];
      if (interfaceLength <= 0.0) continue;

      const double start = std::max(0.0, boundaries[interface] - 0.5 * interfaceLength);
      const double end = std::min(columnLength, boundaries[interface] + 0.5 * interfaceLength);
      if (z < start || z > end || end <= start) continue;

      const double rightFraction = (z - start) / (end - start);
      fractionOfAdsorbent[base + interface] = 1.0 - rightFraction;
      fractionOfAdsorbent[base + interface + 1] = rightFraction;
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        if (ads != interface && ads != interface + 1)
        {
          fractionOfAdsorbent[base + ads] = 0.0;
        }
      }
    }

    totalVoidFraction[grid] = 1.0;
    particleDensity[grid] = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[base + ads];
      hasAdsorbentOfType[base + ads] = fraction > 1.0e-4;
      totalVoidFraction[grid] -= fraction * (1.0 - adsorbentVoidFractions[ads]);
      particleDensity[grid] += fraction * particleDensities[ads];
    }

    if (totalVoidFraction[grid] <= 0.0)
    {
      throw std::runtime_error("Error: total void fraction must be positive");
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (particleDiameters[ads] <= 0.0)
      {
        throw std::runtime_error("Error: ParticleDiameters values must be positive");
      }
      adsorbentScaledVoidFraction[base + ads] =
          (1.0 - adsorbentVoidFractions[ads]) / (totalVoidFraction[grid] * particleDiameters[ads]);
    }
  }

  for (size_t j = 0; j < numberOfComponents; ++j)
  {
    prefactorMassTransfer[j] = components[j].massTransferCoefficient;
  }

  std::fill(partialPressure.begin(), partialPressure.end(), 0.0);
  std::fill(physisorption.begin(), physisorption.end(), 0.0);
  std::fill(concentration.begin(), concentration.end(), 0.0);
  std::fill(moleFraction.begin(), moleFraction.end(), 0.0);
  std::fill(stateDot.begin(), stateDot.end(), 0.0);

  std::vector<double> initialPressure(numberOfGridPoints + 1, 0.0);

  auto gridRatio = [&](size_t i) -> double
  { return columnLength <= 0.0 ? 0.0 : columnDistances[i] / columnLength; };

  auto fillPressure = [&](double pressure) { std::fill(initialPressure.begin(), initialPressure.end(), pressure); };

  auto fillVelocity = [&](double velocity)
  { std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), velocity); };

  switch (boundaryCondition)
  {
    case BoundaryCondition::InletPressureInletVelocity:
    {
      fillPressure(inletPressure);
      fillVelocity(columnEntranceVelocity);
      break;
    }
    case BoundaryCondition::InletPressureOutletPressure:
    {
      for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
      {
        initialPressure[i] = inletPressure + gridRatio(i) * (outletPressure - inletPressure);
      }
      fillVelocity(columnEntranceVelocity);
      break;
    }
    case BoundaryCondition::InletVelocityOutletPressure:
    {
      fillPressure(outletPressure);
      fillVelocity(columnEntranceVelocity);
      break;
    }
    case BoundaryCondition::FixedVelocity:
    {
      fillPressure(inletPressure > 0.0 ? inletPressure : outletPressure);
      fillVelocity(columnEntranceVelocity);
      break;
    }
    case BoundaryCondition::FixedPressureInletVelocity:
    {
      for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
      {
        initialPressure[i] = inletPressure + pressureGradient * gridRatio(i) / columnLength;
      }
      fillVelocity(columnEntranceVelocity);
      break;
    }
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    if (initialPressure[grid] <= 0.0)
    {
      throw std::runtime_error("Error: initialized pressure must be positive");
    }
  }

  // Internal nodes initially contain carrier gas; inlet node contains the feed mixture.
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      partialPressure[grid * numberOfComponents + comp] = components[comp].isCarrierGas ? initialPressure[grid] : 0.0;
      moleFraction[grid * numberOfComponents + comp] = components[comp].isCarrierGas ? 1.0 : 0.0;
    }
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    partialPressure[0 * numberOfComponents + comp] = initialPressure[0] * components[comp].initialGasMoleFraction;
    moleFraction[0 * numberOfComponents + comp] = components[comp].initialGasMoleFraction;
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    totalPressure[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      totalPressure[grid] += std::max(0.0, partialPressure[grid * numberOfComponents + comp]);
    }
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentration[grid * numberOfComponents + comp] = partialPressure[grid * numberOfComponents + comp] / (R * externalTemperature);
    }
  }

  const double gasDensityTemperature = std::abs(influxTemperature) < 1e-10 ? externalTemperature : influxTemperature;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double molecularWeight = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      molecularWeight += moleFraction[grid * numberOfComponents + comp] * components[comp].molecularWeight;
    }

    gasDensity[grid] = molecularWeight > 0.0 && gasDensityTemperature > 0.0
                           ? totalPressure[grid] * molecularWeight / (R * gasDensityTemperature)
                           : 0.0;
  }

  if (std::abs(influxTemperature) < 1e-10) influxTemperature = externalTemperature;
  std::fill(gasTemperature.begin(), gasTemperature.end(), influxTemperature);
  std::fill(solidTemperature.begin(), solidTemperature.end(), influxTemperature);
  std::fill(wallTemperature.begin(), wallTemperature.end(), influxTemperature);

  for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
  {
    double sum = 0.0;
    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      idealGasMolFractions[j] = std::max(partialPressure[i * numberOfComponents + j] / initialPressure[i], 0.0);
      sum += idealGasMolFractions[j];
    }

    if (sum <= 0.0)
    {
      throw std::runtime_error("Error: initialized gas mol-fraction sum must be positive");
    }

    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      idealGasMolFractions[j] /= sum;
    }

    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      equilibriumAdsorption[i * numberOfComponents + j] = 0.0;
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (!hasAdsorbentOfType[i * numberOfAdsorbents + ads]) continue;

      iastPerformance += mixture[ads].predictMixture(
          idealGasMolFractions, initialPressure[i], adsorbedMolFractions, numberOfMolecules,
          std::span<double>(&cachedPressure[(i * numberOfAdsorbents + ads) * numberOfComponents * maxIsothermTerms],
                            numberOfComponents * maxIsothermTerms),
          std::span<double>(&cachedGrandPotential[(i * numberOfAdsorbents + ads) * maxIsothermTerms],
                            maxIsothermTerms),
          gasTemperature[i]);

      for (size_t j = 0; j < numberOfComponents; ++j)
      {
        equilibriumAdsorption[i * numberOfComponents + j] +=
            fractionOfAdsorbent[i * numberOfAdsorbents + ads] * numberOfMolecules[j];
      }
    }
  }
}

void ColumnMultibed::setTemperature(double temperature)
{
  influxTemperature = temperature;
  externalTemperature = temperature;
  if (!energyBalance)
  {
    std::fill(gasTemperature.begin(), gasTemperature.end(), temperature);
    std::fill(solidTemperature.begin(), solidTemperature.end(), temperature);
    std::fill(wallTemperature.begin(), wallTemperature.end(), temperature);
  }
}

void ColumnMultibed::writeOutputHeader(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream) const
{
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    std::print(componentStreams[i], "# component index: {}\n", i);
    std::print(componentStreams[i], "# component name: {}\n", components[i].name);
    std::print(componentStreams[i], "# column 1: Dimensionless time, τ = tv/L [-]\n");
    std::print(componentStreams[i], "# column 2: Time, t [min]\n");
    std::print(componentStreams[i], "# column 3: Column position, z [m]\n");
    std::print(componentStreams[i], "# column 4: Concentration, c_i [mol/m^3]\n");
    std::print(componentStreams[i], "# column 5: Concentration time derivative, dc_i/dt [mol/m^3/s]\n");
    std::print(componentStreams[i], "# column 6: Adsorption, q_i [mol/kg]\n");
    std::print(componentStreams[i], "# column 7: Adsorption time derivative, dq_i/dt [mol/kg/s]\n");
    std::print(componentStreams[i], "# column 8: Partial pressure, p_i [Pa]\n");
    std::print(componentStreams[i], "# column 9: Equilibrium physisorption, q_i^* [mol/kg]\n");
    std::print(componentStreams[i], "# column 10: Normalized partial pressure, p_i / (p_t y_i,0) [-]\n");
  }

  std::print(columnStream, "# column 1: Dimensionless time, τ = tv/L [-]\n");
  std::print(columnStream, "# column 2: Time, t [min]\n");
  std::print(columnStream, "# column 3: Column position, z [m]\n");
  std::print(columnStream, "# column 4: Interstitial gas velocity, v [m/s]\n");
  std::print(columnStream, "# column 5: Total pressure, p_t [Pa]\n");
  std::print(columnStream, "# column 6: Gas temperature, T_g [K]\n");
  std::print(columnStream, "# column 7: Gas temperature time derivative, dT_g/dt [K/s]\n");
  std::print(columnStream, "# column 8: Solid temperature, T_s [K]\n");
  std::print(columnStream, "# column 9: Solid temperature time derivative, dT_s/dt [K/s]\n");
  std::print(columnStream, "# column 10: Wall temperature, T_w [K]\n");
  std::print(columnStream, "# column 11: Wall temperature time derivative, dT_w/dt [K/s]\n");
  std::print(columnStream, "# column 12: Gas density, rho_g [kg/m^3]\n");
}

void ColumnMultibed::writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream, double time) const
{
  // column.data
  // column 1: dimensionless time [-]
  // column 2: time [min]
  // column 3: column position [m]
  // column 4: interstitial gas velocity [m/s]
  // column 5: total pressure [Pa]
  // column 6: gas temperature [K]
  // column 7: gas temperature time derivative [K/s]
  // column 8: solid temperature [K]
  // column 9: solid temperature time derivative [K/s]
  // column 10: wall temperature [K]
  // column 11: wall temperature time derivative [K/s]
  // column 12: gas density [kg/m^3]
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    std::print(columnStream, "{} {} {} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor, time / 60.0,
               columnDistances[grid], interstitialGasVelocity[grid], totalPressure[grid],
               gasTemperature[grid], gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid],
               wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid]);
  }
  std::print(columnStream, "\n\n");

  // Per-component files
  // column 1: dimensionless time [-]
  // column 2: time [min]
  // column 3: column position [m]
  // column 4: concentration [mol/m^3]
  // column 5: concentration time derivative [mol/m^3/s]
  // column 6: physisorption [mol/kg]
  // column 7: physisorption time derivative [mol/kg/s]
  // column 8: partial pressure [Pa]
  // column 9: equilibrium physisorption [mol/kg]
  // column 10: normalized partial pressure [-]
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      const size_t index = grid * numberOfComponents + comp;

      double normalizedPressure = 0.0;
      if (components[comp].initialGasMoleFraction > 0.0 && totalPressure[grid] > 0.0)
      {
        normalizedPressure = partialPressure[index] / (totalPressure[grid] * components[comp].initialGasMoleFraction);
      }

      std::print(componentStreams[comp], "{} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor, time / 60.0,
                 columnDistances[grid], concentration[index], concentrationDot[index],
                 physisorption[index], physisorptionDot[index], partialPressure[index], equilibriumAdsorption[index],
                 normalizedPressure);
    }
    std::print(componentStreams[comp], "\n\n");
  }
}

std::string ColumnMultibed::repr() const
{
  return std::format(
      "Temperature:                           {} [K]\n"
      "Column length:                         {} [m]\n"
      "Adsorbent lengths:                     {} [m]\n"
      "Adsorbent interface lengths:           {} [m]\n"
      "Adsorbent void-fractions:              {} [-]\n"
      "Particle densities:                    {} [kg/m^3]\n"
      "Inlet pressure:                        {} [Pa]\n"
      "Outlet pressure:                       {} [Pa]\n"
      "Pressure gradient:                     {} [Pa]\n"
      "Column entrance interstitial velocity: {} [m/s]\n"
      "\n\n",
      externalTemperature, columnLength, vectorToString(adsorbentLengths), vectorToString(adsorbentInterfaceLengths),
      vectorToString(adsorbentVoidFractions), vectorToString(particleDensities), inletPressure, outletPressure,
      pressureGradient, columnEntranceVelocity);
}

void ColumnMultibed::writeJSON(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("ColumnMultibed::writeJSON: cannot open file '" + filename + "'");

  nlohmann::json j;
  j["numberOfGridPoints"] = numberOfGridPoints;
  j["numberOfComponents"] = numberOfComponents;
  j["maxIsothermTerms"] = maxIsothermTerms;
  j["columnDistances"] = columnDistances;

  j["prefactorMassTransfer"] = prefactorMassTransfer;
  j["idealGasMolFractions"] = idealGasMolFractions;
  j["adsorbedMolFractions"] = adsorbedMolFractions;
  j["numberOfMolecules"] = numberOfMolecules;

  j["interstitialGasVelocity"] = interstitialGasVelocity;
  j["gasDensity"] = gasDensity;
  j["totalVoidFraction"] = totalVoidFraction;
  j["particleDensity"] = particleDensity;
  j["totalPressure"] = totalPressure;
  j["gasTemperature"] = toVector(gasTemperature);
  j["gasTemperatureDot"] = toVector(gasTemperatureDot);
  j["solidTemperature"] = toVector(solidTemperature);
  j["solidTemperatureDot"] = toVector(solidTemperatureDot);
  j["wallTemperature"] = toVector(wallTemperature);
  j["wallTemperatureDot"] = toVector(wallTemperatureDot);

  j["concentration"] = toVector(concentration);
  j["concentrationDot"] = toVector(concentrationDot);
  j["physisorption"] = toVector(physisorption);
  j["physisorptionDot"] = toVector(physisorptionDot);
  j["partialPressure"] = partialPressure;
  j["equilibriumAdsorption"] = equilibriumAdsorption;
  j["moleFraction"] = moleFraction;
  j["fractionOfAdsorbent"] = fractionOfAdsorbent;
  j["adsorbentScaledVoidFraction"] = adsorbentScaledVoidFraction;

  j["cachedPressure"] = cachedPressure;
  j["cachedGrandPotential"] = cachedGrandPotential;

  out << j.dump(4) << "\n";
}

void ColumnMultibed::readJSON(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("ColumnMultibed::readJSONFile: cannot open file '" + filename + "'");

  nlohmann::json j;
  in >> j;

  auto requireSizeT = [&](const char* key) -> size_t
  {
    if (!j.contains(key)) throw std::runtime_error(std::string("ColumnMultibed::readJSON: missing required key '") + key + "'");
    return j.at(key).get<size_t>();
  };

  const size_t fileNgrid = requireSizeT("numberOfGridPoints");
  const size_t fileNcomp = requireSizeT("numberOfComponents");
  const size_t fileMaxIsothermTerms = requireSizeT("maxIsothermTerms");

  if (fileNgrid != numberOfGridPoints)
    throw std::runtime_error("ColumnMultibed::readJSON: numberOfGridPoints mismatch (file " + std::to_string(fileNgrid) + ", column " +
                             std::to_string(numberOfGridPoints) + ")");
  if (fileNcomp != numberOfComponents)
    throw std::runtime_error("ColumnMultibed::readJSON: numberOfComponents mismatch (file " + std::to_string(fileNcomp) + ", column " +
                             std::to_string(numberOfComponents) + ")");
  if (fileMaxIsothermTerms != maxIsothermTerms)
    throw std::runtime_error("ColumnMultibed::readJSON: maxIsothermTerms mismatch (file " +
                             std::to_string(fileMaxIsothermTerms) + ", column " + std::to_string(maxIsothermTerms) +
                             ")");

  auto loadVectorChecked = [&](const char* key, std::vector<double>& dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array()) throw std::runtime_error(std::string("ColumnMultibed::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("ColumnMultibed::readJSON: size mismatch for '" + std::string(key) + "'");

    dst = a.get<std::vector<double>>();
  };

  auto loadSpanChecked = [&](const char* key, std::span<double> dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array()) throw std::runtime_error(std::string("ColumnMultibed::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("ColumnMultibed::readJSON: size mismatch for '" + std::string(key) + "'");

    const std::vector<double> tmp = a.get<std::vector<double>>();
    std::copy(tmp.begin(), tmp.end(), dst.begin());
  };

  loadVectorChecked("prefactorMassTransfer", prefactorMassTransfer);
  loadVectorChecked("idealGasMolFractions", idealGasMolFractions);
  loadVectorChecked("adsorbedMolFractions", adsorbedMolFractions);
  loadVectorChecked("numberOfMolecules", numberOfMolecules);

  loadVectorChecked("interstitialGasVelocity", interstitialGasVelocity);
  loadVectorChecked("gasDensity", gasDensity);
  loadVectorChecked("totalVoidFraction", totalVoidFraction);
  loadVectorChecked("particleDensity", particleDensity);
  loadVectorChecked("totalPressure", totalPressure);

  loadSpanChecked("gasTemperature", gasTemperature);
  loadSpanChecked("gasTemperatureDot", gasTemperatureDot);
  loadSpanChecked("solidTemperature", solidTemperature);
  loadSpanChecked("solidTemperatureDot", solidTemperatureDot);
  loadSpanChecked("wallTemperature", wallTemperature);
  loadSpanChecked("wallTemperatureDot", wallTemperatureDot);

  loadSpanChecked("concentration", concentration);
  loadSpanChecked("concentrationDot", concentrationDot);
  loadSpanChecked("physisorption", physisorption);
  loadSpanChecked("physisorptionDot", physisorptionDot);

  loadVectorChecked("partialPressure", partialPressure);
  loadVectorChecked("equilibriumAdsorption", equilibriumAdsorption);
  loadVectorChecked("moleFraction", moleFraction);
  loadVectorChecked("fractionOfAdsorbent", fractionOfAdsorbent);
  loadVectorChecked("adsorbentScaledVoidFraction", adsorbentScaledVoidFraction);

  loadVectorChecked("cachedPressure", cachedPressure);
  loadVectorChecked("cachedGrandPotential", cachedGrandPotential);
  loadVectorChecked("columnDistances", columnDistances);
  validateColumnDistances(columnDistances, numberOfGridPoints, columnLength);
}
