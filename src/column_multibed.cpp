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

#include "breakthrough.h"
#include "column.h"
#include "component.h"
#include "inputreader.h"
#include "json.h"
#include "mixture_prediction.h"
#include "utils.h"

namespace
{
std::vector<double> toVector(std::span<const double> s) { return std::vector<double>(s.begin(), s.end()); }
}  // namespace

size_t Column::stateSize() const noexcept { return (2 * numberOfComponents + 3) * (numberOfGridPoints + 1); }

void Column::bindStateViews() noexcept
{
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;

  double* base = state.data();
  double* baseDot = stateDot.data();

  moleFraction = std::span<double>(base + 0 * big, big);
  adsorption = std::span<double>(base + 1 * big, big);
  gasTemperature = std::span<double>(base + 2 * big + 0 * small, small);
  solidTemperature = std::span<double>(base + 2 * big + 1 * small, small);
  wallTemperature = std::span<double>(base + 2 * big + 2 * small, small);

  moleFractionDot = std::span<double>(baseDot + 0 * big, big);
  adsorptionDot = std::span<double>(baseDot + 1 * big, big);
  gasTemperatureDot = std::span<double>(baseDot + 2 * big + 0 * small, small);
  solidTemperatureDot = std::span<double>(baseDot + 2 * big + 1 * small, small);
  wallTemperatureDot = std::span<double>(baseDot + 2 * big + 2 * small, small);
}

Column::Column(const Column& other)
    : mixture(other.mixture),
      components(other.components),
      velocityProfile(other.velocityProfile),
      boundaryCondition(other.boundaryCondition),
      energyBalance(other.energyBalance),
      numberOfGridPoints(other.numberOfGridPoints),
      numberOfComponents(other.numberOfComponents),
      maxIsothermTerms(other.maxIsothermTerms),
      numberOfCalls(other.numberOfCalls),
      carrierGasComponent(other.carrierGasComponent),
      externalTemperature(other.externalTemperature),
      inletPressure(other.inletPressure),
      outletPressure(other.outletPressure),
      pressureGradient(other.pressureGradient),
      voidFraction(other.voidFraction),
      particleDensity(other.particleDensity),
      columnEntranceVelocity(other.columnEntranceVelocity),
      columnLength(other.columnLength),
      dynamicViscosity(other.dynamicViscosity),
      particleDiameter(other.particleDiameter),
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
      concentration(other.concentration),
      partialPressure(other.partialPressure),
      equilibriumAdsorption(other.equilibriumAdsorption),
      cachedPressure(other.cachedPressure),
      cachedGrandPotential(other.cachedGrandPotential),
      coeffGasGas(other.coeffGasGas),
      coeffGasSolid(other.coeffGasSolid),
      coeffGasWall(other.coeffGasWall),
      coeffDiffusion(other.coeffDiffusion),
      facePressures(other.facePressures),
      massFlux(other.massFlux),
      state(other.state),
      stateDot(other.stateDot)
{
  bindStateViews();
}

Column& Column::operator=(const Column& other)
{
  if (this == &other) return *this;

  mixture = other.mixture;
  components = other.components;
  velocityProfile = other.velocityProfile;
  boundaryCondition = other.boundaryCondition;
  energyBalance = other.energyBalance;
  numberOfGridPoints = other.numberOfGridPoints;
  numberOfComponents = other.numberOfComponents;
  maxIsothermTerms = other.maxIsothermTerms;
  numberOfCalls = other.numberOfCalls;
  carrierGasComponent = other.carrierGasComponent;
  externalTemperature = other.externalTemperature;
  inletPressure = other.inletPressure;
  outletPressure = other.outletPressure;
  pressureGradient = other.pressureGradient;
  voidFraction = other.voidFraction;
  particleDensity = other.particleDensity;
  columnEntranceVelocity = other.columnEntranceVelocity;
  columnLength = other.columnLength;
  dynamicViscosity = other.dynamicViscosity;
  particleDiameter = other.particleDiameter;
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
  concentration = other.concentration;
  partialPressure = other.partialPressure;
  equilibriumAdsorption = other.equilibriumAdsorption;
  cachedPressure = other.cachedPressure;
  cachedGrandPotential = other.cachedGrandPotential;
  coeffGasGas = other.coeffGasGas;
  coeffGasSolid = other.coeffGasSolid;
  coeffGasWall = other.coeffGasWall;
  coeffDiffusion = other.coeffDiffusion;
  facePressures = other.facePressures;
  massFlux = other.massFlux;
  state = other.state;
  stateDot = other.stateDot;

  bindStateViews();
  return *this;
}

void Column::initialize()
{
  for (size_t j = 0; j < numberOfComponents; ++j)
  {
    prefactorMassTransfer[j] =
        ((1.0 - voidFraction) / voidFraction) * particleDensity * components[j].massTransferCoefficient;
  }

  std::fill(partialPressure.begin(), partialPressure.end(), 0.0);
  std::fill(adsorption.begin(), adsorption.end(), 0.0);
  std::fill(moleFraction.begin(), moleFraction.end(), 0.0);
  std::fill(concentration.begin(), concentration.end(), 0.0);
  std::fill(stateDot.begin(), stateDot.end(), 0.0);

  std::vector<double> initialPressure(numberOfGridPoints + 1, 0.0);

  auto gridRatio = [&](size_t i) -> double
  { return numberOfGridPoints == 0 ? 0.0 : static_cast<double>(i) / static_cast<double>(numberOfGridPoints); };

  auto fillPressure = [&](double pressure) { std::fill(initialPressure.begin(), initialPressure.end(), pressure); };

  auto fillVelocity = [&](double velocity)
  { std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), velocity); };

  bool setErgunVelocitiesFromPressureProfile = false;

  switch (velocityProfile)
  {
    case VelocityProfile::FixedVelocity:
    {
      switch (boundaryCondition)
      {
        case BoundaryCondition::InletPressureInletVelocity:
          fillPressure(inletPressure);
          break;
        case BoundaryCondition::InletPressureOutletPressure:
          fillPressure(outletPressure);
          break;
        case BoundaryCondition::InletVelocityOutletPressure:
          fillPressure(outletPressure);
          break;
      }

      fillVelocity(columnEntranceVelocity);
      break;
    }

    case VelocityProfile::FixedPressureGradient:
    {
      if (boundaryCondition != BoundaryCondition::InletPressureInletVelocity)
      {
        throw std::runtime_error("Error: FixedPressureGradient is only allowed with InletPressureInletVelocity");
      }

      for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
      {
        initialPressure[i] = inletPressure + gridRatio(i) * pressureGradient;
      }

      fillVelocity(columnEntranceVelocity);
      break;
    }

    case VelocityProfile::Ergun:
    {
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

          setErgunVelocitiesFromPressureProfile = true;
          break;
        }

        case BoundaryCondition::InletVelocityOutletPressure:
        {
          fillPressure(outletPressure);
          fillVelocity(columnEntranceVelocity);
          break;
        }
      }

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
      totalConcentration[grid] = totalPressure[grid] / (R * externalTemperature);
      concentration[grid * numberOfComponents + comp] =
          moleFraction[grid * numberOfComponents + comp] * totalConcentration[grid];
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

  if (setErgunVelocitiesFromPressureProfile)
  {
    if (columnLength <= 0.0)
    {
      throw std::runtime_error("Error: Ergun velocity calculation requires ColumnLength > 0");
    }

    if (particleDiameter <= 0.0)
    {
      throw std::runtime_error("Error: Ergun velocity calculation requires ParticleDiameter > 0");
    }

    if (voidFraction <= 0.0)
    {
      throw std::runtime_error("Error: Ergun velocity calculation requires ColumnVoidFraction > 0");
    }

    const double eps3 = voidFraction * voidFraction * voidFraction;
    const double viscousCoefficient = 150.0 * dynamicViscosity * (1.0 - voidFraction) * (1.0 - voidFraction) /
                                      (particleDiameter * particleDiameter * eps3);

    auto pressureDerivative = [&](size_t grid) -> double
    {
      if (grid == 0)
      {
        return (initialPressure[1] - initialPressure[0]) / resolution;
      }

      if (grid == numberOfGridPoints)
      {
        return (initialPressure[numberOfGridPoints] - initialPressure[numberOfGridPoints - 1]) / resolution;
      }

      return (initialPressure[grid + 1] - initialPressure[grid - 1]) / (2.0 * resolution);
    };

    auto ergunInterstitialVelocity = [&](double dPdz, double density) -> double
    {
      const double drivingForce = -dPdz;

      if (std::abs(drivingForce) < 1e-30)
      {
        return 0.0;
      }

      const double sign = drivingForce < 0.0 ? -1.0 : 1.0;
      const double rhs = std::abs(drivingForce);

      const double inertialCoefficient =
          density > 0.0 ? 1.75 * density * (1.0 - voidFraction) / (particleDiameter * eps3) : 0.0;

      double superficialVelocity = 0.0;

      if (inertialCoefficient > 0.0)
      {
        superficialVelocity = (-viscousCoefficient +
                               std::sqrt(viscousCoefficient * viscousCoefficient + 4.0 * inertialCoefficient * rhs)) /
                              (2.0 * inertialCoefficient);
      }
      else if (viscousCoefficient > 0.0)
      {
        superficialVelocity = rhs / viscousCoefficient;
      }
      else
      {
        throw std::runtime_error("Error: invalid Ergun coefficients");
      }

      return sign * superficialVelocity / voidFraction;
    };

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      interstitialGasVelocity[grid] = ergunInterstitialVelocity(pressureDerivative(grid), gasDensity[grid]);
    }
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

    iastPerformance +=
        mixture.predictMixture(idealGasMolFractions, initialPressure[i], adsorbedMolFractions, numberOfMolecules,
                               &cachedPressure[i * numberOfComponents * maxIsothermTerms],
                               &cachedGrandPotential[i * maxIsothermTerms], gasTemperature[i]);

    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      equilibriumAdsorption[i * numberOfComponents + j] = numberOfMolecules[j];
    }
  }
}

void Column::setTemperature(double temperature)
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

void Column::writeOutputHeader(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream) const
{
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    std::print(componentStreams[i], "# component index: {}\n", i);
    std::print(componentStreams[i], "# component name: {}\n", components[i].name);
    std::print(componentStreams[i], "# column 1: Dimensionless time, τ = tv/L [-]\n");
    std::print(componentStreams[i], "# column 2: Time, t [min]\n");
    std::print(componentStreams[i], "# column 3: Column position, z [m]\n");
    std::print(componentStreams[i], "# column 4: Mole fraction, y_i [-]\n");
    std::print(componentStreams[i], "# column 5: Mole fraction time derivative, dy_i/dt [1/s]\n");
    std::print(componentStreams[i], "# column 6: Concentration, c_i [mol/m^3]\n");
    std::print(componentStreams[i], "# column 7: Adsorption, q_i [mol/kg]\n");
    std::print(componentStreams[i], "# column 8: Adsorption time derivative, dq_i/dt [mol/kg/s]\n");
    std::print(componentStreams[i], "# column 9: Partial pressure, p_i [Pa]\n");
    std::print(componentStreams[i], "# column 10: Equilibrium adsorption, q_i^* [mol/kg]\n");
    std::print(componentStreams[i], "# column 11: Normalized partial pressure, p_i / (p_t y_i,0) [-]\n");
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

void Column::writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream, double time) const
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
               static_cast<double>(grid) * resolution, interstitialGasVelocity[grid], totalPressure[grid],
               gasTemperature[grid], gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid],
               wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid]);
  }
  std::print(columnStream, "\n\n");

  // Per-component files
  // column 1: dimensionless time [-]
  // column 2: time [min]
  // column 3: column position [m]
  // column 4: mole fraction [-]
  // column 5: mole fraction time derivative [1/s]
  // column 6: concentration [mol/m^3]
  // column 7: adsorption [mol/kg]
  // column 8: adsorption time derivative [mol/kg/s]
  // column 9: partial pressure [Pa]
  // column 10: equilibrium adsorption [mol/kg]
  // column 11: normalized partial pressure [-]
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

      std::print(componentStreams[comp], "{} {} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor,
                 time / 60.0, static_cast<double>(grid) * resolution, moleFraction[index], moleFractionDot[index],
                 concentration[index], adsorption[index], adsorptionDot[index], partialPressure[index],
                 equilibriumAdsorption[index], normalizedPressure);
    }
    std::print(componentStreams[comp], "\n\n");
  }
}

std::string Column::repr() const
{
  return std::format(
      "Temperature:                           {} [K]\n"
      "Column length:                         {} [m]\n"
      "Column void-fraction:                  {} [-]\n"
      "Particle density:                      {} [kg/m^3]\n"
      "Inlet pressure:                        {} [Pa]\n"
      "Outlet pressure:                       {} [Pa]\n"
      "Pressure gradient:                     {} [Pa]\n"
      "Column entrance interstitial velocity: {} [m/s]\n"
      "\n\n",
      externalTemperature, columnLength, voidFraction, particleDensity, inletPressure, outletPressure, pressureGradient,
      columnEntranceVelocity);
}

void Column::writeJSON(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Column::writeJSON: cannot open file '" + filename + "'");

  nlohmann::json j;
  j["numberOfGridPoints"] = numberOfGridPoints;
  j["numberOfComponents"] = numberOfComponents;
  j["maxIsothermTerms"] = maxIsothermTerms;

  j["prefactorMassTransfer"] = prefactorMassTransfer;
  j["idealGasMolFractions"] = idealGasMolFractions;
  j["adsorbedMolFractions"] = adsorbedMolFractions;
  j["numberOfMolecules"] = numberOfMolecules;

  j["interstitialGasVelocity"] = interstitialGasVelocity;
  j["gasDensity"] = gasDensity;
  j["totalConcentration"] = totalConcentration;
  j["totalPressure"] = totalPressure;
  j["gasTemperature"] = toVector(gasTemperature);
  j["gasTemperatureDot"] = toVector(gasTemperatureDot);
  j["solidTemperature"] = toVector(solidTemperature);
  j["solidTemperatureDot"] = toVector(solidTemperatureDot);
  j["wallTemperature"] = toVector(wallTemperature);
  j["wallTemperatureDot"] = toVector(wallTemperatureDot);

  j["concentration"] = concentration;
  j["adsorption"] = toVector(adsorption);
  j["adsorptionDot"] = toVector(adsorptionDot);
  j["partialPressure"] = partialPressure;
  j["equilibriumAdsorption"] = equilibriumAdsorption;
  j["moleFraction"] = toVector(moleFraction);
  j["moleFractionDot"] = toVector(moleFractionDot);

  j["cachedPressure"] = cachedPressure;
  j["cachedGrandPotential"] = cachedGrandPotential;

  std::print(out, "{}\n", j.dump(4));
}

void Column::readJSON(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Column::readJSONFile: cannot open file '" + filename + "'");

  nlohmann::json j;
  in >> j;

  auto requireSizeT = [&](const char* key) -> size_t
  {
    if (!j.contains(key)) throw std::runtime_error(std::string("Column::readJSON: missing required key '") + key + "'");
    return j.at(key).get<size_t>();
  };

  const size_t fileNgrid = requireSizeT("numberOfGridPoints");
  const size_t fileNcomp = requireSizeT("numberOfComponents");
  const size_t fileMaxIsothermTerms = requireSizeT("maxIsothermTerms");

  if (fileNgrid != numberOfGridPoints)
    throw std::runtime_error("Column::readJSON: numberOfGridPoints mismatch (file " + std::to_string(fileNgrid) +
                             ", column " + std::to_string(numberOfGridPoints) + ")");
  if (fileNcomp != numberOfComponents)
    throw std::runtime_error("Column::readJSON: numberOfComponents mismatch (file " + std::to_string(fileNcomp) +
                             ", column " + std::to_string(numberOfComponents) + ")");
  if (fileMaxIsothermTerms != maxIsothermTerms)
    throw std::runtime_error("Column::readJSON: maxIsothermTerms mismatch (file " +
                             std::to_string(fileMaxIsothermTerms) + ", column " + std::to_string(maxIsothermTerms) +
                             ")");

  auto loadVectorChecked = [&](const char* key, std::vector<double>& dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array()) throw std::runtime_error(std::string("Column::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("Column::readJSON: size mismatch for '" + std::string(key) + "'");

    dst = a.get<std::vector<double>>();
  };

  auto loadSpanChecked = [&](const char* key, std::span<double> dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array()) throw std::runtime_error(std::string("Column::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("Column::readJSON: size mismatch for '" + std::string(key) + "'");

    const std::vector<double> tmp = a.get<std::vector<double>>();
    std::copy(tmp.begin(), tmp.end(), dst.begin());
  };

  loadVectorChecked("prefactorMassTransfer", prefactorMassTransfer);
  loadVectorChecked("idealGasMolFractions", idealGasMolFractions);
  loadVectorChecked("adsorbedMolFractions", adsorbedMolFractions);
  loadVectorChecked("numberOfMolecules", numberOfMolecules);

  loadVectorChecked("interstitialGasVelocity", interstitialGasVelocity);
  loadVectorChecked("gasDensity", gasDensity);
  loadVectorChecked("totalConcentration", totalConcentration);
  loadVectorChecked("totalPressure", totalPressure);

  loadSpanChecked("gasTemperature", gasTemperature);
  loadSpanChecked("gasTemperatureDot", gasTemperatureDot);
  loadSpanChecked("solidTemperature", solidTemperature);
  loadSpanChecked("solidTemperatureDot", solidTemperatureDot);
  loadSpanChecked("wallTemperature", wallTemperature);
  loadSpanChecked("wallTemperatureDot", wallTemperatureDot);

  loadVectorChecked("concentration", concentration);
  loadSpanChecked("adsorption", adsorption);
  loadSpanChecked("adsorptionDot", adsorptionDot);

  loadVectorChecked("partialPressure", partialPressure);
  loadVectorChecked("equilibriumAdsorption", equilibriumAdsorption);
  loadSpanChecked("moleFraction", moleFraction);
  loadSpanChecked("moleFractionDot", moleFractionDot);

  loadVectorChecked("cachedPressure", cachedPressure);
  loadVectorChecked("cachedGrandPotential", cachedGrandPotential);
}
