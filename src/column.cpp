#include "column.h"

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
#include "component.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "integrators/sorption.h"
#include "integrators/transport.h"
#include "json.h"
#include "mixture_prediction.h"
#include "utils.h"

namespace
{
std::vector<double> toVector(std::span<const double> s) { return std::vector<double>(s.begin(), s.end()); }

double sumChemisorptionSites(std::span<const double> values, size_t componentBlockSize, size_t numberOfSites,
                             size_t index)
{
  double sum = 0.0;
  for (size_t site = 0; site < numberOfSites; ++site)
  {
    sum += values[site * componentBlockSize + index];
  }
  return sum;
}
}  // namespace

size_t Column::stateSize() const noexcept { return stateLayout().stateSize(); }

bool Column::requiresSurfacePoreTransport(const std::vector<Component>& components) noexcept
{
  return std::any_of(components.begin(), components.end(),
                     [](const Component& component) { return component.chemisorption.usesSurfacePoreTransport(); });
}

size_t Column::maximumChemisorptionSites(const std::vector<Component>& components) noexcept
{
  size_t maximum = 1;
  for (const Component& component : components)
  {
    maximum = std::max(maximum, component.chemisorption.numberOfSites);
  }
  return maximum;
}

MixturePrediction Column::makeChemisorptionMixture(const MixturePrediction& physisorptionMixture,
                                                   const std::vector<Component>& components)
{
  return MixturePrediction::makeChemisorptionPrediction(physisorptionMixture, components);
}

ColumnStateLayout Column::stateLayout() const noexcept
{
  return ColumnStateLayout{numberOfGridPoints, numberOfComponents, maxChemisorptionSites, surfacePoreTransportEnabled};
}

void Column::bindStateViews() noexcept
{
  const ColumnStateLayout layout = stateLayout();
  double* base = state.data();
  double* baseDot = stateDot.data();

  concentration = layout.concentration(base);
  physisorption = layout.physisorption(base);
  chemisorption = layout.chemisorption(base);
  surfaceConcentration = layout.surfaceConcentration(base);
  poreConcentration = layout.poreConcentration(base);
  gasTemperature = layout.gasTemperature(base);
  solidTemperature = layout.solidTemperature(base);
  wallTemperature = layout.wallTemperature(base);

  concentrationDot = layout.concentration(baseDot);
  physisorptionDot = layout.physisorption(baseDot);
  chemisorptionDot = layout.chemisorption(baseDot);
  surfaceConcentrationDot = layout.surfaceConcentration(baseDot);
  poreConcentrationDot = layout.poreConcentration(baseDot);
  gasTemperatureDot = layout.gasTemperature(baseDot);
  solidTemperatureDot = layout.solidTemperature(baseDot);
  wallTemperatureDot = layout.wallTemperature(baseDot);
}

Column::Column(const Column& other)
    : physisorptionMixture(other.physisorptionMixture),
      chemisorptionMixture(other.chemisorptionMixture),
      components(other.components),
      boundaryCondition(other.boundaryCondition),
      fluidPhase(other.fluidPhase),
      pHMode(other.pHMode),
      energyBalance(other.energyBalance),
      geometry(other.geometry),
      reactions(other.reactions),
      numberOfGridPoints(other.numberOfGridPoints),
      numberOfComponents(other.numberOfComponents),
      maxIsothermTerms(other.maxIsothermTerms),
      maxChemisorptionSites(other.maxChemisorptionSites),
      numberOfCalls(other.numberOfCalls),
      carrierGasComponent(other.carrierGasComponent),
      pHComponent(other.pHComponent),
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
      liquidDensity(other.liquidDensity),
      pHValue(other.pHValue),
      pKw(other.pKw),
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
      surfacePoreTransportEnabled(other.surfacePoreTransportEnabled),
      iastPerformance(other.iastPerformance),
      prefactorMassTransfer(other.prefactorMassTransfer),
      idealGasMolFractions(other.idealGasMolFractions),
      adsorbedMolFractions(other.adsorbedMolFractions),
      numberOfMolecules(other.numberOfMolecules),
      interstitialGasVelocity(other.interstitialGasVelocity),
      gasDensity(other.gasDensity),
      totalConcentration(other.totalConcentration),
      totalPressure(other.totalPressure),
      pH(other.pH),
      moleFraction(other.moleFraction),
      partialPressure(other.partialPressure),
      equilibriumPhysisorption(other.equilibriumPhysisorption),
      equilibriumChemisorption(other.equilibriumChemisorption),
      cachedPressure(other.cachedPressure),
      cachedGrandPotential(other.cachedGrandPotential),
      cachedChemisorptionPressure(other.cachedChemisorptionPressure),
      cachedChemisorptionGrandPotential(other.cachedChemisorptionGrandPotential),
      coeffDiffusion(other.coeffDiffusion),
      facePressures(other.facePressures),
      massFlux(other.massFlux),
      bulkSpeciesSink(other.bulkSpeciesSink),
      reactionPhysisorptionSource(other.reactionPhysisorptionSource),
      reactionChemisorptionSource(other.reactionChemisorptionSource),
      reactionPoreConcentrationSource(other.reactionPoreConcentrationSource),
      reactionHeat(other.reactionHeat),
      state(other.state),
      stateDot(other.stateDot)
{
  bindStateViews();
}

Column& Column::operator=(const Column& other)
{
  if (this == &other) return *this;

  physisorptionMixture = other.physisorptionMixture;
  chemisorptionMixture = other.chemisorptionMixture;
  components = other.components;
  boundaryCondition = other.boundaryCondition;
  fluidPhase = other.fluidPhase;
  pHMode = other.pHMode;
  energyBalance = other.energyBalance;
  geometry = other.geometry;
  reactions = other.reactions;
  numberOfGridPoints = other.numberOfGridPoints;
  numberOfComponents = other.numberOfComponents;
  maxIsothermTerms = other.maxIsothermTerms;
  maxChemisorptionSites = other.maxChemisorptionSites;
  numberOfCalls = other.numberOfCalls;
  carrierGasComponent = other.carrierGasComponent;
  pHComponent = other.pHComponent;
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
  liquidDensity = other.liquidDensity;
  pHValue = other.pHValue;
  pKw = other.pKw;
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
  surfacePoreTransportEnabled = other.surfacePoreTransportEnabled;
  iastPerformance = other.iastPerformance;
  prefactorMassTransfer = other.prefactorMassTransfer;
  idealGasMolFractions = other.idealGasMolFractions;
  adsorbedMolFractions = other.adsorbedMolFractions;
  numberOfMolecules = other.numberOfMolecules;
  interstitialGasVelocity = other.interstitialGasVelocity;
  gasDensity = other.gasDensity;
  totalConcentration = other.totalConcentration;
  totalPressure = other.totalPressure;
  pH = other.pH;
  moleFraction = other.moleFraction;
  partialPressure = other.partialPressure;
  equilibriumPhysisorption = other.equilibriumPhysisorption;
  equilibriumChemisorption = other.equilibriumChemisorption;
  cachedPressure = other.cachedPressure;
  cachedGrandPotential = other.cachedGrandPotential;
  cachedChemisorptionPressure = other.cachedChemisorptionPressure;
  cachedChemisorptionGrandPotential = other.cachedChemisorptionGrandPotential;
  coeffDiffusion = other.coeffDiffusion;
  facePressures = other.facePressures;
  massFlux = other.massFlux;
  bulkSpeciesSink = other.bulkSpeciesSink;
  reactionPhysisorptionSource = other.reactionPhysisorptionSource;
  reactionChemisorptionSource = other.reactionChemisorptionSource;
  reactionPoreConcentrationSource = other.reactionPoreConcentrationSource;
  reactionHeat = other.reactionHeat;
  state = other.state;
  stateDot = other.stateDot;

  bindStateViews();
  return *this;
}

void Column::initialize()
{
  for (size_t j = 0; j < numberOfComponents; ++j)
  {
    prefactorMassTransfer[j] = geometry.loadingPrefactor(particleDensity) * components[j].massTransferCoefficient;
  }

  std::fill(partialPressure.begin(), partialPressure.end(), 0.0);
  std::fill(physisorption.begin(), physisorption.end(), 0.0);
  std::fill(chemisorption.begin(), chemisorption.end(), 0.0);
  std::fill(surfaceConcentration.begin(), surfaceConcentration.end(), 0.0);
  std::fill(poreConcentration.begin(), poreConcentration.end(), 0.0);
  std::fill(concentration.begin(), concentration.end(), 0.0);
  std::fill(moleFraction.begin(), moleFraction.end(), 0.0);
  std::fill(bulkSpeciesSink.begin(), bulkSpeciesSink.end(), 0.0);
  std::fill(reactionPhysisorptionSource.begin(), reactionPhysisorptionSource.end(), 0.0);
  std::fill(reactionChemisorptionSource.begin(), reactionChemisorptionSource.end(), 0.0);
  std::fill(reactionPoreConcentrationSource.begin(), reactionPoreConcentrationSource.end(), 0.0);
  std::fill(reactionHeat.begin(), reactionHeat.end(), 0.0);
  std::fill(stateDot.begin(), stateDot.end(), 0.0);

  std::vector<double> initialPressure(numberOfGridPoints + 1, 0.0);

  auto gridRatio = [&](size_t i) -> double
  { return numberOfGridPoints == 0 ? 0.0 : static_cast<double>(i) / static_cast<double>(numberOfGridPoints); };

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
      if (inletPressure > 0.0)
      {
        fillPressure(inletPressure);
      }
      else if (outletPressure > 0.0)
      {
        fillPressure(outletPressure);
      }
      else
      {
        throw std::runtime_error("Error: FixedVelocity requires InletPressure or OutletPressure");
      }

      fillVelocity(columnEntranceVelocity);
      break;
    }

    case BoundaryCondition::FixedPressureInletVelocity:
    {
      for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
      {
        initialPressure[i] = inletPressure + pressureGradient * columnLength * gridRatio(i);
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

  if (fluidPhase == FluidPhase::Gas)
  {
    // Internal nodes initially contain carrier gas; inlet node contains the feed mixture.
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        partialPressure[grid * numberOfComponents + comp] =
            components[comp].isCarrierGas ? initialPressure[grid] : 0.0;
        moleFraction[grid * numberOfComponents + comp] = components[comp].isCarrierGas ? 1.0 : 0.0;
      }
    }

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      partialPressure[comp] = initialPressure[0] * components[comp].initialGasMoleFraction;
      moleFraction[comp] = components[comp].initialGasMoleFraction;
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
      totalConcentration[grid] = totalPressure[grid] / (R * externalTemperature);
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        concentration[grid * numberOfComponents + comp] =
            moleFraction[grid * numberOfComponents + comp] * totalConcentration[grid];
      }
    }
  }
  else
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      totalPressure[grid] = initialPressure[grid];
      totalConcentration[grid] = 0.0;
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        const size_t index = grid * numberOfComponents + comp;
        concentration[index] = grid == 0 ? components[comp].inletLiquidConcentration
                                         : components[comp].initialLiquidConcentration;
        partialPressure[index] = concentration[index];
        totalConcentration[grid] += concentration[index];
      }
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        const size_t index = grid * numberOfComponents + comp;
        moleFraction[index] = totalConcentration[grid] > 0.0 ? concentration[index] / totalConcentration[grid] : 0.0;
      }
    }
  }

  if (surfacePoreTransportEnabled)
  {
    const size_t componentBlockSize = concentration.size();
    for (size_t site = 0; site < maxChemisorptionSites; ++site)
    {
      std::copy(concentration.begin(), concentration.end(),
                surfaceConcentration.begin() + static_cast<std::ptrdiff_t>(site * componentBlockSize));
      std::copy(concentration.begin(), concentration.end(),
                poreConcentration.begin() + static_cast<std::ptrdiff_t>(site * componentBlockSize));
    }
  }

  const double gasDensityTemperature = std::abs(influxTemperature) < 1e-10 ? externalTemperature : influxTemperature;

  if (fluidPhase == FluidPhase::Liquid)
  {
    std::fill(gasDensity.begin(), gasDensity.end(), liquidDensity);
  }
  else for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
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

  computeBulkSpeciesSink(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites, geometry,
                         particleDensity, concentration, physisorptionDot, chemisorptionDot, surfaceConcentration,
                         bulkSpeciesSink, reactionPhysisorptionSource, reactionChemisorptionSource);
  updateVelocityAndPressure(components, boundaryCondition, numberOfGridPoints, numberOfComponents, inletPressure,
                            outletPressure, pressureGradient, columnLength, geometry, columnEntranceVelocity,
                            dynamicViscosity, resolution, interstitialGasVelocity, gasDensity, totalConcentration,
                            totalPressure, concentration, partialPressure, moleFraction, bulkSpeciesSink,
                            gasTemperature, fluidPhase, liquidDensity, pHMode, pHValue, pKw, pHComponent, pH);
  computePhysisorptionEquilibriumLoadings(physisorptionMixture, numberOfGridPoints, numberOfComponents,
                                          maxIsothermTerms, iastPerformance, idealGasMolFractions, adsorbedMolFractions,
                                          numberOfMolecules, totalPressure, equilibriumPhysisorption, cachedPressure,
                                          cachedGrandPotential, moleFraction, gasTemperature,
                                          fluidPhase == FluidPhase::Gas
                                              ? MixturePrediction::DrivingForceInput::MoleFraction
                                              : MixturePrediction::DrivingForceInput::Concentration,
                                          concentration, pH);
  computeChemisorptionEquilibriumLoadings(
      chemisorptionMixture, numberOfGridPoints, numberOfComponents, maxChemisorptionSites, iastPerformance,
      idealGasMolFractions, adsorbedMolFractions, numberOfMolecules, totalPressure, equilibriumChemisorption,
      cachedChemisorptionPressure, cachedChemisorptionGrandPotential, moleFraction, gasTemperature,
      fluidPhase == FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                    : MixturePrediction::DrivingForceInput::Concentration,
      concentration, pH);
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
    std::print(componentStreams[i], "# column 4: Concentration, c_i [mol/m^3]\n");
    std::print(componentStreams[i], "# column 5: Concentration time derivative, dc_i/dt [mol/m^3/s]\n");
    std::print(componentStreams[i], "# column 6: Mole fraction, y_i [-]\n");
    std::print(componentStreams[i], "# column 7: Physisorption, q_phy_i [mol/kg]\n");
    std::print(componentStreams[i], "# column 8: Physisorption time derivative, dq_phy_i/dt [mol/kg/s]\n");
    std::print(componentStreams[i], "# column 9: Chemisorption, q_chem_i [mol/kg]\n");
    std::print(componentStreams[i], "# column 10: Chemisorption time derivative, dq_chem_i/dt [mol/kg/s]\n");
    if (fluidPhase == FluidPhase::Gas)
      std::print(componentStreams[i], "# column 11: Partial pressure, p_i [Pa]\n");
    else
      std::print(componentStreams[i], "# column 11: Liquid concentration driving force, c_i [mol/m^3]\n");
    std::print(componentStreams[i], "# column 12: Equilibrium physisorption, q_phy_i^* [mol/kg]\n");
    if (fluidPhase == FluidPhase::Gas)
      std::print(componentStreams[i], "# column 13: Normalized partial pressure, p_i / (p_t y_i,0) [-]\n");
    else
      std::print(componentStreams[i], "# column 13: Normalized liquid concentration, c_i / c_i,feed [-]\n");
    std::print(componentStreams[i], "# column 14: Equilibrium chemisorption, q_chem_i^* [mol/kg]\n");
    size_t outputColumn = 15;
    for (size_t site = 0; site < components[i].chemisorption.numberOfSites; ++site)
    {
      if (!components[i].chemisorption.sites[site].usesSurfacePoreTransport()) continue;
      std::print(componentStreams[i], "# column {}: Chemisorption site {} surface concentration, c_s_i [mol/m^3]\n",
                 outputColumn++, site);
      std::print(componentStreams[i], "# column {}: Chemisorption site {} pore concentration, c_p_i [mol/m^3]\n",
                 outputColumn++, site);
    }
  }

  std::print(columnStream, "# column 1: Dimensionless time, τ = tv/L [-]\n");
  std::print(columnStream, "# column 2: Time, t [min]\n");
  std::print(columnStream, "# column 3: Column position, z [m]\n");
  if (fluidPhase == FluidPhase::Gas)
    std::print(columnStream, "# column 4: Interstitial gas velocity, v [m/s]\n");
  else
    std::print(columnStream, "# column 4: Interstitial liquid velocity, v [m/s]\n");
  std::print(columnStream, "# column 5: Total pressure, p_t [Pa]\n");
  std::print(columnStream, "# column 6: Gas temperature, T_g [K]\n");
  std::print(columnStream, "# column 7: Gas temperature time derivative, dT_g/dt [K/s]\n");
  std::print(columnStream, "# column 8: Solid temperature, T_s [K]\n");
  std::print(columnStream, "# column 9: Solid temperature time derivative, dT_s/dt [K/s]\n");
  std::print(columnStream, "# column 10: Wall temperature, T_w [K]\n");
  std::print(columnStream, "# column 11: Wall temperature time derivative, dT_w/dt [K/s]\n");
  if (fluidPhase == FluidPhase::Gas)
    std::print(columnStream, "# column 12: Gas density, rho_g [kg/m^3]\n");
  else
  {
    std::print(columnStream, "# column 12: Liquid density, rho_l [kg/m^3]\n");
    std::print(columnStream, "# column 13: pH [-]\n");
  }
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
    if (fluidPhase == FluidPhase::Gas)
      std::print(columnStream, "{} {} {} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor, time / 60.0,
                 static_cast<double>(grid) * resolution, interstitialGasVelocity[grid], totalPressure[grid],
                 gasTemperature[grid], gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid],
                 wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid]);
    else
      std::print(columnStream, "{} {} {} {} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor,
                 time / 60.0, static_cast<double>(grid) * resolution, interstitialGasVelocity[grid],
                 totalPressure[grid], gasTemperature[grid], gasTemperatureDot[grid], solidTemperature[grid],
                 solidTemperatureDot[grid], wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid],
                 pH[grid]);
  }
  std::print(columnStream, "\n\n");

  // Per-component files
  // column 1: dimensionless time [-]
  // column 2: time [min]
  // column 3: column position [m]
  // column 4: concentration [mol/m^3]
  // column 5: concentration time derivative [mol/m^3/s]
  // column 6: mole fraction [-]
  // column 7: physisorption [mol/kg]
  // column 8: physisorption time derivative [mol/kg/s]
  // column 9: chemisorption [mol/kg]
  // column 10: chemisorption time derivative [mol/kg/s]
  // column 11: partial pressure [Pa]
  // column 12: equilibrium adsorption [mol/kg]
  // column 13: normalized partial pressure [-]
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      const size_t index = grid * numberOfComponents + comp;
      const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
      const size_t numberOfChemisorptionSites = components[comp].chemisorption.numberOfSites;
      const double totalChemisorption =
          sumChemisorptionSites(chemisorption, componentBlockSize, numberOfChemisorptionSites, index);
      const double totalChemisorptionDot =
          sumChemisorptionSites(chemisorptionDot, componentBlockSize, numberOfChemisorptionSites, index);

      double normalizedDrivingForce = 0.0;
      if (fluidPhase == FluidPhase::Liquid && components[comp].inletLiquidConcentration > 0.0)
      {
        normalizedDrivingForce = concentration[index] / components[comp].inletLiquidConcentration;
      }
      else if (components[comp].initialGasMoleFraction > 0.0 && totalPressure[grid] > 0.0)
      {
        normalizedDrivingForce =
            partialPressure[index] / (totalPressure[grid] * components[comp].initialGasMoleFraction);
      }

      const double totalEquilibriumChemisorption =
          sumChemisorptionSites(equilibriumChemisorption, componentBlockSize, numberOfChemisorptionSites, index);

      std::print(componentStreams[comp], "{} {} {} {} {} {} {} {} {} {} {} {} {} {}", time * timeNormalizationFactor,
                 time / 60.0, static_cast<double>(grid) * resolution, concentration[index], concentrationDot[index],
                 moleFraction[index], physisorption[index], physisorptionDot[index], totalChemisorption,
                 totalChemisorptionDot, partialPressure[index], equilibriumPhysisorption[index], normalizedDrivingForce,
                 totalEquilibriumChemisorption);
      for (size_t site = 0; site < numberOfChemisorptionSites; ++site)
      {
        if (!components[comp].chemisorption.sites[site].usesSurfacePoreTransport()) continue;
        const size_t siteIndex = site * componentBlockSize + index;
        std::print(componentStreams[comp], " {} {}", surfaceConcentration[siteIndex], poreConcentration[siteIndex]);
      }
      std::print(componentStreams[comp], "\n");
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
  j["maxChemisorptionSites"] = maxChemisorptionSites;
  j["surfacePoreTransportEnabled"] = surfacePoreTransportEnabled;

  j["prefactorMassTransfer"] = prefactorMassTransfer;
  j["idealGasMolFractions"] = idealGasMolFractions;
  j["adsorbedMolFractions"] = adsorbedMolFractions;
  j["numberOfMolecules"] = numberOfMolecules;

  j["interstitialGasVelocity"] = interstitialGasVelocity;
  j["gasDensity"] = gasDensity;
  j["totalConcentration"] = totalConcentration;
  j["totalPressure"] = totalPressure;
  j["pH"] = pH;
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
  j["chemisorption"] = toVector(chemisorption);
  j["chemisorptionDot"] = toVector(chemisorptionDot);
  j["bulkSpeciesSink"] = bulkSpeciesSink;
  if (surfacePoreTransportEnabled)
  {
    j["surfaceConcentration"] = toVector(surfaceConcentration);
    j["surfaceConcentrationDot"] = toVector(surfaceConcentrationDot);
    j["poreConcentration"] = toVector(poreConcentration);
    j["poreConcentrationDot"] = toVector(poreConcentrationDot);
  }
  j["partialPressure"] = partialPressure;
  j["equilibriumPhysisorption"] = equilibriumPhysisorption;
  j["equilibriumChemisorption"] = equilibriumChemisorption;
  j["moleFraction"] = toVector(moleFraction);

  j["cachedPressure"] = cachedPressure;
  j["cachedGrandPotential"] = cachedGrandPotential;
  j["cachedChemisorptionPressure"] = cachedChemisorptionPressure;
  j["cachedChemisorptionGrandPotential"] = cachedChemisorptionGrandPotential;

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
  const size_t fileMaxChemisorptionSites =
      j.contains("maxChemisorptionSites") ? j.at("maxChemisorptionSites").get<size_t>() : 1;

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
  if (fileMaxChemisorptionSites != maxChemisorptionSites)
    throw std::runtime_error("Column::readJSON: maxChemisorptionSites mismatch (file " +
                             std::to_string(fileMaxChemisorptionSites) + ", column " +
                             std::to_string(maxChemisorptionSites) + ")");

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
  loadVectorChecked("pH", pH);

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
  loadSpanChecked("chemisorption", chemisorption);
  loadSpanChecked("chemisorptionDot", chemisorptionDot);
  loadVectorChecked("bulkSpeciesSink", bulkSpeciesSink);
  loadSpanChecked("surfaceConcentration", surfaceConcentration);
  loadSpanChecked("surfaceConcentrationDot", surfaceConcentrationDot);
  loadSpanChecked("poreConcentration", poreConcentration);
  loadSpanChecked("poreConcentrationDot", poreConcentrationDot);

  loadVectorChecked("partialPressure", partialPressure);
  loadVectorChecked("equilibriumPhysisorption", equilibriumPhysisorption);
  loadVectorChecked("equilibriumChemisorption", equilibriumChemisorption);
  loadSpanChecked("moleFraction", moleFraction);

  loadVectorChecked("cachedPressure", cachedPressure);
  loadVectorChecked("cachedGrandPotential", cachedGrandPotential);
  loadVectorChecked("cachedChemisorptionPressure", cachedChemisorptionPressure);
  loadVectorChecked("cachedChemisorptionGrandPotential", cachedChemisorptionGrandPotential);
}
