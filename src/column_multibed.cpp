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

#include "column.h"
#include "component.h"
#include "inputreader.h"
#include "integrators/compute_multibed.h"
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

size_t componentChemisorptionSites(const MultibedColumn& column, size_t comp)
{
  size_t maximum = 0;
  for (const MixturePrediction& mixture : column.physisorptionMixtures)
  {
    maximum = std::max(maximum, mixture.components[comp].chemisorption.numberOfSites);
  }
  return maximum;
}

bool componentSiteUsesSurfacePoreTransport(const MultibedColumn& column, size_t comp, size_t site)
{
  for (const MixturePrediction& mixture : column.physisorptionMixtures)
  {
    const MultiSiteChemisorption& multisite = mixture.components[comp].chemisorption;
    if (site < multisite.numberOfSites && multisite.sites[site].usesSurfacePoreTransport()) return true;
  }
  return false;
}

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

void validateAdsorbentVectors(const MultibedColumn& column)
{
  if (column.numberOfAdsorbents == 0)
  {
    throw std::runtime_error("Error: multibed column requires at least one adsorbent");
  }
  for (const MixturePrediction& mixture : column.physisorptionMixtures)
  {
    if (mixture.components.size() != column.numberOfComponents)
    {
      throw std::runtime_error("Error: every adsorbent must define the same component set");
    }
  }
  if (column.adsorbentMixFractions.empty() && column.adsorbentLengths.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: AdsorbentLengths size must equal numberOfAdsorbents");
  }
  if (column.adsorbentMixFractions.empty() &&
      column.adsorbentInterfaceLengths.size() != column.numberOfAdsorbents - 1)
  {
    throw std::runtime_error("Error: AdsorbentInterfaceLengths size must equal numberOfAdsorbents - 1");
  }
  if (!column.adsorbentMixFractions.empty())
  {
    if (!column.adsorbentLengths.empty())
    {
      throw std::runtime_error("Error: AdsorbentMixFractions cannot be combined with adsorbent lengths");
    }
    if (column.adsorbentMixFractions.size() != column.numberOfAdsorbents)
    {
      throw std::runtime_error("Error: AdsorbentMixFractions size must equal numberOfAdsorbents");
    }
    const double sum = std::reduce(column.adsorbentMixFractions.begin(), column.adsorbentMixFractions.end(), 0.0);
    if (std::abs(sum - 1.0) > 1.0e-12 ||
        std::ranges::any_of(column.adsorbentMixFractions,
                            [](double fraction)
                            { return !std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0; }))
    {
      throw std::runtime_error("Error: AdsorbentMixFractions values must be between 0 and 1 and sum to 1");
    }
    if (!column.adsorbentInterfaceLengths.empty())
    {
      throw std::runtime_error("Error: AdsorbentMixFractions cannot be combined with interface lengths");
    }
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
  if (column.geometries.size() != column.numberOfAdsorbents)
  {
    throw std::runtime_error("Error: Geometries size must equal numberOfAdsorbents");
  }
  if (column.columnLength <= 0.0)
  {
    throw std::runtime_error("Error: column length must be positive");
  }
}
}  // namespace

size_t MultibedColumn::stateSize() const noexcept { return stateLayout().stateSize(); }

bool MultibedColumn::requiresSurfacePoreTransport(const std::vector<MixturePrediction>& mixtures) noexcept
{
  return std::ranges::any_of(mixtures, [](const MixturePrediction& mixture)
                             { return Column::requiresSurfacePoreTransport(mixture.components); });
}

size_t MultibedColumn::maximumChemisorptionSites(const std::vector<MixturePrediction>& mixtures) noexcept
{
  size_t maximum = 1;
  for (const MixturePrediction& mixture : mixtures)
  {
    maximum = std::max(maximum, Column::maximumChemisorptionSites(mixture.components));
  }
  return maximum;
}

std::vector<Geometry> MultibedColumn::makePackedBedGeometries(std::span<const double> voidFractions,
                                                               std::span<const double> particleDiameters,
                                                               double internalDiameter, double outerDiameter)
{
  if (voidFractions.size() != particleDiameters.size())
  {
    throw std::runtime_error("Error: geometry input sizes must match");
  }
  std::vector<Geometry> result;
  result.reserve(voidFractions.size());
  for (size_t ads = 0; ads < voidFractions.size(); ++ads)
  {
    result.push_back(makeGeometry(PackedBedTubeSpec{.voidFraction = voidFractions[ads],
                                                    .particleDiameter = particleDiameters[ads],
                                                    .internalDiameter = internalDiameter,
                                                    .outerDiameter = outerDiameter}));
  }
  return result;
}

std::vector<MixturePrediction> MultibedColumn::makeChemisorptionMixtures(
    const std::vector<MixturePrediction>& physisorptionMixtures)
{
  std::vector<MixturePrediction> mixtures;
  mixtures.reserve(physisorptionMixtures.size());
  for (const MixturePrediction& physisorptionMixture : physisorptionMixtures)
  {
    mixtures.push_back(Column::makeChemisorptionMixture(physisorptionMixture, physisorptionMixture.components));
  }
  return mixtures;
}

MultibedColumnStateLayout MultibedColumn::stateLayout() const noexcept
{
  return MultibedColumnStateLayout{numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                                   surfacePoreTransportEnabled};
}

void MultibedColumn::bindStateViews() noexcept
{
  const MultibedColumnStateLayout layout = stateLayout();
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

MultibedColumn::MultibedColumn(const MultibedColumn& other)
    : physisorptionMixtures(other.physisorptionMixtures),
      chemisorptionMixtures(other.chemisorptionMixtures),
      components(other.components),
      boundaryCondition(other.boundaryCondition),
      fluidPhase(other.fluidPhase),
      pHMode(other.pHMode),
      energyBalance(other.energyBalance),
      reactions(other.reactions),
      numberOfGridPoints(other.numberOfGridPoints),
      numberOfComponents(other.numberOfComponents),
      numberOfAdsorbents(other.numberOfAdsorbents),
      maxIsothermTerms(other.maxIsothermTerms),
      maxChemisorptionSites(other.maxChemisorptionSites),
      numberOfCalls(other.numberOfCalls),
      carrierGasComponent(other.carrierGasComponent),
      pHComponent(other.pHComponent),
      adsorbentLengths(other.adsorbentLengths),
      adsorbentInterfaceLengths(other.adsorbentInterfaceLengths),
      adsorbentMixFractions(other.adsorbentMixFractions),
      adsorbentGridPoints(other.adsorbentGridPoints),
      adsorbentVoidFractions(other.adsorbentVoidFractions),
      particleDensities(other.particleDensities),
      particleDiameters(other.particleDiameters),
      geometries(other.geometries),
      externalTemperature(other.externalTemperature),
      inletPressure(other.inletPressure),
      outletPressure(other.outletPressure),
      pressureGradient(other.pressureGradient),
      columnEntranceVelocity(other.columnEntranceVelocity),
      columnLength(other.columnLength),
      dynamicViscosity(other.dynamicViscosity),
      liquidDensity(other.liquidDensity),
      pHValue(other.pHValue),
      pKw(other.pKw),
      columnDistances(other.columnDistances),
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
      totalVoidFraction(other.totalVoidFraction),
      particleDensity(other.particleDensity),
      moleFraction(other.moleFraction),
      partialPressure(other.partialPressure),
      equilibriumPhysisorption(other.equilibriumPhysisorption),
      equilibriumChemisorption(other.equilibriumChemisorption),
      fractionOfAdsorbent(other.fractionOfAdsorbent),
      hasAdsorbentOfType(other.hasAdsorbentOfType),
      adsorbentScaledVoidFraction(other.adsorbentScaledVoidFraction),
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

MultibedColumn& MultibedColumn::operator=(const MultibedColumn& other)
{
  if (this == &other) return *this;

  physisorptionMixtures = other.physisorptionMixtures;
  chemisorptionMixtures = other.chemisorptionMixtures;
  components = other.components;
  boundaryCondition = other.boundaryCondition;
  fluidPhase = other.fluidPhase;
  pHMode = other.pHMode;
  energyBalance = other.energyBalance;
  reactions = other.reactions;
  numberOfGridPoints = other.numberOfGridPoints;
  numberOfComponents = other.numberOfComponents;
  numberOfAdsorbents = other.numberOfAdsorbents;
  maxIsothermTerms = other.maxIsothermTerms;
  maxChemisorptionSites = other.maxChemisorptionSites;
  numberOfCalls = other.numberOfCalls;
  carrierGasComponent = other.carrierGasComponent;
  pHComponent = other.pHComponent;
  adsorbentLengths = other.adsorbentLengths;
  adsorbentInterfaceLengths = other.adsorbentInterfaceLengths;
  adsorbentMixFractions = other.adsorbentMixFractions;
  adsorbentGridPoints = other.adsorbentGridPoints;
  adsorbentVoidFractions = other.adsorbentVoidFractions;
  particleDensities = other.particleDensities;
  particleDiameters = other.particleDiameters;
  geometries = other.geometries;
  externalTemperature = other.externalTemperature;
  inletPressure = other.inletPressure;
  outletPressure = other.outletPressure;
  pressureGradient = other.pressureGradient;
  columnEntranceVelocity = other.columnEntranceVelocity;
  columnLength = other.columnLength;
  dynamicViscosity = other.dynamicViscosity;
  liquidDensity = other.liquidDensity;
  pHValue = other.pHValue;
  pKw = other.pKw;
  columnDistances = other.columnDistances;
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
  totalVoidFraction = other.totalVoidFraction;
  particleDensity = other.particleDensity;
  moleFraction = other.moleFraction;
  partialPressure = other.partialPressure;
  equilibriumPhysisorption = other.equilibriumPhysisorption;
  equilibriumChemisorption = other.equilibriumChemisorption;
  fractionOfAdsorbent = other.fractionOfAdsorbent;
  hasAdsorbentOfType = other.hasAdsorbentOfType;
  adsorbentScaledVoidFraction = other.adsorbentScaledVoidFraction;
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

void MultibedColumn::initialize()
{
  validateAdsorbentVectors(*this);

  std::vector<double> boundaries;
  if (adsorbentMixFractions.empty())
  {
    boundaries.assign(numberOfAdsorbents > 1 ? numberOfAdsorbents - 1 : 0, 0.0);
    double cumulativeLength = 0.0;
    for (size_t ads = 0; ads + 1 < numberOfAdsorbents; ++ads)
    {
      cumulativeLength += adsorbentLengths[ads];
      boundaries[ads] = cumulativeLength;
    }
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const double z = columnDistances[grid];
    const size_t base = grid * numberOfAdsorbents;

    if (!adsorbentMixFractions.empty())
    {
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        fractionOfAdsorbent[base + ads] = adsorbentMixFractions[ads];
      }
    }
    else
    {
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
    }

    totalVoidFraction[grid] = 1.0;
    particleDensity[grid] = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[base + ads];
      hasAdsorbentOfType[base + ads] = fraction > 0.0;
      totalVoidFraction[grid] -= fraction * (1.0 - geometries[ads].voidFraction);
      particleDensity[grid] += fraction * particleDensities[ads];
    }

    if (totalVoidFraction[grid] <= 0.0)
    {
      throw std::runtime_error("Error: total void fraction must be positive");
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (geometries[ads].kind == GeometryKind::PackedBed && particleDiameters[ads] <= 0.0)
      {
        throw std::runtime_error("Error: ParticleDiameters values must be positive");
      }
      adsorbentScaledVoidFraction[base + ads] = geometries[ads].kind == GeometryKind::PackedBed
                                                    ? (1.0 - geometries[ads].voidFraction) /
                                                          (totalVoidFraction[grid] * particleDiameters[ads])
                                                    : 0.0;
    }
  }

  for (size_t j = 0; j < numberOfComponents; ++j)
  {
    prefactorMassTransfer[j] = components[j].massTransferCoefficient;
  }

  std::fill(partialPressure.begin(), partialPressure.end(), 0.0);
  std::fill(physisorption.begin(), physisorption.end(), 0.0);
  std::fill(chemisorption.begin(), chemisorption.end(), 0.0);
  std::fill(concentration.begin(), concentration.end(), 0.0);
  std::fill(moleFraction.begin(), moleFraction.end(), 0.0);
  std::fill(bulkSpeciesSink.begin(), bulkSpeciesSink.end(), 0.0);
  std::fill(stateDot.begin(), stateDot.end(), 0.0);

  std::vector<double> initialPressure(numberOfGridPoints + 1, 0.0);

  auto gridRatio = [&](size_t i) -> double { return columnLength <= 0.0 ? 0.0 : columnDistances[i] / columnLength; };

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
        initialPressure[i] = inletPressure + pressureGradient * columnDistances[i];
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

  computeBulkSpeciesSink(physisorptionMixtures, numberOfGridPoints, numberOfComponents, numberOfAdsorbents,
                         maxChemisorptionSites, geometries, adsorbentVoidFractions, particleDensities,
                         particleDiameters,
                         fractionOfAdsorbent, totalVoidFraction, concentration, physisorptionDot, chemisorptionDot,
                         surfaceConcentration, bulkSpeciesSink, reactionPhysisorptionSource,
                         reactionChemisorptionSource);
  updateVelocityAndPressure(components, boundaryCondition, numberOfGridPoints, numberOfComponents, inletPressure,
                            outletPressure, pressureGradient, columnLength, numberOfAdsorbents, columnEntranceVelocity,
                            dynamicViscosity, columnDistances, fractionOfAdsorbent, geometries,
                            adsorbentScaledVoidFraction, totalVoidFraction,
                            interstitialGasVelocity, gasDensity, totalConcentration, totalPressure, concentration,
                            partialPressure, moleFraction, bulkSpeciesSink, gasTemperature, fluidPhase, liquidDensity,
                            pHMode, pHValue, pKw, pHComponent, pH);
  computePhysisorptionEquilibriumLoadings(physisorptionMixtures, numberOfGridPoints, numberOfComponents,
                                          numberOfAdsorbents, fractionOfAdsorbent, hasAdsorbentOfType, maxIsothermTerms,
                                          iastPerformance, idealGasMolFractions, adsorbedMolFractions,
                                          numberOfMolecules, totalPressure, equilibriumPhysisorption, cachedPressure,
                                          cachedGrandPotential, moleFraction, gasTemperature,
                                          fluidPhase == FluidPhase::Gas
                                              ? MixturePrediction::DrivingForceInput::MoleFraction
                                              : MixturePrediction::DrivingForceInput::Concentration,
                                          concentration, pH);
  computeChemisorptionEquilibriumLoadings(
      chemisorptionMixtures, numberOfGridPoints, numberOfComponents, numberOfAdsorbents, fractionOfAdsorbent,
      hasAdsorbentOfType, maxChemisorptionSites, iastPerformance, idealGasMolFractions, adsorbedMolFractions,
      numberOfMolecules, totalPressure, equilibriumChemisorption, cachedChemisorptionPressure,
      cachedChemisorptionGrandPotential, moleFraction, gasTemperature,
      fluidPhase == FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                    : MixturePrediction::DrivingForceInput::Concentration,
      concentration, pH);
}

void MultibedColumn::setTemperature(double temperature)
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

void MultibedColumn::writeOutputHeader(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream) const
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
    for (size_t site = 0; site < componentChemisorptionSites(*this, i); ++site)
    {
      if (!componentSiteUsesSurfacePoreTransport(*this, i, site)) continue;
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

void MultibedColumn::writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream,
                                 double time) const
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
                 columnDistances[grid], interstitialGasVelocity[grid], totalPressure[grid], gasTemperature[grid],
                 gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid], wallTemperature[grid],
                 wallTemperatureDot[grid], gasDensity[grid]);
    else
      std::print(columnStream, "{} {} {} {} {} {} {} {} {} {} {} {} {}\n", time * timeNormalizationFactor,
                 time / 60.0, columnDistances[grid], interstitialGasVelocity[grid], totalPressure[grid],
                 gasTemperature[grid], gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid],
                 wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid], pH[grid]);
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
  // column 12: equilibrium physisorption [mol/kg]
  // column 13: normalized partial pressure [-]
  // column 14: equilibrium chemisorption [mol/kg]
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      const size_t index = grid * numberOfComponents + comp;
      const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
      const size_t numberOfChemisorptionSites = componentChemisorptionSites(*this, comp);
      const double totalChemisorption =
          sumChemisorptionSites(chemisorption, componentBlockSize, numberOfChemisorptionSites, index);
      const double totalChemisorptionDot =
          sumChemisorptionSites(chemisorptionDot, componentBlockSize, numberOfChemisorptionSites, index);
      const double totalEquilibriumChemisorption =
          sumChemisorptionSites(equilibriumChemisorption, componentBlockSize, numberOfChemisorptionSites, index);

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

      std::print(componentStreams[comp], "{} {} {} {} {} {} {} {} {} {} {} {} {} {}", time * timeNormalizationFactor,
                 time / 60.0, columnDistances[grid], concentration[index], concentrationDot[index], moleFraction[index],
                 physisorption[index], physisorptionDot[index], totalChemisorption, totalChemisorptionDot,
                 partialPressure[index], equilibriumPhysisorption[index], normalizedDrivingForce,
                 totalEquilibriumChemisorption);
      for (size_t site = 0; site < numberOfChemisorptionSites; ++site)
      {
        if (!componentSiteUsesSurfacePoreTransport(*this, comp, site)) continue;
        const size_t siteIndex = site * componentBlockSize + index;
        std::print(componentStreams[comp], " {} {}", surfaceConcentration[siteIndex], poreConcentration[siteIndex]);
      }
      std::print(componentStreams[comp], "\n");
    }
    std::print(componentStreams[comp], "\n\n");
  }
}

std::string MultibedColumn::repr() const
{
  return std::format(
      "Temperature:                           {} [K]\n"
      "Column length:                         {} [m]\n"
      "Adsorbent lengths:                     {} [m]\n"
      "Adsorbent interface lengths:           {} [m]\n"
      "Uniform adsorbent mix fractions:       {} [-]\n"
      "Adsorbent void-fractions:              {} [-]\n"
      "Particle densities:                    {} [kg/m^3]\n"
      "Inlet pressure:                        {} [Pa]\n"
      "Outlet pressure:                       {} [Pa]\n"
      "Pressure gradient:                     {} [Pa]\n"
      "Column entrance interstitial velocity: {} [m/s]\n"
      "\n\n",
      externalTemperature, columnLength, vectorToString(adsorbentLengths), vectorToString(adsorbentInterfaceLengths),
      vectorToString(adsorbentMixFractions), vectorToString(adsorbentVoidFractions),
      vectorToString(particleDensities), inletPressure, outletPressure, pressureGradient, columnEntranceVelocity);
}

void MultibedColumn::writeJSON(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("MultibedColumn::writeJSON: cannot open file '" + filename + "'");

  nlohmann::json j;
  j["numberOfGridPoints"] = numberOfGridPoints;
  j["numberOfComponents"] = numberOfComponents;
  j["maxIsothermTerms"] = maxIsothermTerms;
  j["maxChemisorptionSites"] = maxChemisorptionSites;
  j["surfacePoreTransportEnabled"] = surfacePoreTransportEnabled;
  j["columnDistances"] = columnDistances;

  j["prefactorMassTransfer"] = prefactorMassTransfer;
  j["idealGasMolFractions"] = idealGasMolFractions;
  j["adsorbedMolFractions"] = adsorbedMolFractions;
  j["numberOfMolecules"] = numberOfMolecules;

  j["interstitialGasVelocity"] = interstitialGasVelocity;
  j["gasDensity"] = gasDensity;
  j["totalConcentration"] = totalConcentration;
  j["totalVoidFraction"] = totalVoidFraction;
  j["particleDensity"] = particleDensity;
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
  j["moleFraction"] = moleFraction;
  j["fractionOfAdsorbent"] = fractionOfAdsorbent;
  j["adsorbentScaledVoidFraction"] = adsorbentScaledVoidFraction;

  j["cachedPressure"] = cachedPressure;
  j["cachedGrandPotential"] = cachedGrandPotential;
  j["cachedChemisorptionPressure"] = cachedChemisorptionPressure;
  j["cachedChemisorptionGrandPotential"] = cachedChemisorptionGrandPotential;

  out << j.dump(4) << "\n";
}

void MultibedColumn::readJSON(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("MultibedColumn::readJSONFile: cannot open file '" + filename + "'");

  nlohmann::json j;
  in >> j;

  auto requireSizeT = [&](const char* key) -> size_t
  {
    if (!j.contains(key))
      throw std::runtime_error(std::string("MultibedColumn::readJSON: missing required key '") + key + "'");
    return j.at(key).get<size_t>();
  };

  const size_t fileNgrid = requireSizeT("numberOfGridPoints");
  const size_t fileNcomp = requireSizeT("numberOfComponents");
  const size_t fileMaxIsothermTerms = requireSizeT("maxIsothermTerms");
  const size_t fileMaxChemisorptionSites =
      j.contains("maxChemisorptionSites") ? j.at("maxChemisorptionSites").get<size_t>() : 1;

  if (fileNgrid != numberOfGridPoints)
    throw std::runtime_error("MultibedColumn::readJSON: numberOfGridPoints mismatch (file " +
                             std::to_string(fileNgrid) + ", column " + std::to_string(numberOfGridPoints) + ")");
  if (fileNcomp != numberOfComponents)
    throw std::runtime_error("MultibedColumn::readJSON: numberOfComponents mismatch (file " +
                             std::to_string(fileNcomp) + ", column " + std::to_string(numberOfComponents) + ")");
  if (fileMaxIsothermTerms != maxIsothermTerms)
    throw std::runtime_error("MultibedColumn::readJSON: maxIsothermTerms mismatch (file " +
                             std::to_string(fileMaxIsothermTerms) + ", column " + std::to_string(maxIsothermTerms) +
                             ")");
  if (fileMaxChemisorptionSites != maxChemisorptionSites)
    throw std::runtime_error("MultibedColumn::readJSON: maxChemisorptionSites mismatch (file " +
                             std::to_string(fileMaxChemisorptionSites) + ", column " +
                             std::to_string(maxChemisorptionSites) + ")");

  auto loadVectorChecked = [&](const char* key, std::vector<double>& dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array())
      throw std::runtime_error(std::string("MultibedColumn::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("MultibedColumn::readJSON: size mismatch for '" + std::string(key) + "'");

    dst = a.get<std::vector<double>>();
  };

  auto loadSpanChecked = [&](const char* key, std::span<double> dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array())
      throw std::runtime_error(std::string("MultibedColumn::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("MultibedColumn::readJSON: size mismatch for '" + std::string(key) + "'");

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
  loadVectorChecked("totalVoidFraction", totalVoidFraction);
  loadVectorChecked("particleDensity", particleDensity);
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
  loadVectorChecked("moleFraction", moleFraction);
  loadVectorChecked("fractionOfAdsorbent", fractionOfAdsorbent);
  loadVectorChecked("adsorbentScaledVoidFraction", adsorbentScaledVoidFraction);

  loadVectorChecked("cachedPressure", cachedPressure);
  loadVectorChecked("cachedGrandPotential", cachedGrandPotential);
  loadVectorChecked("cachedChemisorptionPressure", cachedChemisorptionPressure);
  loadVectorChecked("cachedChemisorptionGrandPotential", cachedChemisorptionGrandPotential);
  loadVectorChecked("columnDistances", columnDistances);
  validateColumnDistances(columnDistances, numberOfGridPoints, columnLength);
}
