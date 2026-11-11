#include "column.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

size_t componentChemisorptionSites(const Column& column, size_t comp)
{
  size_t maximum = 0;
  for (const MixturePrediction& mixture : column.physisorptionMixtures)
  {
    maximum = std::max(maximum, mixture.components[comp].chemisorption.numberOfSites);
  }
  return maximum;
}

bool componentSiteUsesSurfacePoreTransport(const Column& column, size_t comp, size_t site)
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

void validateAdsorbentVectors(const Column& column)
{
  if (column.numberOfAdsorbents == 0)
  {
    throw std::runtime_error("Error: a column requires at least one adsorbent");
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
  if (column.adsorbentMixFractions.empty() &&
      column.adsorbentInterfaceGridPoints.size() != column.numberOfAdsorbents - 1)
  {
    throw std::runtime_error("Error: AdsorbentInterfaceGridPoints size must equal numberOfAdsorbents - 1");
  }
  if (std::ranges::any_of(column.adsorbentInterfaceLengths, [](double length) { return length < 0.0; }))
  {
    throw std::runtime_error("Error: interface lengths must be non-negative");
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

size_t Column::stateSize() const noexcept { return stateLayout().stateSize(); }

bool ColumnData::requiresSurfacePoreTransport(const std::vector<MixturePrediction>& mixtures) noexcept
{
  return std::ranges::any_of(mixtures,
                             [](const MixturePrediction& mixture)
                             {
                               return std::ranges::any_of(mixture.components, [](const Component& component)
                                                          { return component.chemisorption.usesSurfacePoreTransport(); });
                             });
}

size_t ColumnData::maximumChemisorptionSites(const std::vector<MixturePrediction>& mixtures) noexcept
{
  size_t maximum = 1;
  for (const MixturePrediction& mixture : mixtures)
  {
    for (const Component& component : mixture.components)
    {
      maximum = std::max(maximum, component.chemisorption.numberOfSites);
    }
  }
  return maximum;
}

std::vector<Geometry> ColumnData::makePackedBedGeometries(std::span<const double> voidFractions,
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

std::vector<MixturePrediction> ColumnData::makeChemisorptionMixtures(
    const std::vector<MixturePrediction>& physisorptionMixtures)
{
  std::vector<MixturePrediction> mixtures;
  mixtures.reserve(physisorptionMixtures.size());
  for (const MixturePrediction& physisorptionMixture : physisorptionMixtures)
  {
    mixtures.push_back(
        MixturePrediction::makeChemisorptionPrediction(physisorptionMixture, physisorptionMixture.components));
  }
  return mixtures;
}

ColumnStateLayout Column::stateLayout() const noexcept
{
  return ColumnStateLayout{numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                                   surfacePoreTransportEnabled};
}

void Column::bindStateViews() noexcept
{
  static_cast<ColumnStateViews&>(*this) = stateLayout().views(state.data(), stateDot.data());
}

void Column::initialize()
{
  validateAdsorbentVectors(*this);

  std::vector<double> interfaceStarts;
  if (adsorbentMixFractions.empty())
  {
    interfaceStarts.resize(adsorbentInterfaceLengths.size());
    double cumulativeLength = 0.0;
    for (size_t ads = 0; ads + 1 < numberOfAdsorbents; ++ads)
    {
      cumulativeLength += adsorbentLengths[ads];
      interfaceStarts[ads] = cumulativeLength;
      cumulativeLength += adsorbentInterfaceLengths[ads];
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
      while (region + 1 < numberOfAdsorbents && z > interfaceStarts[region] + adsorbentInterfaceLengths[region])
      {
        ++region;
      }

      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        fractionOfAdsorbent[base + ads] = ads == region ? 1.0 : 0.0;
      }

      if (region + 1 < numberOfAdsorbents && z > interfaceStarts[region] && adsorbentInterfaceLengths[region] > 0.0)
      {
        const double rightFraction =
            z == interfaceStarts[region] + adsorbentInterfaceLengths[region]
                ? 1.0
                : std::clamp((z - interfaceStarts[region]) / adsorbentInterfaceLengths[region], 0.0, 1.0);
        fractionOfAdsorbent[base + region] = 1.0 - rightFraction;
        fractionOfAdsorbent[base + region + 1] = rightFraction;
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
    }
  }
  computeAdsorbentBlends();

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

  // In a reverse step the feed enters at z = L: the velocity points towards z = 0, and the inlet pressure applies at
  // z = L.
  const bool reverse = flowDirection == FlowDirection::Reverse;
  auto fillVelocity = [&](double velocity)
  { std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), reverse ? -velocity : velocity); };

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
        const double fromInlet = reverse ? 1.0 - gridRatio(i) : gridRatio(i);
        initialPressure[i] = inletPressure + fromInlet * (outletPressure - inletPressure);
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
        initialPressure[i] =
            inletPressure + pressureGradient * (reverse ? columnLength - columnDistances[i] : columnDistances[i]);
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

  // The column starts at the feed temperature, which sets the initial gas, solid, and wall temperatures below, and
  // with them the initial and feed concentrations P/(R T) (C-BC-2).
  const double initialGasTemperature = std::abs(influxTemperature) < 1e-10 ? externalTemperature : influxTemperature;

  if (fluidPhase == FluidPhase::Gas)
  {
    // Internal nodes initially contain carrier gas; inlet node contains the feed mixture, which a reverse step takes
    // from the product-end gas at the pressure of z = L instead.
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        partialPressure[grid * numberOfComponents + comp] =
            components[comp].isCarrierGas ? initialPressure[grid] : 0.0;
        moleFraction[grid * numberOfComponents + comp] = components[comp].isCarrierGas ? 1.0 : 0.0;
      }
    }

    productEndConcentration.clear();
    for (size_t comp = 0; comp < numberOfComponents && !reverse; ++comp)
    {
      partialPressure[comp] = initialPressure[0] * components[comp].initialGasMoleFraction;
      moleFraction[comp] = components[comp].initialGasMoleFraction;
    }
    for (size_t comp = 0; comp < numberOfComponents && reverse; ++comp)
    {
      productEndConcentration.push_back(initialPressure[numberOfGridPoints] *
                                        components[comp].initialGasMoleFraction / (R * initialGasTemperature));
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
      totalConcentration[grid] = totalPressure[grid] / (R * initialGasTemperature);
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        concentration[grid * numberOfComponents + comp] =
            moleFraction[grid * numberOfComponents + comp] * totalConcentration[grid];
      }
    }
  }
  else
  {
    productEndConcentration.clear();
    for (size_t comp = 0; comp < numberOfComponents && reverse; ++comp)
    {
      productEndConcentration.push_back(components[comp].inletLiquidConcentration);
    }
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      totalPressure[grid] = initialPressure[grid];
      totalConcentration[grid] = 0.0;
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        const size_t index = grid * numberOfComponents + comp;
        concentration[index] = grid == 0 && !reverse ? components[comp].inletLiquidConcentration
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

    gasDensity[grid] = molecularWeight > 0.0 && initialGasTemperature > 0.0
                           ? totalPressure[grid] * molecularWeight / (R * initialGasTemperature)
                           : 0.0;
  }

  if (std::abs(influxTemperature) < 1e-10) influxTemperature = externalTemperature;
  productEndTemperature = reverse ? influxTemperature : 0.0;
  std::fill(gasTemperature.begin(), gasTemperature.end(), influxTemperature);
  std::fill(solidTemperature.begin(), solidTemperature.end(), influxTemperature);
  std::fill(wallTemperature.begin(), wallTemperature.end(), influxTemperature);

  computeBulkSpeciesSink(physisorptionMixtures, numberOfGridPoints, numberOfComponents, numberOfAdsorbents,
                         maxChemisorptionSites, geometries, fractionOfAdsorbent, loadingPrefactor, concentration,
                         physisorptionDot, chemisorptionDot, surfaceConcentration, bulkSpeciesSink,
                         reactionPhysisorptionSource, reactionChemisorptionSource, reactionBulkSource);
  updateVelocityAndPressure(components, boundaryCondition, numberOfGridPoints, numberOfComponents, inletPressure,
                            outletPressure, pressureGradient, columnLength, numberOfAdsorbents, columnEntranceVelocity,
                            dynamicViscosity, columnDistances, fractionOfAdsorbent, geometries, pressureDropLaw,
                            totalVoidFraction, interstitialGasVelocity, gasDensity, totalConcentration, totalPressure,
                            concentration, partialPressure, moleFraction, bulkSpeciesSink, gasTemperature, fluidPhase,
                            liquidDensity, pHMode, pHValue, pKw, pHComponent, pH, molecularDispersionCoefficient,
                            dispersionLaw, dispersionReferenceTemperature, dispersionReferencePressure,
                            axialDispersionCoefficient, convectionScheme, inletCondition, nullptr, flowEnds());
  computePhysisorptionEquilibriumLoadings(
      physisorptionMixtures, numberOfGridPoints, numberOfComponents, numberOfAdsorbents, fractionOfAdsorbent,
      hasAdsorbentOfType, maxIsothermTerms, iastPerformance, idealGasMolFractions, adsorbedMolFractions,
      numberOfMolecules, totalPressure, equilibriumPhysisorption, cachedPressure, cachedGrandPotential, moleFraction,
      gasTemperature,
      fluidPhase == FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
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

void Column::computeAdsorbentBlends()
{
  const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const size_t base = grid * numberOfAdsorbents;
    // The surface and pore reservoirs hold the solid volume, so the pore reservoir loses the solid-volume-weighted
    // particle density times dq/dt (rho_p in a pure section).
    double solidLoadingDensity = 0.0;
    double solidVolumeDensity = 0.0;
    PressureDropLaw law{};
    DispersionLaw dispersion{};
    solidFluidContactArea[grid] = 0.0;
    fluidSolidContactArea[grid] = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[base + ads];
      const double solidVolume = fraction * geometries[ads].solidToFluidVolumeRatio * geometries[ads].voidFraction;
      solidLoadingDensity += solidVolume * particleDensities[ads];
      solidVolumeDensity += solidVolume;
      solidFluidContactArea[grid] += fraction * geometries[ads].contactAreas.solidFluidPerSolidVolume;
      fluidSolidContactArea[grid] += fraction * geometries[ads].contactAreas.fluidSolidPerFluidVolume;
      law.viscousCoefficient += fraction * geometries[ads].pressureDrop.viscousCoefficient;
      law.inertialCoefficient += fraction * geometries[ads].pressureDrop.inertialCoefficient;
      dispersion.dispersivity += fraction * geometries[ads].dispersion.dispersivity;
      dispersion.taylorArisCoefficient += fraction * geometries[ads].dispersion.taylorArisCoefficient;
    }
    loadingPrefactor[grid] = solidLoadingDensity / std::max(1e-10, totalVoidFraction[grid]);
    poreLoadingConversion[grid] = solidLoadingDensity / std::max(1.0e-30, solidVolumeDensity);
    pressureDropLaw[grid] = law;
    dispersionLaw[grid] = dispersion;

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      massTransferCoefficient[index] = 0.0;
      molecularDispersionCoefficient[index] = 0.0;
      heatOfAdsorption[index] = 0.0;
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        const double fraction = fractionOfAdsorbent[base + ads];
        const Component& component = physisorptionMixtures[ads].components[comp];
        massTransferCoefficient[index] += fraction * component.massTransferCoefficient;
        molecularDispersionCoefficient[index] += fraction * component.axialDispersionCoefficient;
        heatOfAdsorption[index] += fraction * component.heatOfAdsorption;
      }
      // updateVelocityAndPressure applies the dispersion law; until then the dispersion is its molecular part.
      axialDispersionCoefficient[index] = molecularDispersionCoefficient[index];
      for (size_t site = 0; site < maxChemisorptionSites; ++site)
      {
        double heat = 0.0;
        for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
        {
          const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
          if (site < multisite.numberOfSites)
          {
            heat += fractionOfAdsorbent[base + ads] * multisite.sites[site].heatOfChemisorption;
          }
        }
        heatOfChemisorption[site * componentBlockSize + index] = heat;
      }
    }
  }
}

double Column::idealGasDeviation() const
{
  if (fluidPhase != FluidPhase::Gas || boundaryCondition == BoundaryCondition::FixedVelocity ||
      boundaryCondition == BoundaryCondition::FixedPressureInletVelocity)
  {
    return 0.0;
  }
  double deviation = 0.0;
  for (size_t grid = flowDirection == FlowDirection::Reverse ? 1 : 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += concentration[grid * numberOfComponents + comp];
    }
    deviation = std::max(deviation, std::abs(R * gasTemperature[grid] * concentrationSum / totalPressure[grid] - 1.0));
  }
  return deviation;
}

std::vector<double> Column::gasInventory() const
{
  std::vector<double> amount(numberOfComponents, 0.0);
  for (size_t grid = 1; grid <= numberOfGridPoints; ++grid)
  {
    const double fluidVolume = totalVoidFraction[grid] * gridSpacing(columnDistances, grid);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      amount[comp] += fluidVolume * concentration[grid * numberOfComponents + comp];
    }
  }
  return amount;
}

std::vector<double> Column::adsorbedInventory() const
{
  std::vector<double> amount(numberOfComponents, 0.0);
  const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
  for (size_t grid = 1; grid <= numberOfGridPoints; ++grid)
  {
    // loadingPrefactor is the solid mass per fluid volume.
    const double solidMass = loadingPrefactor[grid] * totalVoidFraction[grid] * gridSpacing(columnDistances, grid);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      double loading = physisorption[index];
      for (size_t site = 0; site < maxChemisorptionSites; ++site)
      {
        loading += chemisorption[site * componentBlockSize + index];
      }
      amount[comp] += solidMass * loading;
    }
  }
  return amount;
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
    for (const size_t site : surfacePoreTransportSites(i))
    {
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

std::vector<size_t> Column::surfacePoreTransportSites(size_t comp) const
{
  std::vector<size_t> sites;
  for (size_t site = 0; site < componentChemisorptionSites(*this, comp); ++site)
  {
    if (componentSiteUsesSurfacePoreTransport(*this, comp, site)) sites.push_back(site);
  }
  return sites;
}

void Column::appendColumnOutput(std::vector<double>& row, size_t grid, double time) const
{
  // column 1: dimensionless time [-]
  // column 2: time [min]
  // column 3: column position [m]
  // column 4: interstitial velocity [m/s]
  // column 5: total pressure [Pa]
  // column 6: gas temperature [K]
  // column 7: gas temperature time derivative [K/s]
  // column 8: solid temperature [K]
  // column 9: solid temperature time derivative [K/s]
  // column 10: wall temperature [K]
  // column 11: wall temperature time derivative [K/s]
  // column 12: gas (liquid) density [kg/m^3]
  // column 13: pH [-], liquid only
  row.insert(row.end(), {time * timeNormalizationFactor, time / 60.0, columnDistances[grid],
                         interstitialGasVelocity[grid], totalPressure[grid], gasTemperature[grid],
                         gasTemperatureDot[grid], solidTemperature[grid], solidTemperatureDot[grid],
                         wallTemperature[grid], wallTemperatureDot[grid], gasDensity[grid]});
  if (fluidPhase == FluidPhase::Liquid) row.push_back(pH[grid]);
}

void Column::appendComponentOutput(std::vector<double>& row, size_t comp, size_t grid, double time) const
{
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
  // then the surface and pore concentrations [mol/m^3] of each site of surfacePoreTransportSites
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
    normalizedDrivingForce = partialPressure[index] / (totalPressure[grid] * components[comp].initialGasMoleFraction);
  }

  row.insert(row.end(), {time * timeNormalizationFactor, time / 60.0, columnDistances[grid], concentration[index],
                         concentrationDot[index], moleFraction[index], physisorption[index], physisorptionDot[index],
                         totalChemisorption, totalChemisorptionDot, partialPressure[index],
                         equilibriumPhysisorption[index], normalizedDrivingForce, totalEquilibriumChemisorption});
  for (size_t site = 0; site < numberOfChemisorptionSites; ++site)
  {
    if (!componentSiteUsesSurfacePoreTransport(*this, comp, site)) continue;
    const size_t siteIndex = site * componentBlockSize + index;
    row.insert(row.end(), {surfaceConcentration[siteIndex], poreConcentration[siteIndex]});
  }
}

void Column::writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream,
                                 double time) const
{
  // Each file gets one formatted buffer per snapshot: std::print on a stream checks whether the stream is a
  // terminal (an isatty system call) on every call, which cost more than the formatting itself.
  std::string buffer;
  std::vector<double> row;
  auto appendRow = [&]()
  {
    for (size_t value = 0; value < row.size(); ++value)
    {
      if (value > 0) buffer += ' ';
      std::format_to(std::back_inserter(buffer), "{}", row[value]);
    }
    buffer += '\n';
    row.clear();
  };

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    appendColumnOutput(row, grid, time);
    appendRow();
  }
  buffer += "\n\n";
  columnStream.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    buffer.clear();
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      appendComponentOutput(row, comp, grid, time);
      appendRow();
    }
    buffer += "\n\n";
    componentStreams[comp].write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  }
}

std::string Column::repr() const
{
  std::vector<double> dispersivities;
  std::vector<double> taylorArisCoefficients;
  for (const Geometry& geometry : geometries)
  {
    dispersivities.push_back(geometry.dispersion.dispersivity);
    taylorArisCoefficients.push_back(geometry.dispersion.taylorArisCoefficient);
  }
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
      "Axial dispersivities:                  {} [m]\n"
      "Taylor-Aris coefficients d_h^2/192:    {} [m^2]\n"
      "Axial-dispersion reference state:      {}\n"
      "\n\n",
      externalTemperature, columnLength, vectorToString(adsorbentLengths), vectorToString(adsorbentInterfaceLengths),
      vectorToString(adsorbentMixFractions), vectorToString(adsorbentVoidFractions),
      vectorToString(particleDensities), inletPressure, outletPressure, pressureGradient, columnEntranceVelocity,
      vectorToString(dispersivities), vectorToString(taylorArisCoefficients),
      dispersionReferenceTemperature > 0.0
          ? std::format("{} [K], {} [Pa]", dispersionReferenceTemperature, dispersionReferencePressure)
          : std::string{"none (constant molecular part)"});
}

void Column::writeJSON(const std::string& filename) const
{
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Column::writeJSON: cannot open file '" + filename + "'");

  // The state only: the derived fields (velocity, pressure, equilibrium loadings, rates, mixture-prediction caches)
  // follow from it and are recomputed after reading (Breakthrough constructor).
  nlohmann::json j;
  j["numberOfGridPoints"] = numberOfGridPoints;
  j["numberOfComponents"] = numberOfComponents;
  j["maxIsothermTerms"] = maxIsothermTerms;
  j["maxChemisorptionSites"] = maxChemisorptionSites;
  j["surfacePoreTransportEnabled"] = surfacePoreTransportEnabled;

  j["gasTemperature"] = toVector(gasTemperature);
  j["solidTemperature"] = toVector(solidTemperature);
  j["wallTemperature"] = toVector(wallTemperature);

  j["concentration"] = toVector(concentration);
  j["physisorption"] = toVector(physisorption);
  j["chemisorption"] = toVector(chemisorption);
  if (surfacePoreTransportEnabled)
  {
    j["surfaceConcentration"] = toVector(surfaceConcentration);
    j["poreConcentration"] = toVector(poreConcentration);
  }

  out << j.dump(4) << "\n";
}

void Column::readJSON(const std::string& filename)
{
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Column::readJSONFile: cannot open file '" + filename + "'");

  nlohmann::json j;
  in >> j;

  auto requireSizeT = [&](const char* key) -> size_t
  {
    if (!j.contains(key))
      throw std::runtime_error(std::string("Column::readJSON: missing required key '") + key + "'");
    return j.at(key).get<size_t>();
  };

  const size_t fileNgrid = requireSizeT("numberOfGridPoints");
  const size_t fileNcomp = requireSizeT("numberOfComponents");
  const size_t fileMaxIsothermTerms = requireSizeT("maxIsothermTerms");
  const size_t fileMaxChemisorptionSites =
      j.contains("maxChemisorptionSites") ? j.at("maxChemisorptionSites").get<size_t>() : 1;

  if (fileNgrid != numberOfGridPoints)
    throw std::runtime_error("Column::readJSON: numberOfGridPoints mismatch (file " +
                             std::to_string(fileNgrid) + ", column " + std::to_string(numberOfGridPoints) + ")");
  if (fileNcomp != numberOfComponents)
    throw std::runtime_error("Column::readJSON: numberOfComponents mismatch (file " +
                             std::to_string(fileNcomp) + ", column " + std::to_string(numberOfComponents) + ")");
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
    if (!a.is_array())
      throw std::runtime_error(std::string("Column::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("Column::readJSON: size mismatch for '" + std::string(key) + "'");

    dst = a.get<std::vector<double>>();
  };

  auto loadSpanChecked = [&](const char* key, std::span<double> dst)
  {
    if (!j.contains(key)) return;
    const auto& a = j.at(key);
    if (!a.is_array())
      throw std::runtime_error(std::string("Column::readJSON: key '") + key + "' is not an array");

    if (a.size() != dst.size())
      throw std::runtime_error("Column::readJSON: size mismatch for '" + std::string(key) + "'");

    const std::vector<double> tmp = a.get<std::vector<double>>();
    std::copy(tmp.begin(), tmp.end(), dst.begin());
  };

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

  loadVectorChecked("cachedPressure", cachedPressure);
  loadVectorChecked("cachedGrandPotential", cachedGrandPotential);
  loadVectorChecked("cachedChemisorptionPressure", cachedChemisorptionPressure);
  loadVectorChecked("cachedChemisorptionGrandPotential", cachedChemisorptionGrandPotential);
  loadVectorChecked("columnDistances", columnDistances);
  validateColumnDistances(columnDistances, numberOfGridPoints, columnLength);
  computeAdsorbentBlends();
}
