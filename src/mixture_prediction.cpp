#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <print>
#include <string>
#include <vector>
#if __cplusplus >= 201703L && __has_include(<filesystem>)
#include <filesystem>
#elif __cplusplus >= 201703L && __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#else
#include <sys/stat.h>
#endif

#include "mixture_prediction.h"
#include "utils.h"

namespace
{
constexpr double isothermValueCeiling = 1.0e300;
constexpr double isothermExponentClip = 700.0;
constexpr double isothermDenominatorFloor = 1.0e-300;

double finiteNonnegative(double value, double ceiling = isothermValueCeiling)
{
  if (std::isnan(value) || value <= 0.0) return 0.0;
  if (!std::isfinite(value)) return ceiling;
  return std::min(value, ceiling);
}

double safeMultiply(double lhs, double rhs, double ceiling = isothermValueCeiling)
{
  const double a = finiteNonnegative(lhs, ceiling);
  const double b = finiteNonnegative(rhs, ceiling);
  if (a == 0.0 || b == 0.0) return 0.0;

  return std::exp(std::min(std::log(a) + std::log(b), std::log(ceiling)));
}

double safePositivePower(double base, double exponent)
{
  const double b = finiteNonnegative(base);
  if (b == 0.0) return 0.0;

  double e = exponent;
  if (!std::isfinite(e)) e = 1.0;
  e = std::clamp(e, 0.0, 1.0e6);
  return std::exp(std::clamp(e * std::log(b), -isothermExponentClip, isothermExponentClip));
}

double safeDenominator(double value, double floor = isothermDenominatorFloor)
{
  if (std::isnan(value) || value < floor) return floor;
  return value;
}

double sanitizeUnboundedLoading(double loading)
{
  if (std::isnan(loading) || loading <= 0.0) return 0.0;
  if (!std::isfinite(loading)) return isothermValueCeiling;
  return std::min(loading, isothermValueCeiling);
}

double sanitizeBoundedLoading(double loading, double saturationLoading)
{
  const double maximum = finiteNonnegative(saturationLoading);
  if (std::isnan(loading) || loading <= 0.0) return 0.0;
  if (!std::isfinite(loading)) return maximum;
  return std::clamp(loading, 0.0, maximum);
}

double boundedRatioLoading(double saturationLoading, double numerator, double denominator)
{
  const double loading = safeMultiply(saturationLoading, numerator) / safeDenominator(denominator);
  return sanitizeBoundedLoading(loading, saturationLoading);
}

double safeActivitySum(std::span<const double> activities)
{
  double sum = 0.0;
  for (double activity : activities)
  {
    const double value = finiteNonnegative(activity);
    if (value >= isothermValueCeiling - sum) return isothermValueCeiling;
    sum += value;
  }
  return sum;
}

double safePureSiteLoading(const Isotherm& isotherm, double partialPressure, double scale, double pH)
{
  if (!isotherm.enabled()) return 0.0;

  const double pressure = finiteNonnegative(partialPressure);
  const std::vector<double>& parameters = isotherm.parameters;

  switch (isotherm.type)
  {
    case Isotherm::Type::Langmuir:
    {
      const double activity = safeMultiply(safeMultiply(scale, parameters[1]), pressure);
      return boundedRatioLoading(parameters[0], activity, 1.0 + activity);
    }
    case Isotherm::Type::Langmuir_pH:
    {
      const double pHFactor = std::pow(10.0, std::clamp(parameters[2] - pH, -300.0, 300.0));
      const double activity = safeMultiply(safeMultiply(scale, parameters[1]), pressure) / (1.0 + pHFactor);
      return boundedRatioLoading(parameters[0], activity, 1.0 + activity);
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      const double activity = safeMultiply(parameters[1], pressure);
      return sanitizeUnboundedLoading(safeMultiply(parameters[0], pressure) / safeDenominator(1.0 - activity, 1.0e-12));
    }
    case Isotherm::Type::Henry:
      return sanitizeUnboundedLoading(safeMultiply(parameters[0], pressure));
    case Isotherm::Type::Freundlich:
    {
      const double exponent = 1.0 / std::max(finiteNonnegative(parameters[1]), isothermDenominatorFloor);
      return sanitizeUnboundedLoading(safeMultiply(parameters[0], safePositivePower(pressure, exponent)));
    }
    case Isotherm::Type::Sips:
    {
      const double activity = safeMultiply(safeMultiply(scale, parameters[1]), pressure);
      const double exponent = 1.0 / std::max(finiteNonnegative(parameters[2]), isothermDenominatorFloor);
      const double term = safePositivePower(activity, exponent);
      return boundedRatioLoading(parameters[0], term, 1.0 + term);
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      const double pressurePower = safePositivePower(pressure, parameters[2]);
      const double term = safeMultiply(safeMultiply(scale, parameters[1]), pressurePower);
      return boundedRatioLoading(parameters[0], term, 1.0 + term);
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      const double numerator = safeMultiply(parameters[0], pressure);
      const double denominatorTerm = safeMultiply(parameters[1], safePositivePower(pressure, parameters[2]));
      return sanitizeUnboundedLoading(numerator / safeDenominator(1.0 + denominatorTerm));
    }
    case Isotherm::Type::Toth:
    {
      const double activity = safeMultiply(parameters[1], pressure);
      const double exponent = std::max(finiteNonnegative(parameters[2]), isothermDenominatorFloor);
      const double denominator =
          safeDenominator(safePositivePower(1.0 + safePositivePower(activity, exponent), 1.0 / exponent), 1.0e-12);
      return boundedRatioLoading(parameters[0], activity, denominator);
    }
    default:
      return sanitizeUnboundedLoading(isotherm.value(pressure, scale, pH));
  }
}
}  // namespace

bool LangmuirLoadingSorter(Component const& lhs, Component const& rhs)
{
  const bool lhsEnabled = !lhs.isCarrierGas && lhs.isotherm.enabled();
  const bool rhsEnabled = !rhs.isCarrierGas && rhs.isotherm.enabled();
  if (!lhsEnabled) return false;
  if (!rhsEnabled) return true;
  return lhs.isotherm.sites[0].parameters[0] < rhs.isotherm.sites[0].parameters[0];
}

MixturePrediction::MixturePrediction(const InputReader& inputreader)
    : displayName(inputreader.displayName),
      components(inputreader.adsorbentComponents.empty() ? inputreader.components : inputreader.adsorbentComponents[0]),
      sortedComponents(components),
      numberOfComponents(components.size()),
      numberOfSortedComponents(components.size() - inputreader.numberOfCarrierGases),
      numberOfCarrierGases(inputreader.numberOfCarrierGases),
      carrierGasComponent(inputreader.carrierGasComponent),
      predictionMethod(PredictionMethod(inputreader.mixturePredictionMethod)),
      iastMethod(IASTMethod(inputreader.IASTMethod)),
      maxIsothermTerms(inputreader.maxIsothermTerms),
      segregatedSortedComponents(maxIsothermTerms, std::vector<Component>(components)),
      segregatedNumberOfSortedComponents(maxIsothermTerms, 0),
      equilibriumSiteLoadings(maxIsothermTerms * numberOfComponents, 0.0),
      firstExplicitIsothermAlpha(numberOfComponents),
      secondExplicitIsothermAlpha(numberOfComponents),
      explicitIsothermAlphaProduct(numberOfComponents),
      adsorbedMoleFractionsScratch(numberOfComponents),
      hypotheticalPressure(numberOfSortedComponents),
      reducedGrandPotential(numberOfSortedComponents),
      residualVector(numberOfSortedComponents),
      correctionVector(numberOfSortedComponents),
      jacobianMatrix(numberOfSortedComponents * numberOfSortedComponents),
      temperature(inputreader.temperature),
      pressureStart(inputreader.pressureStart),
      pressureEnd(inputreader.pressureEnd),
      numberOfPressurePoints(inputreader.numberOfPressurePoints),
      pressureScale(PressureScale(inputreader.pressureScale))
{
  if (predictionMethod == PredictionMethod::MPD)
  {
    if (!inputreader.mpdSettings.has_value())
    {
      throw std::runtime_error("Error: MPD mixture prediction requires MPDSettings");
    }
    macrostateParticleDistribution.emplace(*inputreader.mpdSettings);
    for (const Component& component : components)
    {
      if (!component.isCarrierGas) mpdComponentIds.push_back(component.id);
    }
    if (mpdComponentIds.size() != macrostateParticleDistribution->rank())
    {
      throw std::runtime_error("Error: MPD rank does not match the number of non-carrier components");
    }
  }
  sortComponents();
}

MixturePrediction::MixturePrediction(std::string _displayName, std::vector<Component> _components,
                                     size_t _numberOfCarrierGases, size_t _carrierGasComponent, double _temperature,
                                     double _pressureStart, double _pressureEnd, size_t _numberOfPressurePoints,
                                     size_t _pressureScale, size_t _predictionMethod, size_t _iastMethod)
    : displayName(_displayName),
      components(_components),
      sortedComponents(components),
      numberOfComponents(components.size()),
      numberOfSortedComponents(components.size() - _numberOfCarrierGases),
      numberOfCarrierGases(_numberOfCarrierGases),
      carrierGasComponent(_carrierGasComponent),
      predictionMethod(PredictionMethod(_predictionMethod)),
      iastMethod(IASTMethod(_iastMethod)),
      firstExplicitIsothermAlpha(numberOfComponents),
      secondExplicitIsothermAlpha(numberOfComponents),
      explicitIsothermAlphaProduct(numberOfComponents),
      adsorbedMoleFractionsScratch(numberOfComponents),
      hypotheticalPressure(numberOfSortedComponents),
      reducedGrandPotential(numberOfSortedComponents),
      residualVector(numberOfSortedComponents),
      correctionVector(numberOfSortedComponents),
      jacobianMatrix(numberOfSortedComponents * numberOfSortedComponents),
      temperature(_temperature),
      pressureStart(_pressureStart),
      pressureEnd(_pressureEnd),
      numberOfPressurePoints(_numberOfPressurePoints),
      pressureScale(PressureScale(_pressureScale))
{
  maxIsothermTerms = 0;
  if (!components.empty())
  {
    std::vector<Component>::iterator maxIsothermTermsIterator =
        std::max_element(_components.begin(), _components.end(), [](Component& lhs, Component& rhs)
                         { return lhs.isotherm.sites.size() < rhs.isotherm.sites.size(); });
    maxIsothermTerms = maxIsothermTermsIterator->isotherm.sites.size();
  }
  segregatedSortedComponents =
      std::vector<std::vector<Component>>(maxIsothermTerms, std::vector<Component>(components));
  segregatedNumberOfSortedComponents.assign(maxIsothermTerms, 0);
  equilibriumSiteLoadings.assign(maxIsothermTerms * numberOfComponents, 0.0);

  sortComponents();
}

std::pair<size_t, size_t> MixturePrediction::predictMixture(std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules,
                                                            std::span<double> cachedPressure,
                                                            std::span<double> cachedGrandPotential,
                                                            double& gasTemperature, double pH,
                                                            DrivingForceInput input)
{
  const double tiny = 1.0e-10;
  std::fill(equilibriumSiteLoadings.begin(), equilibriumSiteLoadings.end(), 0.0);

  if (externalPressure < 0.0)
  {
    printErrorStatus(0.0, 0.0, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
    throw std::runtime_error("Error (IAST): negative total pressure\n");
  }

  double sumYi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    if (!std::isfinite(idealGasMolFractions[i]) ||
        (input == DrivingForceInput::Concentration && idealGasMolFractions[i] < 0.0))
    {
      throw std::runtime_error("Error: mixture-prediction driving forces must be finite and non-negative");
    }
    sumYi += idealGasMolFractions[i];
  }
  if (input == DrivingForceInput::MoleFraction && std::abs(sumYi - 1.0) > 1e-15)
  {
    printErrorStatus(0.0, sumYi, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
    throw std::runtime_error("Error (IAST): sum idealGasMolFractions at IAST start not unity\n");
  }

  if (predictionMethod == PredictionMethod::MPD)
  {
    return computeMPD(idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                      gasTemperature);
  }

  double adsorbingGasFraction = 0.0;
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    adsorbingGasFraction += idealGasMolFractions[sortedComponents[i].id];
  }

  // No component represented by this competitive phase is present.
  if (adsorbingGasFraction < tiny)
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      adsorbedMolFractions[i] = 0.0;
      numberOfMolecules[i] = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  if (numberOfSortedComponents == 1)
  {
    std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
    std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

    const Component& component = sortedComponents.front();
    const size_t comp = component.id;
    const double partialPressure = idealGasMolFractions[comp] * externalPressure;
    numberOfMolecules[comp] = component.isotherm.value(partialPressure, component.scale(gasTemperature), pH);
    adsorbedMolFractions[comp] = numberOfMolecules[comp] > tiny ? 1.0 : 0.0;
    for (size_t site = 0; site < component.isotherm.sites.size(); ++site)
    {
      equilibriumSiteLoadings[site * numberOfComponents + comp] =
          component.isotherm.value(site, partialPressure, component.scale(gasTemperature), pH);
    }
    return std::make_pair(0, 1);
  }

  std::pair<size_t, size_t> result;
  switch (predictionMethod)
  {
    case PredictionMethod::IAST:
    default:
      switch (iastMethod)
      {
        case IASTMethod::FastIAST:
        default:
          result = computeFastIAST(idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                                   cachedPressure, cachedGrandPotential, gasTemperature);
          break;
        case IASTMethod::NestedLoopBisection:
          result = computeIASTNestedLoopBisection(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                                  numberOfMolecules, cachedPressure, cachedGrandPotential,
                                                  gasTemperature);
          break;
      }
      for (size_t i = 0; i < numberOfSortedComponents; ++i)
      {
        const Component& component = sortedComponents[i];
        const size_t comp = component.id;
        const double scale = component.scale(gasTemperature);
        const double pureTotal = component.isotherm.value(cachedPressure[comp], scale);
        if (pureTotal <= tiny) continue;
        for (size_t site = 0; site < component.isotherm.sites.size(); ++site)
        {
          equilibriumSiteLoadings[site * numberOfComponents + comp] =
              numberOfMolecules[comp] * component.isotherm.value(site, cachedPressure[comp], scale) / pureTotal;
        }
      }
      return result;
    case PredictionMethod::SIAST:
      switch (iastMethod)
      {
        case IASTMethod::FastIAST:
        default:
          return computeFastSIAST(idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                                  cachedPressure, cachedGrandPotential, gasTemperature);
        case IASTMethod::NestedLoopBisection:
          return computeSIASTNestedLoopBisection(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                                 numberOfMolecules, cachedPressure, cachedGrandPotential,
                                                 gasTemperature);
      }
    case PredictionMethod::EI:
      result = computeExplicitIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                       numberOfMolecules, gasTemperature);
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumSiteLoadings[comp] = numberOfMolecules[comp];
      }
      return result;
    case PredictionMethod::SEI:
      return computeSegratedExplicitIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                             numberOfMolecules, gasTemperature);
    case PredictionMethod::SCI:
      return computeSegregatedCompetitiveIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                                  numberOfMolecules, gasTemperature);
    case PredictionMethod::SPI:
      return computeSegregatedPureIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                           numberOfMolecules, gasTemperature, pH);
  }
}

void MixturePrediction::predictPureComponentLoadings(double externalPressure,
                                                     std::span<double> pureComponentLoadings,
                                                     double gasTemperature)
{
  if (pureComponentLoadings.size() != numberOfComponents)
  {
    throw std::runtime_error("Error: pure-component loading output size does not match the component count");
  }

  std::fill(pureComponentLoadings.begin(), pureComponentLoadings.end(), 0.0);
  if (predictionMethod != PredictionMethod::MPD)
  {
    for (const Component& component : components)
    {
      if (component.id >= pureComponentLoadings.size())
      {
        throw std::runtime_error("Error: component index is outside the pure-component loading vector");
      }
      pureComponentLoadings[component.id] =
          component.isotherm.value(externalPressure, component.scale(gasTemperature));
    }
    return;
  }

  std::vector<double> oneHotMoleFractions(numberOfComponents, 0.0);
  std::vector<double> adsorbedMolFractions(numberOfComponents, 0.0);
  std::vector<double> numberOfMolecules(numberOfComponents, 0.0);
  std::vector<double> cachedPressure(numberOfComponents * maxIsothermTerms, 0.0);
  std::vector<double> cachedGrandPotential(maxIsothermTerms, 0.0);

  for (const Component& component : components)
  {
    if (component.isCarrierGas) continue;
    if (component.id >= numberOfComponents)
    {
      throw std::runtime_error("Error: component index is outside the one-hot gas-composition vector");
    }

    oneHotMoleFractions[component.id] = 1.0;
    double pureTemperature = gasTemperature;
    predictMixture(oneHotMoleFractions, externalPressure, adsorbedMolFractions, numberOfMolecules, cachedPressure,
                   cachedGrandPotential, pureTemperature);
    pureComponentLoadings[component.id] = numberOfMolecules[component.id];
    oneHotMoleFractions[component.id] = 0.0;
  }
}

std::pair<size_t, size_t> MixturePrediction::computeMPD(std::span<const double> idealGasMolFractions,
                                                         double fugacity,
                                                         std::span<double> adsorbedMolFractions,
                                                         std::span<double> numberOfMolecules,
                                                         double gasTemperature)
{
  if (!macrostateParticleDistribution.has_value())
  {
    throw std::runtime_error("Error: MPD mixture prediction has no loaded particle distribution");
  }

  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);
  std::fill(equilibriumSiteLoadings.begin(), equilibriumSiteLoadings.end(), 0.0);

  std::vector<double> mpdGasMoleFractions;
  mpdGasMoleFractions.reserve(mpdComponentIds.size());
  for (size_t component : mpdComponentIds)
  {
    if (component >= idealGasMolFractions.size())
    {
      throw std::runtime_error("Error: MPD component index is outside the gas-composition vector");
    }
    mpdGasMoleFractions.push_back(idealGasMolFractions[component]);
  }

  const std::vector<double> meanParticleNumbers =
      macrostateParticleDistribution->meanParticleNumbers(mpdGasMoleFractions, fugacity, gasTemperature);
  constexpr double avogadroConstant = 6.02214076e23;  // mol^-1, exact SI definition
  const double loadingConversion =
      1.0 / (avogadroConstant * macrostateParticleDistribution->referenceFrameworkMass());

  double totalLoading = 0.0;
  for (size_t dimension = 0; dimension < mpdComponentIds.size(); ++dimension)
  {
    const size_t component = mpdComponentIds[dimension];
    const double loading = meanParticleNumbers[dimension] * loadingConversion;
    numberOfMolecules[component] = loading;
    if (component < equilibriumSiteLoadings.size()) equilibriumSiteLoadings[component] = loading;
    totalLoading += loading;
  }
  if (totalLoading > 0.0)
  {
    for (size_t component : mpdComponentIds)
    {
      adsorbedMolFractions[component] = numberOfMolecules[component] / totalLoading;
    }
  }

  return {0, 1};
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastIAST(std::span<const double> idealGasMolFractions,
                                                             const double& externalPressure,
                                                             std::span<double> adsorbedMolFractions,
                                                             std::span<double> numberOfMolecules,
                                                             std::span<double> cachedPressure,
                                                             std::span<double> cachedGrandPotential,
                                                             double& gasTemperature)
{
  const double tiny = 1.0e-13;

  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  size_t numberOfIASTSteps = 0;

  std::fill(hypotheticalPressure.begin(), hypotheticalPressure.end(), 0.0);
  std::fill(residualVector.begin(), residualVector.end(), 0.0);
  std::fill(correctionVector.begin(), correctionVector.end(), 0.0);
  std::fill(jacobianMatrix.begin(), jacobianMatrix.end(), 0.0);

  std::vector<double> componentScale(numberOfComponents);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    componentScale[i] = sortedComponents[i].scale(gasTemperature);
  }

  if (cachedGrandPotential[0] > 0.0)
  {
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      hypotheticalPressure[i] = cachedPressure[sortedComponents[i].id];
    }
  }
  else
  {
    double initial_psi = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      double temp_psi = idealGasMolFractions[sortedComponents[i].id] *
                        sortedComponents[i].isotherm.psiForPressure(externalPressure, componentScale[i]);
      initial_psi += temp_psi;
    }
    cachedGrandPotential[0] = initial_psi;

    double cachevalue = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      hypotheticalPressure[i] =
          1.0 / sortedComponents[i].isotherm.inversePressureForPsi(initial_psi, cachevalue, componentScale[i]);
    }
  }

  double error = 1.0;
  double sum_xi = 0.0;
  do
  {
    // compute residualVector
    for (size_t i = 0; i < numberOfSortedComponents - 1; ++i)
    {
      residualVector[i] =
          sortedComponents[i].isotherm.psiForPressure(hypotheticalPressure[i], componentScale[i]) -
          sortedComponents[numberOfSortedComponents - 1].isotherm.psiForPressure(
              hypotheticalPressure[numberOfSortedComponents - 1], componentScale[numberOfSortedComponents - 1]);
    }

    residualVector[numberOfSortedComponents - 1] = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; i++)
    {
      residualVector[numberOfSortedComponents - 1] +=
          idealGasMolFractions[sortedComponents[i].id] * externalPressure / hypotheticalPressure[i];
    }
    residualVector[numberOfSortedComponents - 1] -= 1.0;

    // compute Jacobian matrix jacobianMatrix
    for (size_t i = 0; i < numberOfSortedComponents - 1; i++)
    {
      jacobianMatrix[i + i * numberOfSortedComponents] =
          sortedComponents[i].isotherm.value(hypotheticalPressure[i], componentScale[i]) / hypotheticalPressure[i];
    }
    for (size_t i = 0; i < numberOfSortedComponents - 1; i++)
    {
      jacobianMatrix[i + (numberOfSortedComponents - 1) * numberOfSortedComponents] =
          -sortedComponents[numberOfSortedComponents - 1].isotherm.value(
              hypotheticalPressure[numberOfSortedComponents - 1], componentScale[numberOfSortedComponents - 1]) /
          hypotheticalPressure[numberOfSortedComponents - 1];
    }
    for (size_t i = 0; i < numberOfSortedComponents; i++)
    {
      jacobianMatrix[(numberOfSortedComponents - 1) + i * numberOfSortedComponents] =
          -idealGasMolFractions[sortedComponents[i].id] * externalPressure /
          (hypotheticalPressure[i] * hypotheticalPressure[i]);
    }

    // corrections
    for (size_t i = 0; i < numberOfSortedComponents - 1; i++)
    {
      jacobianMatrix[(numberOfSortedComponents - 1) + (numberOfSortedComponents - 1) * numberOfSortedComponents] -=
          jacobianMatrix[(numberOfSortedComponents - 1) + i * numberOfSortedComponents] *
          jacobianMatrix[i + (numberOfSortedComponents - 1) * numberOfSortedComponents] /
          jacobianMatrix[i + i * numberOfSortedComponents];
      residualVector[numberOfSortedComponents - 1] -=
          jacobianMatrix[(numberOfSortedComponents - 1) + i * numberOfSortedComponents] * residualVector[i] /
          jacobianMatrix[i + i * numberOfSortedComponents];
    }

    // compute correctionVector
    correctionVector[numberOfSortedComponents - 1] =
        residualVector[numberOfSortedComponents - 1] /
        jacobianMatrix[(numberOfSortedComponents - 1) + (numberOfSortedComponents - 1) * numberOfSortedComponents];

    // trick to loop downward from numberOfSortedComponents - 2 to and including zero (still using size_t as index)
    for (size_t i = numberOfSortedComponents - 1; i-- != 0;)
    {
      correctionVector[i] =
          (residualVector[i] - correctionVector[numberOfSortedComponents - 1] *
                                   jacobianMatrix[i + (numberOfSortedComponents - 1) * numberOfSortedComponents]) /
          jacobianMatrix[i + i * numberOfSortedComponents];
    }

    // update hypotheticalPressure
    for (size_t i = 0; i < numberOfSortedComponents; i++)
    {
      double newvalue = hypotheticalPressure[i] - correctionVector[i];
      if (newvalue > 0.0)
        hypotheticalPressure[i] = newvalue;
      else
      {
        hypotheticalPressure[i] = 0.5 * hypotheticalPressure[i];
      }
    }

    // compute error in reducedGrandPotential's
    for (size_t i = 0; i < numberOfSortedComponents; i++)
    {
      reducedGrandPotential[i] =
          sortedComponents[i].isotherm.psiForPressure(hypotheticalPressure[i], componentScale[i]);
    }

    sum_xi = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      sum_xi +=
          idealGasMolFractions[sortedComponents[i].id] * externalPressure / std::max(hypotheticalPressure[i], 1e-15);
    }

    double avg = std::accumulate(std::begin(reducedGrandPotential), std::end(reducedGrandPotential), 0.0) /
                 static_cast<double>(reducedGrandPotential.size());

    double accum = 0.0;
    std::for_each(std::begin(reducedGrandPotential), std::end(reducedGrandPotential),
                  [&](const double d) { accum += (d - avg) * (d - avg); });

    error = std::sqrt(accum / static_cast<double>(reducedGrandPotential.size() - 1));

    numberOfIASTSteps++;
  } while (!(((error < tiny) && (std::fabs(sum_xi - 1.0) < 1e-10)) || (numberOfIASTSteps >= 50)));

  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    cachedPressure[sortedComponents[i].id] = hypotheticalPressure[i];
  }

  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    adsorbedMolFractions[sortedComponents[i].id] =
        idealGasMolFractions[sortedComponents[i].id] * externalPressure / std::max(hypotheticalPressure[i], 1e-15);
  }
  double sum = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sum += adsorbedMolFractions[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] /= sum;
  }

  double inverse_q_total = 0.0;
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    inverse_q_total += adsorbedMolFractions[sortedComponents[i].id] /
                       sortedComponents[i].isotherm.value(hypotheticalPressure[i], componentScale[i]);
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    numberOfMolecules[i] = adsorbedMolFractions[i] / inverse_q_total;
  }
  return std::make_pair(numberOfIASTSteps, 1);
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastSIAST(std::span<const double> idealGasMolFractions,
                                                              const double& externalPressure,
                                                              std::span<double> adsorbedMolFractions,
                                                              std::span<double> numberOfMolecules,
                                                              std::span<double> cachedPressure,
                                                              std::span<double> cachedGrandPotential,
                                                              double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  std::vector<double> previous(numberOfComponents, 0.0);
  for (size_t site = 0; site < maxIsothermTerms; ++site)
  {
    std::copy(numberOfMolecules.begin(), numberOfMolecules.end(), previous.begin());
    const size_t activeComponents = segregatedNumberOfSortedComponents[site];
    if (activeComponents == 1)
    {
      const Component& component = segregatedSortedComponents[site][0];
      const double partialPressure = idealGasMolFractions[component.id] * externalPressure;
      numberOfMolecules[component.id] +=
          component.isotherm.value(partialPressure, component.scale(gasTemperature));
      cachedPressure[site * numberOfComponents + component.id] = partialPressure;
      acc += std::make_pair<size_t, size_t>(0, 1);
    }
    else if (activeComponents > 1)
    {
      acc += computeFastSIAST(site, idealGasMolFractions, externalPressure, adsorbedMolFractions,
                              numberOfMolecules, cachedPressure, cachedGrandPotential, gasTemperature);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumSiteLoadings[site * numberOfComponents + comp] = numberOfMolecules[comp] - previous[comp];
    }
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = N > 0.0 ? numberOfMolecules[i] / N : 0.0;
  }

  return acc;
}

// computes IAST per term
// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastSIAST(size_t site, std::span<const double> idealGasMolFractions,
                                                              const double& externalPressure,
                                                              std::span<double> adsorbedMolFractions,
                                                              std::span<double> numberOfMolecules,
                                                              std::span<double> cachedPressure,
                                                              std::span<double> cachedGrandPotential,
                                                              double& gasTemperature)
{
  const double tiny = 1.0e-13;
  const std::vector<Component>& siteComponents = segregatedSortedComponents[site];
  const size_t activeComponents = segregatedNumberOfSortedComponents[site];

  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);

  size_t numberOfIASTSteps = 0;

  std::fill(hypotheticalPressure.begin(), hypotheticalPressure.end(), 0.0);
  std::fill(residualVector.begin(), residualVector.end(), 0.0);
  std::fill(correctionVector.begin(), correctionVector.end(), 0.0);
  std::fill(jacobianMatrix.begin(), jacobianMatrix.end(), 0.0);

  std::vector<double> componentScale(numberOfComponents);
  for (size_t i = 0; i < activeComponents; ++i)
  {
    componentScale[i] = siteComponents[i].scale(gasTemperature);
  }

  if (cachedGrandPotential[site] > tiny)
  {
    for (size_t i = 0; i < activeComponents; ++i)
    {
      hypotheticalPressure[i] = cachedPressure[siteComponents[i].id + site * numberOfComponents];
    }
  }
  else
  {
    double initial_psi = 0.0;
    for (size_t i = 0; i < activeComponents; ++i)
    {
      double temp_psi = idealGasMolFractions[siteComponents[i].id] *
                        siteComponents[i].isotherm.psiForPressure(externalPressure, componentScale[i]);
      initial_psi += temp_psi;
    }
    cachedGrandPotential[site] = initial_psi;

    double cachevalue = 0.0;
    for (size_t i = 0; i < activeComponents; ++i)
    {
      hypotheticalPressure[i] =
          1.0 / siteComponents[i].isotherm.inversePressureForPsi(initial_psi, cachevalue, componentScale[i]);
    }
  }

  double error = 1.0;
  double sum_xi = 1.0;
  do
  {
    // compute residualVector
    for (size_t i = 0; i < activeComponents - 1; ++i)
    {
      residualVector[i] =
          siteComponents[i].isotherm.psiForPressure(hypotheticalPressure[i], componentScale[i]) -
          siteComponents[activeComponents - 1].isotherm.psiForPressure(
              hypotheticalPressure[activeComponents - 1], componentScale[activeComponents - 1]);
    }

    residualVector[activeComponents - 1] = 0.0;
    for (size_t i = 0; i < activeComponents; i++)
    {
      residualVector[activeComponents - 1] +=
          idealGasMolFractions[siteComponents[i].id] * externalPressure / hypotheticalPressure[i];
    }
    residualVector[activeComponents - 1] -= 1.0;

    // compute Jacobian matrix jacobianMatrix
    for (size_t i = 0; i < activeComponents - 1; i++)
    {
      jacobianMatrix[i + i * activeComponents] =
          siteComponents[i].isotherm.value(hypotheticalPressure[i], componentScale[i]) /
          hypotheticalPressure[i];
    }
    for (size_t i = 0; i < activeComponents - 1; i++)
    {
      jacobianMatrix[i + (activeComponents - 1) * activeComponents] =
          -siteComponents[activeComponents - 1].isotherm.value(
              hypotheticalPressure[activeComponents - 1], componentScale[activeComponents - 1]) /
          hypotheticalPressure[activeComponents - 1];
    }
    for (size_t i = 0; i < activeComponents; i++)
    {
      jacobianMatrix[(activeComponents - 1) + i * activeComponents] =
          -idealGasMolFractions[siteComponents[i].id] * externalPressure /
          (hypotheticalPressure[i] * hypotheticalPressure[i]);
    }

    // corrections
    for (size_t i = 0; i < activeComponents - 1; i++)
    {
      jacobianMatrix[(activeComponents - 1) + (activeComponents - 1) * activeComponents] -=
          jacobianMatrix[(activeComponents - 1) + i * activeComponents] *
          jacobianMatrix[i + (activeComponents - 1) * activeComponents] /
          jacobianMatrix[i + i * activeComponents];
      residualVector[activeComponents - 1] -=
          jacobianMatrix[(activeComponents - 1) + i * activeComponents] * residualVector[i] /
          jacobianMatrix[i + i * activeComponents];
    }

    // compute correctionVector
    correctionVector[activeComponents - 1] =
        residualVector[activeComponents - 1] /
        jacobianMatrix[(activeComponents - 1) + (activeComponents - 1) * activeComponents];

    // trick to loop downward from numberOfSortedComponents - 2 to and including zero (still using size_t as index)
    for (size_t i = activeComponents - 1; i-- != 0;)
    {
      correctionVector[i] =
          (residualVector[i] - correctionVector[activeComponents - 1] *
                                   jacobianMatrix[i + (activeComponents - 1) * activeComponents]) /
          jacobianMatrix[i + i * activeComponents];
    }

    // update hypotheticalPressure
    for (size_t i = 0; i < activeComponents; i++)
    {
      double newvalue = hypotheticalPressure[i] - correctionVector[i];
      if (newvalue > 0.0)
        hypotheticalPressure[i] = newvalue;
      else
      {
        hypotheticalPressure[i] = 0.5 * hypotheticalPressure[i];
      }
    }

    // compute error in reducedGrandPotential's
    for (size_t i = 0; i < activeComponents; i++)
    {
      reducedGrandPotential[i] = siteComponents[i].isotherm.psiForPressure(hypotheticalPressure[i], componentScale[i]);
    }

    sum_xi = 0.0;
    for (size_t i = 0; i < activeComponents; ++i)
    {
      sum_xi +=
          idealGasMolFractions[siteComponents[i].id] * externalPressure / std::max(hypotheticalPressure[i], 1e-15);
    }

    double avg = std::accumulate(reducedGrandPotential.begin(),
                                 reducedGrandPotential.begin() + static_cast<std::ptrdiff_t>(activeComponents), 0.0) /
                 static_cast<double>(activeComponents);

    double accum = 0.0;
    std::for_each(reducedGrandPotential.begin(),
                  reducedGrandPotential.begin() + static_cast<std::ptrdiff_t>(activeComponents),
                  [&](const double d) { accum += (d - avg) * (d - avg); });

    error = std::sqrt(accum / static_cast<double>(activeComponents - 1));

    numberOfIASTSteps++;
  } while (!(((error < tiny) && (std::fabs(sum_xi - 1.0) < 1e-10)) || (numberOfIASTSteps >= 50)));

  for (size_t i = 0; i < activeComponents; ++i)
  {
    cachedPressure[siteComponents[i].id + site * numberOfComponents] = hypotheticalPressure[i];
  }

  for (size_t i = 0; i < activeComponents; ++i)
  {
    adsorbedMolFractions[siteComponents[i].id] =
        idealGasMolFractions[siteComponents[i].id] * externalPressure / std::max(hypotheticalPressure[i], 1e-15);
  }
  double sum = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sum += adsorbedMolFractions[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] /= sum;
  }

  double inverse_q_total = 0.0;
  for (size_t i = 0; i < activeComponents; ++i)
  {
    inverse_q_total += adsorbedMolFractions[siteComponents[i].id] /
                       siteComponents[i].isotherm.value(hypotheticalPressure[i], componentScale[i]);
  }
  for (size_t i = 0; i < activeComponents; ++i)
  {
    const size_t comp = siteComponents[i].id;
    numberOfMolecules[comp] += adsorbedMolFractions[comp] / inverse_q_total;
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeIASTNestedLoopBisection(
    std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, double& gasTemperature)
{
  const double tiny = 1.0e-15;

  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::vector<double> componentScale(numberOfSortedComponents);
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    componentScale[i] = sortedComponents[i].scale(gasTemperature);
  }

  double initial_psi = 0.0;
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    initial_psi += idealGasMolFractions[sortedComponents[i].id] *
                   sortedComponents[i].isotherm.psiForPressure(externalPressure, componentScale[i]);
  }

  if (initial_psi < tiny)
  {
    // nothing is adsorbing
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      adsorbedMolFractions[i] = 0.0;
      numberOfMolecules[i] = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedGrandPotential[0] > tiny)
  {
    initial_psi = cachedGrandPotential[0];
  }
  auto sumAdsorbedFractions = [&](double psi)
  {
    double sum = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      const size_t comp = sortedComponents[i].id;
      sum += idealGasMolFractions[comp] * externalPressure *
             sortedComponents[i].isotherm.inversePressureForPsi(psi, cachedPressure[comp], componentScale[i]);
    }
    return sum;
  };

  double sumXi = sumAdsorbedFractions(initial_psi);

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = sumAdsorbedFractions(right_bracket);
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (sumXi > 1.0);
  }
  else
  {
    // Make an initial estimate for the reduced grandpotential when the
    // sum of the molefractions is larger than 1
    do
    {
      left_bracket *= 0.5;

      sumXi = sumAdsorbedFractions(left_bracket);
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (sumXi < 1.0);
  }

  // bisection algorithm
  size_t numberOfIASTSteps = 0;
  do
  {
    psi_value = 0.5 * (left_bracket + right_bracket);

    sumXi = sumAdsorbedFractions(psi_value);

    if (sumXi > 1.0)
    {
      left_bracket = psi_value;
    }
    else
    {
      right_bracket = psi_value;
    }

    ++numberOfIASTSteps;
    if (numberOfIASTSteps > 100000)
    {
      throw std::runtime_error("Error (IAST bisection): NO convergence\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);  // convergence test

  psi_value = 0.5 * (left_bracket + right_bracket);

  sumXi = sumAdsorbedFractions(psi_value);

  // cache the value of reducedGrandPotential for subsequent use
  cachedGrandPotential[0] = psi_value;

  double inverse_q_total = 0.0;
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    const size_t comp = sortedComponents[i].id;
    double ip = sortedComponents[i].isotherm.inversePressureForPsi(
        psi_value, cachedPressure[comp], componentScale[i]);
    cachedPressure[comp] = 1.0 / ip;
    adsorbedMolFractions[comp] = idealGasMolFractions[comp] * externalPressure * ip / sumXi;

    if (adsorbedMolFractions[comp] > tiny)
    {
      inverse_q_total += adsorbedMolFractions[comp] /
                         sortedComponents[i].isotherm.value(1.0 / ip, componentScale[i]);
    }
  }

  if (inverse_q_total > 0.0)
  {
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      const size_t comp = sortedComponents[i].id;
      numberOfMolecules[comp] = adsorbedMolFractions[comp] / inverse_q_total;
    }
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeSIASTNestedLoopBisection(
    std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  std::vector<double> previous(numberOfComponents, 0.0);
  for (size_t site = 0; site < maxIsothermTerms; ++site)
  {
    std::copy(numberOfMolecules.begin(), numberOfMolecules.end(), previous.begin());
    const size_t activeComponents = segregatedNumberOfSortedComponents[site];
    if (activeComponents == 1)
    {
      const Component& component = segregatedSortedComponents[site][0];
      const double partialPressure = idealGasMolFractions[component.id] * externalPressure;
      numberOfMolecules[component.id] +=
          component.isotherm.value(partialPressure, component.scale(gasTemperature));
      cachedPressure[site * numberOfComponents + component.id] = partialPressure;
      acc += std::make_pair<size_t, size_t>(0, 1);
    }
    else if (activeComponents > 1)
    {
      acc += computeSIASTNestedLoopBisection(site, idealGasMolFractions, externalPressure,
                                             adsorbedMolFractions, numberOfMolecules, cachedPressure,
                                             cachedGrandPotential, gasTemperature);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumSiteLoadings[site * numberOfComponents + comp] = numberOfMolecules[comp] - previous[comp];
    }
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = N > 0.0 ? numberOfMolecules[i] / N : 0.0;
  }

  return acc;
}

// computes IAST per term
// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeSIASTNestedLoopBisection(
    size_t site, std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, double& gasTemperature)
{
  const double tiny = 1.0e-15;
  const std::vector<Component>& siteComponents = segregatedSortedComponents[site];
  const size_t activeComponents = segregatedNumberOfSortedComponents[site];

  std::vector<double> componentScale(activeComponents);
  for (size_t i = 0; i < activeComponents; ++i)
  {
    componentScale[i] = siteComponents[i].scale(gasTemperature);
  }

  double initial_psi = 0.0;
  for (size_t i = 0; i < activeComponents; ++i)
  {
    initial_psi += idealGasMolFractions[siteComponents[i].id] *
                   siteComponents[i].isotherm.psiForPressure(externalPressure, componentScale[i]);
  }

  if (initial_psi < tiny)
  {
    // nothing is adsorbing
    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedGrandPotential[site] > tiny)
  {
    initial_psi = cachedGrandPotential[site];
  }
  auto sumAdsorbedFractions = [&](double psi)
  {
    double sum = 0.0;
    for (size_t i = 0; i < activeComponents; ++i)
    {
      const size_t comp = siteComponents[i].id;
      sum += idealGasMolFractions[comp] * externalPressure *
             siteComponents[i].isotherm.inversePressureForPsi(
                 psi, cachedPressure[comp + numberOfComponents * site], componentScale[i]);
    }
    return sum;
  };

  double sumXi = sumAdsorbedFractions(initial_psi);

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = sumAdsorbedFractions(right_bracket);
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        printErrorStatus(0.0, sumXi, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (sumXi > 1.0);
  }
  else
  {
    // Make an initial estimate for the reduced grandpotential when the
    // sum of the molefractions is larger than 1
    do
    {
      left_bracket *= 0.5;

      sumXi = sumAdsorbedFractions(left_bracket);
      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        printErrorStatus(0.0, sumXi, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
        throw std::runtime_error("Error (IAST bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (sumXi < 1.0);
  }

  // bisection algorithm
  size_t numberOfIASTSteps = 0;
  do
  {
    psi_value = 0.5 * (left_bracket + right_bracket);

    sumXi = sumAdsorbedFractions(psi_value);

    if (sumXi > 1.0)
    {
      left_bracket = psi_value;
    }
    else
    {
      right_bracket = psi_value;
    }

    ++numberOfIASTSteps;
    if (numberOfIASTSteps > 100000)
    {
      throw std::runtime_error("Error (IAST bisection): NO convergence\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);  // convergence test

  psi_value = 0.5 * (left_bracket + right_bracket);

  // cache the value of reducedGrandPotential for subsequent use
  cachedGrandPotential[site] = psi_value;

  sumXi = sumAdsorbedFractions(psi_value);
  double inverse_q_total = 0.0;
  for (size_t i = 0; i < activeComponents; ++i)
  {
    const size_t comp = siteComponents[i].id;
    double ip = siteComponents[i].isotherm.inversePressureForPsi(
        psi_value, cachedPressure[comp + numberOfComponents * site], componentScale[i]);
    cachedPressure[comp + numberOfComponents * site] = 1.0 / ip;
    adsorbedMolFractions[comp] = idealGasMolFractions[comp] * externalPressure * ip / sumXi;

    if (adsorbedMolFractions[comp] > tiny)
    {
      inverse_q_total +=
          adsorbedMolFractions[comp] / siteComponents[i].isotherm.value(1.0 / ip, componentScale[i]);
    }
  }

  if (inverse_q_total > 0.0)
  {
    for (size_t i = 0; i < activeComponents; ++i)
    {
      const size_t comp = siteComponents[i].id;
      numberOfMolecules[comp] += adsorbedMolFractions[comp] / inverse_q_total;
    }
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// solve the mixed-langmuir equations derived by Assche et al.
// T. R. Van Assche, residualVector.V. Baron, and J. F. Denayer
// An explicit multicomponent adsorption isotherm model:
// Accounting for the size-effect for components with Langmuir adsorption behavior.
// Adsorption, 24(6), 517-530 (2018)

// An explicit multicomponent adsorption isotherm model: accounting for the
// size-effect for components with Langmuir adsorption behavior

// In the input file molecules must be added in the following order:
// Largest molecule should be the first component or the component with
// smallest saturation(Nimax) loading should be the first component
// Last component is the carrier gas

// At present, only single site isotherms are considered for pure components

std::pair<size_t, size_t> MixturePrediction::computeExplicitIsotherm(std::span<const double> idealGasMolFractions,
                                                                     const double& externalPressure,
                                                                     std::span<double> adsorbedMolFractions,
                                                                     std::span<double> numberOfMolecules,
                                                                     double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::vector<double> componentScale(numberOfSortedComponents);
  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    componentScale[i] = sortedComponents[i].scale(gasTemperature);
  }

  adsorbedMoleFractionsScratch[0] = 1.0;
  for (size_t i = 1; i < numberOfSortedComponents; ++i)
  {
    adsorbedMoleFractionsScratch[i] =
        sortedComponents[i].isotherm.sites[0].parameters[0] / sortedComponents[i - 1].isotherm.sites[0].parameters[0];
  }

  const size_t last = numberOfSortedComponents - 1;
  double b = componentScale[last] * sortedComponents[last].isotherm.sites[0].parameters[1];
  firstExplicitIsothermAlpha[last] =
      std::pow((1.0 + b * idealGasMolFractions[sortedComponents[last].id] * externalPressure),
               adsorbedMoleFractionsScratch[last]);
  secondExplicitIsothermAlpha[last] = 1.0 + b * idealGasMolFractions[sortedComponents[last].id] * externalPressure;
  for (size_t i = numberOfSortedComponents - 2; i > 0; i--)
  {
    b = componentScale[i] * sortedComponents[i].isotherm.sites[0].parameters[1];
    firstExplicitIsothermAlpha[i] = std::pow(
        (firstExplicitIsothermAlpha[i + 1] + b * idealGasMolFractions[sortedComponents[i].id] * externalPressure),
        adsorbedMoleFractionsScratch[i]);
    secondExplicitIsothermAlpha[i] =
        firstExplicitIsothermAlpha[i + 1] + b * idealGasMolFractions[sortedComponents[i].id] * externalPressure;
  }

  b = componentScale[0] * sortedComponents[0].isotherm.sites[0].parameters[1];
  firstExplicitIsothermAlpha[0] =
      firstExplicitIsothermAlpha[1] + b * idealGasMolFractions[sortedComponents[0].id] * externalPressure;
  secondExplicitIsothermAlpha[0] =
      firstExplicitIsothermAlpha[1] + b * idealGasMolFractions[sortedComponents[0].id] * externalPressure;

  double beta = secondExplicitIsothermAlpha[0];

  explicitIsothermAlphaProduct[0] = 1.0;
  for (size_t i = 1; i < numberOfSortedComponents; ++i)
  {
    explicitIsothermAlphaProduct[i] =
        (firstExplicitIsothermAlpha[i] / secondExplicitIsothermAlpha[i]) * explicitIsothermAlphaProduct[i - 1];
  }

  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    size_t index = sortedComponents[i].id;
    b = componentScale[i] * sortedComponents[i].isotherm.sites[0].parameters[1];
    numberOfMolecules[index] = sortedComponents[i].isotherm.sites[0].parameters[0] * b * idealGasMolFractions[index] *
                               externalPressure * explicitIsothermAlphaProduct[i] / beta;
  }
  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = N > 0.0 ? numberOfMolecules[i] / N : 0.0;
  }

  return std::make_pair(1, 1);
}

std::pair<size_t, size_t> MixturePrediction::computeSegratedExplicitIsotherm(
    std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  std::vector<double> previous(numberOfComponents, 0.0);
  for (size_t site = 0; site < maxIsothermTerms; ++site)
  {
    std::copy(numberOfMolecules.begin(), numberOfMolecules.end(), previous.begin());
    const size_t activeComponents = segregatedNumberOfSortedComponents[site];
    if (activeComponents == 1)
    {
      const Component& component = segregatedSortedComponents[site][0];
      numberOfMolecules[component.id] += component.isotherm.value(
          idealGasMolFractions[component.id] * externalPressure, component.scale(gasTemperature));
      acc += std::make_pair<size_t, size_t>(0, 1);
    }
    else if (activeComponents > 1)
    {
      acc += computeSegratedExplicitIsotherm(site, idealGasMolFractions, externalPressure,
                                             adsorbedMolFractions, numberOfMolecules, gasTemperature);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumSiteLoadings[site * numberOfComponents + comp] = numberOfMolecules[comp] - previous[comp];
    }
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = N > 0.0 ? numberOfMolecules[i] / N : 0.0;
  }

  return acc;
}

std::pair<size_t, size_t> MixturePrediction::computeSegratedExplicitIsotherm(
    size_t site, std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> /*adsorbedMolFractions*/, std::span<double> numberOfMolecules, double& gasTemperature)
{
  const std::vector<Component>& siteComponents = segregatedSortedComponents[site];
  const size_t activeComponents = segregatedNumberOfSortedComponents[site];
  std::vector<double> componentScale(activeComponents);
  for (size_t i = 0; i < activeComponents; ++i)
  {
    componentScale[i] = siteComponents[i].scale(gasTemperature);
  }

  adsorbedMoleFractionsScratch[0] = 1.0;
  for (size_t i = 1; i < activeComponents; ++i)
  {
    adsorbedMoleFractionsScratch[i] = siteComponents[i].isotherm.sites[0].parameters[0] /
                                      siteComponents[i - 1].isotherm.sites[0].parameters[0];
  }

  const size_t last = activeComponents - 1;
  double b = componentScale[last] * siteComponents[last].isotherm.sites[0].parameters[1];
  firstExplicitIsothermAlpha[last] =
      std::pow(1.0 + b * idealGasMolFractions[siteComponents[last].id] * externalPressure,
               adsorbedMoleFractionsScratch[last]);
  secondExplicitIsothermAlpha[last] =
      1.0 + b * idealGasMolFractions[siteComponents[last].id] * externalPressure;
  for (size_t i = activeComponents - 2; i > 0; --i)
  {
    b = componentScale[i] * siteComponents[i].isotherm.sites[0].parameters[1];
    firstExplicitIsothermAlpha[i] =
        std::pow((firstExplicitIsothermAlpha[i + 1] +
                  b * idealGasMolFractions[siteComponents[i].id] * externalPressure),
                 adsorbedMoleFractionsScratch[i]);
    secondExplicitIsothermAlpha[i] =
        firstExplicitIsothermAlpha[i + 1] +
        b * idealGasMolFractions[siteComponents[i].id] * externalPressure;
  }

  b = componentScale[0] * siteComponents[0].isotherm.sites[0].parameters[1];
  firstExplicitIsothermAlpha[0] = firstExplicitIsothermAlpha[1] +
                                  b * idealGasMolFractions[siteComponents[0].id] * externalPressure;
  secondExplicitIsothermAlpha[0] = firstExplicitIsothermAlpha[1] +
                                   b * idealGasMolFractions[siteComponents[0].id] * externalPressure;

  double beta = secondExplicitIsothermAlpha[0];

  explicitIsothermAlphaProduct[0] = 1.0;
  for (size_t i = 1; i < activeComponents; ++i)
  {
    explicitIsothermAlphaProduct[i] =
        (firstExplicitIsothermAlpha[i] / secondExplicitIsothermAlpha[i]) * explicitIsothermAlphaProduct[i - 1];
  }

  for (size_t i = 0; i < activeComponents; ++i)
  {
    size_t index = siteComponents[i].id;
    b = componentScale[i] * siteComponents[i].isotherm.sites[0].parameters[1];
    numberOfMolecules[index] += siteComponents[i].isotherm.sites[0].parameters[0] * b *
                                idealGasMolFractions[index] * externalPressure * explicitIsothermAlphaProduct[i] / beta;
  }

  return std::make_pair(1, 1);
}

std::pair<size_t, size_t> MixturePrediction::computeSegregatedCompetitiveIsotherm(
    std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  std::vector<double> previous(numberOfComponents, 0.0);
  for (size_t site = 0; site < maxIsothermTerms; ++site)
  {
    if (segregatedNumberOfSortedComponents[site] == 0) continue;

    std::copy(numberOfMolecules.begin(), numberOfMolecules.end(), previous.begin());
    acc += computeSegregatedCompetitiveIsotherm(site, idealGasMolFractions, externalPressure, numberOfMolecules,
                                                gasTemperature);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumSiteLoadings[site * numberOfComponents + comp] = numberOfMolecules[comp] - previous[comp];
    }
  }

  const double totalLoading = std::accumulate(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    adsorbedMolFractions[comp] = totalLoading > 0.0 ? numberOfMolecules[comp] / totalLoading : 0.0;
  }

  return acc;
}

std::pair<size_t, size_t> MixturePrediction::computeSegregatedCompetitiveIsotherm(
    size_t site, std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> numberOfMolecules, double& gasTemperature)
{
  const std::vector<Component>& siteComponents = segregatedSortedComponents[site];
  const size_t activeComponents = segregatedNumberOfSortedComponents[site];
  const Isotherm::Type model = siteComponents.front().isotherm.sites.front().type;

  std::vector<double> denominatorTerms(activeComponents, 0.0);
  std::vector<double> numeratorTerms(activeComponents, 0.0);
  std::vector<double> exponents(activeComponents, 1.0);

  for (size_t i = 0; i < activeComponents; ++i)
  {
    const Component& component = siteComponents[i];
    const Isotherm& isotherm = component.isotherm.sites.front();
    if (isotherm.type != model)
    {
      throw std::runtime_error("Error: SCI requires one isotherm model per segregated site");
    }

    const std::vector<double>& parameters = isotherm.parameters;
    const double pressure = finiteNonnegative(idealGasMolFractions[component.id] * externalPressure);
    const double scale = component.scale(gasTemperature);

    switch (model)
    {
      case Isotherm::Type::Langmuir:
      {
        denominatorTerms[i] = safeMultiply(safeMultiply(scale, parameters[1]), pressure);
        numeratorTerms[i] = denominatorTerms[i];
        break;
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        denominatorTerms[i] =
            safeMultiply(safeMultiply(scale, parameters[1]), safePositivePower(pressure, parameters[2]));
        numeratorTerms[i] = denominatorTerms[i];
        break;
      }
      case Isotherm::Type::Sips:
      {
        const double activity = safeMultiply(safeMultiply(scale, parameters[1]), pressure);
        exponents[i] = 1.0 / std::max(finiteNonnegative(parameters[2]), isothermDenominatorFloor);
        denominatorTerms[i] = safePositivePower(activity, exponents[i]);
        numeratorTerms[i] = denominatorTerms[i];
        break;
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        denominatorTerms[i] = safeMultiply(parameters[1], pressure);
        numeratorTerms[i] = safeMultiply(parameters[0], pressure);
        break;
      }
      case Isotherm::Type::Toth:
      {
        numeratorTerms[i] = safeMultiply(parameters[1], pressure);
        exponents[i] = std::max(finiteNonnegative(parameters[2]), isothermDenominatorFloor);
        denominatorTerms[i] = safePositivePower(numeratorTerms[i], exponents[i]);
        break;
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        numeratorTerms[i] = safeMultiply(parameters[0], pressure);
        denominatorTerms[i] = safeMultiply(parameters[1], safePositivePower(pressure, parameters[2]));
        break;
      }
      default:
        throw std::runtime_error("Error: unsupported isotherm model for SCI mixture prediction");
    }
  }

  const double sum = safeActivitySum(denominatorTerms);
  for (size_t i = 0; i < activeComponents; ++i)
  {
    const Component& component = siteComponents[i];
    const Isotherm& isotherm = component.isotherm.sites.front();
    double loading = 0.0;

    switch (model)
    {
      case Isotherm::Type::Langmuir:
      case Isotherm::Type::Langmuir_Freundlich:
      case Isotherm::Type::Sips:
        loading = boundedRatioLoading(isotherm.parameters[0], numeratorTerms[i], 1.0 + sum);
        break;
      case Isotherm::Type::Anti_Langmuir:
        loading = sanitizeUnboundedLoading(numeratorTerms[i] / safeDenominator(1.0 - sum, 1.0e-12));
        break;
      case Isotherm::Type::Toth:
      {
        const double denominator = safeDenominator(safePositivePower(1.0 + sum, 1.0 / exponents[i]), 1.0e-12);
        loading = boundedRatioLoading(isotherm.parameters[0], numeratorTerms[i], denominator);
        break;
      }
      case Isotherm::Type::Redlich_Peterson:
        loading = sanitizeUnboundedLoading(numeratorTerms[i] / safeDenominator(1.0 + sum));
        break;
      default:
        break;
    }

    numberOfMolecules[component.id] += loading;
  }

  return std::make_pair(1, 1);
}

std::pair<size_t, size_t> MixturePrediction::computeSegregatedPureIsotherm(std::span<const double> idealGasMolFractions,
                                                                           const double& externalPressure,
                                                                           std::span<double> adsorbedMolFractions,
                                                                           std::span<double> numberOfMolecules,
                                                                           double& gasTemperature, double pH)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t site = 0; site < maxIsothermTerms; ++site)
  {
    const std::vector<Component>& siteComponents = segregatedSortedComponents[site];
    const size_t activeComponents = segregatedNumberOfSortedComponents[site];
    if (activeComponents == 0) continue;

    for (size_t i = 0; i < activeComponents; ++i)
    {
      const Component& component = siteComponents[i];
      const Isotherm& isotherm = component.isotherm.sites.front();
      const double partialPressure = idealGasMolFractions[component.id] * externalPressure;
      const double loading = safePureSiteLoading(isotherm, partialPressure, component.scale(gasTemperature), pH);
      numberOfMolecules[component.id] += loading;
      equilibriumSiteLoadings[site * numberOfComponents + component.id] = loading;
    }
    acc += std::make_pair<size_t, size_t>(1, 1);
  }

  const double totalLoading = std::accumulate(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    adsorbedMolFractions[comp] = totalLoading > 0.0 ? numberOfMolecules[comp] / totalLoading : 0.0;
  }

  return acc;
}

void MixturePrediction::print() const { std::print("{}", repr()); }

std::string MixturePrediction::repr() const
{
  std::string s;
  s += "Component data\n";
  s += "=======================================================\n";
  s += "maximum isotherm terms:        " + std::to_string(maxIsothermTerms) + "\n";
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    s += sortedComponents[i].repr();
    s += "\n";
  }
  return s;
}

void MixturePrediction::run()
{
  std::vector<double> idealGasMolFractions(numberOfComponents);
  std::vector<double> adsorbedMolFractions(numberOfComponents);
  std::vector<double> numberOfMolecules(numberOfComponents);
  std::vector<double> pureComponentLoadings(numberOfComponents);
  std::vector<double> cachedPressure(numberOfComponents * maxIsothermTerms);
  std::vector<double> cachedGrandPotential(maxIsothermTerms);

  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    idealGasMolFractions[i] = components[i].initialGasMoleFraction;
  }

  std::vector<double> pressures = initPressures();

  // create the output files
  std::vector<std::ofstream> streams;
  for (size_t i = 0; i < numberOfComponents; i++)
  {
    std::string fileName = "component_" + std::to_string(i) + "_" + components[i].name + ".data";
    streams.emplace_back(std::ofstream{fileName});
  }

  for (size_t i = 0; i < numberOfComponents; i++)
  {
    if (predictionMethod == PredictionMethod::MPD)
    {
      std::print(streams[i], "# column 1: target fugacity [Pa]\n");
      std::print(streams[i], "# column 2: pure-component MPD loading (one-hot gas composition)\n");
    }
    else
    {
      std::print(streams[i], "# column 1: total pressure [Pa]\n");
      std::print(streams[i], "# column 2: pure component isotherm value\n");
    }
    std::print(streams[i], "# column 3: mixture component isotherm value\n");
    std::print(streams[i], "# column 4: gas-phase mol-fraction y_i\n");
    std::print(streams[i], "# column 5: adsorbed phase mol-fraction x_i\n");
    if (predictionMethod == PredictionMethod::MPD)
    {
      std::print(streams[i], "# column 6: not applicable for MPD (zero)\n");
      std::print(streams[i], "# column 7: not applicable for MPD (zero)\n");
    }
    else
    {
      std::print(streams[i], "# column 6: hypothetical pressure p_i^*\n");
      std::print(streams[i], "# column 7: reduced grand potential psi_i\n");
    }
  }

  for (size_t i = 0; i < numberOfPressurePoints; ++i)
  {
    predictPureComponentLoadings(pressures[i], pureComponentLoadings, temperature);
    std::pair<double, double> performance =
        predictMixture(idealGasMolFractions, pressures[i], adsorbedMolFractions, numberOfMolecules, cachedPressure,
                       cachedGrandPotential, temperature);
    std::print("Pressure: {} iterations: {}\n", pressures[i], performance.first);

    for (size_t j = 0; j < numberOfComponents; j++)
    {
      const bool hasHypotheticalPressure =
          predictionMethod != PredictionMethod::MPD && adsorbedMolFractions[j] > 0.0;
      const double p_star = hasHypotheticalPressure
                                ? idealGasMolFractions[j] * pressures[i] / adsorbedMolFractions[j]
                                : 0.0;
      const double scale = components[j].scale(temperature);
      const double pureLoading = pureComponentLoadings[j];
      const double grandPotential =
          hasHypotheticalPressure ? components[j].isotherm.psiForPressure(p_star, scale) : 0.0;
      std::print(streams[j], "{:.14g} {:.14g} {:.14g} {:.14g} {:.14g} {:.14g} {:.14g}\n", pressures[i],
                 pureLoading, numberOfMolecules[j], idealGasMolFractions[j], adsorbedMolFractions[j], p_star,
                 grandPotential);
    }
  }
}

std::vector<double> MixturePrediction::initPressures()
{
  std::vector<double> pressures(numberOfPressurePoints);
  if (numberOfPressurePoints > 1)
  {
    switch (pressureScale)
    {
      case PressureScale::Log:
      default:
        for (size_t i = 0; i < numberOfPressurePoints; ++i)
        {
          pressures[i] = std::pow(10, std::log10(pressureStart) +
                                          ((std::log10(pressureEnd) - log10(pressureStart)) *
                                           (static_cast<double>(i) / static_cast<double>(numberOfPressurePoints - 1))));
        }
        break;
      case PressureScale::Linear:
        for (size_t i = 0; i < numberOfPressurePoints; ++i)
        {
          pressures[i] = pressureStart + (pressureEnd - pressureStart) *
                                             (static_cast<double>(i) / static_cast<double>(numberOfPressurePoints - 1));
        }
        break;
    }
  }
  else
  {
    pressures[0] = pressureStart;
  }
  return pressures;
}

void MixturePrediction::printErrorStatus(double psi_value, double sum, double externalPressure,
                                         std::span<const double> idealGasMolFractions,
                                         std::span<double> cachedPressure, double gasTemperature)
{
  std::print("reducedGrandPotential: {}\n", psi_value);
  std::print("sum: {}\n", sum);
  std::print("T: {}\n", gasTemperature);
  for (size_t i = 0; i < numberOfComponents; ++i) std::print("cachedPressure: {}\n", cachedPressure[i]);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    if (components[i].isCarrierGas) continue;
    double value = components[i].isotherm.inversePressureForPsi(psi_value, cachedPressure[i], gasTemperature);
    std::print("inversePressure: {}\n", value);
  }
  std::print("externalPressure: {}\n", externalPressure);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    std::print("idealGasMolFractions[i] {} {}\n", i, idealGasMolFractions[i]);
  }
}

void MixturePrediction::sortComponents()
{
  const auto isActiveAdsorbate = [this](const Component& component)
  {
    return !component.isCarrierGas &&
           (predictionMethod == PredictionMethod::MPD || component.isotherm.enabled());
  };

  if (predictionMethod == PredictionMethod::EI)
  {
    std::sort(sortedComponents.begin(), sortedComponents.end(), &LangmuirLoadingSorter);
  }
  else if (predictionMethod == PredictionMethod::SIAST || predictionMethod == PredictionMethod::SEI ||
           predictionMethod == PredictionMethod::SCI || predictionMethod == PredictionMethod::SPI)
  {
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      segregatedSortedComponents[i] = components;
      size_t activeComponents = 0;
      for (size_t j = 0; j < numberOfComponents; ++j)
      {
        if (!components[j].isCarrierGas && i < components[j].isotherm.sites.size())
        {
          segregatedSortedComponents[i][j].isotherm = MultiSiteIsotherm({components[j].isotherm.sites[i]});
          const bool siteEnabled = components[j].isotherm.sites[i].enabled();
          segregatedSortedComponents[i][j].isCarrierGas = !siteEnabled;
          if (siteEnabled) ++activeComponents;
        }
        else
        {
          segregatedSortedComponents[i][j].isotherm = MultiSiteIsotherm{};
          segregatedSortedComponents[i][j].isCarrierGas = true;
        }
      }
      segregatedNumberOfSortedComponents[i] = activeComponents;
    }
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      if (predictionMethod == PredictionMethod::SEI)
      {
        std::sort(segregatedSortedComponents[i].begin(), segregatedSortedComponents[i].end(), &LangmuirLoadingSorter);
      }
      else
      {
        std::stable_partition(segregatedSortedComponents[i].begin(), segregatedSortedComponents[i].end(),
                              [](const Component& component) { return !component.isCarrierGas; });
      }
    }
    std::stable_partition(sortedComponents.begin(), sortedComponents.end(), isActiveAdsorbate);
  }
  else
  {
    std::stable_partition(sortedComponents.begin(), sortedComponents.end(), isActiveAdsorbate);
  }

  numberOfSortedComponents =
      static_cast<size_t>(std::count_if(sortedComponents.begin(), sortedComponents.end(), isActiveAdsorbate));
  hypotheticalPressure.resize(numberOfSortedComponents);
  reducedGrandPotential.resize(numberOfSortedComponents);
  residualVector.resize(numberOfSortedComponents);
  correctionVector.resize(numberOfSortedComponents);
  jacobianMatrix.resize(numberOfSortedComponents * numberOfSortedComponents);
}
