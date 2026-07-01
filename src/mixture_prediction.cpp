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

bool LangmuirLoadingSorter(Component const& lhs, Component const& rhs)
{
  if (lhs.isCarrierGas) return false;
  if (rhs.isCarrierGas) return true;
  return lhs.isotherm.sites[0].parameters[0] < rhs.isotherm.sites[0].parameters[0];
}

MixturePrediction::MixturePrediction(const InputReader& inputreader)
    : displayName(inputreader.displayName),
      components(inputreader.components),
      sortedComponents(components),
      numberOfComponents(components.size()),
      numberOfSortedComponents(components.size() - inputreader.numberOfCarrierGases),
      numberOfCarrierGases(inputreader.numberOfCarrierGases),
      carrierGasComponent(inputreader.carrierGasComponent),
      predictionMethod(PredictionMethod(inputreader.mixturePredictionMethod)),
      iastMethod(IASTMethod(inputreader.IASTMethod)),
      maxIsothermTerms(inputreader.maxIsothermTerms),
      segregatedSortedComponents(maxIsothermTerms, std::vector<Component>(components)),
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
                         { return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites; });
    maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
  }
  segregatedSortedComponents =
      std::vector<std::vector<Component>>(maxIsothermTerms, std::vector<Component>(components));

  sortComponents();
}

std::pair<size_t, size_t> MixturePrediction::predictMixture(std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules, double* cachedPressure,
                                                            double* cachedGrandPotential, double& gasTemperature)
{
  const double tiny = 1.0e-10;

  if (externalPressure < 0.0)
  {
    printErrorStatus(0.0, 0.0, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
    throw std::runtime_error("Error (IAST): negative total pressure\n");
  }

  double sumYi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sumYi += idealGasMolFractions[i];
  }
  if (std::abs(sumYi - 1.0) > 1e-15)
  {
    printErrorStatus(0.0, sumYi, externalPressure, idealGasMolFractions, cachedPressure, gasTemperature);
    throw std::runtime_error("Error (IAST): sum idealGasMolFractions at IAST start not unity\n");
  }

  // if only an inert component present
  // this happens at the beginning of the simulation when the whole column is filled with the carrier gas
  if (std::abs(idealGasMolFractions[carrierGasComponent] - 1.0) < tiny)
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      adsorbedMolFractions[i] = 0.0;
      numberOfMolecules[i] = 0.0;
    }

    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  switch (predictionMethod)
  {
    case PredictionMethod::IAST:
    default:
      switch (iastMethod)
      {
        case IASTMethod::FastIAST:
        default:
          return computeFastIAST(idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                                 cachedPressure, cachedGrandPotential, gasTemperature);
        case IASTMethod::NestedLoopBisection:
          return computeIASTNestedLoopBisection(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                                numberOfMolecules, cachedPressure, cachedGrandPotential,
                                                gasTemperature);
      }
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
      return computeExplicitIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                                     gasTemperature);
    case PredictionMethod::SEI:
      return computeSegratedExplicitIsotherm(idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                             numberOfMolecules, gasTemperature);
  }
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeFastIAST(std::span<const double> idealGasMolFractions,
                                                             const double& externalPressure,
                                                             std::span<double> adsorbedMolFractions,
                                                             std::span<double> numberOfMolecules,
                                                             double* cachedPressure, double* cachedGrandPotential,
                                                             double& gasTemperature)
{
  const double tiny = 1.0e-13;

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
  if (numberOfCarrierGases > 0)
  {
    adsorbedMolFractions[carrierGasComponent] = 0.0;
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
  if (numberOfCarrierGases > 0)
  {
    numberOfMolecules[carrierGasComponent] = 0.0;
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
                                                              double* cachedPressure, double* cachedGrandPotential,
                                                              double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeFastSIAST(i, idealGasMolFractions, externalPressure, adsorbedMolFractions, numberOfMolecules,
                            cachedPressure, cachedGrandPotential, gasTemperature);
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = numberOfMolecules[i] / N;
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
                                                              double* cachedPressure, double* cachedGrandPotential,
                                                              double& gasTemperature)
{
  const double tiny = 1.0e-13;

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

  if (cachedGrandPotential[site] > tiny)
  {
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      hypotheticalPressure[i] = cachedPressure[sortedComponents[i].id + site * numberOfComponents];
    }
  }
  else
  {
    double initial_psi = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      double temp_psi = idealGasMolFractions[sortedComponents[i].id] *
                        sortedComponents[i].isotherm.psiForPressure(site, externalPressure, componentScale[i]);
      initial_psi += temp_psi;
    }
    cachedGrandPotential[site] = initial_psi;

    double cachevalue = 0.0;
    for (size_t i = 0; i < numberOfSortedComponents; ++i)
    {
      hypotheticalPressure[i] =
          1.0 / sortedComponents[i].isotherm.inversePressureForPsi(site, initial_psi, cachevalue, componentScale[i]);
    }
  }

  double error = 1.0;
  double sum_xi = 1.0;
  do
  {
    // compute residualVector
    for (size_t i = 0; i < numberOfSortedComponents - 1; ++i)
    {
      residualVector[i] =
          sortedComponents[i].isotherm.psiForPressure(site, hypotheticalPressure[i], componentScale[i]) -
          sortedComponents[numberOfSortedComponents - 1].isotherm.psiForPressure(
              site, hypotheticalPressure[numberOfSortedComponents - 1], componentScale[numberOfSortedComponents - 1]);
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
          sortedComponents[i].isotherm.value(site, hypotheticalPressure[i], componentScale[i]) /
          hypotheticalPressure[i];
    }
    for (size_t i = 0; i < numberOfSortedComponents - 1; i++)
    {
      jacobianMatrix[i + (numberOfSortedComponents - 1) * numberOfSortedComponents] =
          -sortedComponents[numberOfSortedComponents - 1].isotherm.value(
              site, hypotheticalPressure[numberOfSortedComponents - 1], componentScale[numberOfSortedComponents - 1]) /
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
          sortedComponents[i].isotherm.psiForPressure(site, hypotheticalPressure[i], componentScale[i]);
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
    cachedPressure[sortedComponents[i].id + site * numberOfComponents] = hypotheticalPressure[i];
  }

  for (size_t i = 0; i < numberOfSortedComponents; ++i)
  {
    adsorbedMolFractions[sortedComponents[i].id] =
        idealGasMolFractions[sortedComponents[i].id] * externalPressure / std::max(hypotheticalPressure[i], 1e-15);
  }
  if (numberOfCarrierGases > 0)
  {
    adsorbedMolFractions[carrierGasComponent] = 0.0;
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
                       sortedComponents[i].isotherm.value(site, hypotheticalPressure[i], componentScale[i]);
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    numberOfMolecules[i] += adsorbedMolFractions[i] / inverse_q_total;
  }
  if (numberOfCarrierGases > 0)
  {
    numberOfMolecules[carrierGasComponent] = 0.0;
  }

  return std::make_pair(numberOfIASTSteps, 1);
}

// idealGasMolFractions  = gas phase molefraction
// externalPressure   = total pressure
// adsorbedMolFractions  = adsorbed phase molefraction
// numberOfMolecules  = number of adsorbed molecules of component i
std::pair<size_t, size_t> MixturePrediction::computeIASTNestedLoopBisection(
    std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double* cachedPressure,
    double* cachedGrandPotential, double& gasTemperature)
{
  const double tiny = 1.0e-15;

  std::vector<double> componentScale(numberOfComponents);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    componentScale[i] = components[i].scale(gasTemperature);
  }

  double initial_psi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    initial_psi += idealGasMolFractions[i] * components[i].isotherm.psiForPressure(externalPressure, componentScale[i]);
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

  // condition 1: same reduced grand potential for all components (done by using a single variable)
  // condition 2: mol-fractions add up to unity

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedGrandPotential[0] > tiny)
  {
    initial_psi = cachedGrandPotential[0];
  }
  // for this initial estimate 'initial_psi' compute the sum of mol-fractions
  double sumXi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sumXi += idealGasMolFractions[i] * externalPressure *
             components[i].isotherm.inversePressureForPsi(initial_psi, cachedPressure[i], componentScale[i]);
  }

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = 0.0;
      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        sumXi += idealGasMolFractions[i] * externalPressure *
                 components[i].isotherm.inversePressureForPsi(right_bracket, cachedPressure[i], componentScale[i]);
      }
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

      sumXi = 0.0;
      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        sumXi += idealGasMolFractions[i] * externalPressure *
                 components[i].isotherm.inversePressureForPsi(left_bracket, cachedPressure[i], componentScale[i]);
      }
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

    sumXi = 0.0;
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      sumXi += idealGasMolFractions[i] * externalPressure *
               components[i].isotherm.inversePressureForPsi(psi_value, cachedPressure[i], componentScale[i]);
    }

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

  sumXi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sumXi += idealGasMolFractions[i] * externalPressure *
             components[i].isotherm.inversePressureForPsi(psi_value, cachedPressure[i], componentScale[i]);
  }

  // cache the value of reducedGrandPotential for subsequent use
  cachedGrandPotential[0] = psi_value;

  // calculate mol-fractions in adsorbed phase and total loading
  double inverse_q_total = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    double ip = components[i].isotherm.inversePressureForPsi(psi_value, cachedPressure[i], componentScale[i]);
    adsorbedMolFractions[i] = idealGasMolFractions[i] * externalPressure * ip / sumXi;

    if (adsorbedMolFractions[i] > tiny)
    {
      inverse_q_total += adsorbedMolFractions[i] / components[i].isotherm.value(1.0 / ip, componentScale[i]);
    }
    else
    {
      adsorbedMolFractions[i] = 0.0;
    }
  }

  // calculate loading for all of the components
  if (inverse_q_total == 0.0)
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      numberOfMolecules[i] = 0.0;
    }
  }
  else
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      numberOfMolecules[i] = adsorbedMolFractions[i] / inverse_q_total;
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
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double* cachedPressure,
    double* cachedGrandPotential, double& gasTemperature)
{
  std::fill(adsorbedMolFractions.begin(), adsorbedMolFractions.end(), 0.0);
  std::fill(numberOfMolecules.begin(), numberOfMolecules.end(), 0.0);

  std::pair<size_t, size_t> acc;
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeSIASTNestedLoopBisection(i, idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                           numberOfMolecules, cachedPressure, cachedGrandPotential, gasTemperature);
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = numberOfMolecules[i] / N;
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
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double* cachedPressure,
    double* cachedGrandPotential, double& gasTemperature)
{
  const double tiny = 1.0e-15;

  std::vector<double> componentScale(numberOfComponents);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    componentScale[i] = components[i].scale(gasTemperature);
  }

  double initial_psi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    initial_psi +=
        idealGasMolFractions[i] * components[i].isotherm.psiForPressure(site, externalPressure, componentScale[i]);
  }

  if (initial_psi < tiny)
  {
    // nothing is adsorbing
    // do not count it for the IAST statistics
    return std::make_pair(0, 0);
  }

  // condition 1: same reduced grand potential for all components (done by using a single variable)
  // condition 2: mol-fractions add up to unity

  double psi_value = 0.0;
  size_t nr_steps = 0;
  if (cachedGrandPotential[site] > tiny)
  {
    initial_psi = cachedGrandPotential[site];
  }
  // for this initial estimate 'initial_psi' compute the sum of mol-fractions
  double sumXi = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    sumXi += idealGasMolFractions[i] * externalPressure *
             components[i].isotherm.inversePressureForPsi(
                 site, initial_psi, cachedPressure[i + numberOfComponents * site], componentScale[i]);
  }

  // initialize the bisection algorithm
  double left_bracket = initial_psi;
  double right_bracket = initial_psi;
  if (sumXi > 1.0)
  {
    do
    {
      right_bracket *= 2.0;

      sumXi = 0.0;
      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        sumXi += idealGasMolFractions[i] * externalPressure *
                 components[i].isotherm.inversePressureForPsi(
                     site, right_bracket, cachedPressure[i + numberOfComponents * site], componentScale[i]);
      }
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

      sumXi = 0.0;
      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        sumXi += idealGasMolFractions[i] * externalPressure *
                 components[i].isotherm.inversePressureForPsi(
                     site, left_bracket, cachedPressure[i + numberOfComponents * site], componentScale[i]);
      }
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

    sumXi = 0.0;
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      sumXi += idealGasMolFractions[i] * externalPressure *
               components[i].isotherm.inversePressureForPsi(
                   site, psi_value, cachedPressure[i + numberOfComponents * site], componentScale[i]);
    }

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

  // calculate mol-fractions in adsorbed phase and total loading
  double inverse_q_total = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    double ip = components[i].isotherm.inversePressureForPsi(
        site, psi_value, cachedPressure[i + numberOfComponents * site], componentScale[i]);
    adsorbedMolFractions[i] = idealGasMolFractions[i] * externalPressure * ip;

    if (adsorbedMolFractions[i] > tiny)
    {
      inverse_q_total += adsorbedMolFractions[i] / components[i].isotherm.value(site, 1.0 / ip, componentScale[i]);
    }
  }

  // calculate loading for all of the components
  if (inverse_q_total > 0.0)
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      numberOfMolecules[i] += adsorbedMolFractions[i] / inverse_q_total;
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
  std::vector<double> componentScale(numberOfComponents);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    componentScale[i] = sortedComponents[i].scale(gasTemperature);
  }

  adsorbedMoleFractionsScratch[0] = 1.0;
  for (size_t i = 1; i < numberOfComponents; ++i)
  {
    adsorbedMoleFractionsScratch[i] =
        sortedComponents[i].isotherm.sites[0].parameters[0] / sortedComponents[i - 1].isotherm.sites[0].parameters[0];
  }

  double b =
      componentScale[numberOfComponents - 1] * sortedComponents[numberOfComponents - 1].isotherm.sites[0].parameters[1];
  firstExplicitIsothermAlpha[numberOfComponents - 1] =
      std::pow((1.0 + b * idealGasMolFractions[sortedComponents[numberOfComponents - 1].id] * externalPressure),
               adsorbedMoleFractionsScratch[numberOfComponents - 1]);
  secondExplicitIsothermAlpha[numberOfComponents - 1] =
      1.0 + b * idealGasMolFractions[sortedComponents[numberOfComponents - 1].id] * externalPressure;
  for (size_t i = numberOfComponents - 2; i > 0; i--)
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
  for (size_t i = 1; i < numberOfComponents; ++i)
  {
    explicitIsothermAlphaProduct[i] =
        (firstExplicitIsothermAlpha[i] / secondExplicitIsothermAlpha[i]) * explicitIsothermAlphaProduct[i - 1];
  }

  for (size_t i = 0; i < numberOfComponents; ++i)
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
    adsorbedMolFractions[i] = numberOfMolecules[i] / N;
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
  for (size_t i = 0; i < maxIsothermTerms; ++i)
  {
    acc += computeSegratedExplicitIsotherm(i, idealGasMolFractions, externalPressure, adsorbedMolFractions,
                                           numberOfMolecules, gasTemperature);
  }

  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = numberOfMolecules[i] / N;
  }

  return acc;
}

std::pair<size_t, size_t> MixturePrediction::computeSegratedExplicitIsotherm(
    size_t site, std::span<const double> idealGasMolFractions, const double& externalPressure,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, double& gasTemperature)
{
  std::vector<std::vector<double>> componentScale(maxIsothermTerms, std::vector<double>(numberOfComponents));
  for (size_t j = 0; j < maxIsothermTerms; ++j)
  {
    for (size_t i = 0; i < numberOfComponents; ++i)
    {
      componentScale[j][i] = sortedComponents[i].scale(gasTemperature);
    }
  }

  adsorbedMoleFractionsScratch[0] = 1.0;
  for (size_t i = 1; i < numberOfComponents; ++i)
  {
    adsorbedMoleFractionsScratch[i] = segregatedSortedComponents[site][i].isotherm.sites[0].parameters[0] /
                                      segregatedSortedComponents[site][i - 1].isotherm.sites[0].parameters[0];
  }

  double b = componentScale[site][numberOfComponents - 1] *
             segregatedSortedComponents[site][numberOfComponents - 1].isotherm.sites[site].parameters[1];
  firstExplicitIsothermAlpha[numberOfComponents - 1] = std::pow(
      (1.0 + b * idealGasMolFractions[segregatedSortedComponents[site][numberOfComponents - 1].id] * externalPressure),
      adsorbedMoleFractionsScratch[numberOfComponents - 1]);
  secondExplicitIsothermAlpha[numberOfComponents - 1] =
      1.0 + b * idealGasMolFractions[segregatedSortedComponents[site][numberOfComponents - 1].id] * externalPressure;
  for (size_t i = numberOfComponents - 2; i > 0; i--)
  {
    b = componentScale[site][i] * segregatedSortedComponents[site][i].isotherm.sites[site].parameters[1];
    firstExplicitIsothermAlpha[i] =
        std::pow((firstExplicitIsothermAlpha[i + 1] +
                  b * idealGasMolFractions[segregatedSortedComponents[site][i].id] * externalPressure),
                 adsorbedMoleFractionsScratch[i]);
    secondExplicitIsothermAlpha[i] =
        firstExplicitIsothermAlpha[i + 1] +
        b * idealGasMolFractions[segregatedSortedComponents[site][i].id] * externalPressure;
  }

  b = componentScale[site][0] * segregatedSortedComponents[site][0].isotherm.sites[site].parameters[1];
  firstExplicitIsothermAlpha[0] = firstExplicitIsothermAlpha[1] +
                                  b * idealGasMolFractions[segregatedSortedComponents[site][0].id] * externalPressure;
  secondExplicitIsothermAlpha[0] = firstExplicitIsothermAlpha[1] +
                                   b * idealGasMolFractions[segregatedSortedComponents[site][0].id] * externalPressure;

  double beta = secondExplicitIsothermAlpha[0];

  explicitIsothermAlphaProduct[0] = 1.0;
  for (size_t i = 1; i < numberOfComponents; ++i)
  {
    explicitIsothermAlphaProduct[i] =
        (firstExplicitIsothermAlpha[i] / secondExplicitIsothermAlpha[i]) * explicitIsothermAlphaProduct[i - 1];
  }

  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    size_t index = segregatedSortedComponents[site][i].id;
    b = componentScale[site][i] * segregatedSortedComponents[site][i].isotherm.sites[site].parameters[1];
    numberOfMolecules[index] += segregatedSortedComponents[site][i].isotherm.sites[0].parameters[0] * b *
                                idealGasMolFractions[index] * externalPressure * explicitIsothermAlphaProduct[i] / beta;
  }
  double N = 0.0;
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    N += numberOfMolecules[i];
  }
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
    adsorbedMolFractions[i] = numberOfMolecules[i] / N;
  }

  return std::make_pair(1, 1);
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
    std::print(streams[i], "# column 1: total pressure [Pa]\n");
    std::print(streams[i], "# column 2: pure component isotherm value\n");
    std::print(streams[i], "# column 3: mixture component isotherm value\n");
    std::print(streams[i], "# column 4: gas-phase mol-fraction y_i\n");
    std::print(streams[i], "# column 5: adsorbed phase mol-fraction x_i\n");
    std::print(streams[i], "# column 6: hypothetical pressure p_i^*\n");
    std::print(streams[i], "# column 7: reduced grand potential psi_i\n");
  }

  for (size_t i = 0; i < numberOfPressurePoints; ++i)
  {
    std::pair<double, double> performance =
        predictMixture(idealGasMolFractions, pressures[i], adsorbedMolFractions, numberOfMolecules, &cachedPressure[0],
                       &cachedGrandPotential[0], temperature);
    std::print("Pressure: {} iterations: {}\n", pressures[i], performance.first);

    for (size_t j = 0; j < numberOfComponents; j++)
    {
      double p_star = idealGasMolFractions[j] * pressures[i] / adsorbedMolFractions[j];
      double scale = components[j].scale(temperature);
      std::print(streams[j], "{:.14g} {:.14g} {:.14g} {:.14g} {:.14g} {:.14g}\n", pressures[i],
                 components[j].isotherm.value(pressures[i], scale), numberOfMolecules[j], idealGasMolFractions[j],
                 adsorbedMolFractions[j], components[j].isotherm.psiForPressure(p_star, scale));
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
      case PressureScale::Normal:
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
                                         std::span<const double> idealGasMolFractions, double cachedPressure[],
                                         double gasTemperature)
{
  std::print("reducedGrandPotential: {}\n", psi_value);
  std::print("sum: {}\n", sum);
  std::print("T: {}\n", gasTemperature);
  for (size_t i = 0; i < numberOfComponents; ++i) std::print("cachedPressure: {}\n", cachedPressure[i]);
  for (size_t i = 0; i < numberOfComponents; ++i)
  {
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
  if (predictionMethod == PredictionMethod::EI)
  {
    std::sort(sortedComponents.begin(), sortedComponents.end(), &LangmuirLoadingSorter);
  }
  else if (predictionMethod == PredictionMethod::SEI)
  {
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      for (size_t j = 0; j < numberOfComponents; ++j)
      {
        if (j != carrierGasComponent)
        {
          segregatedSortedComponents[i][j].isotherm.sites[0] = components[j].isotherm.sites[i];
          segregatedSortedComponents[i][j].isotherm.numberOfSites = 1;
        }
      }
    }
    for (size_t i = 0; i < maxIsothermTerms; ++i)
    {
      std::sort(segregatedSortedComponents[i].begin(), segregatedSortedComponents[i].end(), &LangmuirLoadingSorter);
    }
  }
  else
  {
    auto it = sortedComponents.begin() + static_cast<std::ptrdiff_t>(carrierGasComponent);
    std::rotate(it, it + 1, sortedComponents.end());
  }
}
