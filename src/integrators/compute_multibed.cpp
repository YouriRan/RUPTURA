#include "compute_multibed.h"

#include <algorithm>
#include <mdspan>

#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;

void updateVelocityAndPressure(
    const std::vector<Component>& components, const ColumnMultibed::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const double> adsorbentScaledVoidFraction, std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature)
{
  auto refreshNode = [&](size_t grid)
  {
    totalConcentration[grid] = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid]));
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }

    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = moleFraction[index] * totalPressure[grid];
      gasDensity[grid] += moleFraction[index] * totalConcentration[grid] * components[comp].molecularWeight;
    }
  };

  auto ergunGrad = [&](size_t grid)
  {
    double grad = 0.0;
    double visc = 150.0 * dynamicViscosity * interstitialGasVelocity[grid];
    double iner = 1.75 * gasDensity[grid] * interstitialGasVelocity[grid] * interstitialGasVelocity[grid];
    for (size_t ads = 0; ads < numberOfAdsorbents; ads++)
    {
      double prefactor = adsorbentScaledVoidFraction[grid * numberOfAdsorbents + ads];
      grad += visc * prefactor * prefactor * fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      grad += iner * prefactor * fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
    }
    return grad;
  };

  auto sinkTerm = [&](size_t grid)
  {
    const auto begin = static_cast<std::ptrdiff_t>(grid * numberOfComponents);
    const auto end = static_cast<std::ptrdiff_t>((grid + 1) * numberOfComponents);
    return gridSpacing(columnDistances, grid) *
           std::reduce(bulkSpeciesSink.begin() + begin, bulkSpeciesSink.begin() + end);
  };

  auto gridRatio = [&](size_t grid) -> double
  { return columnLength <= 0.0 ? 0.0 : columnDistances[grid] / columnLength; };

  auto bisection = [&](auto&& func)
  {
    double a = 1e-7;
    double b = std::max(10.0 * std::abs(columnEntranceVelocity), 1e-6);
    double tolerance = 1e-6;
    double c = a;
    double fa = func(a);
    double fb = func(b);

    for (size_t expansion = 0; fa * fb > 0.0 && expansion < 20; expansion++)
    {
      b *= 2.0;
      fb = func(b);
    }

    if (fa * fb > 0.0)
    {
      throw std::runtime_error("Bounds for bisection method to solve for velocity improperly set.\n");
    }

    while ((b - a) > tolerance)
    {
      c = 0.5 * (a + b);
      double fc = func(c);

      if (fc == 0.0) break;

      if (fa * fc < 0.0)
      {
        b = c;
        fb = fc;
      }
      else
      {
        a = c;
        fa = fc;
      }
    }

    return c;
  };

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }

    totalConcentration[grid] = concentrationSum;
    totalPressure[grid] = totalConcentration[grid] * R * std::max(1e-10, gasTemperature[grid]);
    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = concentration[index] * R * std::max(1e-10, gasTemperature[grid]);
      gasDensity[grid] += concentration[index] * components[comp].molecularWeight;
    }
  }

  if (boundaryCondition == ColumnMultibed::BoundaryCondition::InletPressureInletVelocity)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[0] = inletPressure;
    refreshNode(0);

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
      refreshNode(grid);
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
    }
  }
  else if (boundaryCondition == ColumnMultibed::BoundaryCondition::InletPressureOutletPressure)
  {
    auto shoot = [&](double v0)
    {
      interstitialGasVelocity[0] = v0;
      totalPressure[0] = inletPressure;
      refreshNode(0);

      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        const double cprev = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid - 1]));
        const double cnow = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid]));
        interstitialGasVelocity[grid] =
            (interstitialGasVelocity[grid - 1] * cprev - sinkTerm(grid)) / std::max(1e-10, cnow);
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
        if (totalPressure[grid] <= 0.0)
        {
          return -outletPressure;
        }
        refreshNode(grid);
      }

      return totalPressure[numberOfGridPoints] - outletPressure;
    };

    columnEntranceVelocity = bisection(shoot);
    shoot(columnEntranceVelocity);
  }
  else if (boundaryCondition == ColumnMultibed::BoundaryCondition::InletVelocityOutletPressure)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[numberOfGridPoints] = outletPressure;
    refreshNode(numberOfGridPoints);

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
    }
    for (size_t grid = numberOfGridPoints; grid > 0; --grid)
    {
      const size_t current = grid - 1;
      totalPressure[current] =
          totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
      refreshNode(current);
    }
  }
  else if (boundaryCondition == ColumnMultibed::BoundaryCondition::FixedVelocity)
  {
    std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), columnEntranceVelocity);

    if (inletPressure > 0.0)
    {
      totalPressure[0] = inletPressure;
      refreshNode(0);
      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
        refreshNode(grid);
      }
    }
    else
    {
      totalPressure[numberOfGridPoints] = outletPressure;
      refreshNode(numberOfGridPoints);
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        const size_t current = grid - 1;
        totalPressure[current] =
            totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
        refreshNode(current);
      }
    }
  }
  else if (boundaryCondition == ColumnMultibed::BoundaryCondition::FixedPressureInletVelocity)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      totalPressure[grid] = inletPressure + pressureGradient * gridRatio(grid) / columnLength;
      refreshNode(grid);
    }

    for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
    }
  }

  if (totalPressure[numberOfGridPoints] <= 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computePhysisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(0.0, moleFraction[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    if (sum > 0.0)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] /= sum;
      }
    }
    else
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] = 1.0 / static_cast<double>(numberOfComponents);
      }
    }

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumPhysisorption[grid * numberOfComponents + comp] = 0.0;
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (!hasAdsorbentOfType[grid * numberOfAdsorbents + ads]) continue;

      // use Yi and Pt[i] to compute the loadings in the physisorption mixture via mixture prediction
      std::span<double> spanCachedPressure =
          cachedPressure.subspan((grid * numberOfAdsorbents + ads) * numberOfComponents * maxIsothermTerms,
                                 numberOfComponents * maxIsothermTerms);
      std::span<double> spanCachedGrandPotential =
          cachedGrandPotential.subspan((grid * numberOfAdsorbents + ads) * maxIsothermTerms, maxIsothermTerms);
      iastPerformance += physisorptionMixtures[ads].predictMixture(
          idealGasMolFractions, totalPressure[grid], adsorbedMolFractions, numberOfMolecules, spanCachedPressure,
          spanCachedGrandPotential, gasTemperature[grid]);

      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumPhysisorption[grid * numberOfComponents + comp] +=
            fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * numberOfMolecules[comp];
      }
    }
  }
}

void computePhysisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents,
                          std::span<const double> fractionOfAdsorbent, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double massTransferCoefficient = 0.0;
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        massTransferCoefficient += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                                   physisorptionMixtures[ads].components[comp].massTransferCoefficient;
      }

      const size_t index = grid * numberOfComponents + comp;
      physisorptionDot[index] = massTransferCoefficient * (equilibriumPhysisorption[index] - physisorption[index]);
    }
  }
}

void computeBulkSpeciesSink(size_t numberOfGridPoints, size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                            std::span<const double> fractionOfAdsorbent, std::span<const double> totalVoidFraction,
                            std::span<const double> physisorptionDot, std::span<double> bulkSpeciesSink)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double solidLoadingDensity = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      solidLoadingDensity += fraction * (1.0 - adsorbentVoidFractions[ads]) * particleDensities[ads];
    }

    const double loadingPrefactor = solidLoadingDensity / std::max(1e-10, totalVoidFraction[grid]);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      bulkSpeciesSink[index] = loadingPrefactor * physisorptionDot[index];
    }
  }
}

void computeMassDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  auto axialDispersion = [&](size_t grid, size_t comp)
  {
    double coefficient = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      coefficient += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                     physisorptionMixtures[ads].components[comp].axialDispersionCoefficient;
    }
    return coefficient;
  };

  auto secondDerivative = [&](size_t grid, size_t comp)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 *
           (hm * spanConcentration[grid + 1, comp] - (hm + hp) * spanConcentration[grid, comp] +
            hp * spanConcentration[grid - 1, comp]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = 0.0;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    const double hm = gridSpacing(columnDistances, grid);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double dvcDz = (interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
                            interstitialGasVelocity[grid - 1] * spanConcentration[grid - 1, comp]) /
                           std::max(1e-30, hm);
      spanConcentrationDot[grid, comp] =
          -dvcDz + axialDispersion(grid, comp) * secondDerivative(grid, comp) - spanBulkSpeciesSink[grid, comp];
    }
  }

  const size_t outlet = numberOfGridPoints;
  const double hm = gridSpacing(columnDistances, outlet);
  const double invHm = 1.0 / std::max(1e-30, hm);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double dvcDz = (interstitialGasVelocity[outlet] * spanConcentration[outlet, comp] -
                          interstitialGasVelocity[outlet - 1] * spanConcentration[outlet - 1, comp]) *
                         invHm;
    const double d2cDz2 = (spanConcentration[outlet - 1, comp] - spanConcentration[outlet, comp]) * invHm * invHm;
    spanConcentrationDot[outlet, comp] =
        -dvcDz + axialDispersion(outlet, comp) * d2cDz2 - spanBulkSpeciesSink[outlet, comp];
  }
}

void computeEnergyDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                              size_t numberOfComponents, size_t numberOfAdsorbents, double externalTemperature,
                              std::span<const double> totalVoidFraction, std::span<const double> particleDensities,
                              std::span<const double> particleDiameters, std::span<const double> fractionOfAdsorbent,
                              double internalDiameter, double outerDiameter, double wallDensity,
                              double gasThermalConductivity, double wallThermalConductivity,
                              double heatTransferGasSolid, double heatTransferGasWall, double heatTransferWallExternal,
                              double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
                              std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
                              std::span<const double> gasDensity, std::span<double> coeffDiffusion,
                              std::span<const double> physisorptionDot, std::span<const double> gasTemperature,
                              std::span<double> gasTemperatureDot, std::span<const double> solidTemperature,
                              std::span<double> solidTemperatureDot, std::span<const double> wallTemperature,
                              std::span<double> wallTemperatureDot)
{
  mdspan2d_const spanAdsorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  auto idx = [&](size_t grid) { return 1.0 / std::max(1e-30, gridSpacing(columnDistances, grid)); };
  auto d2Temperature = [&](std::span<const double> values, size_t grid)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 * (hm * values[grid + 1] - (hm + hp) * values[grid] + hp * values[grid - 1]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(numberOfGridPoints + 1);
  // computeWENO(gasTemperature, gasTemperatureFlux);

  // is this extra 1/eps necessary? It's in python not in eqs
  double prefactorGasWall = 4.0 * heatTransferGasWall / (heatCapacityGas * internalDiameter);

  auto accessibleSurface = [&](size_t grid)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * 6.0 / particleDiameters[ads];
    }
    return value;
  };

  auto solidDensity = [&](size_t grid)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * particleDensities[ads];
    }
    return value;
  };

  auto coeffSolidGasAt = [&](size_t grid)
  {
    return accessibleSurface(grid) * heatTransferGasSolid / (heatCapacitySolid * std::max(1e-10, solidDensity(grid)));
  };

  // prefactor 4 in python, 2 in eqs
  double internalArea = 4.0 * internalDiameter / (outerDiameter * outerDiameter - internalDiameter * internalDiameter);
  double externalArea = 4.0 * outerDiameter / (outerDiameter * outerDiameter - internalDiameter * internalDiameter);
  double invHeatDensityWall = 1.0 / (heatCapacityWall * wallDensity);
  double coeffWallGas = heatTransferGasWall * internalArea * invHeatDensityWall;
  double coeffWallExternal = heatTransferWallExternal * externalArea * invHeatDensityWall;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  auto gasHeatExchange = [&](size_t grid)
  {
    const double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    const double relativeVolume = (1.0 - totalVoidFraction[grid]) / std::max(1e-10, totalVoidFraction[grid]);
    const double gasSolidExchange =
        relativeVolume * accessibleSurface(grid) * heatTransferGasSolid * invGasDensity / heatCapacityGas;
    const double gasWallExchange = prefactorGasWall * invGasDensity;
    return gasSolidExchange * (solidTemperature[grid] - gasTemperature[grid]) +
           gasWallExchange * (wallTemperature[grid] - gasTemperature[grid]);
  };

  auto solidHeatExchange = [&](size_t grid)
  { return coeffSolidGasAt(grid) * (gasTemperature[grid] - solidTemperature[grid]); };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  auto heatOfAdsorption = [&](size_t grid, size_t comp)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
               physisorptionMixtures[ads].components[comp].heatOfAdsorption;
    }
    return value;
  };

  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = solidHeatExchange(0);
  wallTemperatureDot[0] = wallHeatExchange(0);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[0] += heatOfAdsorption(0, comp) * spanAdsorptionDot[0, comp] / heatCapacitySolid;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    gasTemperatureDot[grid] = gasHeatExchange(grid);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    gasTemperatureDot[grid] -=
        interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx(grid);
    gasTemperatureDot[grid] += coeffDiffusion[grid] * d2Temperature(gasTemperature, grid);

    wallTemperatureDot[grid] += (wallThermalConductivity * invHeatDensityWall) * d2Temperature(wallTemperature, grid);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] += heatOfAdsorption(grid, comp) * spanAdsorptionDot[grid, comp] / heatCapacitySolid;
    }
  }

  gasTemperatureDot[numberOfGridPoints] = gasHeatExchange(numberOfGridPoints);
  solidTemperatureDot[numberOfGridPoints] = solidHeatExchange(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] = wallHeatExchange(numberOfGridPoints);

  gasTemperatureDot[numberOfGridPoints] -=
      interstitialGasVelocity[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints] - gasTemperature[numberOfGridPoints - 1]) * idx(numberOfGridPoints);
  const double outletTemperatureIdx2 = idx(numberOfGridPoints) * idx(numberOfGridPoints);
  gasTemperatureDot[numberOfGridPoints] +=
      coeffDiffusion[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints - 1] - gasTemperature[numberOfGridPoints]) * outletTemperatureIdx2;

  wallTemperatureDot[numberOfGridPoints] +=
      outletTemperatureIdx2 * (wallThermalConductivity * invHeatDensityWall) *
      (wallTemperature[numberOfGridPoints - 1] - wallTemperature[numberOfGridPoints]);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[numberOfGridPoints] +=
        heatOfAdsorption(numberOfGridPoints, comp) * spanAdsorptionDot[numberOfGridPoints, comp] / heatCapacitySolid;
  }
}

// void computeTVD(std::span<double> input, std::span<double> output, bool clamp)
// {
//   double tol = 1e-10;
//   double r_value, flux_limiter;

//   size_t size = input.size();
//   if (size < 3)
//   {
//     throw std::runtime_error("Unable to call TVD with numberOfGridPoints smaller than 3, increase number of grid
//     points.");
//   }

//   if (clamp && input[size - 1] >= 1.0) output[size - 1] = 1.0;

//   // For right wall of 1st Node, r_value is calculated using half-cell approximation
//   r_value = (2.0 * (input[1] - input[0]) + tol) / (input[2] - input[1] + tol);
//   flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
//   output[1] += 0.5 * flux_limiter * (input[2] - input[1]);

//   // For right walls of 2nd to numberOfGridPoints-1 node
//   for (size_t i = 2; i < size - 1; ++i)
//   {
//     r_value = ((input[i] - input[i - 1]) + tol) / ((input[i + 1] - input[i]) + tol);
//     flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
//     output[i] = input[i] + 0.5 * flux_limiter * (input[i + 1] - input[i]);
//   }
// }
