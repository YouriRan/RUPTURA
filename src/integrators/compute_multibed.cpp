#include <algorithm>
#include <mdspan>

#include "compute_multibed.h"
#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;

void updateVelocityAndPressure(const std::vector<Component>& components,
                               const ColumnMultibed::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, size_t numberOfAdsorbents,
                               std::span<const double> particleDensities, double& columnEntranceVelocity,
                               double dynamicViscosity, std::span<const double> columnDistances,
                               std::span<const double> fractionOfAdsorbent,
                               std::span<const double> adsorbentScaledVoidFraction,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> physisorptionDot,
                               std::span<const double> gasTemperature)
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
    double adsorbentSinkPrefactor = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      adsorbentSinkPrefactor += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                                adsorbentScaledVoidFraction[grid * numberOfAdsorbents + ads] *
                                particleDensities[ads];
    }

    return gridSpacing(columnDistances, grid) * adsorbentSinkPrefactor *
           std::reduce(physisorptionDot.begin() + begin, physisorptionDot.begin() + end);
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
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
      totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
      refreshNode(grid);
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
      totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
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
        totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
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

void updateVelocityAndPressure(ColumnMultibed& column)
{
  updateVelocityAndPressure(column.components, column.boundaryCondition, column.numberOfGridPoints,
                            column.numberOfComponents, column.inletPressure, column.outletPressure,
                            column.pressureGradient, column.columnLength, column.numberOfAdsorbents,
                            column.particleDensities, column.columnEntranceVelocity, column.dynamicViscosity,
                            column.columnDistances, column.fractionOfAdsorbent, column.adsorbentScaledVoidFraction,
                            column.interstitialGasVelocity,
                            column.gasDensity, column.totalConcentration, column.totalPressure, column.concentration,
                            column.partialPressure, column.moleFraction, column.physisorptionDot, column.gasTemperature);
}

void computeEquilibriumLoadings(ColumnMultibed& column)
{
  computeEquilibriumLoadings(column.mixture, column.numberOfGridPoints, column.numberOfComponents,
                             column.numberOfAdsorbents, column.fractionOfAdsorbent, column.hasAdsorbentOfType,
                             column.maxIsothermTerms, column.iastPerformance, column.idealGasMolFractions,
                             column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
                             column.equilibriumAdsorption, column.cachedPressure, column.cachedGrandPotential,
                             column.moleFraction, column.gasTemperature);
}

void computeEquilibriumLoadings(std::vector<MixturePrediction>& mixture, size_t numberOfGridPoints,
                                size_t numberOfComponents, size_t numberOfAdsorbents,
                                std::span<const double> fractionOfAdsorbent,
                                const std::vector<bool>& hasAdsorbentOfType, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<const double> totalPressure,
                                std::span<double> equilibriumAdsorption, std::span<double> cachedPressure,
                                std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
                                std::span<double> gasTemperature)
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
      equilibriumAdsorption[grid * numberOfComponents + comp] = 0.0;
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
      iastPerformance +=
          mixture[ads].predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions,
                                      numberOfMolecules, spanCachedPressure, spanCachedGrandPotential,
                                      gasTemperature[grid]);

      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumAdsorption[grid * numberOfComponents + comp] +=
            fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * numberOfMolecules[comp];
      }
    }
  }
}

void computeDerivatives(ColumnMultibed& column)
{
  computeMassDerivatives(column);
  if (column.energyBalance)
  {
    computeEnergyDerivatives(column);
  }
}

void computeMassDerivatives(ColumnMultibed& column)
{
  computeMassDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.columnDistances, column.totalVoidFraction, column.particleDensity,
                         column.interstitialGasVelocity,
                         column.concentration, column.concentrationDot, column.physisorptionDot);
}

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents,
                            std::span<const double> columnDistances,
                            std::span<const double> totalVoidFraction,
                            std::span<const double> particleDensity,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanAdsorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);

  auto physisorptionPrefactor = [&](size_t grid)
  {
    return particleDensity[grid] * (1.0 - totalVoidFraction[grid]) / std::max(1e-10, totalVoidFraction[grid]);
  };

  auto secondDerivative = [&](size_t grid, size_t comp)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 * (hm * spanConcentration[grid + 1, comp] - (hm + hp) * spanConcentration[grid, comp] +
                  hp * spanConcentration[grid - 1, comp]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = 0.0;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double hm = gridSpacing(columnDistances, grid);
      const double dvcDz = (interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
                           interstitialGasVelocity[grid - 1] * spanConcentration[grid - 1, comp]) /
                           std::max(1e-30, hm);

      spanConcentrationDot[grid, comp] =
          -dvcDz + components[comp].axialDispersionCoefficient * secondDerivative(grid, comp) -
          physisorptionPrefactor(grid) * spanAdsorptionDot[grid, comp];
    }
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double hm = gridSpacing(columnDistances, numberOfGridPoints);
    const double invHm = 1.0 / std::max(1e-30, hm);
    const double dvcDz = (interstitialGasVelocity[numberOfGridPoints] *
                              spanConcentration[numberOfGridPoints, comp] -
                          interstitialGasVelocity[numberOfGridPoints - 1] *
                              spanConcentration[numberOfGridPoints - 1, comp]) *
                         invHm;
    const double d2cDz2 =
        (spanConcentration[numberOfGridPoints - 1, comp] - spanConcentration[numberOfGridPoints, comp]) * invHm * invHm;

    spanConcentrationDot[numberOfGridPoints, comp] =
        -dvcDz + components[comp].axialDispersionCoefficient * d2cDz2 -
        physisorptionPrefactor(numberOfGridPoints) * spanAdsorptionDot[numberOfGridPoints, comp];
  }
}

void computeEnergyDerivatives(ColumnMultibed& column)
{
  computeEnergyDerivatives(
      column.components, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.externalTemperature, column.totalVoidFraction, column.particleDensities, column.particleDiameters,
      column.fractionOfAdsorbent, column.internalDiameter, column.outerDiameter, column.wallDensity,
      column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
      column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
      column.heatCapacityWall, column.columnDistances, column.interstitialGasVelocity, column.gasDensity,
      column.coeffDiffusion, column.physisorptionDot, column.gasTemperature, column.gasTemperatureDot,
      column.solidTemperature, column.solidTemperatureDot, column.wallTemperature, column.wallTemperatureDot);
}

void computeEnergyDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, double externalTemperature, std::span<const double> totalVoidFraction,
    std::span<const double> particleDensities, std::span<const double> particleDiameters,
    std::span<const double> fractionOfAdsorbent, double internalDiameter, double outerDiameter, double wallDensity,
    double gasThermalConductivity, double wallThermalConductivity, double heatTransferGasSolid,
    double heatTransferGasWall, double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
    double heatCapacityWall, std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
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
    return accessibleSurface(grid) * heatTransferGasSolid /
           (heatCapacitySolid * std::max(1e-10, solidDensity(grid)));
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
  {
    return coeffSolidGasAt(grid) * (gasTemperature[grid] - solidTemperature[grid]);
  };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = solidHeatExchange(0);
  wallTemperatureDot[0] = wallHeatExchange(0);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[0] += components[comp].heatOfAdsorption * spanAdsorptionDot[0, comp] / heatCapacitySolid;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    gasTemperatureDot[grid] = gasHeatExchange(grid);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    gasTemperatureDot[grid] -=
        interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx(grid);
    gasTemperatureDot[grid] += coeffDiffusion[grid] * d2Temperature(gasTemperature, grid);

    wallTemperatureDot[grid] +=
        (wallThermalConductivity * invHeatDensityWall) * d2Temperature(wallTemperature, grid);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] +=
          components[comp].heatOfAdsorption * spanAdsorptionDot[grid, comp] / heatCapacitySolid;
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
        components[comp].heatOfAdsorption * spanAdsorptionDot[numberOfGridPoints, comp] / heatCapacitySolid;
  }
}

void computeWENO(std::span<const double> input, std::span<double> output)
{
  double tol = 1e-10;
  double df0, df1, alpha_0, alpha_1, beta_0, beta_1, first_term, second_term;

  size_t size = input.size();
  if (size < 3)
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 3.");
  }

  // inlet boundary flux: prescribed from Dirichlet inflow state
  output[0] = input[0];

  // first interior interface, special one-sided closure
  df0 = input[2] - input[1];
  df1 = input[1] - input[0];
  beta_0 = df0 * df0;
  beta_1 = df1 * df1;

  alpha_0 = (2.0 / 3.0) / ((beta_0 + tol) * (beta_0 + tol));
  alpha_1 = (1.0 / 3.0) / (16.0 * (beta_1 + tol) * (beta_1 + tol));

  first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (input[2] + input[1]);
  second_term = (alpha_1 / (alpha_0 + alpha_1)) * (2.0 * input[1] - input[0]);
  output[1] = first_term + second_term;

  // interior interfaces
  for (size_t i = 2; i < size - 1; ++i)
  {
    df0 = input[i + 1] - input[i];
    df1 = input[i] - input[i - 1];
    beta_0 = df0 * df0;
    beta_1 = df1 * df1;

    alpha_0 = (2.0 / 3.0) / ((beta_0 + tol) * (beta_0 + tol));
    alpha_1 = (1.0 / 3.0) / ((beta_1 + tol) * (beta_1 + tol));

    first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (input[i + 1] + input[i]);
    second_term = (alpha_1 / (alpha_0 + alpha_1)) * ((3.0 / 2.0) * input[i] - (1.0 / 2.0) * input[i - 1]);
    output[i] = first_term + second_term;
  }

  // outlet boundary flux: outflow closure
  output[size - 1] = input[size - 1];  // simplest option
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
