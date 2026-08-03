#include "compute.h"

#include <algorithm>
#include <mdspan>

#include "transport.h"
#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;

void updateVelocityAndPressure(const std::vector<Component>& components,
                               const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, const Geometry& geometry,
                               double& columnEntranceVelocity, double dynamicViscosity, double resolution,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> bulkSpeciesSink,
                               std::span<const double> gasTemperature)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; grid++)
  {
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; comp++)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }

    totalConcentration[grid] = concentrationSum;
    totalPressure[grid] = totalConcentration[grid] * R * std::max(1e-10, gasTemperature[grid]);
    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; comp++)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = concentration[index] * R * std::max(1e-10, gasTemperature[grid]);
      gasDensity[grid] += concentration[index] * components[comp].molecularWeight;
    }
  }

  auto ergunGrad = [&](size_t grid)
  {
    return geometry.pressureGradient(dynamicViscosity, gasDensity[grid], interstitialGasVelocity[grid]);
  };

  auto sinkTerm = [&](size_t grid)
  {
    const auto begin = static_cast<std::ptrdiff_t>(grid * numberOfComponents);
    const auto end = static_cast<std::ptrdiff_t>((grid + 1) * numberOfComponents);
    return resolution * std::reduce(bulkSpeciesSink.begin() + begin, bulkSpeciesSink.begin() + end);
  };

  auto gridRatio = [&](size_t grid) -> double
  { return numberOfGridPoints == 0 ? 0.0 : static_cast<double>(grid) / static_cast<double>(numberOfGridPoints); };

  auto refreshNode = [&](size_t grid)
  {
    totalConcentration[grid] = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid]));
    gasDensity[grid] = 0.0;
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = moleFraction[index] * totalPressure[grid];
      gasDensity[grid] += moleFraction[index] * totalConcentration[grid] * components[comp].molecularWeight;
    }
  };

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

  if (boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[0] = inletPressure;

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
      totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
  {
    auto shoot = [&](double v0)
    {
      interstitialGasVelocity[0] = v0;
      totalPressure[0] = inletPressure;

      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        const double cprev = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid - 1]));
        const double cnow = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid]));
        interstitialGasVelocity[grid] =
            (interstitialGasVelocity[grid - 1] * cprev - sinkTerm(grid)) / std::max(1e-10, cnow);
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
        if (totalPressure[grid] <= 0.0)
        {
          return -outletPressure;
        }
      }

      return totalPressure[numberOfGridPoints] - outletPressure;
    };

    columnEntranceVelocity = bisection(shoot);
    shoot(columnEntranceVelocity);
  }
  else if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[numberOfGridPoints] = outletPressure;

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
    }
    for (size_t grid = numberOfGridPoints; grid > 0; --grid)
    {
      const size_t current = grid - 1;
      totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * resolution;
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedVelocity)
  {
    std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), columnEntranceVelocity);

    if (inletPressure > 0.0)
    {
      totalPressure[0] = inletPressure;
      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
      }
    }
    else
    {
      totalPressure[numberOfGridPoints] = outletPressure;
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        const size_t current = grid - 1;
        totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * resolution;
      }
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
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

void updateVelocityAndPressure(const std::vector<Component>& components,
                               const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, double voidFraction,
                               double /*particleDensity*/, double& columnEntranceVelocity, double dynamicViscosity,
                               double particleDiameter, double resolution,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> bulkSpeciesSink,
                               std::span<const double> gasTemperature)
{
  const Geometry geometry{HollowTube{voidFraction, particleDiameter}};
  updateVelocityAndPressure(components, boundaryCondition, numberOfGridPoints, numberOfComponents,
                            inletPressure, outletPressure, pressureGradient, columnLength, geometry,
                            columnEntranceVelocity, dynamicViscosity, resolution, interstitialGasVelocity,
                            gasDensity, totalConcentration, totalPressure, concentration, partialPressure,
                            moleFraction, bulkSpeciesSink, gasTemperature);
}

void updateVelocityAndPressure(Column& column)
{
  computeBulkSpeciesSink(column);
  updateVelocityAndPressure(column.components, column.boundaryCondition, column.numberOfGridPoints,
                            column.numberOfComponents, column.inletPressure, column.outletPressure,
                            column.pressureGradient, column.columnLength, column.geometry,
                            column.columnEntranceVelocity, column.dynamicViscosity, column.resolution,
                            column.interstitialGasVelocity, column.gasDensity,
                            column.totalConcentration, column.totalPressure, column.concentration,
                            column.partialPressure, column.moleFraction, column.bulkSpeciesSink,
                            column.gasTemperature);
}

void computeEquilibriumLoadings(Column& column)
{
  computeEquilibriumLoadings(column.physisorptionMixture, column.numberOfGridPoints, column.numberOfComponents,
                             column.maxIsothermTerms, column.iastPerformance, column.idealGasMolFractions,
                             column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
                             column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
                             column.moleFraction, column.gasTemperature);

  computeChemisorptionEquilibriumLoadings(
      column.chemisorptionMixture, column.numberOfGridPoints, column.numberOfComponents,
      column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
      column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
      column.equilibriumChemisorption, column.cachedChemisorptionPressure,
      column.cachedChemisorptionGrandPotential, column.moleFraction, column.gasTemperature);
}

void computeEquilibriumLoadings(MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
                                size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<const double> totalPressure,
                                std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
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

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    std::span<double> spanCachedPressure =
        cachedPressure.subspan(grid * numberOfComponents * maxIsothermTerms, numberOfComponents * maxIsothermTerms);
    std::span<double> spanCachedGrandPotential =
        cachedGrandPotential.subspan(grid * maxIsothermTerms, maxIsothermTerms);
    iastPerformance +=
        mixture.predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions, numberOfMolecules,
                               spanCachedPressure, spanCachedGrandPotential, gasTemperature[grid]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumPhysisorption[grid * numberOfComponents + comp] = numberOfMolecules[comp];
    }
  }
}

void computeChemisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, std::pair<size_t, size_t>& iastPerformance,
    std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
    std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
    std::span<double> gasTemperature)
{
  std::fill(equilibriumChemisorption.begin(), equilibriumChemisorption.end(), 0.0);
  const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(0.0, moleFraction[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = sum > 0.0
                                       ? idealGasMolFractions[comp] / sum
                                       : 1.0 / static_cast<double>(numberOfComponents);
    }

    std::span<double> spanCachedPressure = cachedPressure.subspan(
        grid * numberOfComponents * maxChemisorptionSites,
        numberOfComponents * maxChemisorptionSites);
    std::span<double> spanCachedGrandPotential = cachedGrandPotential.subspan(
        grid * maxChemisorptionSites, maxChemisorptionSites);
    iastPerformance += mixture.predictMixture(
        idealGasMolFractions, totalPressure[grid], adsorbedMolFractions, numberOfMolecules,
        spanCachedPressure, spanCachedGrandPotential, gasTemperature[grid]);

    for (size_t site = 0; site < mixture.maxIsothermTerms; ++site)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumChemisorption[site * componentBlockSize + grid * numberOfComponents + comp] =
            mixture.equilibriumSiteLoadings[site * numberOfComponents + comp];
      }
    }
  }
}

void computeDerivatives(Column& column)
{
  computeMassDerivatives(column);
  if (column.energyBalance)
  {
    computeEnergyDerivatives(column);
  }
}

void computeMassDerivatives(Column& column)
{
  computeBulkSpeciesSink(column);
  computeMassDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents, column.resolution,
                         column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                         column.bulkSpeciesSink);
}

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double resolution, const ShapeParameters& geometry, double particleDensity,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot)
{
  const double prefactor = geometry.loadingPrefactor(particleDensity);
  std::vector<double> bulkSpeciesSink((numberOfGridPoints + 1) * numberOfComponents, 0.0);
  for (size_t i = 0; i < bulkSpeciesSink.size(); ++i)
  {
    double totalChemisorptionDot = 0.0;
    for (size_t site = 0; site < maxChemisorptionSites; ++site)
    {
      totalChemisorptionDot += chemisorptionDot[site * bulkSpeciesSink.size() + i];
    }
    bulkSpeciesSink[i] = prefactor * (physisorptionDot[i] + totalChemisorptionDot);
  }
  computeMassDerivatives(components, numberOfGridPoints, numberOfComponents, resolution, interstitialGasVelocity,
                         concentration, concentrationDot, bulkSpeciesSink);
}

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double resolution, double voidFraction,
                            double particleDensity, std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot)
{
  const Geometry geometry{HollowTube{voidFraction}};
  computeMassDerivatives(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                         resolution, geometry.shapeParameters(), particleDensity, interstitialGasVelocity,
                         concentration, concentrationDot, physisorptionDot, chemisorptionDot);
}

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, double resolution,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> bulkSpeciesSink)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = 0.0;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double dvcDz = (interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
                           interstitialGasVelocity[grid - 1] * spanConcentration[grid - 1, comp]) *
                           idx;
      const double d2cDz2 =
          (spanConcentration[grid + 1, comp] - 2.0 * spanConcentration[grid, comp] +
           spanConcentration[grid - 1, comp]) *
          idx2;

      spanConcentrationDot[grid, comp] = -dvcDz + components[comp].axialDispersionCoefficient * d2cDz2 -
                                         spanBulkSpeciesSink[grid, comp];
    }
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double dvcDz = (interstitialGasVelocity[numberOfGridPoints] *
                              spanConcentration[numberOfGridPoints, comp] -
                         interstitialGasVelocity[numberOfGridPoints - 1] *
                             spanConcentration[numberOfGridPoints - 1, comp]) *
                         idx;
    const double d2cDz2 =
        (spanConcentration[numberOfGridPoints - 1, comp] - spanConcentration[numberOfGridPoints, comp]) * idx2;

    spanConcentrationDot[numberOfGridPoints, comp] =
        -dvcDz + components[comp].axialDispersionCoefficient * d2cDz2 -
        spanBulkSpeciesSink[numberOfGridPoints, comp];
  }
}

void computeEnergyDerivatives(Column& column)
{
  computeEnergyDerivatives(
      column.components, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
      column.externalTemperature, column.geometry.shapeParameters(), column.particleDensity,
      column.wallDensity, column.gasThermalConductivity,
      column.wallThermalConductivity, column.heatTransferGasSolid, column.heatTransferGasWall,
      column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid, column.heatCapacityWall,
      column.resolution, column.interstitialGasVelocity, column.gasDensity, column.coeffDiffusion,
      column.physisorptionDot, column.chemisorptionDot, column.gasTemperature, column.gasTemperatureDot,
      column.solidTemperature, column.solidTemperatureDot, column.wallTemperature, column.wallTemperatureDot);
}

void computeEnergyDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, const ShapeParameters& geometry,
    double particleDensity, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> interstitialGasVelocity, std::span<const double> gasDensity,
    std::span<double> coeffDiffusion, std::span<const double> physisorptionDot,
    std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot)
{
  mdspan2d_const spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  std::mdspan<const double, std::dextents<size_t, 3>> spanChemisorptionDot(
      chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1, numberOfComponents);
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(numberOfGridPoints + 1);
  // computeWENO(gasTemperature, gasTemperatureFlux);

  const double prefactorGasSolid =
      geometry.fluidSolidContactAreaPerFluidVolume * heatTransferGasSolid / heatCapacityGas;

  // is this extra 1/eps necessary? It's in python not in eqs
  const double prefactorGasWall =
      geometry.fluidWallContactAreaPerFluidVolume * heatTransferGasWall / heatCapacityGas;

  const double coeffSolidGas =
      geometry.solidFluidContactAreaPerSolidVolume * heatTransferGasSolid /
      (heatCapacitySolid * particleDensity);

  // prefactor 4 in python, 2 in eqs
  double invHeatDensityWall = 1.0 / (heatCapacityWall * wallDensity);
  const double coeffWallGas = heatTransferGasWall * geometry.wallInnerContactAreaPerWallVolume * invHeatDensityWall;
  const double coeffWallExternal =
      heatTransferWallExternal * geometry.wallOuterContactAreaPerWallVolume * invHeatDensityWall;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  auto gasHeatExchange = [&](size_t grid)
  {
    const double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    const double gasSolidExchange = prefactorGasSolid * invGasDensity;
    const double gasWallExchange = prefactorGasWall * invGasDensity;
    return gasSolidExchange * (solidTemperature[grid] - gasTemperature[grid]) +
           gasWallExchange * (wallTemperature[grid] - gasTemperature[grid]);
  };

  auto solidHeatExchange = [&](size_t grid)
  {
    return coeffSolidGas * (gasTemperature[grid] - solidTemperature[grid]);
  };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  auto chemisorptionHeat = [&](size_t grid, size_t comp)
  {
    double heat = 0.0;
    const MultiSiteChemisorption& multisite = components[comp].chemisorption;
    for (size_t site = 0; site < multisite.numberOfSites; ++site)
    {
      heat += multisite.sites[site].heatOfChemisorption * spanChemisorptionDot[site, grid, comp];
    }
    return heat;
  };

  // first grid point
  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = solidHeatExchange(0);
  wallTemperatureDot[0] = wallHeatExchange(0);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[0] +=
        (components[comp].heatOfAdsorption * spanPhysisorptionDot[0, comp] +
         chemisorptionHeat(0, comp)) /
        heatCapacitySolid;
  }

  // middle grid points
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    // heat flux transfer
    gasTemperatureDot[grid] = gasHeatExchange(grid);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    // flux from gas diffusion and advection
    gasTemperatureDot[grid] -= interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx;
    gasTemperatureDot[grid] += coeffDiffusion[grid] *
                               (gasTemperature[grid - 1] - 2.0 * gasTemperature[grid] + gasTemperature[grid + 1]) *
                               idx2;

    // flux from heat diffusion in wall
    wallTemperatureDot[grid] += idx2 * (wallThermalConductivity * invHeatDensityWall) *
                                (wallTemperature[grid - 1] - 2.0 * wallTemperature[grid] + wallTemperature[grid + 1]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] +=
          (components[comp].heatOfAdsorption * spanPhysisorptionDot[grid, comp] +
           chemisorptionHeat(grid, comp)) /
          heatCapacitySolid;
    }
  }

  // last gridpoint
  // heat flux transfer
  gasTemperatureDot[numberOfGridPoints] = gasHeatExchange(numberOfGridPoints);
  solidTemperatureDot[numberOfGridPoints] = solidHeatExchange(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] = wallHeatExchange(numberOfGridPoints);

  // flux from gas diffusion and advection
  gasTemperatureDot[numberOfGridPoints] -=
      interstitialGasVelocity[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints] - gasTemperature[numberOfGridPoints - 1]) * idx;
  gasTemperatureDot[numberOfGridPoints] +=
      coeffDiffusion[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints - 1] - gasTemperature[numberOfGridPoints]) * idx2;

  // flux from heat diffusion in wall
  wallTemperatureDot[numberOfGridPoints] +=
      idx2 * (wallThermalConductivity * invHeatDensityWall) *
      (wallTemperature[numberOfGridPoints - 1] - wallTemperature[numberOfGridPoints]);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[numberOfGridPoints] +=
        (components[comp].heatOfAdsorption * spanPhysisorptionDot[numberOfGridPoints, comp] +
         chemisorptionHeat(numberOfGridPoints, comp)) /
        heatCapacitySolid;
  }
}

void computeEnergyDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, double voidFraction,
    double particleDensity, double particleDiameter,
    double internalDiameter, double outerDiameter, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> interstitialGasVelocity, std::span<const double> gasDensity,
    std::span<double> coeffDiffusion, std::span<const double> physisorptionDot,
    std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot)
{
  const Geometry geometry{HollowTube{voidFraction, particleDiameter, internalDiameter, outerDiameter}};
  computeEnergyDerivatives(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                           externalTemperature, geometry.shapeParameters(), particleDensity, wallDensity,
                           gasThermalConductivity, wallThermalConductivity, heatTransferGasSolid,
                           heatTransferGasWall, heatTransferWallExternal, heatCapacityGas, heatCapacitySolid,
                           heatCapacityWall, resolution, interstitialGasVelocity, gasDensity, coeffDiffusion,
                           physisorptionDot, chemisorptionDot, gasTemperature, gasTemperatureDot,
                           solidTemperature, solidTemperatureDot, wallTemperature, wallTemperatureDot);
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
