#include "compute.h"

#include <mdspan>

#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;

void computePressure(Column& column)
{
  computePressure(column.components, column.velocityProfile, column.boundaryCondition, column.numberOfGridPoints,
                  column.numberOfComponents, column.inletPressure, column.outletPressure, column.voidFraction,
                  column.dynamicViscosity, column.particleDiameter, column.resolution, column.interstitialGasVelocity,
                  column.gasDensity, column.totalConcentration, column.totalPressure, column.concentration,
                  column.partialPressure, column.moleFraction, column.gasTemperature);
}

void computePressure(const std::vector<Component>& components, const Column::VelocityProfile& velocityProfile,
                     const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                     size_t numberOfComponents, double inletPressure, double outletPressure, double voidFraction,
                     double dynamicViscosity, double particleDiameter, double resolution,
                     std::span<const double> interstitialGasVelocity, std::span<double> gasDensity,
                     std::span<double> totalConcentration, std::span<double> totalPressure,
                     std::span<double> concentration, std::span<double> partialPressure,
                     std::span<const double> moleFraction, std::span<const double> gasTemperature)
{
  double prefactor = (1 - voidFraction) / (voidFraction * particleDiameter);

  auto ergunGrad = [prefactor, dynamicViscosity, interstitialGasVelocity, gasDensity](size_t grid)
  {
    double visc = 150.0 * prefactor * prefactor * dynamicViscosity * interstitialGasVelocity[grid];
    double iner = 1.75 * prefactor * gasDensity[grid] * interstitialGasVelocity[grid] * interstitialGasVelocity[grid];
    return visc + iner;
  };

  if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; grid++)
    {
      totalPressure[grid] =
          inletPressure + grid * (outletPressure - inletPressure) / static_cast<double>(numberOfGridPoints);
    }
  }
  else
  {
    if (velocityProfile == Column::VelocityProfile::Ergun)
    {
      if (boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
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
    else if (boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
    {
      totalPressure[0] = inletPressure;
    }
    else if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure)
    {
      totalPressure[numberOfGridPoints] = outletPressure;
    }
  }

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const double invTemperature = 1.0 / std::max(1e-10, gasTemperature[grid]);
    totalConcentration[grid] = totalPressure[grid] / R * invTemperature;

    double moleFractionSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      moleFractionSum += std::max(0.0, moleFraction[grid * numberOfComponents + comp]);
    }
    const double invMoleFractionSum = moleFractionSum > 0.0 ? 1.0 / moleFractionSum : 0.0;

    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double yi = moleFractionSum > 0.0
                            ? std::max(0.0, moleFraction[grid * numberOfComponents + comp]) * invMoleFractionSum
                            : 1.0 / static_cast<double>(numberOfComponents);
      partialPressure[grid * numberOfComponents + comp] = yi * totalPressure[grid];
      concentration[grid * numberOfComponents + comp] = yi * totalConcentration[grid];
      gasDensity[grid] += concentration[grid * numberOfComponents + comp] * components[comp].molecularWeight;
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (totalPressure[numberOfGridPoints] < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeEquilibriumLoadings(Column& column)
{
  computeEquilibriumLoadings(column.mixture, column.numberOfGridPoints, column.numberOfComponents,
                             column.maxIsothermTerms, column.iastPerformance, column.idealGasMolFractions,
                             column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
                             column.equilibriumAdsorption, column.cachedPressure, column.cachedGrandPotential,
                             column.moleFraction, column.gasTemperature);
}

void computeEquilibriumLoadings(MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
                                size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance,
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

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    iastPerformance +=
        mixture.predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions, numberOfMolecules,
                               &cachedPressure[grid * numberOfComponents * maxIsothermTerms],
                               &cachedGrandPotential[grid * maxIsothermTerms], gasTemperature[grid]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumAdsorption[grid * numberOfComponents + comp] = numberOfMolecules[comp];
    }
  }
}

void computeVelocity(Column& column)
{
  switch (column.velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column.boundaryCondition, column.components, column.numberOfGridPoints,
                                   column.numberOfComponents, column.pressureGradient, column.columnEntranceVelocity,
                                   column.resolution, column.prefactorMassTransfer, column.interstitialGasVelocity,
                                   column.totalConcentration, column.totalPressure, column.equilibriumAdsorption,
                                   column.moleFraction, column.adsorption);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column.boundaryCondition, column.numberOfGridPoints, column.voidFraction,
                           column.columnEntranceVelocity, column.columnLength, column.dynamicViscosity,
                           column.particleDiameter, column.resolution, column.interstitialGasVelocity,
                           column.gasDensity, column.totalPressure);
      break;
    }
    case Column::VelocityProfile::FixedVelocity:
      break;
    default:
      break;
  }
}

void computeVelocityFixedGradient(const Column::BoundaryCondition& boundaryCondition,
                                  const std::vector<Component>& components, size_t numberOfGridPoints,
                                  size_t numberOfComponents, double pressureGradient, double columnEntranceVelocity,
                                  double resolution, std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalConcentration,
                                  std::span<const double> totalPressure, std::span<const double> equilibriumAdsorption,
                                  std::span<const double> moleFraction, std::span<const double> adsorption)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  if (boundaryCondition != Column::BoundaryCondition::InletPressureInletVelocity)
  {
    throw std::logic_error(
        "Fixed pressure gradient velocity profile only allowed to use in "
        "combination with input pressure boundary condition.");
  }

  // first grid point
  interstitialGasVelocity[0] = columnEntranceVelocity;

  // middle grid points
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    double sum = 0.0;

    auto concentrationSecondDerivative = [&](size_t grid, size_t comp) -> double
    {
      const double cLeft = totalConcentration[grid - 1] * moleFraction[(grid - 1) * numberOfComponents + comp];
      const double cCenter = totalConcentration[grid] * moleFraction[grid * numberOfComponents + comp];
      const double cRight = totalConcentration[grid + 1] * moleFraction[(grid + 1) * numberOfComponents + comp];
      return (cLeft - 2.0 * cCenter + cRight) * idx2;
    };

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      // mass transfer term
      sum -= prefactorMassTransfer[comp] *
             (equilibriumAdsorption[grid * numberOfComponents + comp] - adsorption[grid * numberOfComponents + comp]);

      sum += components[comp].axialDispersionCoefficient * concentrationSecondDerivative(grid, comp);
    }

    double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[grid]);

    // explicit update
    interstitialGasVelocity[grid] =
        interstitialGasVelocity[grid - 1] +
        resolution *
            (sum * invTotalConcentration - interstitialGasVelocity[grid - 1] * pressureGradient / totalPressure[grid]);
  }

  // last grid point
  double sum = 0.0;
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    sum -= prefactorMassTransfer[comp] * (equilibriumAdsorption[numberOfGridPoints * numberOfComponents + comp] -
                                          adsorption[numberOfGridPoints * numberOfComponents + comp]);

    const double cLeft =
        totalConcentration[numberOfGridPoints - 1] * moleFraction[(numberOfGridPoints - 1) * numberOfComponents + comp];
    const double cOutlet =
        totalConcentration[numberOfGridPoints] * moleFraction[numberOfGridPoints * numberOfComponents + comp];
    sum += components[comp].axialDispersionCoefficient * (cLeft - cOutlet) * idx2;
  }
  double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[numberOfGridPoints]);

  interstitialGasVelocity[numberOfGridPoints] =
      interstitialGasVelocity[numberOfGridPoints - 1] +
      resolution * (sum * invTotalConcentration - interstitialGasVelocity[numberOfGridPoints - 1] * pressureGradient /
                                                      totalPressure[numberOfGridPoints]);
}

void computeVelocityErgun(const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                          double voidFraction, double columnEntranceVelocity, double columnLength,
                          double dynamicViscosity, double particleDiameter, double resolution,
                          std::span<double> interstitialGasVelocity, std::span<const double> gasDensity,
                          std::span<const double> totalPressure)
{
  const double idx = 1.0 / resolution;

  const double voidFractionPrefactorA = (1.0 - voidFraction) / voidFraction;
  const double voidFractionPrefactorB = (1.0 - voidFraction) * (1.0 - voidFraction) / (voidFraction * voidFraction);

  const double prefactorA = 1.75 * voidFractionPrefactorA / particleDiameter;
  const double B = 150.0 * voidFractionPrefactorB * dynamicViscosity / (particleDiameter * particleDiameter);

  auto computedPdz = [&](size_t grid) -> double
  {
    if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
    {
      return (totalPressure[numberOfGridPoints] - totalPressure[0]) / columnLength;
    }

    return (totalPressure[grid] - totalPressure[grid - 1]) * idx;
  };

  auto solveVelocity = [&](size_t grid, double C) -> double
  {
    const double A = prefactorA * gasDensity[grid];

    if (std::abs(A) < 1e-10) return -C / B;

    const double discriminant = std::max(0.0, B * B - 4.0 * A * C);
    return (-B + std::sqrt(discriminant)) / (2.0 * A);
  };

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    interstitialGasVelocity[grid] = solveVelocity(grid, computedPdz(grid));
  }

  if (boundaryCondition != Column::BoundaryCondition::InletPressureOutletPressure)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
  }
}
void computeDerivatives(Column& column)
{
  if (column.energyBalance)
  {
    computeDerivativesEnergyBalance(
        column.components, column.numberOfGridPoints, column.numberOfComponents, column.externalTemperature,
        column.voidFraction, column.particleDensity, column.particleDiameter, column.internalDiameter,
        column.outerDiameter, column.wallDensity, column.gasThermalConductivity, column.wallThermalConductivity,
        column.heatTransferGasSolid, column.heatTransferGasWall, column.heatTransferWallExternal,
        column.heatCapacityGas, column.heatCapacitySolid, column.heatCapacityWall, column.resolution,
        column.prefactorMassTransfer, column.interstitialGasVelocity, column.gasDensity, column.totalConcentration,
        column.equilibriumAdsorption, column.coeffGasGas, column.coeffGasSolid, column.coeffGasWall,
        column.coeffDiffusion, column.moleFraction, column.moleFractionDot, column.adsorption, column.adsorptionDot,
        column.gasTemperature, column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot,
        column.wallTemperature, column.wallTemperatureDot);
  }
  else
  {
    computeDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents, column.resolution,
                       column.prefactorMassTransfer, column.interstitialGasVelocity, column.totalConcentration,
                       column.equilibriumAdsorption, column.moleFraction, column.moleFractionDot, column.adsorption,
                       column.adsorptionDot);
  }
}

void computeDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
                        double resolution, std::span<const double> prefactorMassTransfer,
                        std::span<const double> interstitialGasVelocity, std::span<const double> totalConcentration,
                        std::span<const double> equilibriumAdsorption, std::span<const double> moleFraction,
                        std::span<double> moleFractionDot, std::span<const double> adsorption,
                        std::span<double> adsorptionDot)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanMoleFraction(moleFraction.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanMoleFractionDot(moleFractionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanAdsorption(adsorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanAdsorptionDot(adsorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);

  // first gridpoint
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    spanAdsorptionDot[0, comp] = components[comp].massTransferCoefficient * diffAdsorption;
    spanMoleFractionDot[0, comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      spanAdsorptionDot[grid, comp] = components[comp].massTransferCoefficient * diffAdsorption;

      const double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[grid]);
      const double dctDz = (totalConcentration[grid] - totalConcentration[grid - 1]) * idx;
      const double d2ctDz2 =
          (totalConcentration[grid + 1] - 2.0 * totalConcentration[grid] + totalConcentration[grid - 1]) * idx2;
      const double dyDz = (spanMoleFraction[grid, comp] - spanMoleFraction[grid - 1, comp]) * idx;
      const double d2yDz2 =
          (spanMoleFraction[grid + 1, comp] - 2.0 * spanMoleFraction[grid, comp] + spanMoleFraction[grid - 1, comp]) *
          idx2;
      const double adsorptionSource = -prefactorMassTransfer[comp] * diffAdsorption;

      spanMoleFractionDot[grid, comp] =
          -interstitialGasVelocity[grid] * dyDz +
          components[comp].axialDispersionCoefficient * invTotalConcentration *
              (spanMoleFraction[grid, comp] * d2ctDz2 + 2.0 * dctDz * dyDz + totalConcentration[grid] * d2yDz2) +
          adsorptionSource * invTotalConcentration;
    }
  }

  // last gridpoint
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    double diffAdsorption =
        spanEquilibriumAdsorption[numberOfGridPoints, comp] - spanAdsorption[numberOfGridPoints, comp];

    spanAdsorptionDot[numberOfGridPoints, comp] = components[comp].massTransferCoefficient * diffAdsorption;

    const double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[numberOfGridPoints]);
    const double dctDz = (totalConcentration[numberOfGridPoints] - totalConcentration[numberOfGridPoints - 1]) * idx;
    const double d2ctDz2 = (totalConcentration[numberOfGridPoints - 1] - totalConcentration[numberOfGridPoints]) * idx2;
    const double dyDz =
        (spanMoleFraction[numberOfGridPoints, comp] - spanMoleFraction[numberOfGridPoints - 1, comp]) * idx;
    const double d2yDz2 =
        (spanMoleFraction[numberOfGridPoints - 1, comp] - spanMoleFraction[numberOfGridPoints, comp]) * idx2;
    const double adsorptionSource = -prefactorMassTransfer[comp] * diffAdsorption;

    spanMoleFractionDot[numberOfGridPoints, comp] =
        -interstitialGasVelocity[numberOfGridPoints] * dyDz +
        components[comp].axialDispersionCoefficient * invTotalConcentration *
            (spanMoleFraction[numberOfGridPoints, comp] * d2ctDz2 + 2.0 * dctDz * dyDz +
             totalConcentration[numberOfGridPoints] * d2yDz2) +
        adsorptionSource * invTotalConcentration;
  }
}

void computeDerivativesEnergyBalance(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    double externalTemperature, double voidFraction, double particleDensity, double particleDiameter,
    double internalDiameter, double outerDiameter, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> prefactorMassTransfer, std::span<const double> interstitialGasVelocity,
    std::span<const double> gasDensity, std::span<const double> totalConcentration,
    std::span<const double> equilibriumAdsorption, std::span<double> coeffGasGas, std::span<double> coeffGasSolid,
    std::span<double> coeffGasWall, std::span<double> coeffDiffusion, std::span<const double> moleFraction,
    std::span<double> moleFractionDot, std::span<const double> adsorption, std::span<double> adsorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanAdsorption(adsorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanAdsorptionDot(adsorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanMoleFraction(moleFraction.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanMoleFractionDot(moleFractionDot.data(), numberOfGridPoints + 1, numberOfComponents);

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(numberOfGridPoints + 1);
  // computeWENO(gasTemperature, gasTemperatureFlux);

  double relativeVolume = ((1 - voidFraction) / voidFraction);
  double accessibleSurface = 6.0 / particleDiameter;  // prefactor 6.0 in python and 2.0 in eqs

  double prefactorGasSolid = relativeVolume * accessibleSurface * heatTransferGasSolid / heatCapacityGas;

  // is this extra 1/eps necessary? It's in python not in eqs
  double prefactorGasWall = 4.0 * heatTransferGasWall / (heatCapacityGas * internalDiameter);
  double prefactorGasGas = -(prefactorGasSolid + prefactorGasWall);

  double coeffSolidGas = accessibleSurface * heatTransferGasSolid / (heatCapacitySolid * particleDensity);
  double coeffSolidSolid = -coeffSolidGas;

  // prefactor 4 in python, 2 in eqs
  double internalArea = 4.0 * internalDiameter / (outerDiameter * outerDiameter - internalDiameter * internalDiameter);
  double externalArea = 4.0 * outerDiameter / (outerDiameter * outerDiameter - internalDiameter * internalDiameter);
  double invHeatDensityWall = 1.0 / (heatCapacityWall * wallDensity);
  double coeffWallGas = heatTransferGasWall * internalArea * invHeatDensityWall;
  double coeffWallWall = -coeffWallGas - heatTransferWallExternal * externalArea * invHeatDensityWall;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffGasGas[grid] = prefactorGasGas * invGasDensity;
    coeffGasSolid[grid] = prefactorGasSolid * invGasDensity;
    coeffGasWall[grid] = prefactorGasWall * invGasDensity;
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  // first grid point
  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = coeffSolidGas * gasTemperature[0] + coeffSolidSolid * solidTemperature[0];

  wallTemperatureDot[0] = coeffWallGas * gasTemperature[0] + coeffWallWall * wallTemperature[0];

  // add effect of ambient temperature to wall
  wallTemperatureDot[0] += heatTransferWallExternal * externalArea * externalTemperature * invHeatDensityWall;

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    spanAdsorptionDot[0, comp] = components[comp].massTransferCoefficient * diffAdsorption;
    spanMoleFractionDot[0, comp] = 0.0;

    solidTemperatureDot[0] += components[comp].heatOfAdsorption * spanAdsorptionDot[0, comp] / heatCapacitySolid;
  }

  // middle grid points
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    // heat flux transfer
    gasTemperatureDot[grid] = coeffGasGas[grid] * gasTemperature[grid] + coeffGasSolid[grid] * solidTemperature[grid] +
                              coeffGasWall[grid] * wallTemperature[grid];
    solidTemperatureDot[grid] = coeffSolidGas * gasTemperature[grid] + coeffSolidSolid * solidTemperature[grid];

    wallTemperatureDot[grid] = coeffWallGas * gasTemperature[grid] + coeffWallWall * wallTemperature[grid];

    // flux from gas diffusion and advection
    gasTemperatureDot[grid] -= interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx;
    gasTemperatureDot[grid] += coeffDiffusion[grid] *
                               (gasTemperature[grid - 1] - 2.0 * gasTemperature[grid] + gasTemperature[grid + 1]) *
                               idx2;

    // flux from heat diffusion in wall
    wallTemperatureDot[grid] += idx2 * (wallThermalConductivity * invHeatDensityWall) *
                                (wallTemperature[grid - 1] - 2.0 * wallTemperature[grid] + wallTemperature[grid + 1]);

    // add effect of ambient temperature to wall
    wallTemperatureDot[grid] += heatTransferWallExternal * externalArea * externalTemperature * invHeatDensityWall;

    double invGasTemperature = 1.0 / std::max(1e-10, gasTemperature[grid]);
    double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[grid]);
    double dctDz = (totalConcentration[grid] - totalConcentration[grid - 1]) * idx;
    double dTgDz = (gasTemperature[grid] - gasTemperature[grid - 1]) * idx;

    // mass balance
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      spanAdsorptionDot[grid, comp] = components[comp].massTransferCoefficient * diffAdsorption;

      // now add heat exchange from heat of adsorption
      solidTemperatureDot[grid] +=
          components[comp].heatOfAdsorption * spanAdsorptionDot[grid, comp] / heatCapacitySolid;

      const double dyDz = (spanMoleFraction[grid, comp] - spanMoleFraction[grid - 1, comp]) * idx;
      const double d2yDz2 =
          (spanMoleFraction[grid + 1, comp] - 2.0 * spanMoleFraction[grid, comp] + spanMoleFraction[grid - 1, comp]) *
          idx2;
      const double adsorptionSource = -prefactorMassTransfer[comp] * diffAdsorption;

      spanMoleFractionDot[grid, comp] =
          -interstitialGasVelocity[grid] * dyDz +
          components[comp].axialDispersionCoefficient *
              (invTotalConcentration * dctDz * dyDz + invGasTemperature * dTgDz * dyDz + d2yDz2) +
          adsorptionSource * invTotalConcentration;
    }
  }

  // last gridpoint
  // heat flux transfer
  gasTemperatureDot[numberOfGridPoints] = coeffGasGas[numberOfGridPoints] * gasTemperature[numberOfGridPoints] +
                                          coeffGasSolid[numberOfGridPoints] * solidTemperature[numberOfGridPoints] +
                                          coeffGasWall[numberOfGridPoints] * wallTemperature[numberOfGridPoints];
  solidTemperatureDot[numberOfGridPoints] =
      coeffSolidGas * gasTemperature[numberOfGridPoints] + coeffSolidSolid * solidTemperature[numberOfGridPoints];

  wallTemperatureDot[numberOfGridPoints] =
      coeffWallGas * gasTemperature[numberOfGridPoints] + coeffWallWall * wallTemperature[numberOfGridPoints];

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

  // add effect of ambient temperature to wall
  wallTemperatureDot[numberOfGridPoints] +=
      heatTransferWallExternal * externalArea * externalTemperature * invHeatDensityWall;

  double invGasTemperature = 1.0 / std::max(1e-10, gasTemperature[numberOfGridPoints]);
  double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[numberOfGridPoints]);
  double dctDz = (totalConcentration[numberOfGridPoints] - totalConcentration[numberOfGridPoints - 1]) * idx;
  double dTgDz = (gasTemperature[numberOfGridPoints] - gasTemperature[numberOfGridPoints - 1]) * idx;
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    double diffAdsorption =
        spanEquilibriumAdsorption[numberOfGridPoints, comp] - spanAdsorption[numberOfGridPoints, comp];

    spanAdsorptionDot[numberOfGridPoints, comp] = components[comp].massTransferCoefficient * diffAdsorption;

    // now add heat exchange from heat of adsorption
    solidTemperatureDot[numberOfGridPoints] +=
        components[comp].heatOfAdsorption * spanAdsorptionDot[numberOfGridPoints, comp] / heatCapacitySolid;

    const double dyDz =
        (spanMoleFraction[numberOfGridPoints, comp] - spanMoleFraction[numberOfGridPoints - 1, comp]) * idx;
    const double d2yDz2 =
        (spanMoleFraction[numberOfGridPoints - 1, comp] - spanMoleFraction[numberOfGridPoints, comp]) * idx2;
    const double adsorptionSource = -prefactorMassTransfer[comp] * diffAdsorption;

    spanMoleFractionDot[numberOfGridPoints, comp] =
        -interstitialGasVelocity[numberOfGridPoints] * dyDz +
        components[comp].axialDispersionCoefficient *
            (invTotalConcentration * dctDz * dyDz + invGasTemperature * dTgDz * dyDz + d2yDz2) +
        adsorptionSource * invTotalConcentration;
  }
}

void enforceBoundaryCondition(Column& column)
{
  switch (column.boundaryCondition)
  {
    case Column::BoundaryCondition::InletPressureInletVelocity:
    {
      for (size_t j = 0; j < column.numberOfComponents; ++j)
      {
        column.moleFraction[0 * column.numberOfComponents + j] = column.components[j].initialGasMoleFraction;
        column.totalPressure[0] = column.inletPressure;
        column.totalConcentration[0] = column.totalPressure[0] / (R * std::max(1e-10, column.gasTemperature[0]));
        column.partialPressure[0 * column.numberOfComponents + j] =
            column.components[j].initialGasMoleFraction * column.totalPressure[0];
        column.concentration[0 * column.numberOfComponents + j] =
            column.moleFraction[0 * column.numberOfComponents + j] * column.totalConcentration[0];
        column.moleFractionDot[0 * column.numberOfComponents + j] = 0.0;
      }
      break;
    }
    default:
      break;
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
