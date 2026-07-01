#include <mdspan>

#include "compute.h"
#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;

void computePressure(Column& column)
{
  computePressure(column.velocityProfile, column.boundaryCondition, column.components, column.Ngrid, column.Ncomp,
                  column.inletPressure, column.outletPressure, column.voidFraction, column.dynamicViscosity,
                  column.particleDiameter, column.resolution, column.interstitialGasVelocity, column.gasDensity,
                  column.totalConcentration, column.totalPressure, column.gasTemperature, column.concentration,
                  column.partialPressure, column.moleFraction);
}

void computePressure(Column::VelocityProfile& velocityProfile, Column::BoundaryCondition& boundaryCondition,
                     const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, size_t Nads,
                     double inletPressure, double outletPressure, std::span<const double> fractionOfAdsorbent,
                     std::span<const double> adsorbentScaledVoidFractions, double dynamicViscosity,
                     std::span<const double> particleDiameter, double resolution,
                     std::span<const double> interstitialGasVelocity, std::span<double> gasDensity,
                     std::span<double> totalConcentration, std::span<double> totalPressure,
                     std::span<double> gasTemperature, std::span<const double> concentration,
                     std::span<double> partialPressure, std::span<double> moleFraction)
{
  std::vector<double> prefactors((Ngrid + 1) * Nads);
  for (std::size_t grid = 0; grid < Ngrid + 1; grid++)
  {
    for (std::size_t ads = 0; ads < Nads; ads++)
    {
      prefactors[grid * Nads + ads] = adsorbentScaledVoidFractions[grid * Nads + ads] / particleDiameter[ads];
    }
  }

  auto ergunGrad = [Nads, prefactors, dynamicViscosity, interstitialGasVelocity, gasDensity](size_t grid)
  {
    double visc = 150.0 * dynamicViscosity * dynamicViscosity * interstitialGasVelocity[grid] *
                  std::transform_reduce(prefactors.begin() + grid * Nads, prefactors.begin() + (grid + 1) * Nads, 0.0,
                                        std::plus<>(), [](double acc, double x) { return acc + x * x; });
    double iner = 1.75 * gasDensity[grid] * interstitialGasVelocity[grid] * interstitialGasVelocity[grid] *
                  std::reduce(prefactors.begin() + grid * Nads, prefactors.begin() + (grid + 1) * Nads, 0.0);
    return visc + iner;
  };

  if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
  {
    for (size_t grid = 0; grid < Ngrid + 1; grid++)
    {
      totalPressure[grid] = inletPressure + grid * (outletPressure - inletPressure) / static_cast<double>(Ngrid);
    }
  }
  else
  {
    if (velocityProfile == Column::VelocityProfile::Ergun)
    {
      if (boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
      {
        totalPressure[0] = inletPressure;
        for (size_t grid = 1; grid < Ngrid + 1; ++grid)
        {
          totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
        }
      }
      else
      {
        totalPressure[Ngrid] = outletPressure;
        for (int grid = Ngrid - 1; grid >= 0; grid--)
        {
          totalPressure[grid] = totalPressure[grid + 1] + ergunGrad(static_cast<size_t>(grid)) * resolution;
        }
      }
    }
    else
    {
      totalPressure[0] = inletPressure;
      for (size_t grid = 1; grid < Ngrid + 1; ++grid)
      {
        totalPressure[grid] = 0.0;
        for (size_t comp = 0; comp < Ncomp; ++comp)
        {
          totalPressure[grid] += concentration[grid * Ncomp + comp] * R * gasTemperature[grid];
        }
      }
    }
  }

  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    totalConcentration[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      totalConcentration[grid] += std::max(0.0, concentration[grid * Ncomp + comp]);
    }

    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      gasDensity[grid] += concentration[grid * Ncomp + comp] * components[comp].molecularWeight;
      partialPressure[grid * Ncomp + comp] = moleFraction[grid * Ncomp + comp] * totalPressure[grid];
      if (totalConcentration[grid] > 0.0)
      {
        moleFraction[grid * Ncomp + comp] =
            std::max(0.0, concentration[grid * Ncomp + comp]) / totalConcentration[grid];
      }
      else
      {
        moleFraction[grid * Ncomp + comp] = 1.0 / Ncomp;
      }
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (totalPressure[Ngrid] < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeEquilibriumLoadings(Column& column)
{
  computeEquilibriumLoadings(
      column.mixture, column.Ngrid, column.Ncomp, column.Nads, column.maxIsothermTerms, column.iastPerformance, column.totalPressure,
      column.gasTemperature, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
      column.equilibriumAdsorption, column.moleFraction, column.cachedPressure, column.cachedGrandPotential);
}

void computeEquilibriumLoadings(MixturePrediction& mixture, size_t Ngrid, size_t Ncomp, size_t Nads, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance, std::span<double> totalPressure,
                                std::span<double> gasTemperature, std::span<double> idealGasMolFractions,
                                std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules,
                                std::span<double> equilibriumAdsorption, std::span<double> moleFraction,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential)
{
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      idealGasMolFractions[comp] = moleFraction[grid * Ncomp + comp];
      equilibriumAdsorption[grid * Ncomp + comp] = 0.0;
    }

    for (size_t ads = 0; ads < Nads; ++ads)
    {
      std::span<double> spanCachedPressure{cachedPressure.data(), (grid * Nads + ads) * Ncomp * maxIsothermTerms,
                                           Ncomp * maxIsothermTerms};
      std::span<double> spanGrandPotential{cachedGrandPotential.data(), (grid * Nads + ads) * maxIsothermTerms,
                                           maxIsothermTerms};
      // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
      iastPerformance +=
          mixture[ads].predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions,
                                      numberOfMolecules, spanCachedPressure, spanGrandPotential, gasTemperature[grid]);

      for (size_t comp = 0; comp < Ncomp; ++comp)
      {
        equilibriumAdsorption[grid * Ncomp + comp] += numberOfMolecules[comp];
      }
    }
  }
}

void computeVelocity(Column& column)
{
  switch (column.velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column.boundaryCondition, column.components, column.Ngrid, column.Ncomp,
                                   column.pressureGradient, column.columnEntranceVelocity, column.resolution,
                                   column.prefactorMassTransfer, column.interstitialGasVelocity, column.totalPressure,
                                   column.totalConcentration, column.adsorption, column.equilibriumAdsorption,
                                   column.concentration);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column.boundaryCondition, column.Ngrid, column.voidFraction, column.columnEntranceVelocity,
                           column.columnLength, column.dynamicViscosity, column.particleDiameter, column.resolution,
                           column.interstitialGasVelocity, column.gasDensity, column.totalPressure);
      break;
    }
    case Column::VelocityProfile::FixedVelocity:
      break;
    default:
      break;
  }
}

void computeVelocityFixedGradient(Column::BoundaryCondition& boundaryCondition,
                                  const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                                  double pressureGradient, double columnEntranceVelocity, double resolution,
                                  std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                                  std::span<const double> totalConcentration, std::span<const double> adsorption,
                                  std::span<const double> equilibriumAdsorption, std::span<const double> concentration)
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
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    double sum = 0.0;

    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      // mass transfer term
      sum -=
          prefactorMassTransfer[comp] * (equilibriumAdsorption[grid * Ncomp + comp] - adsorption[grid * Ncomp + comp]);

      // diffusion term in concentration form
      sum += components[comp].D *
             (concentration[(grid - 1) * Ncomp + comp] - 2.0 * concentration[grid * Ncomp + comp] +
              concentration[(grid + 1) * Ncomp + comp]) *
             idx2;
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
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    sum -=
        prefactorMassTransfer[comp] * (equilibriumAdsorption[Ngrid * Ncomp + comp] - adsorption[Ngrid * Ncomp + comp]);

    sum +=
        components[comp].D * (concentration[(Ngrid - 1) * Ncomp + comp] - concentration[Ngrid * Ncomp + comp]) * idx2;
  }
  double invTotalConcentration = 1.0 / std::max(1e-10, totalConcentration[Ngrid]);

  interstitialGasVelocity[Ngrid] =
      interstitialGasVelocity[Ngrid - 1] +
      resolution *
          (sum * invTotalConcentration - interstitialGasVelocity[Ngrid - 1] * pressureGradient / totalPressure[Ngrid]);
}

void computeVelocityErgun(Column::BoundaryCondition& boundaryCondition, size_t Ngrid, double voidFraction,
                          double columnEntranceVelocity, double columnLength, double dynamicViscosity,
                          double particleDiameter, double resolution, std::span<double> interstitialGasVelocity,
                          std::span<double> gasDensity, std::span<double> totalPressure)
{
  const double idx = 1.0 / resolution;


  auto computedPdz = [&](size_t grid) -> double
  {
    if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
    {
      return (totalPressure[Ngrid] - totalPressure[0]) / columnLength;
    }

    return (totalPressure[grid] - totalPressure[grid - 1]) * idx;
  };

  auto solveVelocity = [&](size_t grid, double C) -> double
  {
    const double A = gasDensity[grid] * 1.75 * adsorbentVoidFractions / particleDiameter;
    const double B = 150.0 * adsorbentVoidFractions * adsorbentVoidFractions * dynamicViscosity / (particleDiameter * particleDiameter);

    if (std::abs(A) < 1e-10) return -C / B;

    const double discriminant = std::max(0.0, B * B - 4.0 * A * C);
    return (-B + std::sqrt(discriminant)) / (2.0 * A);
  };

  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
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
        column.components, column.Ngrid, column.Ncomp, column.externalTemperature, column.voidFraction,
        column.particleDensity, column.particleDiameter, column.influxTemperature, column.internalDiameter,
        column.outerDiameter, column.wallDensity, column.gasThermalConductivity, column.wallThermalConductivity,
        column.heatTransferGasSolid, column.heatTransferGasWall, column.heatTransferWallExternal,
        column.heatCapacityGas, column.heatCapacitySolid, column.heatCapacityWall, column.resolution,
        column.prefactorMassTransfer, column.interstitialGasVelocity, column.totalPressure, column.gasTemperature,
        column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot, column.wallTemperature,
        column.wallTemperatureDot, column.concentration, column.concentrationDot, column.adsorption,
        column.adsorptionDot, column.equilibriumAdsorption, column.moleFraction, column.gasDensity, column.coeffGasGas,
        column.coeffGasSolid, column.coeffGasWall, column.coeffDiffusion, column.facePressures, column.massFlux);
  }
  else
  {
    computeDerivatives(column.components, column.Ngrid, column.Ncomp, column.resolution, column.prefactorMassTransfer,
                       column.interstitialGasVelocity, column.concentration, column.concentrationDot, column.adsorption,
                       column.adsorptionDot, column.equilibriumAdsorption, column.massFlux);
  }
}

void computeDerivatives(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double resolution,
                        std::span<const double> prefactorMassTransfer, std::span<const double> interstitialGasVelocity,
                        std::span<const double> concentration, std::span<double> concentrationDot,
                        std::span<const double> adsorption, std::span<double> adsorptionDot,
                        std::span<const double> equilibriumAdsorption, std::span<double> massFlux)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanConcentration(concentration.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), Ngrid + 1, Ncomp);
  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), Ngrid + 1, Ncomp);
  mdspan2d_const spanAdsorption(adsorption.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanAdsorptionDot(adsorptionDot.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanMassFlux(massFlux.data(), Ngrid + 1, Ncomp);

  // commented out the parts for weno, seems to be unstable
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    // std::vector<double> compVelocityPressure(Ngrid + 1);
    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      // compVelocityPressure[grid] = interstitialGasVelocity[grid] * spanPartialPressure[grid, comp];
      spanMassFlux[grid, comp] = interstitialGasVelocity[grid] * spanConcentration[grid, comp];
    }
  }

  // first gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    spanAdsorptionDot[0, comp] = components[comp].Kl * diffAdsorption;
    spanConcentrationDot[0, comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      spanAdsorptionDot[grid, comp] = components[comp].Kl * diffAdsorption;

      spanConcentrationDot[grid, comp] = idx * (spanMassFlux[grid - 1, comp] - spanMassFlux[grid, comp]) +
                                         components[comp].D *
                                             (spanConcentration[grid + 1, comp] - 2.0 * spanConcentration[grid, comp] +
                                              spanConcentration[grid - 1, comp]) *
                                             idx2 -
                                         prefactorMassTransfer[comp] * diffAdsorption;
    }
  }

  // last gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    spanAdsorptionDot[Ngrid, comp] = components[comp].Kl * diffAdsorption;

    spanConcentrationDot[Ngrid, comp] =
        (spanMassFlux[Ngrid - 1, comp] - spanMassFlux[Ngrid, comp]) * idx +
        components[comp].D * (spanConcentration[Ngrid - 1, comp] - spanConcentration[Ngrid, comp]) * idx2 -
        prefactorMassTransfer[comp] * diffAdsorption;
  }
}

void computeDerivativesEnergyBalance(
    const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double externalTemperature,
    double voidFraction, double particleDensity, double particleDiameter, double influxTemperature,
    double internalDiameter, double outerDiameter, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> prefactorMassTransfer, std::span<const double> interstitialGasVelocity,
    std::span<const double> totalPressure, std::span<double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<double> solidTemperature, std::span<double> solidTemperatureDot, std::span<double> wallTemperature,
    std::span<double> wallTemperatureDot, std::span<double> concentration, std::span<double> concentrationDot,
    std::span<double> adsorption, std::span<double> adsorptionDot, std::span<double> equilibriumAdsorption,
    std::span<double> moleFraction, std::span<double> gasDensity, std::span<double> coeffGasGas,
    std::span<double> coeffGasSolid, std::span<double> coeffGasWall, std::span<double> coeffDiffusion,
    std::span<double> facePressures, std::span<double> massFlux)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), Ngrid + 1, Ncomp);
  mdspan2d_const spanAdsorption(adsorption.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanAdsorptionDot(adsorptionDot.data(), Ngrid + 1, Ncomp);
  mdspan2d_const spanConcentration(concentration.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanMoleFraction(moleFraction.data(), Ngrid + 1, Ncomp);
  mdspan2d_mut spanMassFlux(massFlux.data(), Ngrid + 1, Ncomp);

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(Ngrid + 1);
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

  for (size_t grid = 0; grid < Ngrid; ++grid)
  {
    facePressures[grid] = 0.5 * (totalPressure[grid] + totalPressure[grid + 1]) / R;
  }
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      spanMassFlux[grid, comp] = spanConcentration[grid, comp] * interstitialGasVelocity[grid];
    }
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

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    spanAdsorptionDot[0, comp] = components[comp].Kl * diffAdsorption;
    spanConcentrationDot[0, comp] = 0.0;

    solidTemperatureDot[0] += components[comp].heatOfAdsorption * spanAdsorptionDot[0, comp] / heatCapacitySolid;
  }

  // middle grid points
  for (size_t grid = 1; grid < Ngrid; ++grid)
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

    // mass balance
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      spanAdsorptionDot[grid, comp] = components[comp].Kl * diffAdsorption;

      // now add heat exchange from heat of adsorption
      solidTemperatureDot[grid] +=
          components[comp].heatOfAdsorption * spanAdsorptionDot[grid, comp] / heatCapacitySolid;

      // flow due to velocity
      double dcdt1 = -(spanMassFlux[grid, comp] - spanMassFlux[grid - 1, comp]) * idx;

      // diffusion due to pressure
      double dcdt2 = components[comp].D * idx2 * invGasTemperature *
                     (facePressures[grid] * (spanMoleFraction[grid + 1, comp] - spanMoleFraction[grid, comp]) -
                      facePressures[grid - 1] * (spanMoleFraction[grid, comp] - spanMoleFraction[grid - 1, comp]));

      // pressure change due to adsorption
      double dcdt3 = -prefactorMassTransfer[comp] * diffAdsorption;

      spanConcentrationDot[grid, comp] = dcdt1 + dcdt2 + dcdt3;
    }
  }

  // last gridpoint
  // heat flux transfer
  gasTemperatureDot[Ngrid] = coeffGasGas[Ngrid] * gasTemperature[Ngrid] +
                             coeffGasSolid[Ngrid] * solidTemperature[Ngrid] +
                             coeffGasWall[Ngrid] * wallTemperature[Ngrid];
  solidTemperatureDot[Ngrid] = coeffSolidGas * gasTemperature[Ngrid] + coeffSolidSolid * solidTemperature[Ngrid];

  wallTemperatureDot[Ngrid] = coeffWallGas * gasTemperature[Ngrid] + coeffWallWall * wallTemperature[Ngrid];

  // flux from gas diffusion and advection
  gasTemperatureDot[Ngrid] -=
      interstitialGasVelocity[Ngrid] * (gasTemperature[Ngrid] - gasTemperature[Ngrid - 1]) * idx;
  gasTemperatureDot[Ngrid] += coeffDiffusion[Ngrid] * (gasTemperature[Ngrid - 1] - gasTemperature[Ngrid]) * idx2;

  // flux from heat diffusion in wall
  wallTemperatureDot[Ngrid] +=
      idx2 * (wallThermalConductivity * invHeatDensityWall) * (wallTemperature[Ngrid - 1] - wallTemperature[Ngrid]);

  // add effect of ambient temperature to wall
  wallTemperatureDot[Ngrid] += heatTransferWallExternal * externalArea * externalTemperature * invHeatDensityWall;

  double invGasTemperature = 1.0 / std::max(1e-10, gasTemperature[Ngrid]);
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    spanAdsorptionDot[Ngrid, comp] = components[comp].Kl * diffAdsorption;

    // now add heat exchange from heat of adsorption
    solidTemperatureDot[Ngrid] +=
        components[comp].heatOfAdsorption * spanAdsorptionDot[Ngrid, comp] / heatCapacitySolid;

    // flow due to velocity
    double dcdt1 = -(spanMassFlux[Ngrid, comp] - spanMassFlux[Ngrid - 1, comp]) * idx;

    // diffusion due to pressure
    double dcdt2 = -components[comp].D * idx2 * invGasTemperature * facePressures[Ngrid - 1] *
                   (spanMoleFraction[Ngrid, comp] - spanMoleFraction[Ngrid - 1, comp]);

    // pressure change due to adsorption
    double dcdt3 = -prefactorMassTransfer[comp] * diffAdsorption;

    spanConcentrationDot[Ngrid, comp] = dcdt1 + dcdt2 + dcdt3;
  }
}

void enforceBoundaryCondition(Column& column)
{
  switch (column.boundaryCondition)
  {
    case Column::BoundaryCondition::InletPressureInletVelocity:
    {
      for (size_t j = 0; j < column.Ncomp; ++j)
      {
        column.moleFraction[0 * column.Ncomp + j] = column.components[j].Yi0;
        column.totalPressure[0] = column.inletPressure;
        column.partialPressure[0 * column.Ncomp + j] = column.components[j].Yi0 * column.totalPressure[0];
        column.concentration[0 * column.Ncomp + j] =
            column.partialPressure[0 * column.Ncomp + j] / (R * column.externalTemperature);
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
    throw std::runtime_error("Unable to call WENO with Ngrid smaller than 3.");
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
//     throw std::runtime_error("Unable to call TVD with Ngrid smaller than 3, increase number of grid points.");
//   }

//   if (clamp && input[size - 1] >= 1.0) output[size - 1] = 1.0;

//   // For right wall of 1st Node, r_value is calculated using half-cell approximation
//   r_value = (2.0 * (input[1] - input[0]) + tol) / (input[2] - input[1] + tol);
//   flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
//   output[1] += 0.5 * flux_limiter * (input[2] - input[1]);

//   // For right walls of 2nd to Ngrid-1 node
//   for (size_t i = 2; i < size - 1; ++i)
//   {
//     r_value = ((input[i] - input[i - 1]) + tol) / ((input[i + 1] - input[i]) + tol);
//     flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
//     output[i] = input[i] + 0.5 * flux_limiter * (input[i + 1] - input[i]);
//   }
// }
