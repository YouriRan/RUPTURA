#include "compute.h"

#include <mdspan>

#include "utils.h"

void computeEquilibriumLoadings(Column& column)
{
  computeEquilibriumLoadings(column.mixture, column.Ngrid, column.Ncomp, column.maxIsothermTerms,
                             column.iastPerformance, column.pressureGradient, column.columnLength, column.totalPressure,
                             column.partialPressure, column.idealGasMolFractions, column.adsorbedMolFractions,
                             column.numberOfMolecules, column.equilibriumAdsorption, column.cachedPressure,
                             column.cachedGrandPotential);
}
void computeEquilibriumLoadings(MixturePrediction& mixture, size_t Ngrid, size_t Ncomp, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance, double pressureGradient,
                                double columnLength, std::span<double> totalPressure, std::span<double> partialPressure,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<double> equilibriumAdsorption,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential)
{
  // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    totalPressure[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      totalPressure[grid] += std::max(0.0, partialPressure[grid * Ncomp + comp]);
    }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      idealGasMolFractions[comp] = std::max(partialPressure[grid * Ncomp + comp], 0.0);
      sum += idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      if (sum < 1e-10)
      {
        idealGasMolFractions[comp] = 1.0 / static_cast<double>(Ncomp);
      }
      else
      {
        idealGasMolFractions[comp] /= sum;
      }
    }

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    iastPerformance += mixture.predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions,
                                              numberOfMolecules, &cachedPressure[grid * Ncomp * maxIsothermTerms],
                                              &cachedGrandPotential[grid * maxIsothermTerms]);

    for (size_t comp = 0; comp < Ncomp; ++comp)

    {
      equilibriumAdsorption[grid * Ncomp + comp] = numberOfMolecules[comp];
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (totalPressure[0] + pressureGradient * columnLength < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeVelocity(Column& column)
{
  switch (column.velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column.components, column.Ngrid, column.Ncomp, column.pressureGradient,
                                   column.columnEntranceVelocity, column.resolution, column.prefactorMassTransfer,
                                   column.interstitialGasVelocity, column.totalPressure, column.adsorption,
                                   column.equilibriumAdsorption, column.partialPressure);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column.components, column.Ngrid, column.Ncomp, column.externalTemperature,
                           column.voidFraction, column.columnEntranceVelocity, column.dynamicViscosity,
                           column.particleDiameter, column.resolution, column.interstitialGasVelocity,
                           column.totalPressure, column.partialPressure);
      break;
    }
    case Column::VelocityProfile::FixedVelocity:
      break;
    default:
      break;
  }
}

void computeVelocityFixedGradient(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                                  double pressureGradient, double columnEntranceVelocity, double resolution,
                                  std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                                  std::span<const double> adsorption, std::span<const double> equilibriumAdsorption,
                                  std::span<const double> partialPressure)
{
  double idx2 = 1.0 / (resolution * resolution);

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

      // diffusion term
      sum += components[comp].D *
             (partialPressure[(grid - 1) * Ncomp + comp] - 2.0 * partialPressure[grid * Ncomp + comp] +
              partialPressure[(grid + 1) * Ncomp + comp]) *
             idx2;
    }

    // explicit update
    interstitialGasVelocity[grid] =
        interstitialGasVelocity[grid - 1] +
        resolution * (sum - interstitialGasVelocity[grid - 1] * pressureGradient) / totalPressure[grid];
  }

  // last grid point
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      sum -= prefactorMassTransfer[comp] *
             (equilibriumAdsorption[Ngrid * Ncomp + comp] - adsorption[Ngrid * Ncomp + comp]);

      sum += components[comp].D *
             (partialPressure[(Ngrid - 1) * Ncomp + comp] - partialPressure[Ngrid * Ncomp + comp]) * idx2;
    }

    interstitialGasVelocity[Ngrid] =
        interstitialGasVelocity[Ngrid - 1] +
        resolution * (sum - interstitialGasVelocity[Ngrid - 1] * pressureGradient) / totalPressure[Ngrid];
  }
}

void computeVelocityErgun(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                          double externalTemperature, double voidFraction, double columnEntranceVelocity,
                          double dynamicViscosity, double particleDiameter, double resolution,
                          std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                          std::span<const double> partialPressure)
{
  double idx = 1.0 / (resolution);

  // 1.75 * (L * rho / d) * ((1 - epsilon) / epsilon**3) v^2 +
  // 150 (mu * L / d^2) ((1 - epsilon)^2 / epsilon^3) - (dp/dz)
  // C_quad = d(P_Total) / dz;
  // K_viscous = (150 mu L / d^2) ((1 - epsilon)^2 / epsilon^3)
  // K_inertial = (1.75 L rho / d) ((1 - epsilon) / epsilon^3)
  // B_quad = K_viscous * mu
  // A_quad = K_inertial * rho_g
  // A v^2 + b v + c = 0

  // v = -B + sqrt(B^2 - 4 A C) / (2 * A)

  // mu = dynamic viscosity
  // L = column length
  // d = particle diameter?
  // epsilon = bed porosity
  // rho_f = fluid density
  // v = velocity

  double voidFractionPrefactorA = (1 - voidFraction) / voidFraction;
  double voidFractionPrefactorB = (1 - voidFraction) * (1 - voidFraction) / (voidFraction * voidFraction);

  double prefactorA = 1.75 * voidFractionPrefactorA / particleDiameter;
  double B = 150.0 * voidFractionPrefactorB * dynamicViscosity / (particleDiameter * particleDiameter);

  // first grid point
  interstitialGasVelocity[0] = columnEntranceVelocity;

  for (size_t grid = 1; grid < Ngrid + 1; grid++)
  {
    double C = -(totalPressure[grid] - totalPressure[grid - 1]) * idx;

    double gasDensity = 0.0;
    for (size_t comp = 0; comp < Ncomp; comp++)
    {
      gasDensity += components[comp].molecularWeight * partialPressure[grid * Ncomp + comp] / (R * externalTemperature);
    }
    double A = prefactorA * gasDensity;

    double v = 0.0;
    if (std::abs(A) < 1e-10)
    {
      v = -C / B;
    }
    else
    {
      double ac4 = std::max(0.0, B * B - 4.0 * A * C);
      v = (-B + std::sqrt(ac4)) / (2.0 * A);
    }
    interstitialGasVelocity[grid] = std::max(0.0, v);
  }
}

void computeFirstDerivatives(Column& column)
{
  if (column.energyBalance)
  {
    computeFirstDerivativesEnergyBalance(
        column.components, column.Ngrid, column.Ncomp, column.externalTemperature, column.voidFraction,
        column.particleDensity, column.particleDiameter, column.internalDiameter, column.outerDiameter,
        column.wallDensity, column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
        column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
        column.heatCapacityWall, column.resolution, column.prefactorMassTransfer, column.interstitialGasVelocity,
        column.totalPressure, column.gasTemperature, column.gasTemperatureDot, column.solidTemperature,
        column.solidTemperatureDot, column.wallTemperature, column.wallTemperatureDot, column.partialPressure,
        column.partialPressureDot, column.adsorption, column.adsorptionDot, column.equilibriumAdsorption,
        column.moleFraction, column.gasDensity, column.coeffGasGas, column.coeffGasSolid, column.coeffGasWall,
        column.coeffDiffusion, column.facePressures, column.energyFlow);
  }
  else
  {
    computeFirstDerivatives(column.components, column.Ngrid, column.Ncomp, column.resolution,
                            column.prefactorMassTransfer, column.interstitialGasVelocity, column.partialPressure,
                            column.partialPressureDot, column.adsorption, column.adsorptionDot,
                            column.equilibriumAdsorption, column.totalPressureDot);
  }
}

void computeFirstDerivatives(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double resolution,
                             std::span<const double> prefactorMassTransfer,
                             std::span<const double> interstitialGasVelocity, std::span<const double> partialPressure,
                             std::span<double> partialPressureDot, std::span<const double> adsorption,
                             std::span<double> adsorptionDot, std::span<const double> equilibriumAdsorption,
                             std::span<double> totalPressureDot)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  std::mdspan<const double, std::dextents<size_t, 2>> spanPartialPressure(partialPressure.data(), Ngrid + 1, Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanEquilibriumAdsorption(equilibriumAdsorption.data(), Ngrid + 1,
                                                                                Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanAdsorption(adsorption.data(), Ngrid + 1, Ncomp);

  // first gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    adsorptionDot[comp] = components[comp].Kl * diffAdsorption;
    partialPressureDot[comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    totalPressureDot[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      adsorptionDot[grid * Ncomp + comp] = components[comp].Kl * diffAdsorption;

      partialPressureDot[grid * Ncomp + comp] =
          (interstitialGasVelocity[grid - 1] * spanPartialPressure[grid - 1, comp] -
           interstitialGasVelocity[grid] * spanPartialPressure[grid, comp]) *
              idx +
          components[comp].D *
              (spanPartialPressure[grid + 1, comp] - 2.0 * spanPartialPressure[grid, comp] +
               spanPartialPressure[grid - 1, comp]) *
              idx2 -
          prefactorMassTransfer[comp] * diffAdsorption;
      totalPressureDot[grid] += partialPressureDot[grid * Ncomp + comp];
    }
  }

  // last gridpoint
  totalPressureDot[Ngrid] = 0.0;
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    adsorptionDot[Ngrid * Ncomp + comp] = components[comp].Kl * diffAdsorption;

    partialPressureDot[Ngrid * Ncomp + comp] =
        (interstitialGasVelocity[Ngrid - 1] * spanPartialPressure[Ngrid - 1, comp] -
         interstitialGasVelocity[Ngrid] * spanPartialPressure[Ngrid, comp]) *
            idx +
        components[comp].D * (spanPartialPressure[Ngrid - 1, comp] - spanPartialPressure[Ngrid, comp]) * idx2 -
        prefactorMassTransfer[comp] * diffAdsorption;
    totalPressureDot[Ngrid] += partialPressureDot[Ngrid * Ncomp + comp];
  }
}

void computeFirstDerivativesEnergyBalance(
    const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double externalTemperature,
    double voidFraction, double particleDensity, double particleDiameter, double internalDiameter, double outerDiameter,
    double wallDensity, double gasThermalConductivity, double wallThermalConductivity, double heatTransferGasSolid,
    double heatTransferGasWall, double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
    double heatCapacityWall, double resolution, std::span<const double> prefactorMassTransfer,
    std::span<const double> interstitialGasVelocity, std::span<const double> totalPressure,
    std::span<double> gasTemperature, std::span<double> gasTemperatureDot, std::span<double> solidTemperature,
    std::span<double> solidTemperatureDot, std::span<double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<double> partialPressure, std::span<double> partialPressureDot, std::span<double> adsorption,
    std::span<double> adsorptionDot, std::span<double> equilibriumAdsorption, std::span<double> moleFraction,
    std::span<double> gasDensity, std::span<double> coeffGasGas, std::span<double> coeffGasSolid,
    std::span<double> coeffGasWall, std::span<double> coeffDiffusion, std::span<double> facePressures,
    std::span<double> energyFlow)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  std::mdspan<double, std::dextents<size_t, 2>> spanPartialPressure(partialPressure.data(), Ngrid + 1, Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanEquilibriumAdsorption(equilibriumAdsorption.data(), Ngrid + 1,
                                                                          Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanAdsorption(adsorption.data(), Ngrid + 1, Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanPartialPressureDot(partialPressureDot.data(), Ngrid + 1, Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanAdsorptionDot(adsorptionDot.data(), Ngrid + 1, Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanMoleFraction(moleFraction.data(), Ngrid + 1, Ncomp);
  std::mdspan<double, std::dextents<size_t, 2>> spanEnergyFlow(energyFlow.data(), Ngrid + 1, Ncomp);

  double relativeVolume = ((1 - voidFraction) / voidFraction);
  double accessibleSurface = 6.0 / particleDiameter;  // prefactor 6.0 in python and 2.0 in eqs

  double prefactorGasSolid = relativeVolume * accessibleSurface * heatTransferGasSolid / heatCapacityGas;

  // is this extra 1/eps necessary? It's in python not in eqs
  double prefactorGasWall = 4.0 * heatTransferGasWall / (voidFraction * heatCapacityGas * internalDiameter);
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
    facePressures[grid] = 0.5 * (totalPressure[grid] + totalPressure[grid + 1]);
  }

  double invGasConstant = 1.0 / R;
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    double invTotalPressure = 1.0 / std::max(1e-10, totalPressure[grid]);
    double invGasTemperature = 1.0 / std::max(1e-10, gasTemperature[grid]);
    gasDensity[grid] = 0.0;

    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      spanMoleFraction[grid, comp] = spanPartialPressure[grid, comp] * invTotalPressure;
      gasDensity[grid] += components[comp].molecularWeight * invGasTemperature * invGasConstant *
                          std::max(0.0, spanPartialPressure[grid, comp]);
      spanEnergyFlow[grid, comp] = spanPartialPressure[grid, comp] * interstitialGasVelocity[grid] * invGasTemperature;
    }
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffGasGas[grid] = prefactorGasGas * invGasDensity;
    coeffGasSolid[grid] = prefactorGasSolid * invGasDensity;
    coeffGasWall[grid] = prefactorGasWall * invGasDensity;
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  // first grid point
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    spanAdsorptionDot[0, comp] = components[comp].Kl * diffAdsorption;
    spanPartialPressureDot[0, comp] = 0.0;
  }

  // boundary conditions
  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = 0.0;
  wallTemperatureDot[0] = 0.0;

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
    wallTemperatureDot[grid] += (wallThermalConductivity * invHeatDensityWall) *
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
      // python defines as H and eqs as -H, what do?
      solidTemperatureDot[grid] -=
          components[comp].heatOfAdsorption * spanAdsorptionDot[grid, comp] / heatCapacitySolid;

      // flow due to velocity
      double dpdt1 = -gasTemperature[grid] * ((spanEnergyFlow[grid, comp] - spanEnergyFlow[grid - 1, comp])) * idx;

      // pressure change due to heating
      double dpdt2 = spanPartialPressure[grid, comp] * gasTemperatureDot[grid] * invGasTemperature;

      // diffusion due to pressure
      double dpdt3 = components[comp].D * idx2 *
                     (facePressures[grid] * (spanMoleFraction[grid + 1, comp] - spanMoleFraction[grid, comp]) -
                      facePressures[grid - 1] * (spanMoleFraction[grid, comp] - spanMoleFraction[grid - 1, comp]));

      // pressure correction due to advection cause by heating
      double dpdt4 = -components[comp].D * totalPressure[grid] * invGasTemperature * idx2 *
                     (gasTemperature[grid - 1] - gasTemperature[grid]) *
                     (spanMoleFraction[grid - 1, comp] - spanMoleFraction[grid, comp]);

      // pressure change due to adsorption
      double dpdt5 = -prefactorMassTransfer[comp] * diffAdsorption;

      spanPartialPressureDot[grid, comp] = dpdt1 + dpdt2 + dpdt3 + dpdt4 + dpdt5;
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
  gasTemperatureDot[Ngrid] =
      interstitialGasVelocity[Ngrid] * (gasTemperature[Ngrid] - gasTemperature[Ngrid - 1]) * idx +
      coeffDiffusion[Ngrid] * (gasTemperature[Ngrid - 1] - gasTemperature[Ngrid]) * idx2;

  // flux from heat diffusion in wall
  wallTemperatureDot[Ngrid] +=
      (wallThermalConductivity * invHeatDensityWall) * (wallTemperature[Ngrid - 1] - wallTemperature[Ngrid]);

  // add effect of ambient temperature to wall
  wallTemperatureDot[Ngrid] += heatTransferWallExternal * externalArea * externalTemperature * invHeatDensityWall;

  double invGasTemperature = std::max(1e-10, 1.0 / gasTemperature[Ngrid]);
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    spanAdsorptionDot[Ngrid, comp] = components[comp].Kl * diffAdsorption;

    // flow due to velocity
    double dpdt1 = -gasTemperature[Ngrid] * (spanEnergyFlow[Ngrid, comp] - spanEnergyFlow[Ngrid - 1, comp]) * idx;

    // pressure change due to heating
    double dpdt2 = spanPartialPressure[Ngrid, comp] * gasTemperatureDot[Ngrid] * invGasTemperature;

    // diffusion due to pressure
    double dpdt3 = -components[comp].D * idx2 * facePressures[Ngrid - 1] *
                   (spanMoleFraction[Ngrid, comp] - spanMoleFraction[Ngrid - 1, comp]);

    // pressure correction due to advection cause by heating
    double dpdt4 = -components[comp].D * invGasTemperature * totalPressure[Ngrid] * idx2 *
                   (gasTemperature[Ngrid] - gasTemperature[Ngrid - 1]) *
                   (spanMoleFraction[Ngrid, comp] - spanMoleFraction[Ngrid - 1, comp]);

    // pressure change due to adsorption
    double dpdt5 = -prefactorMassTransfer[comp] * diffAdsorption;

    spanPartialPressureDot[Ngrid, comp] = dpdt1 + dpdt2 + dpdt3 + dpdt4 + dpdt5;
  }
}

std::vector<double> computeWENO(std::vector<double>& vec, bool clamp)
{
  if (vec.empty()) return vec;

  double tol = 1e-10;
  double alpha_0, alpha_1, first_term, second_term;

  size_t size = vec.size();
  std::vector<double> out(vec);

  if (clamp && vec[size - 1] >= 1.0) out[size - 1] = 1.0;

  // For right wall of 1st Node, alpha_1 and second term values are calculated using half-cell approximation
  alpha_0 = (2.0 / 3.0) / std::pow((vec[2] - vec[1]) + tol, 4.0);
  alpha_1 = (1.0 / 3.0) / std::pow(2.0 * (vec[1] - vec[0]) + tol, 4.0);

  first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (vec[2] + vec[1]);
  second_term = (alpha_1 / (alpha_0 + alpha_1)) * (2.0 * vec[1] - vec[0]);
  out[1] = first_term + second_term;

  // For right walls of 2nd to Ngrid-1 node
  for (size_t i = 2; i < size - 1; ++i)
  {
    alpha_0 = (2.0 / 3.0) / std::pow(vec[i + 1] - vec[i] + tol, 4.0);
    alpha_1 = (1.0 / 3.0) / std::pow(vec[i] - vec[i - 1] + tol, 4.0);

    first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (vec[i + 1] + vec[i]);
    second_term = (alpha_1 / (alpha_0 + alpha_1)) * ((3.0 / 2.0) * vec[i] - (1.0 / 2.0) * vec[i - 1]);
    out[i] = first_term + second_term;
  }
  return out;
}

std::vector<double> computeTVD(std::vector<double>& vec, bool clamp)
{
  if (vec.empty()) return vec;

  double tol = 1e-10;
  double r_value, flux_limiter;

  size_t size = vec.size();
  std::vector<double> out(vec);

  if (clamp && vec[size - 1] >= 1.0) out[size - 1] = 1.0;

  // For right wall of 1st Node, r_value is calculated using half-cell approximation
  r_value = (2.0 * (vec[1] - vec[0]) + tol) / (vec[2] - vec[1] + tol);
  flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
  out[1] += 0.5 * flux_limiter * (vec[2] - vec[1]);

  // For right walls of 2nd to Ngrid-1 node
  for (size_t i = 2; i < size - 1; ++i)
  {
    r_value = ((vec[i] - vec[i - 1]) + tol) / ((vec[i + 1] - vec[i]) + tol);
    flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
    out[i] = vec[i] + 0.5 * flux_limiter * (vec[i + 1] - vec[i]);
  }
  return out;
}
