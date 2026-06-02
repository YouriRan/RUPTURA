#pragma once
#include <span>
#include <vector>

#include "column.h"
#include "mixture_prediction.h"

void computePressure(Column& column);
void computePressure(Column::VelocityProfile& velocityProfile, Column::BoundaryCondition& boundaryCondition,
                     const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double inletPressure,
                     double outletPressure, double voidFraction, double dynamicViscosity, double particleDiameter,
                     double resolution, std::span<const double> interstitialGasVelocity, std::span<double> gasDensity,
                     std::span<double> totalConcentration, std::span<double> totalPressure,
                     std::span<double> gasTemperature, std::span<const double> concentration,
                     std::span<double> partialPressure, std::span<double> moleFraction);

void computeEquilibriumLoadings(Column& column);
void computeEquilibriumLoadings(MixturePrediction& mixture, size_t Ngrid, size_t Ncomp, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance, std::span<double> totalPressure,
                                std::span<double> gasTemperature, std::span<double> idealGasMolFractions,
                                std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules,
                                std::span<double> equilibriumAdsorption, std::span<double> moleFraction,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential);

void computeVelocity(Column& column);
void computeVelocityFixedGradient(Column::BoundaryCondition& boundaryCondition,
                                  const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                                  double pressureGradient, double columnEntranceVelocity, double resolution,
                                  std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                                  std::span<const double> totalConcentration, std::span<const double> adsorption,
                                  std::span<const double> equilibriumAdsorption, std::span<const double> concentration);
void computeVelocityErgun(Column::BoundaryCondition& boundaryCondition, size_t Ngrid, double voidFraction,
                          double columnEntranceVelocity, double columnLength, double dynamicViscosity,
                          double particleDiameter, double resolution, std::span<double> interstitialGasVelocity,
                          std::span<double> gasDensity, std::span<double> totalPressure);

void computeFirstDerivatives(Column& column);
void computeFirstDerivatives(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double resolution,
                             std::span<const double> prefactorMassTransfer,
                             std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                             std::span<double> concentrationDot, std::span<const double> adsorption,
                             std::span<double> adsorptionDot, std::span<const double> equilibriumAdsorption,
                             std::span<double> massFlux);

void computeFirstDerivativesEnergyBalance(
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
    std::span<double> facePressures, std::span<double> massFlux);

void computeFirstDerivativesWENO(Column& column);
void enforceBoundaryCondition(Column& column);

void computeWENO(std::span<const double> input, std::span<double> output);
