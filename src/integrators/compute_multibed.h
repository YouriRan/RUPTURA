#pragma once
#include <span>
#include <vector>

#include "column_multibed.h"
#include "mixture_prediction.h"

void updateVelocityAndPressure(ColumnMultibed& column);
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
                               std::span<const double> gasTemperature);
void computeEquilibriumLoadings(ColumnMultibed& column);
void computeEquilibriumLoadings(std::vector<MixturePrediction>& mixture, size_t numberOfGridPoints,
                                size_t numberOfComponents, size_t numberOfAdsorbents,
                                std::span<const double> fractionOfAdsorbent,
                                const std::vector<bool>& hasAdsorbentOfType, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<const double> totalPressure,
                                std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
                                std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
                                std::span<double> gasTemperature);

void computeDerivatives(ColumnMultibed& column);
void computeMassDerivatives(ColumnMultibed& column);
void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents,
                            std::span<const double> columnDistances,
                            std::span<const double> totalVoidFraction,
                            std::span<const double> particleDensity,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot);

void computeEnergyDerivatives(ColumnMultibed& column);
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
    std::span<double> wallTemperatureDot);

void computeDerivativesWENO(ColumnMultibed& column);

void computeWENO(std::span<const double> input, std::span<double> output);
