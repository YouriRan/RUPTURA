#pragma once
#include <span>
#include <vector>

#include "column_multibed.h"
#include "mixture_prediction.h"

/**
 * \brief Updates multibed velocity and pressure fields from explicit arrays.
 */
void updateVelocityAndPressure(
    const std::vector<Component>& components, const ColumnMultibed::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const double> adsorbentScaledVoidFraction, std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature);

/**
 * \brief Computes per-bed physisorption equilibrium and blends interface loadings.
 */
void computePhysisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature);

/**
 * \brief Computes multibed linear-driving-force physisorption derivatives.
 */
void computePhysisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents,
                          std::span<const double> fractionOfAdsorbent, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot);

/**
 * \brief Builds the gas-phase sink from bed-weighted solid uptake.
 */
void computeBulkSpeciesSink(size_t numberOfGridPoints, size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                            std::span<const double> fractionOfAdsorbent, std::span<const double> totalVoidFraction,
                            std::span<const double> physisorptionDot, std::span<double> bulkSpeciesSink);

/**
 * \brief Computes multibed concentration derivatives from explicit arrays.
 */
void computeMassDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink);

/**
 * \brief Computes multibed gas, solid, and wall temperature derivatives.
 */
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
                              std::span<double> wallTemperatureDot);
