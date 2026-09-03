#pragma once
#include <span>
#include <vector>

#include "column_multibed.h"
#include "mixture_prediction.h"

/**
 * \brief Updates multibed velocity and pressure fields from explicit arrays.
 */
void updateVelocityAndPressure(
    const std::vector<Component>& components, const MultibedColumn::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const double> adsorbentScaledVoidFraction,
    std::span<const double> totalVoidFraction,
    std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature,
    MultibedColumn::FluidPhase fluidPhase, double liquidDensity, MultibedColumn::PHMode pHMode,
    double pHValue, double pKw, size_t pHComponent, std::span<double> pH);

/**
 * \brief Computes per-bed physisorption equilibrium and blends interface loadings.
 */
void computePhysisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH);

/**
 * \brief Computes per-bed chemisorption equilibrium and blends interface loadings.
 */
void computeChemisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& chemisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxChemisorptionSites, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH);

/**
 * \brief Computes multibed linear-driving-force physisorption derivatives.
 */
void computePhysisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents,
                          std::span<const double> fractionOfAdsorbent, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot);

/**
 * \brief Computes bed-weighted chemisorption kinetic derivatives.
 */
void computeChemisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                          double externalTemperature, double elapsedTime, std::span<const double> fractionOfAdsorbent,
                          std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                          std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
                          std::span<const double> chemisorption, std::span<double> chemisorptionDot,
                          std::span<const double> poreConcentration, std::span<const double> solidTemperature);

/**
 * \brief Computes bed-weighted film/pore transport derivatives for chemisorption.
 */
void computeChemisorptionTransportDerivatives(
    const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, size_t maxChemisorptionSites, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries,
    std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
    std::span<const double> particleDiameters, std::span<const double> totalVoidFraction,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

/**
 * \brief Builds the gas-phase sink from bed-weighted solid uptake.
 */
void computeBulkSpeciesSink(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                            std::span<const Geometry> geometries,
                            std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                            std::span<const double> particleDiameters, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> totalVoidFraction, std::span<const double> concentration,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration, std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource = {},
                            std::span<const double> reactionChemisorptionSource = {});

/**
 * \brief Computes multibed concentration derivatives from explicit arrays.
 */
void computeMassDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink,
                            MultibedColumn::FluidPhase fluidPhase, std::span<const double> totalVoidFraction);

/**
 * \brief Computes multibed gas, solid, and wall temperature derivatives.
 */
void computeEnergyDerivatives(
    const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, double externalTemperature, std::span<const double> totalVoidFraction,
    std::span<const Geometry> geometries,
    std::span<const double> particleDensities, std::span<const double> particleDiameters,
    std::span<const double> fractionOfAdsorbent, double internalDiameter, double outerDiameter, double wallDensity,
    double gasThermalConductivity, double wallThermalConductivity, double heatTransferGasSolid,
    double heatTransferGasWall, double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
    double heatCapacityWall, std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
    std::span<const double> gasDensity, std::span<double> coeffDiffusion, size_t maxChemisorptionSites,
    std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<const double> reactionPhysisorptionSource = {}, std::span<const double> reactionChemisorptionSource = {},
    std::span<const double> reactionHeat = {});
