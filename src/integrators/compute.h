#pragma once
#include <span>
#include <vector>

#include "column.h"
#include "mixture_prediction.h"

/**
 * \brief Updates velocity, pressure, concentration, partial-pressure, and density fields for a column.
 */
void updateVelocityAndPressure(Column& column);

/**
 * \brief Updates velocity and pressure-related fields from explicit arrays and boundary settings.
 */
void updateVelocityAndPressure(const std::vector<Component>& components,
                               const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, const Geometry& geometry,
                               double& columnEntranceVelocity, double dynamicViscosity, double resolution,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> bulkSpeciesSink,
                               std::span<const double> gasTemperature);

void updateVelocityAndPressure(const std::vector<Component>& components,
                               const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, double voidFraction,
                               double particleDensity, double& columnEntranceVelocity, double dynamicViscosity,
                               double particleDiameter, double resolution,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> bulkSpeciesSink,
                               std::span<const double> gasTemperature);

/**
 * \brief Updates equilibrium adsorbed loadings for every grid node in a column.
 */
void computeEquilibriumLoadings(Column& column);

/**
 * \brief Computes equilibrium adsorbed loadings from explicit arrays and cache storage.
 */
void computeEquilibriumLoadings(MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
                                size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<const double> totalPressure,
                                std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
                                std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
                                std::span<double> gasTemperature);

/**
 * \brief Computes site-major chemisorption equilibrium loadings with one multisite mixture model.
 */
void computeChemisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, std::pair<size_t, size_t>& iastPerformance,
    std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
    std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
    std::span<double> gasTemperature);

/**
 * \brief Updates mass derivatives, and temperature derivatives when energy balance is enabled.
 */
void computeDerivatives(Column& column);

/**
 * \brief Computes isothermal concentration derivatives from explicit arrays.
 */
void computeMassDerivatives(Column& column);

/**
 * \brief Computes concentration derivatives from explicit arrays.
 */
void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double resolution, const ShapeParameters& geometry, double particleDensity,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot);

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double resolution, double voidFraction,
                            double particleDensity, std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot);

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, double resolution,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> bulkSpeciesSink);

/**
 * \brief Computes gas, solid, and wall temperature derivatives from explicit arrays.
 */
void computeEnergyDerivatives(Column& column);

/**
 * \brief Computes gas, solid, and wall temperature derivatives from explicit arrays.
 */
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
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot);

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
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot);

/**
 * \brief Updates derivatives using the WENO advection reconstruction.
 */
void computeDerivativesWENO(Column& column);

/**
 * \brief Reconstructs a one-dimensional signal with the WENO stencil.
 */
void computeWENO(std::span<const double> input, std::span<double> output);
