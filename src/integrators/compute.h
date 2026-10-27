#pragma once
#include <cstddef>
#include <span>
#include <vector>

#include "column.h"

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
                               std::span<const double> gasTemperature, Column::FluidPhase fluidPhase,
                               double liquidDensity, Column::PHMode pHMode, double pHValue, double pKw,
                               size_t pHComponent, std::span<double> pH);

/**
 * \brief Adds reaction source terms to explicit adsorbed/pore state arrays.
 *
 * BulkConcentration reactions write their source, per unit fluid volume, into reactionBulkSource; the caller
 * subtracts it from the bulk species sink. A reaction restricted to Adsorbents is weighted by their local
 * fraction (fractionOfAdsorbent, grid-major, numberOfAdsorbents per node; empty for a single-bed column).
 */
void computeReactionDerivatives(
    const std::vector<Component>& components, const std::vector<Reaction>& reactions,
    size_t numberOfGridPoints, size_t numberOfComponents, size_t maxChemisorptionSites,
    double externalTemperature, std::span<const double> physisorption,
    std::span<double> physisorptionDot, std::span<const double> chemisorption,
    std::span<double> chemisorptionDot, std::span<const double> poreConcentration,
    std::span<double> poreConcentrationDot, std::span<const double> solidTemperature,
    std::span<double> reactionPhysisorptionSource,
    std::span<double> reactionChemisorptionSource,
    std::span<double> reactionPoreConcentrationSource, std::span<double> reactionHeat,
    std::span<const double> concentration, std::span<double> reactionBulkSource,
    std::span<const double> fractionOfAdsorbent, size_t numberOfAdsorbents);

/**
 * \brief Computes concentration derivatives from explicit arrays.
 *
 * With WENO3 the convective face values are reconstructed with computeWENO: in the gas phase the mole fractions,
 * with the total molar flux v c_tot kept upwind (the flux that the velocity closures conserve), in the liquid phase
 * the concentrations.
 */
void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, double resolution,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> bulkSpeciesSink, Column::ConvectionScheme convectionScheme,
                            Column::FluidPhase fluidPhase);

/**
 * \brief Computes gas, solid, and wall temperature derivatives from explicit arrays.
 */
void computeEnergyDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, const Geometry& geometry,
    double particleDensity, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> interstitialGasVelocity, std::span<const double> gasDensity,
    std::span<double> coeffDiffusion, std::span<const double> physisorptionDot,
    std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<const double> reactionPhysisorptionSource = {},
    std::span<const double> reactionChemisorptionSource = {},
    std::span<const double> reactionHeat = {});

/**
 * \brief Third-order WENO face values on a uniform grid (v > 0).
 *
 * values[0] is the inlet value at the inlet face z_0, and values[i], i >= 1, the value of the cell (z_{i-1}, z_i].
 * faces[i] is the value at the face z_i reconstructed from the upstream side: faces[0] = values[0], and the outlet
 * face takes the last cell value (first order).
 */
void computeWENO(std::span<const double> values, std::span<double> faces);

/**
 * \brief Third-order WENO face values on the grid with node positions columnDistances (z_0 ... z_N).
 *
 * Same conventions as the uniform overload, to which it reduces for equal spacings.
 */
void computeWENO(std::span<const double> values, std::span<const double> columnDistances, std::span<double> faces);
