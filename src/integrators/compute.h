#pragma once
#include <span>
#include <vector>

#include "column.h"
#include "mixture_prediction.h"

/**
 * \brief Updates pressure, concentration, partial-pressure, and density fields for a column.
 */
void computePressure(Column& column);

/**
 * \brief Computes pressure-related fields from explicit arrays and model settings.
 */
void computePressure(const std::vector<Component>& components, const Column::VelocityProfile& velocityProfile,
                     const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                     size_t numberOfComponents, double inletPressure, double outletPressure, double voidFraction,
                     double dynamicViscosity, double particleDiameter, double resolution,
                     std::span<const double> interstitialGasVelocity, std::span<double> gasDensity,
                     std::span<double> totalConcentration, std::span<double> totalPressure,
                     std::span<double> concentration, std::span<double> partialPressure,
                     std::span<const double> moleFraction, std::span<const double> gasTemperature);

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
                                std::span<double> equilibriumAdsorption, std::span<double> cachedPressure,
                                std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
                                std::span<double> gasTemperature);

/**
 * \brief Updates the interstitial gas velocity field for a column.
 */
void computeVelocity(Column& column);

/**
 * \brief Computes velocity using the fixed-pressure-gradient formulation.
 */
void computeVelocityFixedGradient(const Column::BoundaryCondition& boundaryCondition,
                                  const std::vector<Component>& components, size_t numberOfGridPoints,
                                  size_t numberOfComponents, double pressureGradient, double columnEntranceVelocity,
                                  double resolution, std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalConcentration,
                                  std::span<const double> totalPressure, std::span<const double> equilibriumAdsorption,
                                  std::span<const double> moleFraction, std::span<const double> adsorption);

/**
 * \brief Computes velocity from the Ergun pressure-drop relation.
 */
void computeVelocityErgun(const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                          double voidFraction, double columnEntranceVelocity, double columnLength,
                          double dynamicViscosity, double particleDiameter, double resolution,
                          std::span<double> interstitialGasVelocity, std::span<const double> gasDensity,
                          std::span<const double> totalPressure);

/**
 * \brief Updates mole-fraction and adsorption derivatives for a column.
 */
void computeDerivatives(Column& column);

/**
 * \brief Computes isothermal mole-fraction and adsorption derivatives from explicit arrays.
 */
void computeDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
                        double resolution, std::span<const double> prefactorMassTransfer,
                        std::span<const double> interstitialGasVelocity, std::span<const double> totalConcentration,
                        std::span<const double> equilibriumAdsorption, std::span<const double> moleFraction,
                        std::span<double> moleFractionDot, std::span<const double> adsorption,
                        std::span<double> adsorptionDot);

/**
 * \brief Computes concentration, adsorption, and temperature derivatives with energy balance enabled.
 */
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
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot);

/**
 * \brief Updates derivatives using the WENO advection reconstruction.
 */
void computeDerivativesWENO(Column& column);

/**
 * \brief Applies the selected inlet/outlet boundary condition to the column state.
 */
void enforceBoundaryCondition(Column& column);

/**
 * \brief Reconstructs a one-dimensional signal with the WENO stencil.
 */
void computeWENO(std::span<const double> input, std::span<double> output);
