#pragma once
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "geometry.h"
#include "mixture_prediction.h"

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
 * \brief Computes linear-driving-force physisorption derivatives.
 */
void computePhysisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents,
                          std::span<const double> fractionOfAdsorbent, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot);

/**
 * \brief Largest chemisorption site count of each component over the adsorbents.
 */
std::vector<size_t> chemisorptionSitesPerComponent(const std::vector<MixturePrediction>& mixtures,
                                                   size_t numberOfComponents);

/**
 * \brief Computes bed-weighted chemisorption kinetic derivatives.
 */
void computeChemisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                          double externalTemperature, double elapsedTime, std::span<const double> fractionOfAdsorbent,
                          std::span<const double> particleDensities,
                          std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
                          std::span<const double> chemisorption, std::span<double> chemisorptionDot,
                          std::span<const double> poreConcentration, std::span<const double> solidTemperature);
