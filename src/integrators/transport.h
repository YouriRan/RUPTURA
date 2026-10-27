#pragma once
#include <cstddef>
#include <span>
#include <vector>

#include "geometry.h"
#include "mixture_prediction.h"

/**
 * \brief Computes bed-weighted film/pore transport derivatives for chemisorption.
 */
void computeChemisorptionTransportDerivatives(
    const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, size_t maxChemisorptionSites, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const double> particleDensities,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

/**
 * \brief Builds the gas-phase sink from bed-weighted solid uptake.
 */
void computeBulkSpeciesSink(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                            std::span<const Geometry> geometries, std::span<const double> particleDensities,
                            std::span<const double> fractionOfAdsorbent, std::span<const double> totalVoidFraction,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration, std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource = {},
                            std::span<const double> reactionChemisorptionSource = {},
                            std::span<const double> reactionBulkSource = {});
