#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "component.h"
#include "geometry.h"

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, const Geometry& geometry, double particleDensity,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double voidFraction, double particleDensity, double particleDiameter,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            const Geometry& geometry, double particleDensity,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource = {},
                            std::span<const double> reactionChemisorptionSource = {});

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double voidFraction, double particleDensity, double particleDiameter,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink);
