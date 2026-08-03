#pragma once

#include <span>
#include <vector>

#include "column.h"

void computeChemisorptionTransportDerivatives(Column& column);

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, const ShapeParameters& geometry, double particleDensity,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double voidFraction, double particleDensity, double particleDiameter,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot);

void computeBulkSpeciesSink(Column& column);

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            const ShapeParameters& geometry, double particleDensity,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink);

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double voidFraction, double particleDensity, double particleDiameter,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink);
