#pragma once
#include <span>
#include <vector>

#include "column.h"

void computeSorptionDerivatives(Column& column);

void computeSorptionDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                                size_t numberOfComponents, size_t maxChemisorptionSites,
                                double externalTemperature, double voidFraction, double particleDensity,
                                double particleDiameter,
                                std::span<const double> equilibriumAdsorption,
                                std::span<const double> concentration,
                                std::span<const double> physisorption, std::span<double> physisorptionDot,
                                std::span<const double> chemisorption, std::span<double> chemisorptionDot,
                                std::span<const double> surfaceConcentration,
                                std::span<double> surfaceConcentrationDot,
                                std::span<const double> poreConcentration, std::span<double> poreConcentrationDot,
                                std::span<const double> solidTemperature, std::span<double> bulkSpeciesSink);

void computePhysisorption(Column& column);

void computePhysisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, std::span<const double> equilibriumAdsorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot);

void computeChemisorption(Column& column);

void computeChemisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t maxChemisorptionSites,
                          double externalTemperature, double voidFraction, double particleDensity,
                          std::span<const double> equilibriumAdsorption,
                          std::span<const double> concentration, std::span<const double> chemisorption,
                          std::span<double> chemisorptionDot, std::span<const double> poreConcentration,
                          std::span<const double> solidTemperature);
