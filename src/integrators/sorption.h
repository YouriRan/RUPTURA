#pragma once
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "component.h"
#include "geometry.h"
#include "mixture_prediction.h"

void computePhysisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents, size_t maxIsothermTerms,
    std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH);

void computeChemisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents, size_t maxChemisorptionSites,
    std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH);

void computePhysisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot);

void computeChemisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t maxChemisorptionSites, double externalTemperature,
                          const Geometry& geometry, double particleDensity,
                          std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
                          std::span<const double> chemisorption, std::span<double> chemisorptionDot,
                          std::span<const double> poreConcentration, std::span<const double> solidTemperature);
